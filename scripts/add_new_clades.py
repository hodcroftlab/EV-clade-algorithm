import json, argparse
from collections import defaultdict
import numpy as np
from BCBio import GFF  # Requires Biopython
import os
import pandas as pd
import yaml


def is_nonterminal(n):
    return "children" in n and len(n["children"])

def is_terminal(n):
    return not is_nonterminal(n)

def prepare_tree(n, max_date, min_date):
    '''
    Assigns a node attribute "alive" to 1 if the node is within the date rate, 0 otherwise.
    '''
    if is_nonterminal(n):
        tmp_alive = False
        new_children = []
        for c in n["children"]:
            prepare_tree(c, max_date, min_date)
            if len(c['branch_attrs']['mutations'].get('nuc',[])) or is_terminal(c):
                new_children.append(c)
            else:
                new_children.extend(c['children'])

            tmp_alive = tmp_alive or c['alive']
        n['alive'] = tmp_alive
        n['children'] = new_children
    else:
        try:
            numdate = n['node_attrs']['num_date']['value']
            n['alive'] = 1 if numdate<max_date and numdate>min_date else 0
        except:
            n['alive'] = 1

    n['clade_break_point'] = False


def label_backbone(tree):
    '''
    This function labels all branches/nodes from the root to existing clade labels
    as 'back-bone'. This can be used to prevent introduction of additional clades
    along the backbone
    '''
    def label_backbone_recursive(n):
        on_backbone = False
        if is_nonterminal(n):
            for c in n["children"]:
                label_backbone_recursive(c)
                if c["backbone"]:
                    on_backbone = True

        if on_backbone or n["clade_break_point"]:
            n["backbone"] = True
        else:
            n["backbone"] = False

    label_backbone_recursive(tree)

def dealias(n, aliases, label):
    if label in n["node_attrs"]:
        levels = n["node_attrs"][label]["value"].split('.')
        if levels[0]=='unassigned':
            full_clade = '*'
        else:
            if levels[0] not in aliases:
                raise KeyError(f"Alias for '{levels[0]}' not found in aliases. Current aliases: {aliases.keys()}. Ensure that '{levels[0]}' is defined in the aliases dictionary.")
            full_clade = tuple(list(aliases[levels[0]]) + [int(x) for x in levels[1:]])
        n["node_attrs"][f"full_{label}"] = {"value":full_clade}
    if "children" in n:
        for c in n["children"]:
            dealias(c, aliases, label)


def get_existing_clade_labels(tree, key):
    '''
    Returns a list of existing clade labels in the tree.
    '''
    def add_clades(n, full_clades, key):
        if f"full_{key}" in n["node_attrs"]:
            full_clades.add(tuple(n["node_attrs"][f"full_{key}"]["value"]))
        if "labels" in n["branch_attrs"] and key in n["branch_attrs"]["labels"]:
            n["clade_break_point"]=True

        if "children" in n:
            for c in n["children"]:
                add_clades(c, full_clades, key)

    full_clades = set()
    add_clades(tree, full_clades, key)
    return full_clades


def assign_divergence(n, genes):
    '''
    Assigns a node attribute "div" to each node that counts the number of mutations
    since the root of the tree
    '''
    if is_nonterminal(n):
        for c in n["children"]:
            c['div'] = n['div']
            for gene in genes:
                bad_states = ['-', 'N'] if gene=='nuc' else ['X', '-']
                c['div'] += len([x for x in c['branch_attrs']['mutations'].get(gene,[])
                                if (x[0] not in bad_states) and (x[-1] not in bad_states)])
            assign_divergence(c, genes)


def calc_phylo_score(n, distance=None, ignore_backbone=False):
    '''
    Assigns a node attribute "bushiness" to each node that counts the number of downstream tips.
    This can be distanced similar to the LBI calculation.
    '''
    if is_nonterminal(n):
        n["bushiness"] = 0
        n["ntips"] = 0
        for c in n["children"]:
            calc_phylo_score(c, distance=distance, ignore_backbone=ignore_backbone)

            if (ignore_backbone and c["backbone"])==False:
                n["bushiness"] += c["bushiness"]*np.exp(-distance(c))
                n["bushiness"] += (1-np.exp(-distance(c))) if c['alive'] else 0

            n["ntips"] += c["ntips"]
        n["children"] = sorted(n["children"], key=lambda x:x['ntips'])
    else:
        n["bushiness"] = 1 if n['alive'] else 0
        n["ntips"] = 1

    n['node_attrs']["bushiness_raw"] = {'value': n["bushiness"]}


def calc_phylo_scale(T):
    '''
    Calculates the phylogenetic scale of the tree as the median bushiness of all internal nodes
    '''
    def collect_recursive(n, values):
        if "children" in n and len(n["children"])>0:
            for c in n["children"]:
                collect_recursive(c, values)
            if n['alive']: # only include alive nodes and skip terminals
                values.append(n["bushiness"])

    values = []
    collect_recursive(T, values)
    from scipy.stats import scoreatpercentile
    return scoreatpercentile(values,80)

def extract_protein_names_from_gff(gff_file):
    proteins = set()
    with open(gff_file) as fh:
        for rec in GFF.parse(fh):
            for feat in rec.features:
                if feat.type == "CDS":
                    proteins.add(feat.qualifiers.get("Name", ["CDS"])[0])
    return sorted(proteins)

def score(n, weights=None, bushiness_scale=1, ignore_backbone=False,
          proteins=None, branch_length_scale=4):
    '''
    Assign a score to each node that combines phylogenetic signal (bushiness) and protein-level mutation weights.
    If proteins are not supplied, they are inferred from a provided GFF3 file.    
    '''
    if weights is None:
        weights = {}

    if proteins is None:
        proteins = []

    score = n["bushiness"]/(n["bushiness"] + bushiness_scale)
    n['node_attrs']['bushiness'] = {'value': score}

    aa_weight = 0
    for cds in n['branch_attrs']['mutations'].keys():  #proteins: 
        if cds=='nuc': continue  # Skip nucleotide mutations
        w = weights.get(cds, {})
        for mut in n['branch_attrs']['mutations'].get(cds,[]):
            if mut[-1] in ['-', 'X']: continue # Skip deletions and unknown AAs
            pos = str(mut[1:-1]) # Extract position number; must be str otherwise not taken into account
            aa_weight += w[pos] if pos in w else w.get("default", 1)

    n['node_attrs']['branch_score'] = {'value': aa_weight/(branch_length_scale + aa_weight)}
    score += aa_weight/(branch_length_scale + aa_weight)

    # return 0 if the node is on the backbone and we are ignoring backbone nodes
    if ignore_backbone and n['backbone']:
        return 0.0
    # return 0 if there are no mutations in the proteins of interest
    # if sum([len(n['branch_attrs']['mutations'].get(cds,[])) for cds in proteins])==0:
    #     return 0.0

    return score


def assign_score(n, score=None, **kwargs):
    '''
    recursively assign the clade demarcation score to each branch
    '''
    if is_nonterminal(n):
        for c in n["children"]:
            assign_score(c, score=score, **kwargs)
    n['node_attrs']["score"] = {'value': score(n, **kwargs)}


def assign_clade(n, clade, key):
    '''Assign a clade to a node and recursively to all its children'''
    if is_nonterminal(n):
        for c in n["children"]:
            assign_clade(c, clade, key)
    n['node_attrs'][key] = {'value': clade}


def assign_new_clades_to_branches(n, hierarchy, new_key, new_clades=None,
                                  cutoff=1.0, divergence_addition=None, divergence_base=0,
                                  divergence_scale=4, min_size=5):
    '''
    walk through the tree in pre-order (by recursively calling this function)
    and call a new clade whenever there is a branch that crosses the threshold
    '''
    if divergence_addition:
        delta_div = n['div']-divergence_base  # calculate the divergence since the parent
        div_score = divergence_addition*delta_div/(delta_div+divergence_scale)
    else: div_score=0
    n["node_attrs"]['div_score'] = {'value': div_score}
    n["node_attrs"]['div_impact'] = {'value': (cutoff - (n["node_attrs"]['score']['value']))}  # add the divergence score to the clade score

    n["node_attrs"]['score']['value'] += div_score  # add the divergence score to the clade score

    # trigger = (n["node_attrs"]['score']['value'] + div_score > cutoff) and (n["ntips"]>min_size)
    trigger = (n["node_attrs"]['score']['value'] > cutoff) and (n["ntips"]>min_size)
    if trigger:
        potential_div = n['div']
        child_triggers = [((c["node_attrs"]['score']['value'] + (c['div']-potential_div)/((c['div']-potential_div) +divergence_scale) > cutoff) and (c["ntips"]>min_size)) for c in n["children"]]
        one_daughter = sum(child_triggers)==1

    if trigger and (not one_daughter):
        if 'labels' not in n['branch_attrs']:
            n['branch_attrs']['labels'] = {}

        # determine parent clade
        parent_clade = tuple(n['node_attrs'][f"full_{new_key}"]["value"])
        if parent_clade in hierarchy:
            # determine the number of existing children of the parent and the index of the new subclade
            if len(hierarchy[parent_clade]):
                new_suffix = max(hierarchy[parent_clade])+1
            else:
                new_suffix = 1
            # print(parent_clade, hierarchy[parent_clade], new_suffix)
            if new_suffix>2: # consistency check
                assert new_suffix==hierarchy[parent_clade][-1]+1

            hierarchy[parent_clade].append(new_suffix)
            new_clade = tuple(list(parent_clade) + [new_suffix])
        else:
            new_clade = parent_clade
            print("new clade found, but no parent to attach it to")

        hierarchy[new_clade] = []
        new_clades[new_clade] = n
        assign_clade(n, new_clade, f"full_{new_key}")
        n["clade_break_point"] = True  # mark as clade break_point

    # reset divergence to clade break point.
    if n['clade_break_point']:
        # print(n['div'] - divergence_base, n['branch_attrs'].get('labels',{}))
        divergence_base=n['div']

    if 'children' in n:
        for c in n["children"]:
            assign_new_clades_to_branches(c, hierarchy, new_key,
                                new_clades=new_clades, cutoff=cutoff,
                                divergence_addition=divergence_addition,
                                divergence_base=divergence_base,
                                divergence_scale=divergence_scale,
                                min_size=min_size)


def copy_over_old_clades(tree, old_key, new_key, branch_label):
    def copy_recursive(n, old_key, new_key):
        n["node_attrs"][new_key] = {k:v for k,v in n["node_attrs"][old_key].items()}
        n["node_attrs"][f"full_{new_key}"] = {k:v for k,v in n["node_attrs"][f"full_{old_key}"].items()}
        n['clade_break_point'] = False
        if branch_label in n["branch_attrs"].get("labels",{}):
            n["branch_attrs"]["labels"][new_key] = n["branch_attrs"]["labels"][branch_label]
            n['clade_break_point'] = True

        if is_nonterminal(n):
            for c in n["children"]:
                copy_recursive(c, old_key, new_key)

    copy_recursive(tree, old_key, new_key)

def full_clade_to_short_name(full_clade, aliases):
    # Track the highest letter/number used for each level
    used_names = {}
    for key in aliases.keys():
        top_level = key[0]  # Get the first element of the tuple
        if top_level not in used_names:
            used_names[top_level] = key[1] if len(key)>1 else 0
        else:
            used_names[top_level] = max(used_names[top_level], key[1] if len(key)>1 else 0)
    # Remove unused clades
    used_names = {key: value for key, value in used_names.items() if key != "Unassigned"}
    max_letter = max(used_names.keys(), key=lambda k: ord(k[0]))

    def get_next_name(clade, parent):
        """Generate the next name for a given level."""
        # ipdb.set_trace()
        base = clade[0]
        level = clade[1]
        if level < 3:
            if base in used_names and not used_names[base] == level+1:
                # If the next level is not already used, increment the number           
                return base+str(level+1)
            else:
                if base not in used_names:
                    return chr(ord(base) + 1)
                else:
                     return chr(ord(max_letter) + 1)
        else:
                if base not in used_names:
                     return  chr(ord(base) + 1)
                else:
                     return chr(ord(max_letter) + 1)            

    for full_parent in sorted(aliases.keys(), key=lambda x:len(x), reverse=True):
        if tuple(full_clade[:len(full_parent)])==full_parent:
            clade_name = aliases[full_parent]
            # used_names[full_clade[0]]=len(full_parent)
            # # If the full clade has more than 4 levels, generate a new clade name
            # if len(full_clade)> 4:
            #     # Generate a new name for the 5th level and beyond
            #     clade_name = get_next_name(full_clade, full_parent)
            #     aliases[tuple(full_clade)] = clade_name  # Update aliases with the new name

            if len(full_parent)<len(full_clade):
                clade_name += '.' + '.'.join([str(x) for x in full_clade[len(full_parent):]])
            return clade_name

    return '.'.join([str(x) for x in full_clade])

def get_clade_map(fname):
    if fname is None:
        return {}, {}

    aliases = {}
    with open(fname) as fh:
        old_to_new_clades = json.load(fh)
    for k,v in old_to_new_clades.items():
        if len(v[0])==1:
            aliases[tuple(v[1])] = v[0]
    short_to_full_clades = {v[0]:v[1] for v in old_to_new_clades.values()}

    return short_to_full_clades, aliases

def add_aliases(fname, aliases):
    with open(fname) as fh:
        new_aliases = json.load(fh)
    for k,v in new_aliases.items():
        aliases[tuple(v)] = k

def get_tree(T, max_date=None, min_date=None, old_key=None, new_key=None, branch_label=None, proteins=['nuc']):
    prepare_tree(T, max_date=max_date, min_date=min_date)
    T['div']=0
    # assign_divergence(T, ['nuc'])
    if not proteins:
        proteins = ['nuc']
    assign_divergence(T, proteins)

    hierarchy = defaultdict(list)
    # print("adding existing clades")
    existing_full_clades = get_existing_clade_labels(T, old_key)
    for full_clade in existing_full_clades:
        if len(full_clade)>1:
            hierarchy[full_clade[:-1]].append(full_clade[-1])
        if full_clade not in hierarchy:
            hierarchy[full_clade] = []
    copy_over_old_clades(T, old_key, new_key, branch_label=branch_label)
    label_backbone(T)

    hierarchy = {k:sorted(hierarchy[k]) for k in sorted(hierarchy.keys())}

    return T, hierarchy

def extract_node_stats(root):
    values = {
        "score": [],
        "bushiness": [],
        "div": [],
        "branch_score": []
    }

    def walk(node):
        attrs = node.get("node_attrs", {})

        if "score" in attrs:
            values["score"].append(attrs["score"]["value"])
        if "bushiness" in attrs:
            values["bushiness"].append(attrs["bushiness"]["value"])
        if "div" in attrs:
            values["div"].append(attrs["div"])
        if "branch_score" in attrs:
            values["branch_score"].append(attrs["branch_score"]["value"])

        for child in node.get("children", []):
            walk(child)

    walk(root)
    return pd.Series({k: sum(v) / len(v) for k, v in values.items() if v})


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Assign clades to a tree")
    parser.add_argument('--tree', type=str, required=True, help="JSON file with tree")
    parser.add_argument('--config', type=str, required=True, help="config.yaml")
    parser.add_argument('--aliases', type=str, help="JSON file with aliases")
    parser.add_argument('--weights', type=str, required=True, help="JSON file with weights")
    parser.add_argument('--output', type=str, required=True, help="Output JSON tree")
    parser.add_argument('--clades', type=str, help="Output TSV with clade counts")
    parser.add_argument('--gff', type=str, default=None, help="GFF3 file for protein extraction")
    parser.add_argument('--clade-key', type=str, default="subclade", help="Key for existing clades in node attributes")
    parser.add_argument('--virus', type=str, help="Virus name for aliasing and config")
        
    # Parameter sweep (optional, for --plots mode)
    parser.add_argument('--plots', type=str, default="False", help="Run parameter sweep")
    parser.add_argument('--cutoff-sweep', nargs='+', type=float, help="Cutoff range for sweep")
    parser.add_argument('--div_add-sweep', nargs='+', type=float, help="Div_add range for sweep")
    parser.add_argument('--div_scale-sweep', nargs='+', type=float, help="Div_scale range for sweep")
    parser.add_argument('--min_size-sweep', nargs='+', type=int, help="Min_size range for sweep")
    parser.add_argument('--bush_scale-sweep', nargs='+', type=float, help="Bush_scale range for sweep")
    parser.add_argument('--bls_range-sweep', nargs='+', type=float, help="Bls_range range for sweep")

    args = parser.parse_args()

    new_clade_key = 'new-clade'
    old_clade_key = args.clade_key
    branch_label = 'clade'
    aliases = {}
    gff = args.gff

    with open("config.yaml") as fh:
        config = yaml.safe_load(fh)

    cutoff = config["viruses"][args.virus]["cutoff"]
    div_add = config["viruses"][args.virus]["divergence_addition"]
    div_scale = config["viruses"][args.virus]["divergence_scale"]
    min_size = config["viruses"][args.virus]["min_size"]
    bush_scale = config["viruses"][args.virus]["bushiness_branch_scale"]
    bls_range = config["viruses"][args.virus]["branch_length_scale"]

    if args.aliases:
        add_aliases(args.aliases, aliases)
    reverse_aliases = {v:k for k,v in aliases.items()}
    

    with open(args.weights) as fh:
        weights = json.load(fh)

    with open(args.tree) as fh:
        data = json.load(fh)

    proteins = config["defaults"]["proteins"]
    if not proteins:
        # ipdb.set_trace()
        if args.gff and os.path.isfile(args.gff):
            proteins = extract_protein_names_from_gff(args.gff)
            print("Proteins from GFF:", proteins)
        else:
            proteins = []

    # ipdb.set_trace()
    T = data["tree"]
    dealias(T, reverse_aliases, old_clade_key)
    T, hierarchy = get_tree(T, max_date=config["defaults"]["max_date"], min_date=config["defaults"]["min_date"],
                            new_key=new_clade_key, old_key=old_clade_key, branch_label=branch_label,
                            proteins=proteins)

    # nucleotide branch length excluding gaps and N
    branch_length_function = lambda x:len([y for y in x['branch_attrs']['mutations'].get('nuc',[])
                                           if y[-1] not in ['N', '-'] and y[0] not in ['N', '-']])/bush_scale
    calc_phylo_score(T, branch_length_function, ignore_backbone=True)
    bushiness_scale = calc_phylo_scale(T)
    print("phylo_score_scale", bushiness_scale)

    # compute aggregate score from branches and phylo/bushiness
    assign_score(T, score, weights=weights,
                 bushiness_scale=bushiness_scale, ignore_backbone=True,
                 proteins=proteins, branch_length_scale=bls_range)

    # assign clades while also taking into account the divergence, modifies hierarchy and new_clades in place
    new_clades = {}
    assign_new_clades_to_branches(T, hierarchy, new_clade_key,
        new_clades=new_clades, cutoff=cutoff, divergence_addition=div_add,
        divergence_base=0.0, divergence_scale=div_scale, min_size=min_size)

    # # copy over old clades to count
    count_old = set()
    count_new = set()
    count_old= count_old.union(aliases.values())

    # process and assign human readable clade names
    for new_clade in new_clades:
        clade_name = full_clade_to_short_name(new_clade, aliases)
        n = new_clades[new_clade]
        if "labels" not in n["branch_attrs"]: n["branch_attrs"]["labels"] = {}
        n["branch_attrs"]["labels"][new_clade_key] = clade_name
        assign_clade(n, clade_name, new_clade_key)
        print("suggested clade:", clade_name,
              {k:v for k, v in n["branch_attrs"]["mutations"].items() if k!='nuc'})

        count_new.add(clade_name)
    print("\nnumber of new clades:", len(count_new))

    # count_old = sorted(list(count_old))
    # count_new = sorted(list(count_new))

    # count_old.extend([np.nan] * (len(count_new)-len(count_old)))

    # count = pd.DataFrame({'old_clade': list(count_old), 'new_clade': list(count_new)})
    # count.to_csv(args.clades, index=False, header=True, sep='\t')

    # export
    data['meta']['colorings'].append({'key':new_clade_key, 'type':'ordinal', 'title':new_clade_key})
    data['meta']['colorings'].append({'key':"score", 'type':'continuous', 'title':"Clade score"})
    data['meta']['colorings'].append({'key':"div_score", 'type':'continuous', 'title':"Divergence score"})
    data['meta']['colorings'].append({'key':"div_impact", 'type':'continuous', 'title':"Divergence impact"})
    data['meta']['colorings'].append({'key':"bushiness_raw", 'type':'continuous', 'title':"Bushiness (raw)"}) #  pylo_score_raw
    data['meta']['colorings'].append({'key':"bushiness", 'type':'continuous', 'title':"Bushiness score"}) #  pylo_score
    data['meta']['colorings'].append({'key':"branch_score", 'type':'continuous', 'title':"Branch score"})
    data['meta']['filters'].append({'key':new_clade_key})

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)

    # --------- EXTRA PARAMETER SWEEP ---------
    if args.plots == "True":
        import itertools, copy
        from itertools import product
        from tqdm import tqdm

        all_dfs = []

        param_combos = list(product(
            args.cutoff_sweep,
            args.div_add_sweep,
            args.div_scale_sweep,
            args.bls_range_sweep,
            args.min_size_sweep,
            args.bush_scale_sweep
        ))
        # ipdb.set_trace()

        for cutoff, div_add, div, bls, size, bush in tqdm(param_combos, desc="Sweeping parameters"):
            cfg = copy.deepcopy(config["viruses"][args.virus])
            cfg["cutoff"] = cutoff
            cfg["divergence_addition"] = div_add
            cfg["divergence_scale"] = div
            cfg["branch_length_scale"] = bls
            cfg["min_size"] = size
            cfg["bushiness_branch_scale"] = bush

            T_tmp = copy.deepcopy(data["tree"])
            dealias(T_tmp, reverse_aliases, old_clade_key)
            T_tmp, hierarchy = get_tree(T_tmp, max_date=cfg["max_date"], min_date=cfg["min_date"],
                                        new_key=new_clade_key, old_key=old_clade_key, branch_label=branch_label)

            # nucleotide branch length excluding gaps and N
            branch_length_function = lambda x:len([y for y in x['branch_attrs']['mutations'].get('nuc',[])
                                                if y[-1] not in ['N', '-'] and y[0] not in ['N', '-']])/cfg["bushiness_branch_scale"]

            calc_phylo_score(T_tmp, branch_length_function, ignore_backbone=True)
            bushiness_scale = calc_phylo_scale(T_tmp)

            assign_score(T_tmp, score, weights=weights,
                         bushiness_scale=bushiness_scale, ignore_backbone=True,
                         proteins=proteins, branch_length_scale=cfg["branch_length_scale"])

            new_clades = {}
            assign_new_clades_to_branches(T_tmp, hierarchy, new_clade_key,
                new_clades=new_clades, cutoff=cfg["cutoff"], divergence_addition=cfg["divergence_addition"],
                divergence_base=0.0, divergence_scale=cfg["divergence_scale"], min_size=cfg["min_size"])

            count_old = len(aliases.values())
            count_new = 0

            for new_clade in new_clades:
                clade_name = full_clade_to_short_name(new_clade, aliases)
                count_new += 1

            # max_len = max(len(count_old), len(count_new))
            # count_old += [np.nan] * (max_len - len(count_old))
            # count_new += [np.nan] * (max_len - len(count_new))

            # ipdb.set_trace()

           # Append to list
            df = pd.DataFrame([{
                'old_clade': count_old,
                'new_clade': count_new,
                'cutoff': cutoff,
                'divergence_addition': div_add,
                'divergence_scale': cfg["divergence_scale"],
                'min_size': cfg["min_size"],
                'bushiness_branch_scale': cfg["bushiness_branch_scale"],
                'branch_length_scale': cfg["branch_length_scale"]                
            }])
            stats = extract_node_stats(T_tmp)
            df = pd.concat([df, stats.to_frame().T], axis=1)
            all_dfs.append(df)

            # print(all_dfs)
            # ipdb.set_trace()

            # print(len(count_new), "new clades")

        final_df = pd.concat(all_dfs, ignore_index=True)
        final_df.to_csv(str.replace(args.clades, ".tsv", "_all.tsv"), sep='\t', index=False)
