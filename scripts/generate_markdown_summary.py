#! python
#%%
import yaml
import glob

def generate_lineage_md(subclade, lineage, gene):
    lines = []
    revoked = subclade.get('revoked', False)
    if revoked:
        lines.append(f"## ~~{subclade['name']}~~ (revoked)")
    else:
        lines.append(f"## {subclade['name']}")

    lines.append(f" * parent: [{subclade['parent']}](#{subclade['parent'].replace('.', '')})")
    if 'comment' in subclade and subclade['comment']:
        lines.append(f" * comment: {subclade['comment']}")
    if subclade.get('defining_mutations'):
        snp_str = ', '.join(f"{x['locus']}:{x['position']}{x['state']}" for x in subclade['defining_mutations'])
    else:
        snp_str = "None listed"    
    lines.append(f" * defining mutations or substitutions: {snp_str}")
    if "clade" in subclade and subclade['clade'] != "none":
        lines.append(f" * clade: {subclade['clade']}")

    ref_seqs = []
    for x in subclade['representatives']:
        nextstrain_link = f"[View on Nextstrain](https://nextstrain.org/groups/hodcroftlab/enterovirus/{lineage}/{gene}?branchLabel=new-clade&c=new-clade&label=new-clade&s={x['isolate']})"
        if x['source']=='genbank' and 'accession' in x:
            accession_link = f"[{x['accession']}](https://www.ncbi.nlm.nih.gov/nuccore/{x['accession']})"
        elif x['source']=='gisaid':
            accession_link = x['accession']

        if 'other_accession' in x:
            other_accession = f", {x['other_accession']}"
        else: other_accession=''
        ref_seqs.append(f"{x.get('isolate', x['accession'])} ({accession_link}{other_accession}) {nextstrain_link}")

    if len(ref_seqs)==1:
        lines.append(f" * representative sequence: {ref_seqs[0]}")
    elif len(ref_seqs)>1:
        lines.append(f" * representative sequences:")
        for r in ref_seqs:
            lines.append(f"   - {r}")
    lines.append(f" * [View on Nextstrain](https://nextstrain.org/groups/hodcroftlab/enterovirus/d68/genome?branchLabel=new-clade&c=new-clade&label=new-clade:{subclade['name']})")
    return '\n'.join(lines) + '\n'

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', required=True)
    parser.add_argument('--lineage')
    parser.add_argument('--file', required=True)
    args = parser.parse_args()

    subclades = []
    # Iterate through all lineage definition files
    for yaml_file in sorted(glob.glob(f"{args.input_dir}/*.yml")):
        with open(yaml_file, 'r') as stream:
            yaml_data = yaml.safe_load(stream)
        subclades.append(yaml_data)

    subclades.sort(key=lambda x:x['name'])
    clade_lineage_map = [(x['parent'], x['name'], x['unaliased_name'])
                         for x in subclades if 'name' in x and x['name'] != 'none' and (not x.get('revoked', False))]
    # Write to json file
    with open(args.file, 'w') as outfile:
        outfile.write("# Summary of designated subclades\n")

        for subclade in subclades:
            print("output clade", subclade['name'])
            outfile.write(generate_lineage_md(subclade, args.lineage, "genome") + '\n')
            # ipdb.set_trace()


        if len(clade_lineage_map):
            # write table of clade -- subclade correspondence
            outfile.write("# Clade -- subclade correspondence\n")
            outfile.write(f"|*Clade*|*Subclade*|*full subclade name*|\n")
            outfile.write(f"|-------------|---------|----------------------|\n")
            for clade, lineage, unaliased_name in clade_lineage_map:
                outfile.write(f"|{clade}|[{lineage}](#{lineage.replace('.','')})|{unaliased_name}|\n")


    with open('subclades.tex', 'w') as latexoutfile:
        latexoutfile.write(f"Subclade & Clade & full subclade name\\\\\\hline\n")

        for subclade in subclades:
            latexoutfile.write(f"{subclade['name']} & {subclade.get('clade','')} & {subclade.get('unaliased_name', '')}\\\\\n")
