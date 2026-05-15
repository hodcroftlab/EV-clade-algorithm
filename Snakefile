"""
Snakefile-style workflow for clade suggestion across multiple viruses.

Suggested structural + commenting changes:
- Put *all* user-editable settings in config.yaml (avoid editing this file).
- Split the workflow into clearly marked sections:
  1) config file
  2) helper functions
  3) target rules (rule all)
  4) tree acquisition
  5) defaults generation (weights, aliases)
  6) analysis steps (optimal scales, suggest clades)
  7) optional reporting / visualization
  8) housekeeping (clean)

"""
import json
import os
import re
from pathlib import Path

############################
# 1) CONFIG / GLOBALS
############################

configfile: "config.yaml"

# Viruses configured in config.yaml under `viruses:`
VIRUS_CONFIG = config["viruses"]
VIRUSES = sorted(VIRUS_CONFIG.keys())

wildcard_constraints:
    virus="|".join(VIRUSES)

############################
# 2) PATH HELPERS
############################

def load_json(path: str | Path) -> dict:
    with open(path) as fh:
        return json.load(fh)

def load_scales(virus: str, filename: str = "optimal_scales.json") -> dict:
    """Load per-virus optimization output from `{virus}/results/optimal_scales.json`."""
    return load_json( Path(virus) / "results" / filename)

def clade_key_arg(virus: str) -> str:
    """
    Influenza-specific clade handling (subclade key).
    For enteroviruses you likely want "".
    """
    key = config["defaults"].get("clade_key", "")
    if re.search(r"(NA|HA)", virus):
        key = "subclade"
    return f"--clade-key {key}"

def gff_arg_if_present(virus: str) -> str:
    """Return '--gff ...' only if annotation exists, else empty string."""
    gff = Path(virus) / "resources" / "genome_annotation.gff3"
    return f"--gff {gff}" if gff.exists() else ""

############################
# 3) WORKFLOW TARGETS
############################

rule all: 
    """
    Final products:
    - original tree (fetched or provided)
    - optimal scales (parameter optimization result)
    - suggested tree (clade suggestions added)
    """
    input:
        expand("auspice/original_{virus}.json", virus=VIRUSES),
        expand("{virus}/results/optimal_scales.json", virus=VIRUSES),
        expand("auspice/suggested_{virus}.json", virus=VIRUSES),

rule viz:
    """
    Visualize Auspice Trees. Link: http://localhost:4000/suggested/{virus}?branchLabel=new-clade&c=new-clade&d=tree&p=full&showBranchLabels=all
    """
    shell:
        "auspice view --datasetDir auspice"

############################
# 4) GET TREE & ANNOTATION
############################

## provide  a tree in json format (previous Nextstrain run neccesary, template here: github.com/hodcroftlab/template_nextstrain)
TREE = None ### format: "auspice/original_{virus}.json"
skip_fetch = bool(TREE or config["defaults"].get("local_tree"))

## provide a gff3 annotation, e.g., with this script (github.com/hodcroftlab/template_nextstrain/blob/master/ingest/bin/generate_from_genbank.py)
ANNOTATION = None ### format: "{virus}/resources/genome_annotation.gff3" 

## if you provide a tree in json format, this step will be skipped.
if not skip_fetch:
    rule fetch:
        output:
            tree = "auspice/original_{virus}.json",
            annotation = "{virus}/resources/genome_annotation.gff3",
        params:
            virus = "{virus}",
            tree_url = lambda w: config["viruses"][w.virus]["tree_url"]
        log: 
            "logs/fetch.{virus}.log"
        shell:
            """
            mkdir -p logs {params.virus}/resources
            # Try Nextclade first
            echo "Attempting to fetch from Nextclade: enpen/enterovirus/{params.virus}"
            if nextclade dataset get --name "enpen/enterovirus/{params.virus}" --output-dir .nextclade_tmp 2>/dev/null; then
                if [[ -f ".nextclade_tmp/tree.json" ]]; then
                    echo "✓ Success: Got tree from Nextclade"
                    cp .nextclade_tmp/tree.json "{output.tree}"
                    if [[ -f ".nextclade_tmp/genome_annotation.gff3" ]]; then
                        cp .nextclade_tmp/genome_annotation.gff3 {output.annotation}
                        echo "✓ Got annotation from Nextclade"
                    fi
                    rm -rf .nextclade_tmp
                    exit 0
                fi
            fi
                
            # Fallback to Nextstrain
            echo "✗ Nextclade not available, falling back to Nextstrain"
            curl "{params.tree_url}" -o "{output.tree}" --max-time 10 -A "Mozilla/5.0"
            if [[ ! -f "{output.tree}" ]]; then
                echo "✗ Failed to fetch tree from Nextstrain"
                exit 1
            fi
            echo "✓ Success: Got tree from Nextstrain"

            # Extract genome_annotations from auspice.json and convert to GFF3
            seqid=$(jq -r '.meta.genome_annotations.nuc.seqid' "{output.tree}" | sed 's|.*/||;s|\.gb||')
            seq_start=$(jq -r '.meta.genome_annotations.nuc.start' "{output.tree}")
            seq_end=$(jq -r '.meta.genome_annotations.nuc.end' "{output.tree}")
            
            (
                echo "##gff-version 3"
                echo "#!gff-spec-version 1.21"
                echo "##sequence-region $seqid $seq_start $seq_end"
                jq -r '.meta.genome_annotations | to_entries[] | 
                "\\(.value.seqid)\\t.\\t\\(.value.type)\\t\\(.value.start)\\t\\(.value.end)\\t.\\t\\(.value.strand)\\t.\\tID=\\(.key);Name=\\(.key)"' \
                "{output.tree}"
            ) > "{output.annotation}"
                
            echo "✓ Generated GFF3 from tree metadata"
            """

rule ensure_clades:
    """
    Ensure the input tree has *some* clade annotation.

    Why:
    - If a user provides their own auspice JSON, it may not have clades/labels yet.
    - `add_new_clades.py` expects an "old clade key" to exist (so it can build hierarchy,
      label backbone, etc.). For EV we typically use `clade_membership` (or whatever your
      script uses as old_clade_key).

    What this does:
    - If the tree already contains clade annotations, do nothing (copy through).
    - Otherwise, assign a single dummy clade ("XXX") to all nodes (and optionally label root).
    """
    input:
        tree="auspice/original_{virus}.json",
    output:
        tree="{virus}/results/original_{virus}.clades.json",
    params:
        clade_key=config["defaults"]["clade_key"],
        dummy="XX0",
    log:
        "logs/ensure_clades.{virus}.log"
    run:
        import json

        def tree_has_clades(n):
            # True if *any* node has node_attrs[clade_key]["value"]
            attrs = n.get("node_attrs", {})
            if params.clade_key in attrs and isinstance(attrs[params.clade_key], dict) and "value" in attrs[params.clade_key]:
                return True
            for c in n.get("children", []):
                if tree_has_clades(c):
                    return True
            return False

        def assign_dummy(n):
            n.setdefault("node_attrs", {})
            n["node_attrs"].setdefault(params.clade_key, {"value": params.dummy})
            for c in n.get("children", []):
                assign_dummy(c)

        with open(input.tree) as fh:
            data = json.load(fh)

        T = data.get("tree", data)  # tolerate both raw tree dict or auspice JSON
        if not tree_has_clades(T):
            assign_dummy(T)

        # if it was a raw tree dict, wrap it back? keep original structure:
        if "tree" in data:
            data["tree"] = T
            out = data
        else:
            out = T

        with open(output.tree, "w") as fh:
            json.dump(out, fh)
        

############################
# 5) DEFAULT INPUTS (WEIGHTS / ALIASES)
############################

rule setup_default_weights:
    """
    Generate simple default weights if user hasn't provided any.

    Convention:
    - weights.json is a dict keyed by gene/protein name, with:
      - "default": weight applied to any position not explicitly listed
      - optionally position-specific weights, e.g. {"VP1": {"default": 1, 95: 3}}
    """
    message: "Set up default weights for amino acid positions in specified proteins. Written to {output}"
    output:
        "{virus}/resources/weights.json"
    params:
        default_weights = lambda w: config["viruses"][w.virus].get("weights", config["defaults"].get("weights", {"VP1": {"default": 1}}))
    shell:
        """
        mkdir -p $(dirname {output}) && python3 -c \"import json; json.dump({params.default_weights}, open('{output}', 'w'))\"
        """

rule setup_aliases:
    """
    Extract current alias mapping from the existing tree annotation. Clade names must be present in the tree metadata for this to work.
    This is useful as a "starting point" if you don't have curated aliases yet.
    """
    message: "Extract aliases from tree and save to {output}"
    input:
        tree = "{virus}/results/original_{virus}.clades.json",
    output:
        aliases = "{virus}/resources/aliases_default.json",
    params:
        clade_key=lambda w: clade_key_arg(w.virus),
    shell:
        """
        python3 scripts/extract_aliases_from_tree.py \
            --tree {input.tree} \
            {params.clade_key} \
            --output {output.aliases}
        """

############################
# 6) PARAMETER OPTIMIZATION
############################

rule calculate_optimal_scales:
    """
    Fit/choose (some) scales used by the algorithm, saved per virus.

    Output is later used to populate sweep defaults:
    - bushiness_branch_scale
    - branch_length_scale
    - divergence_scale
    """
    input: 
        tree = "{virus}/results/original_{virus}.clades.json",                    ## auspice.json tree
        weights = "{virus}/resources/weights.json",       ## weights: if some mutations are more important than others (e.g., epitopes) put them in here
    output:
        scales = "{virus}/results/optimal_scales.json"  ## gives you the optimal parameters to run the algorithm with
    params:
        clade_key=lambda w: clade_key_arg(w.virus),
        show_recommendations = "True",
    log: 
        "logs/optimal_scales.{virus}.log"
    shell:
        """
        mkdir -p logs
        python3 scripts/calculate_optimal_scales.py \
            --tree {input.tree} \
            {params.clade_key} \
            --recommend {params.show_recommendations} \
            --weights {input.weights} \
            --output {output.scales} \
            |& tee -a {log}
        """

############################
# 7) CLADE SUGGESTION
############################
## This rule assumes that clades defined in the tree
rule suggest_new_clades:
    input: # define config paths
        tree = "{virus}/results/original_{virus}.clades.json",
        aliases = lambda w: f"{w.virus}/resources/aliases_{config['defaults']['alias']}.json",
        weights = "{virus}/resources/weights.json",
        optimal = "{virus}/results/optimal_scales.json", 
        config= "config.yaml", 
    output:
        tree = "auspice/suggested_{virus}.json",
        clade_file = "{virus}/results/new-clades.tsv",
    params:
        virus = "{virus}",
        # Optional GFF
        gff=lambda w: gff_arg_if_present(w.virus),

        # Clade key (flu-specific, skip for enteroviruses)
        clade_key=lambda w: clade_key_arg(w.virus),
        defaults = config["defaults"],

        ## params for the plots
        plots = "True",  # "True" to generate plots
        p_cutoff = config["sweep"]["cutoff"],
        p_div_add = config["sweep"]["divergence_addition"],
        p_min_size = config["sweep"]["min_size"],
        p_bush_scale = lambda w: load_scales(w.virus, "optimal_scales.json")["optimal_scales"].get("bushiness_branch_scale", config["sweep"]["bushiness_branch_scale"]), 
        p_bls_range = lambda w: load_scales(w.virus, "optimal_scales.json")["optimal_scales"].get("branch_length_scale", config["sweep"]["branch_length_scale"]), 
        p_div_scale = lambda w: load_scales(w.virus, "optimal_scales.json")["optimal_scales"].get("divergence_scale", config["sweep"]["divergence_scale"]), 
    
    shell:
        """
        if [ -n "{params.gff}" ]; then 
            echo "{params.gff}"; 
        else 
            echo "No GFF file provided - Attention the divergence scale will be calculated based on nucleotides only"; 
        fi
        python3 scripts/add_new_clades.py \
            --virus {params.virus} \
            --tree {input.tree} \
            --config {input.config} \
            --optimal-scales {input.optimal} \
            --aliases {input.aliases} \
            --weights {input.weights} \
            --output {output.tree} \
            {params.clade_key} \
            {params.gff} \
            \
            --plots {params.plots} \
            --clades {output.clade_file} \
            --cutoff-sweep {params.p_cutoff} \
            --div_add-sweep {params.p_div_add} \
            --div_scale-sweep {params.p_div_scale} \
            --min_size-sweep {params.p_min_size} \
            --bush_scale-sweep {params.p_bush_scale} \
            --bls_range-sweep {params.p_bls_range}
        
        """ #

        # visualize: http://localhost:4000/suggested/D68?branchLabel=new-clade&c=new-clade&d=tree&p=full&tl=__strain__

############################
# 8) OPTIONAL POST-PROCESSING
############################

rule visualization:
    input:
        clades = "{virus}/results/new-clades.tsv",
        tree = rules.suggest_new_clades.output.tree,
        config = "config.yaml",
    params:
        color = "True",
        grid = "{virus}/results/plots/clade_table_names.tiff",
        d = "{virus}/results/plots/",
        v = "{virus}"
    output:
        tiff = "{virus}/results/plots/violin_plots.tiff",
    shell:
        """
        Rscript scripts/visualize_clades.R \
          --clades {input.clades} \
          --tree {input.tree} \
          --dir-out {params.d} \
          --config {input.config} \
          --virus {params.v} \
          --out {output.tiff} \
          --grid {params.grid} \
          --color {params.color}        
        """


rule extract_clade_ymls:
    input:
        json = rules.suggest_new_clades.output.tree,
    output:
        clades = directory("{virus}/clades"),
    shell:
        """
        python3 scripts/extract_yml_from_json.py \
            --json {input.json} \
            --outdir {output.clades}
        """

# generate the markdown summary of the clade definitions
# rule summary:
#     input:
#         expand("clades/{file}", file=glob_w("clades/{file}").file),
#         clades = rules.extract_clade_ymls.output.clades,
#     output:
#         md = "{virus}/.auto-generated/clades_{virus}.md", #".auto-generated/subclades_D68.md"
        
#     shell:
#         """
#         python3 scripts/generate_markdown_summary.py \
#             --input-dir {input.clades} --virus {w.virus} --file {output.md}
        # """

rule construct_tsv:
    input:
        "{virus}/clades"
    output:
        tsv = "{virus}/.auto-generated/subclades.tsv",
        ptsv = "{virus}/.auto-generated/subclade-proposals.tsv",
        clades = "{virus}/.auto-generated/clades.tsv",
        clades_long = "{virus}/.auto-generated/clades-long.tsv"
    shell:
        """
        python3 scripts/construct_tsv.py --input-dir {input} --output-tsv {output.tsv}
        python3 scripts/construct_tsv.py --input-dir {input} subclade-proposals --output-tsv {output.ptsv}
        python3 scripts/construct_tsv.py --input-dir {input} --aux-input-dir subclades --flat-output --use-short-name --output-tsv {output.clades}
        python3 scripts/construct_tsv.py --input-dir {input} --aux-input-dir subclades --flat-output --output-tsv {output.clades_long}
        """


############################
# 9) HOUSEKEEPING
############################
rule clean:
    shell:
        """
        rm -r .*/.auto-generated \
        .*/clades/ \
        .*/subclades.tex \
        auspice/*
        """