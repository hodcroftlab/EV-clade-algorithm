import numpy as np
import json
import os
import ipdb

# Load config once
configfile: "config.yaml"

VIRUSES = config["viruses"]

wildcard_constraints:
    virus=VIRUSES


# ! adjust the `tree_url` in suggestion_params.json
# ! define some weights (e.g. BC, DE loop in VP1) in weights.json
# ! define the current nomenclature in aliases.json

rule all: # order: fetch > calculate_optimal_scales > suggest_new_clades
    input:
        expand("auspice/original_{virus}.json", virus=VIRUSES),
        expand("{virus}/results/optimal_scales.json", virus=VIRUSES),
        expand("auspice/suggested_{virus}.json", virus=VIRUSES),
        # expand("{virus}/clades", virus=VIRUSES),
        # expand("{virus}/.auto-generated/clades_{virus}.md", virus=VIRUSES)

ruleorder: fetch > calculate_optimal_scales > suggest_new_clades


## if you provide a tree in json format, this step will be skipped.
rule fetch:
    output:
        tree = "auspice/original_{virus}.json",
        annotation = "{virus}/resources/genome_annotation.gff3",
    params:
        virus = "{virus}",
        tree_url = lambda w: config["viruses"][w.virus]["tree_url"]
    log: "logs/fetch.{virus}.log"
    shell:
        """
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

rule setup_default_weights:
    output:
        "{virus}/resources/weights.json"
    params:
        default_weights = {"VP1": {"default": 1}}
    shell:
        "mkdir -p $(dirname {output}) && python3 -c \"import json; json.dump({params.default_weights}, open('{output}', 'w'))\""

rule calculate_optimal_scales:
    input: 
        tree = "auspice/original_{virus}.json",                    ## auspice.json tree
        weights = "{virus}/resources/weights.json",       ## weights: if some mutations are more important than others (e.g., epitopes) put them in here
    output:
        scales = "{virus}/results/optimal_scales.json"  ## gives you the optimal parameters to run the algorithm with
    params:
        clade_key=lambda w: "--clade-key subclade" if re.search(r"(NA|HA)", w.virus) else ""
    log: "logs/optimal_scales.{virus}.log"
    shell:
        """
        python3 scripts/calculate_optimal_scales.py \
            --tree {input.tree} \
            {params.clade_key} \
            --weights {input.weights} \
            --output {output.scales} >> {log} 2>&1
        """

rule setup_aliases:
    input:
        tree = "auspice/original_{virus}.json",
    output:
        aliases = "{virus}/resources/aliases_default.json",
    params:
        clade_key = lambda w: "--clade-key subclade" if re.search(r"(NA|HA)", w.virus) else "",
    shell:
        """
        python3 scripts/extract_aliases_from_tree.py \
            --tree {input.tree} \
            {params.clade_key} \
            --output {output.aliases}
        """

rule suggest_new_clades:
    input: # define config paths
        tree = rules.fetch.output.tree,
        aliases = lambda w: f"{w.virus}/resources/aliases_{config['defaults']['alias']}.json",
        weights = "{virus}/resources/weights.json",
        params = "{virus}/resources/suggestion_params.json",
        optimal = rules.calculate_optimal_scales.output.scales,
        
    params:
        clade = "{virus}/new-clades.tsv",
        gff = lambda w: f"--gff .{w.virus}/config/genome_annotation.gff3" if os.path.exists(f".{w.virus}/config/genome_annotation.gff3") else "",    
        clade_key=lambda w: "--clade-key subclade" if re.search(r"(NA|HA)", w.virus) else "", # flu has different key for clades
        
        ## params for the plots
        plots = "False",  # "True" to generate plots
        cutoff = [round(x, 1) for x in np.arange(0.5, 1.5, 0.1)],
        div_add = [round(x, 2) for x in np.arange(0.5, 1.2, 0.2)],
        min_size = [round(x, 1) for x in np.arange(10, 40, 10)],
        bush_scale = lambda w: load_config(w.virus, "optimal_scales.json")["optimal_scales"].get("bushiness_branch_scale"), #[5,10,100],
        bls_range = lambda w: load_config(w.virus, "optimal_scales.json")["optimal_scales"].get("branch_length_scale"), #[2,10,20],
        div_scale = lambda w: load_config(w.virus, "optimal_scales.json")["optimal_scales"].get("divergence_scale"), #[0.1, 1, 5, 10, 20],
    
    output:
        tree = "auspice/suggested_{virus}.json"
    shell:
        """
        if [ -n "{params.gff}" ]; then 
            echo "{params.gff}"; 
        else 
            echo "No GFF file provided - Attention the divergence scale will be calculated based on nucleotides only"; 
        fi
        python3 scripts/add_new_clades.py \
            --tree {input.tree} \
            --config {input.params} \
            --aliases {input.aliases} \
            --weights {input.weights} \
            --clades {params.clade} \
            --output {output.tree} \
            {params.gff} \
            {params.clade_key} \
            \
            --plots {params.plots} \
            --cutoff {params.cutoff} \
            --div_add {params.div_add} \
            --div_scale {params.div_scale} \
            --min_size {params.min_size} \
            --bush_scale {params.bush_scale} \
            --bls_range {params.bls_range}
        
        """ #

        # visualize: http://localhost:4000/suggested/D68?branchLabel=new-clade&c=new-clade&d=tree&p=full&tl=__strain__


rule visualization:
    input:
        clades = "{virus}/new-clades_all.tsv",
        tree = rules.suggest_new_clades.output.tree,
        d = "{virus}/plots/",
        params = "{virus}/resources/suggestion_params.json",
    params:
        color = "True",
        grid = "{virus}/plots/clade_table_names.tiff"

    output:
        tiff = "{virus}/plots/violin_plots.tiff",
    shell:
        """
        Rscript scripts/visualize_clades.R --clades {input.clades} --tree {input.tree} \
        --dir-out {input.d} --params {input.params} --out {output.tiff} --grid {params.grid} --color {params.color}
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

rule clean:
    shell:
        """
        rm -r .*/.auto-generated \
        .*/clades/ \
        .*/subclades.tex \
        auspice/*
        """