# EV-clade-algorithm

Suggest candidate clade breakpoints and clade labels for non‑polio enteroviruses (EV‑D68, EV‑A71, CV‑A16, CV‑A10, …) from **Nextstrain Auspice JSON trees**.

This repo adapts the influenza clade-suggestion approach to enteroviruses using a combined score from:

- **tree growth / structure** (“bushiness”, LBI-like),
- **amino-acid change patterns** (weighted substitutions),
- **divergence since last clade breakpoint** (saturating contribution).

> This is a *clade suggestion* tool: it proposes breakpoints and labels for review. It does not claim a globally optimal partition.

---

## Repository layout (current)

The repo is organized with one directory per virus:

```text
EV-clade-algorithm/
├── Snakefile
├── config.yaml
├── auspice/
│   ├── original_{virus}.json
│   └── suggested_{virus}.json
├── scripts/
│   ├── add_new_clades.py        # main algorithm (called by Snakefile)
│   ├── calculate_optimal_scales.py     # scale recommendations
│   ├── extract_aliases_from_tree.py    # optional: derive aliases from an annotated tree
│   ├── extract_yml_from_json.py        # export clade definitions to YAML
│   ├── construct_tsv.py                # build TSVs from exported YAMLs
│   └── visualize_clades.R              # optional plots
├── ev-d68/
│   ├── resources/
│   │   ├── weights.json
│   │   ├── aliases_*.json
│   │   └── genome_annotation.gff3      # optional (if not fetched/generated)
│   └── results/
│       ├── original_ev-d68.clades.json
│       ├── optimal_scales.json
│       ├── new-clades.tsv
│       └── (plots/)
├── ev-a71/
├── cva16/
└── cva10/
```

(Directory names must match the keys under `viruses:` in `config.yaml`.)

---

## Quick start

### 0) Install dependencies

This repo uses a `renv.yaml` conda/ micromamba environment with Python + R dependencies and required CLI tools.

```bash
micromamba env create -f renv.yaml -n clade
micromamba activate clade
```

### 1) Configure viruses + defaults (`config.yaml`)

At minimum:

- list viruses under `viruses:`
- provide `tree_url` (if not using local trees)
- set defaults (`defaults:`) for date range, proteins, and thresholds

### 2) Run the workflow for everything

```bash
snakemake -c 1 all
```

Outputs:

- `auspice/suggested_{virus}.json` (viewable in Auspice)
- `{virus}/results/optimal_scales.json`
- `{virus}/results/new-clades.tsv` (+ optional sweep output)

### 3) View results in Auspice locally

```bash
auspice view --datasetDir auspice
```

Then open something like:

```bash
http://localhost:4000/suggested/ev-d68?branchLabel=new-clade&c=new-clade&d=tree&p=full&showBranchLabels=all
```

---

## What the workflow does (Snakefile)

For each `{virus}`:

1. **Fetch input tree + annotation** (unless local tree provided)
   - Tries **Nextclade dataset**: `enpen/enterovirus/{virus}`
   - Falls back to `curl {tree_url}`
   - If only a tree is fetched, a **GFF3** is generated from `meta.genome_annotations` using `jq`

2. **Ensure existing clade annotations exist**
   - If the input tree has no `node_attrs[clade_key]`, assigns a dummy clade to all nodes.
   - This keeps the algorithm happy because it expects an “old clade key” for hierarchy/backbone logic.

3. **Create default inputs**
   - `resources/weights.json` if missing (simple defaults)
   - optional `resources/aliases_default.json` extracted from the tree (if you want a starting point)

4. **Compute recommended scales** (`scripts/calculate_optimal_scales.py`)
   - writes `{virus}/results/optimal_scales.json`
   - used as defaults for sweep ranges / per-virus settings

5. **Suggest new clades** (`scripts/add_new_clades.py`)
   - writes `auspice/suggested_{virus}.json`
   - writes `{virus}/results/new-clades.tsv`
   - optionally runs a **parameter sweep** (`--plots True`) and writes a TSV of sweep results

---

## Core algorithm (summary)

Implemented in `scripts/add_new_clades.py`.

### Per-node quantities

- **alive**: tip in `[min_date, max_date]` based on `node_attrs.num_date.value`
- **ntips**: number of descendant tips
- **bushiness_raw**: LBI-like score using exponential decay along branches  
  (branch length is computed from nucleotide mutations, scaled by `bushiness_branch_scale`)
- **bushiness**: normalized bushiness `raw/(raw + bushiness_scale)`
- **branch_score**: normalized weighted AA score `aa_weight/(branch_length_scale + aa_weight)`
- **div_score**: divergence contribution since last breakpoint:
  `divergence_addition * delta_div/(delta_div + divergence_scale)`

### Trigger rule (current implementation)

A branch triggers a new clade if:

- `score_total > cutoff` and `ntips > min_size`
- plus a “one-daughter” check to avoid creating a clade when exactly one child would also trigger

**Important implementation detail:** the script *adds* `div_score` directly into `node_attrs["score"]["value"]` before testing the cutoff.

---

## Key configuration concepts

### `clade_key` (“old clade key”)

This is the node attribute key containing **existing clade membership**, used to:

- build the existing hierarchy (`full_{clade_key}`)
- label the backbone (avoid placing new clades on branches leading to existing clade roots)

For EV datasets this is usually something like:

- `clade_membership`
- or dataset-specific

The workflow’s `ensure_clades` rule can assign a dummy clade if missing.

### Proteins / genes

The algorithm uses:

- nucleotide mutations (`nuc`) for bushiness branch length
- amino-acid mutations for branch weighting (`weights.json`)
- divergence computed over `defaults.proteins` (if empty, it falls back to `nuc`)

If a `genome_annotation.gff3` is present, the script can infer CDS names (useful if you don’t want to hardcode proteins).

---

## Parameter sweep (recommended for tuning)

The Snakemake rule runs `add_new_clades.py` with sweep ranges from `config.yaml`:

- `--cutoff-sweep`
- `--div_add-sweep`
- `--div_scale-sweep`
- `--min_size-sweep`
- `--bush_scale-sweep`
- `--bls_range-sweep`

Output: a TSV (same path as `--clades`) with one row per parameter set including:

- counts of old/new clades
- mean node stats (`score`, `bushiness`, `div`, `branch_score`)

**Tip:** start with a small grid (e.g. 5×5 over `cutoff` and `divergence_addition`) before adding more axes.

---

## Outputs (what to look at)

In `auspice/suggested_{virus}.json`:

- `new-clade`: suggested clade name
- continuous colorings:
  - `score` (already includes `div_score` in current implementation)
  - `div_score`
  - `div_impact`
  - `bushiness_raw`, `bushiness`
  - `branch_score`

In `{virus}/results/new-clades.tsv`:

- a simple table listing old clades vs newly suggested clades (mainly for quick counts)

---

## Caveats (enterovirus-specific)

- **Recombination** can break the “single tree = single history” assumption. Consider restricting mutation scoring/divergence to VP1 unless you have a clear reason to include additional proteins.
- Sampling density strongly affects bushiness. Re-tune thresholds when sampling changes substantially.

---

## Useful commands

Clean derived outputs:

```bash
snakemake clean -c 1
```

Export clade definitions to YAML (for downstream Nextclade / documentation):

```bash
snakemake extract_clade_ymls -c 1
```

---

## References

- Upstream inspiration: <https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm>
- Nextstrain: <https://nextstrain.org>
- Nextclade: <https://clades.nextstrain.org>
