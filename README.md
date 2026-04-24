# EV-Clade-Algorithm

Automated clade assignment for enteroviruses (EV-D68, EV-A71, CVA16, etc.) based on phylogenetic structure, mutation patterns, and sequence divergence.

Adapts the [influenza clade-suggestion-algorithm](https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm) for non-polio enteroviruses.

## What it does

- Analyzes Nextstrain phylogenetic trees (JSON format)
- Scores branches combining:
  - **Phylogenetic signal** (bushiness: number of downstream tips with exponential decay)
  - **Amino acid mutations** (weighted by position: epitope sites > other sites)
  - **Sequence divergence** (cumulative AA changes since last clade breakpoint)
- Suggests new clades when combined score exceeds threshold
- Outputs tree with clade assignments ready for Nextstrain visualization

## Structure

```bash
# please make sure the directory structure is as follows:
enterovirus-clade-nomenclature/
├── {virus}/
│   ├── config/
│   │   ├── suggestion_params.json
│   │   ├── weights.json
│   │   ├── aliases.json
│   │   ├── genome_annotation.gff3
│   ├── CHANGELOG.md
│   ├── README.md
├── clade-suggestion-algorithm/
│   ├── auspice/
│   ├── scripts/
│   │   ├── add_new_clades.py
│   │   ├── calculate_optimal_scales.py
│   │   ├── construct_tsv.py
│   │   ├── extract_yml_from_json.py
│   │   ├── generate_markdown_summary.py
│   │   ├── visualize_clades.R
│   ├── Snakefile
│   ├── README.md
```

## Quick start

### Single run

```bash
snakemake clade-suggestion-algorithm/Snakefile all
```

## Key configuration parameters

| Parameter | Meaning | Typical range | Notes |
|-----------|---------|--------------|-------|
| `cutoff` | Score threshold for new clade | 0.5–2.0 | Main tuning lever |
| `min_size` | Min. tips per clade | 3–15 | Prevents spurious clades |
| `divergence_addition` | Weight on divergence score | 0.0–1.0 | 0 = ignore divergence |
| `divergence_scale` | Divergence saturation | 2–10 | Larger = more resistant |
| `bushiness_branch_scale` | Phylo decay rate | 2–10 | Larger = shorter memory |
| `branch_length_scale` | Mutation weight saturation | 2–10 | Larger = less sensitive |
| `proteins` | Protein targets | `["VP1"]` or `["VP1", "2C"]` | Virus-dependent |

## Config files

- **`config/suggestion_params.json`** — Hyperparameters (cutoff, divergence, size thresholds)
- **`config/weights.json`** — Per-position mutation weights (epitope sites typically 2–3, others 1, default 0)
- **`config/aliases.json`** — Clade naming scheme (e.g., `(2020, 1)` → `"A"`)
- **`Snakefile`** — Workflow (fetches tree from Nextstrain, runs algorithm, outputs JSON)

## Adapting to other EV species

Minimal changes needed per virus:

```json
{
  "max_date": 2026.0,
  "min_date": 2020.0,
  "proteins": ["VP1"],           // EV-D68, EV-A71: VP1 only
                                  // CVA16, CVA6: ["VP1", "3D"] (recomb. hotspot)
  "cutoff": 1.0,                 // Tune per virus (0.8–1.4)
  "bushiness_branch_scale": 4,
  "branch_length_scale": 4,
  "divergence_addition": 0.5,    // Adjust for EV evolution rate
  "divergence_scale": 4,
  "min_size": 5
}
```

**Steps:**

1. Copy `config/` → `config_EVA71/` (or similar)
2. Adjust `proteins` list and evolution-rate parameters
3. Update `suggestion_params.json` with virus-specific cutoffs
4. Retune `weights.json` if epitope sites differ from EV-D68
5. Run parameter sweep on your tree to find optimal `cutoff` and `divergence_addition`

## Interpreting results

- **`suggested_tree.json`** — Tree with new clades + scoring attributes:
  - `new-clade`: suggested clade name
  - `score`: combined phylo + mutation score
  - `bushiness`: phylogenetic signal (0–1)
  - `branch_score`: mutation weight score (0–1)
  - `div_score`: divergence contribution (0–1)

- **`*_all.tsv`** (parameter sweep output):
  - `cutoff`, `divergence_addition`, `min_size`, etc. — parameter values
  - `new_clade` — count of suggested clades
  - `old_clade` — count of existing clades (baseline)
  - `mean_score`, `max_score`, `score_std` — node scoring statistics

## Workflow (Snakefile)

1. **Fetch** tree from Nextstrain URL (config)
2. **Suggest** new clades using `add_new_clades.py`
3. **Output** `auspice/suggested_{virus}.json` ready for visualization

Edit `Snakefile` to set `tree_url` and adjust virus names for your viruses.

## Next steps

- [ ] Run parameter sweep on your EV-D68 tree
- [ ] Plot parameter sensitivity (heatmap of cutoff vs. divergence_addition)
- [ ] Validate against manually-curated clades
- [ ] Adapt config for EV-A71, CVA16, etc.
- [ ] Integrate output column into Nextclade pipeline

## References

- Original algorithm: [influenza-clade-nomenclature/clade-suggestion-algorithm](https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm)
- Nextstrain: <https://nextstrain.org>
- Nextclade: <https://clades.nextstrain.org>
