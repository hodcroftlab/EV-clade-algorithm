# Algorithm notes: clade suggestion scoring (enteroviruses)

This document describes the scoring model used in this repository to **suggest candidate clade breakpoints** on a phylogenetic tree, and how those suggestions translate into `new-clade` labels in Auspice/Nextstrain JSON output.

If you only want to run the workflow, start with the main [`README.md`](README.md).
If you are tuning parameters or interpreting score outputs, read this file.

---

## Overview

The algorithm is designed to be run repeatedly as datasets grow. It is **not** meant to find a globally optimal partition of an unlabeled tree, but to identify *new* groups that have become large, diverged, and/or mutation-defined enough to justify a clade designation.

At a high level, for every node/branch we compute three components:

1. **Bushiness / growth signal** (`bushiness_raw`, normalized `bushiness`)
2. **Branch mutation score** (`branch_score`) based on **weighted AA substitutions**
3. **Divergence contribution** (`div_score`) based on divergence since the last clade breakpoint

These are combined into a **total score** (`score`) and compared to a cutoff. If the cutoff is exceeded (and a minimum size is met), the branch is marked as a new clade breakpoint and the subtree is assigned a new hierarchical clade ID which is then converted to a human-friendly name.

---

## Input assumptions

The algorithm operates on a Nextstrain **Auspice JSON** export (typically via `augur export`) and expects:

- `tree` with `children`, `node_attrs`, `branch_attrs`
- `branch_attrs.mutations`:
  - must include `nuc` for branch-length proxy used in bushiness
  - should include one or more translated CDS keys (e.g. `VP1`) to score amino-acid mutations
- `node_attrs.num_date.value` for time-filtering ("alive" tips) if available
- an existing clade label present in `node_attrs[clade_key].value` (or the workflow can insert a dummy clade)

---

## Definitions

### Terminal vs non-terminal nodes

- **Terminal** node: no `children`
- **Non-terminal** node: has one or more `children`

The implementation uses these distinctions for recursive traversal and for defining bushiness base cases.

---

## Step 1 — Alive tips and tree "cleanup"

### Alive tips

A node is considered "alive" if it falls within a date window:

- `min_date < num_date < max_date`

This is used to focus the growth/bushiness score on **recent** lineages (configurable), so very old clusters do not dominate clade designation just because they have many historical tips.

**If `num_date` is missing**, the code currently treats the node as alive.

### Collapsing mutation-less internal nodes (tree smoothing)

During preparation, mutation-less internal nodes may be collapsed upward (implementation detail) to reduce chains of internal nodes with no mutations and simplify the tree for scoring.

This does **not** change topology in a way that affects clade membership, but can change where exactly breakpoints land if a long chain of internal nodes has no mutations.

---

## Step 2 — Backbone labeling (protect existing clades)

Before suggesting new clades, the tree is annotated with a boolean attribute:

- `backbone = True` for nodes on paths from root to existing clade breakpoints
- `backbone = False` otherwise

This allows the algorithm to avoid creating new clades along the "spine" of already-defined clades when `ignore_backbone=True` is used in bushiness/scoring.

**What counts as an "existing breakpoint"?**

- A branch with `branch_attrs.labels[clade_key]` (dataset-dependent)

---

## Step 3 — Divergence since last breakpoint

Each node gets a divergence counter:

- `div`: cumulative count of mutations along the path since the root **(then effectively "reset" at breakpoints during traversal)**

In this repo’s implementation the divergence is computed as **counts of mutations in a configurable set of genes/proteins** (typically `VP1` for EV clades).

### What counts as a divergence event?

For amino acids, a mutation string looks like `A123T` meaning A→T at position 123.

The algorithm ignores:

- deletions (`-`)
- unknown amino acids (`X`)
- for nuc, ambiguous states like `N` and gaps `-` where relevant

---

## Step 4 — Bushiness / phylogenetic growth score

### Motivation

We want to prioritize clades that are:

- not just deep splits, but show **recent downstream growth**
- robust to uneven branch lengths
- comparable across different parts of the tree

We compute a local branching index (LBI)-like measure on the tree.

### Recurrence

For a node \( n \), bushiness is:

```math
\phi_n = \sum_{c \in children(n)} \left[(1 - e^{-l_c/d}) \cdot I(c\ \text{alive}) \;+\; \phi_c \, e^{-l_c/d}\right]
```

Where:

- \( c \) are child nodes of \( n \)
- \( l_c \) is the branch length leading to \( c \)
- \( d \) is a distance scale (in code: you pass a `distance(c)` function; the scale is controlled by `bushiness_branch_scale`)
- \( I(c\ \text{alive}) \) is an indicator (1 if alive, else 0)

**Terminal nodes** have:

- \( \phi = 1 \) if alive else 0

### Branch length proxy (implementation)

Enterovirus trees often don’t have consistent continuous branch lengths in JSON.
Instead, this repo uses a *mutation-count proxy*:

- branch length ≈ number of nucleotide mutations on the branch excluding `N` and `-`
- then divided by `bushiness_branch_scale`

Conceptually:

- smaller `bushiness_branch_scale` → longer effective branch lengths → bushiness decays faster
- larger `bushiness_branch_scale` → shorter effective branch lengths → bushiness persists further

### Normalization

`bushiness_raw` can vary widely with sampling density and tree size. We normalize it using saturation:

```math
B = \frac{\phi}{\phi + s_{\phi}}
```

Where \( s_{\phi} \) is a scale derived from the tree, often an upper percentile (implementation uses the 80th percentile of internal-node bushiness among "alive" nodes).

This yields `bushiness` in `[0, 1)`.

---

## Step 5 — Branch mutation score (AA weights)

### Motivation

We want breakpoints to land on branches with **informative AA changes**, e.g. known antigenic loops in EV-D68 VP1.

### Weighting model

For each non-nucleotide CDS key in `branch_attrs.mutations`:

- look up per-position weights in `weights.json`
- sum the weights across AA substitutions on that branch

Let raw AA weight on the branch be \( w \). We normalize via saturation:

```math
S = \frac{w}{w + s_w}
```

Where:

- \( s_w \) corresponds to `branch_length_scale` in config/code
- interpretation: how many "weight units" should constitute a "strong" branch

This yields `branch_score` in `[0, 1)`.

### Practical advice for `weights.json`

Keep it simple:

- encode a small set of known important positions with higher weights (2–3)
- use a conservative `"default"` (often 0 or 1)
- start with VP1 only unless you have a strong reason to include other genes given recombination

---

## Step 6 — Divergence contribution (saturating)

Divergence since the last breakpoint is used to make new clades easier to trigger as a lineage accumulates changes.

Let:

- \( \Delta = div - div_{breakpoint} \) (divergence since last breakpoint)
- \( s_d \) = `divergence_scale`

The divergence contribution is:

```math
D = a_d \cdot \frac{\Delta}{\Delta + s_d}
```

Where:

- \( a_d \) = `divergence_addition` (max contribution to total score)

So `div_score` ranges from `0` to `divergence_addition`.

### Important implementation note (this repo)

In `tree_clade_assignment.py` the divergence contribution is **added into the stored `score`**:

- `node_attrs["score"]["value"] += div_score`

So the exported `score` is a **total score** already including divergence.

---

## Step 7 — Total score and clade triggering

### Total score

Before divergence is added, the combined "base score" is effectively:

```math
score_{base} = bushiness + branch\_score
```

After adding divergence, the total is:

```math
score_{total} = bushiness + branch\_score + div\_score
```

### Trigger condition

A node triggers a new clade breakpoint if:

- `score_total > cutoff`
- and `ntips > min_size`

There is an additional "one-daughter" safeguard:

- if exactly one child would also trigger under the same thresholding logic, the parent node is not used as a breakpoint
- this reduces redundant breakpoints along a single path

### Divergence reset at breakpoints

When a breakpoint is designated, divergence is conceptually "reset" for downstream evaluation by updating the base divergence reference used for \(\Delta\).

---

## Step 8 — Naming and aliasing

Internally, clades are represented hierarchically as tuples, e.g.:

- `('A', 1, 2)` or similar (exact structure depends on your alias map)

An `aliases.json` file provides mapping between:

- a short top-level label (e.g. `"A"`)
- and the corresponding full hierarchical tuple prefix

When a new breakpoint is triggered:

- the algorithm assigns the next available child index under the parent clade in the hierarchy
- then converts the full hierarchical representation into a short name via `full_clade_to_short_name(...)`

The result is written to:

- `branch_attrs.labels["new-clade"]` at the breakpoint branch
- `node_attrs["new-clade"].value` for all downstream nodes

---

## Outputs and how to validate them in Auspice

The suggested tree JSON includes several helpful continuous attributes:

- `score` (total score used for triggering; includes divergence)
- `div_score` (divergence contribution only)
- `div_impact` (implementation-defined helper; roughly "how much divergence is needed to cross cutoff")
- `bushiness_raw` and `bushiness`
- `branch_score`
- `new-clade` (assigned clade label)

### Practical validation checklist

When inspecting in Auspice:

1. Color by **`new-clade`**
   - do proposed clades correspond to visually coherent clusters?
   - are they too nested / too granular?

2. Color by **`score`**
   - do breakpoint branches sit near the upper end of score distribution?

3. Color by **`branch_score`**
   - do breakpoints align with branches containing expected informative VP1 changes?

4. Color by **`div_score`**
   - is divergence driving clade calls everywhere (too high `divergence_addition` or too low `divergence_scale`)?
   - or is it never contributing (too low `divergence_addition` or too high `divergence_scale`)?

5. Check sizes
   - if many tiny clades appear: increase `min_size` and/or `cutoff`

---

## Parameter interpretation and tuning heuristics

### `cutoff`

- Higher cutoff → fewer clades
- Lower cutoff → more clades
- If you want "only when at least two signals agree", target cutoff around ~1.0 (because components are saturating in ~[0,1] and divergence adds up to `divergence_addition`)

### `min_size`

- Most effective lever to prevent "tip noise"
- If you get many small "microclades", raise `min_size`

### `divergence_addition` and `divergence_scale`

- `divergence_addition` sets the **maximum** divergence contribution (0 disables divergence)
- `divergence_scale` sets how quickly divergence saturates:
  - smaller → divergence rises quickly with few changes
  - larger → divergence needs many changes to matter

### `branch_length_scale`

- Larger → branch_score saturates more slowly (less sensitive to a few high-weight mutations)
- Smaller → a few weighted AA changes push branch_score close to 1 quickly

### `bushiness_branch_scale`

Controls the effective decay length used in bushiness:

- smaller → bushiness is "short-memory" (focuses strongly on very local tip growth)
- larger → bushiness is "longer-memory" (signals propagate further up)

---

## Enterovirus-specific caveats

- **Recombination:** genome-wide phylogeny and VP1 phylogeny can disagree. For nomenclature, VP1-only is often most interpretable.
- **Sampling changes:** bushiness depends on sampling density; expect to re-tune when surveillance intensity changes.
- **Time filtering:** `min_date` and `max_date` strongly affect which parts of the tree contribute to bushiness. Use them deliberately (e.g. focus on recent years when defining new names).

---

## References / provenance

This approach is adapted from the influenza clade suggestion algorithm:

- <https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm>

Related tools:

- Nextstrain: <https://nextstrain.org>
- Nextclade: <https://clades.nextstrain.org>
