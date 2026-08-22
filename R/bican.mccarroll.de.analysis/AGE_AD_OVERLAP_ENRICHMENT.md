# Directional DE Overlap Enrichment Between Aging and AD

## Objective

Test whether genes that are significantly differentially expressed with age
in BICAN are enriched for genes that are significantly differentially
expressed in the same direction with AD progression in the external Gabitto
et al. 2024 dataset (PMID 39402379).

## Data

| Role | Path (relative to the data root) |
|---|---|
| Age (BICAN) | `differential_expression/results/LEVEL_6_sea_ad_mtg_mmc/sex_age/cell_type/<cell_type>__age_DE_results.txt` |
| AD (Gabitto et al. 2024) | `differential_expression/external_comparison_PMID_39402379/voom-like/<cell_type>__MTG__ad_cps_DE_results.txt` |
| Matched cell types + display order | `differential_expression/metadata/cell_types_for_sea_ad_mtg_mmc_plots.txt` (33 SEA-AD supertypes) |

The BICAN side is region-pooled, and the AD side is the `ad_cps`
(continuous pseudo-progression score) contrast. All 33 BICAN supertypes
match a Gabitto cell type by exact string. Both DE result files have gene
as the row name and columns `logFC AveExpr t P.Value adj.P.Val B z.std`;
only `logFC` and `adj.P.Val` are used.

## Gene universe

For each cell type, the universe is genes with a finite `logFC` and a
finite `adj.P.Val` in **both** studies.

## Significance and direction

Within the shared universe, for each cell type:

```
age_up   = adj_p_age <= fdr_cutoff_age & beta_age > 0
age_down = adj_p_age <= fdr_cutoff_age & beta_age < 0
ad_up    = adj_p_ad  <= fdr_cutoff_ad  & beta_ad  > 0
ad_down  = adj_p_ad  <= fdr_cutoff_ad  & beta_ad  < 0
```

`fdr_cutoff_age` and `fdr_cutoff_ad` both default to `0.05`.

## Inclusion filter

A cell type is included only if both `n_age_de` (age-up + age-down) and
`n_ad_de` (AD-up + AD-down) are at least `min_de_genes` (default `10`).
Excluded cell types are dropped from the Fisher results but still appear
in the counts table with `included = FALSE`.

At the defaults, 29 of 33 cell types are included (58 Fisher tests). The
four excluded are `Lamp5_2`, `Lamp5_Lhx6_1`, `Oligo_1`, and
`OPC_2_1-SEAAD`.

## Fisher tests

For each included cell type, two one-sided (`alternative = "greater"`)
`stats::fisher.test()` calls, one per direction, with no pseudocounts even
when a cell is zero:

```
             AD up      Not AD up
Age up         a            b
Not age up     c            d
```

```
             AD down    Not AD down
Age down       a            b
Not age down   c            d
```

- `a` = overlap (both significant, same direction)
- `b` = age-significant only
- `c` = AD-significant only
- `d` = neither

### Worked example: L23_IT_1

```
universe_size = 15036
n_age_up = 417     n_age_down = 617
n_ad_up  = 763     n_ad_down  = 2813
```

**Up table:**

| | AD up | Not AD up |
|---|---|---|
| **Age up** | 80 | 337 |
| **Not age up** | 683 | 13,936 |

```
odds_ratio = (80 * 13936) / (337 * 683) = 4.843
p_value    = 1.173e-25
p_adj      = 9.715e-25
```

**Down table:**

| | AD down | Not AD down |
|---|---|---|
| **Age down** | 181 | 436 |
| **Not age down** | 2,632 | 11,787 |

```
odds_ratio = (181 * 11787) / (436 * 2632) = 1.859
p_value    = 4.253e-11
p_adj      = 9.135e-11
```

## Multiple-testing correction

Benjamini-Hochberg correction is applied to the Fisher `p_value`s of
whatever tests actually ran (`n = 2 * K`, `K` = number of included cell
types), via `stats::p.adjust()`'s default denominator.

## Correlation

`pearson_r` and `spearman_rho` are computed on the subset of genes
significant in **both** datasets: `(age_up | age_down) & (ad_up | ad_down)`.
`spearman_rho` additionally gets a `stats::cor.test()` P value
(`spearman_p_value`), Benjamini-Hochberg corrected across the cell types
with at least 2 genes in that subset (`spearman_p_adj`) -- its own
correction family, separate from the Fisher correction.

### Worked example, continued: L23_IT_1

```
both_sig genes   = 277
pearson_r        = 0.800
spearman_rho     = 0.839
spearman_p_value = 0    (underflows double precision)
spearman_p_adj   = 0
```

## Discordant genes

Genes significant in both datasets with opposite effect signs are listed
(not tested with a formal Fisher test):

```
age_up_ad_down = age_up   & ad_down
age_down_ad_up = age_down & ad_up
```

For `L23_IT_1`: 12 genes are `age_up_ad_down` (e.g. `ALOXE3`, `ATP6V1C2`,
`CTPS2`, `CYFIP1`, `KDM2B-DT`, ...) and 4 are `age_down_ad_up`. Each row
also carries `n_cell_types_same_pattern`: how many included cell types
show that exact gene with that exact discordant pattern.

## Output tables

`compute_de_overlap_enrichment(gene_dt, fdr_cutoff_age, fdr_cutoff_ad,
min_de_genes)` returns:

- **`fisher`**: one row per included cell type per direction --
  `cell_type, direction, universe_size, n_age_de, n_ad_de, n_overlap,
  odds_ratio, p_value, p_adj, pearson_r, spearman_rho, spearman_p_value,
  spearman_p_adj`.
- **`counts`**: one row per cell type, including excluded ones --
  `cell_type, included, universe_size, n_age_up, n_age_down, n_ad_up,
  n_ad_down, n_up_up, n_down_down, n_age_up_ad_down, n_age_down_ad_up,
  pearson_r, spearman_rho, spearman_p_value, spearman_p_adj`.
- **`discordant`**: `cell_type, gene, beta_age, adj_p_age, beta_ad,
  adj_p_ad, pattern, n_cell_types_same_pattern`.

## Plot

`plot_de_overlap_enrichment_heatmap(fisher_dt, cell_type_order,
p_adj_threshold)` draws three single-column heatmaps side by side -- Up
odds ratio, Down odds ratio, and Spearman rho -- each on its own linear,
sequential color scale (no log transform, no shared scale between panels),
modeled on Frohlich et al. 2024's Fig. 6a. Odds-ratio cells are asterisked
using `p_adj`; the rho cell is asterisked using `spearman_p_adj`. Separate
scales per panel matter because the Fisher test is one-sided
(enrichment-only): a shared diverging scale would visually imply a tested
depletion signal that was never tested.

## Implementation

- `bican.mccarroll.de.analysis::build_de_overlap_gene_table()` -- the join
  step (one row per tested `cell_type`/`gene` pair), cached.
- `bican.mccarroll.de.analysis::compute_de_overlap_enrichment()` -- the
  Fisher tests, correlation, FDR corrections, and discordant-gene table.
- `bican.mccarroll.de.analysis::plot_de_overlap_enrichment_heatmap()` --
  the figure.
- `bican.mccarroll.figures::plot_de_overlap_enrichment_bican_sea_ad_vs_pmid_39402379(min_de_genes, fdr_cutoff_age, fdr_cutoff_ad, force_recompute)`
  -- the exported wrapper.

Source: `R/de_overlap_enrichment.R` in this package;
`R/de_overlap_enrichment_plots.R` in `bican.mccarroll.figures`.
