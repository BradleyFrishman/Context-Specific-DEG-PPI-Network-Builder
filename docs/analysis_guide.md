# Analysis guide

Use the analysis names below to find results even if table numbers change during manuscript formatting. `principal` means the union of G1, G2, and G4; `support` means G3. `limma` means the subset matched to the initial Principal Response set size, not a G3-sized limma control.

## Main manuscript

| Reported analysis | Generated file | Fields |
|---|---|---|
| Connectivity retention / paired ΔLCC | `results/network_summary.csv` | Filter `metric=lcc_size`; `observed_loss`, `null_mean_loss`, `p_smaller_loss` |
| Principal Response LCC size | Same file | `analysis=principal_32` or `principal_41`; `observed_restricted`, `restricted_null_mean`, `p_restricted` |
| limma LCC size | Same file | `analysis=limma_32` or `limma_41`; same fields |
| Support Group LCC size and edges | Same file | `analysis=support_32` or `support_41`; distinguish `metric=lcc_size` and `lcc_edges` |
| Unrestricted versus restricted results | Same file | `observed_unrestricted`, `unrestricted_null_mean`, `p_unrestricted`, and corresponding restricted fields |
| Housekeeping allocation and enrichment | `results/housekeeping_fisher.csv` | Explicit HK/non-HK contingency counts, percentages, fold enrichment, Fisher p-values |
| Observed graph membership and edges | `results/observed_network_nodes.csv`, `observed_network_edges.csv` | Retained selected nodes and first-/second-order edge indicators |

## Supplementary information

| Reported analysis | Generated file |
|---|---|
| G3 median-degree tests, all and non-housekeeping genes | `results/g3_degree_summary.csv` |
| Fisher contingency counts and HK percentages | `results/housekeeping_fisher.csv` |
| Condition-specific mean, median, sample SD, and sample size | `results/degree_descriptive_statistics.csv` |
| Degree-distribution KDE and histogram panels | `results/figures/degree_distributions.png` and `.pdf` |
| Iteration sensitivity | `results/iteration_sensitivity.csv` |
| Original degree-shift figure | `archive/degree_shift/32C.png` and `41C.png`; see its README |

`--plots` also redraws the degree-shift bar labels transcribed in `archive/degree_shift/reported_means.csv`. This is a redraw of reported summary values, not regeneration of their unavailable simulation draws.

The full degree summaries use the supplied complete gene lists. The network node lists record the nodes actually entering each graph. Summary tables and network inputs are kept as separate outputs so they can be traced directly.

## Statistical definitions

- Network enrichment: upper-tail empirical p = `(1 + count(null >= observed)) / (N + 1)`.
- Lower-tail restricted comparisons: `(1 + count(null <= observed)) / (N + 1)`.
- Connectivity retention: ΔLCC = unrestricted LCC size − restricted LCC size; p = `(1 + count(random ΔLCC <= observed ΔLCC)) / 1001`.
- G3 degree tests: compare observed median degree against 1,000 medians from uniformly sampled sets of the same size, without replacement. No degree matching is used for these tests.
- Housekeeping tests: one-sided Fisher exact tests (`alternative='greater'`). The background comparison removes G3 from the comparison group. The within-model comparison uses G1+G2+G4 as the other group.
- All empirical comparisons include ties. All p-values are unadjusted. Sample SD uses `ddof=1`.

## Deterministic randomization

The network seed is `20261011`. Each draw uses `random.Random(f'{seed}:{case}:{zero_based_iteration}')`. Bin order is ascending; background pool order follows the frozen cache. Neither may be reordered without changing the saved sequence.

Degree tests use NumPy `default_rng(20261021 + temperature)` for full G3 and `default_rng(20261021 + temperature + 1000)` for non-housekeeping G3. Background degrees are ordered by sorted gene ID before sampling. The 1,000 draws match the manuscript's selected checkpoint.

Iteration sensitivity uses the first 100, 500, 1,000, and 1,500 restricted draws from a single sequence per gene set. These are not four independent runs. All 24 upper-tail size/edge comparisons have the same significance classifications at 500, 1,000, and 1,500 draws. The combined 32°C Principal Response result remains near the threshold (p = 0.0493 at 1,500). G3 32°C edge counts have no exceedances at any checkpoint; the smaller p-values at larger N reflect finer resolution.

## Reference verification

Reference network summaries and paired draws were copied from the finalized saved analysis. The sensitivity table was copied from its saved 1,500-draw analysis. The first 1,000 degree draws were retained from the saved G3 degree analysis. The housekeeping references were checked against the reported contingency counts and Fisher results; descriptive references were checked against the reported condition-specific statistics.

The input SHA-256 manifest is checked before every run. `--verify` compares row counts and every reference field, using absolute/relative tolerance `1e-12` for numerical fields. Output includes software versions and parameters in `run_metadata.json`.
