# Analysis guide

Use the analysis names below to find results even if table numbers change during manuscript formatting. `principal` means the union of G1, G2, and G4; `support` means G3. `limma` means the subset matched to the initial Principal Response set size, while `limma_g3` means the subset matched to the original G3 set size.

## Main manuscript

| Reported analysis | Generated file | Fields |
|---|---|---|
| Connectivity retention / paired ΔLCC | `results/network_summary.csv` | Filter `metric=lcc_size`; `observed_loss`, `null_mean_loss`, `p_smaller_loss` |
| Principal Response LCC size | Same file | `analysis=principal_32` or `principal_41`; `observed_restricted`, `restricted_null_mean`, `p_restricted` |
| Principal-Response-sized limma LCC size | Same file | `analysis=limma_32` or `limma_41`; same fields |
| G3-sized limma LCC size, edges, and retention | Same file | `analysis=limma_g3_32` or `limma_g3_41`; both metrics and all restricted/unrestricted/loss fields |
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
| Original degree-shift scenarios and figure | `archive/degree_shift/degree_shift_original.py`, `32C.png`, and `41C.png`; see its README |

Figure S1's four degree-shifted reference scenarios are implemented in `archive/degree_shift/degree_shift_original.py`. Run both temperatures with `python archive/degree_shift/run_degree_shift.py`; this uses the bundled inputs and 20 draws per scenario, writing fresh estimates and comparisons with S1 to `results/degree_shift/comparison.csv`. The observed LCC sizes reproduce exactly; randomized means vary because the original script does not set a seed. `run.py --plots` renders the Figure S1 values from `archive/degree_shift/reported_means.csv`. See the archive README for scenario definitions and execution details.

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

Iteration sensitivity uses the first 100, 500, 1,000, and 1,500 restricted draws from a single sequence per gene set. These are not four independent runs. All 24 original upper-tail size/edge comparisons have the same significance classifications at 500, 1,000, and 1,500 draws. The combined 32°C Principal Response result remains near the threshold (p = 0.0493 at 1,500). G3 32°C edge counts have no exceedances at any checkpoint; the smaller p-values at larger N reflect finer resolution.

## Reference verification

Reference network summaries and paired draws were copied from the finalized saved analysis. The sensitivity table was copied from its saved 1,500-draw analysis. The first 1,000 degree draws were retained from the saved G3 degree analysis. The housekeeping references were checked against the reported contingency counts and Fisher results; descriptive references were checked against the reported condition-specific statistics.

G3-sized limma summary rows and 8,000 paired size/edge draw rows were appended from the saved `Limma_G3_results_and_code.zip` analysis and independently reproduced with the existing pipeline. Its additional larger-loss p-values are calculated from those same saved draws. The G3-sized iteration-sensitivity rows extend the same deterministic sequence to 1,500 draws.

The input SHA-256 manifest is checked before every run. `--verify` compares row counts and every reference field, using absolute/relative tolerance `1e-12` for numerical fields. Output includes software versions and parameters in `run_metadata.json`.

## G3-sized limma comparisons

The first 67 and 53 IDs in the supplied limma order match the original G3 counts before filtering; 63 and 51 are eligible for networks. No replacement genes are added. These subsets use the same 20-bin matching, full Gonzalez DEG bridge lists, 1,000 paired draws, and plus-one empirical tests as the other analyses. They are comparisons against each subset's own null distribution, not direct tests of limma versus G3.

### Restricted networks

| Condition | Network | LCC size | Null mean size | Size p | LCC edges | Null mean edges | Edge p |
|---|---|---:|---:|---:|---:|---:|---:|
| 32°C | Second-order only | 42 | 43.0 | 0.6723 | 152 | 133.8 | 0.2607 |
| 32°C | First + second order | 43 | 43.6 | 0.6374 | 154 | 135.7 | 0.2597 |
| 41°C | Second-order only | 25 | 28.0 | 0.8462 | 55 | 64.3 | 0.6863 |
| 41°C | First + second order | 25 | 28.5 | 0.8731 | 60 | 65.3 | 0.5914 |

### Unrestricted networks

| Condition | Network | LCC size | Null mean size | Size p | LCC edges | Null mean edges | Edge p |
|---|---|---:|---:|---:|---:|---:|---:|
| 32°C | Second-order only | 58 | 56.0 | 0.2048 | 384 | 298.6 | **0.0180** |
| 32°C | First + second order | 58 | 56.1 | 0.2148 | 384 | 299.5 | **0.0190** |
| 41°C | Second-order only | 45 | 41.6 | 0.0809 | 138 | 154.0 | 0.7303 |
| 41°C | First + second order | 45 | 41.7 | 0.0839 | 140 | 154.5 | 0.7103 |

### Connectivity retention

| Condition | Network | Unrestricted LCC | Restricted LCC | Observed loss | Mean random loss | Smaller-loss p |
|---|---|---:|---:|---:|---:|---:|
| 32°C | Second-order only | 58 | 42 | 16 | 13.0 | 0.8322 |
| 32°C | First + second order | 58 | 43 | 15 | 12.5 | 0.7952 |
| 41°C | Second-order only | 45 | 25 | 20 | 13.6 | 0.9500 |
| 41°C | First + second order | 45 | 25 | 20 | 13.2 | 0.9610 |

Enrichment p-values are upper-tail empirical comparisons. Retention p-values test smaller-than-expected loss. All p-values include ties and the plus-one correction; bold indicates p < 0.05.


No restricted-network or smaller-loss test is significant for these subsets. Unrestricted 32°C LCC edge enrichment is significant for both network definitions. Full-precision values, including lower-tail and larger-loss tests, are in `reference/network_summary.csv`.
