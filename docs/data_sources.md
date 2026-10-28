# Data sources

| Input | Source and role |
|---|---|
| `data/background_interactions.csv.gz` | Exact supplied Human Protein Atlas interaction cache for 9,825 candidate genes; frozen for reproduction |
| `data/gene_sets/principal_*.csv` | Supplied Gonzalez-model G1+G2+G4 lists, preserving input order |
| `data/gene_sets/support_*.csv` | Supplied Gonzalez-model G3 lists, preserving input order |
| `data/gene_sets/limma_32.csv`, `limma_41.csv` | Supplied top-ranked limma subsets, matched to initial Principal Response counts |
| Derived `limma_g3_32` and `limma_g3_41` | First 67 and 53 IDs from the existing ordered limma lists, matching original G3 counts; 63 and 51 remain after network eligibility filtering, without replacement |
| `data/deg_bridges/32.csv`, `41.csv` | Full condition-specific Gonzalez DEG lists: 2,810 and 2,679 IDs |
| `data/housekeeping_genes.csv` | Unique Ensembl gene IDs from the supplied Housekeeping and Reference Transcript Atlas v1.0 annotation file |
| `data/manifest.json` | Byte counts and SHA-256 digests of all computational inputs |

The Human Protein Atlas Interaction section aggregates IntAct, BioGRID, BioPlex, and OpenCell physical interaction evidence, as described in the manuscript. The frozen cache stores one Ensembl gene ID per row and semicolon-separated recorded partner IDs in `KNOWN_INTERACTIONS`. Compression does not change its decompressed bytes. Self-interactions are removed when calculating degree or graph edges.

Gene-set and bridge files contain only the `ENSEMBL_ID` column required by the analyses. Their row order is preserved from the supplied files. Housekeeping IDs are deduplicated at gene level and version suffixes are stripped as in the original degree analysis. There are 2,139 housekeeping IDs within the 9,825-gene background.

No live database download is required. The original cache's collection date was July 2025, so the repository identifies the snapshot by its checksum. Downloading present-day interaction data would constitute a different input snapshot and can change results.

## External resources

- Human Protein Atlas: https://www.proteinatlas.org/
- Housekeeping and Reference Transcript Atlas: https://housekeeping.unicamp.br/
- STRING: https://string-db.org/
- Gonzalez et al. temperature-response model: https://doi.org/10.1186/s12864-026-12823-7

The code license does not replace database licensing or attribution requirements. Consult each database's terms when reusing its content.

The G3-sized limma selections and 1,000-draw reference results were checked against `Limma_G3_results_and_code.zip` (26 October 2026). Its limma subset order, background interaction cache, and full Gonzalez DEG bridge IDs match the inputs already in this repository. The subsets are derived in `src/common.py`, so no duplicate input files are required. The source provenance identifies `paper_update_.zip/humans_limaDEG_32.csv` and `humans_limaDEG_41.csv` as the original ordered limma lists.

## Original degree-shift source

`archive/degree_shift/degree_shift_original.py` is the original degree-shift analysis script supplied by the author on 28 October 2026 as `SHIFTS on Build 2nd order PA sets from Background 4.py`. Its contents are preserved unchanged. It implements the SI's four conservative/aggressive up/down reference scenarios, with 20 iterations per scenario and DEG-restricted second-order-only LCC evaluation. Its supplied configuration selects the 32°C Principal Response input; the corresponding 41°C analysis uses that condition's Principal Response and full DEG lists. See the archive README for the exact input schema and sampling behavior.
