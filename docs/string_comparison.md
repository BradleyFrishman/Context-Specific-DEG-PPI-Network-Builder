# STRING comparison

STRING was evaluated separately through its **version 12.0 website**, not through the offline PPI algorithm. The reported settings were:

| Setting | Value |
|---|---|
| Organism | Homo sapiens |
| Network | Full STRING network: functional and physical associations |
| Minimum required interaction score | Medium confidence, 0.400 |
| Evidence channels | Text mining, experiments, databases, co-expression, neighborhood, gene fusion, co-occurrence |
| Statistical background | Default whole-genome background |
| Reported outcome | STRING PPI enrichment p-value |

For the observed comparisons, upload the complete IDs in `data/gene_sets/principal_32.csv`, `principal_41.csv`, `limma_32.csv`, or `limma_41.csv`, selecting the `ENSEMBL_ID` column. Check the identifier mappings returned by STRING before inspecting enrichment. The 9,825-gene offline background is not used for this website comparison.

The manuscript's reported observed values are preserved in `reference/string_reported.csv`. They were transcribed from the earlier full manuscript's STRING results table; the latest Methods attachment specifies the same version and settings. The original STRING session exports, mapping reports, and random-control exports were not supplied. Consequently, the offline verification command does not assert exact reproduction of these website outputs or the reported variability among random controls. Repeating the website procedure against a changed database version may produce different values.
