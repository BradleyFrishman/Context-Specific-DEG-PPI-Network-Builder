# STRING comparison

STRING was evaluated separately through its **version 12.0 website**. The settings used were:

| Setting | Value |
|---|---|
| Organism | Homo sapiens |
| Network | Full STRING network: functional and physical associations |
| Minimum required interaction score | Medium confidence, 0.400 |
| Evidence channels | Text mining, experiments, databases, co-expression, neighborhood, gene fusion, co-occurrence |
| Statistical background | Default whole-genome background |
| Reported outcome | STRING PPI enrichment p-value |

For the observed comparisons, upload the complete IDs in `data/gene_sets/principal_32.csv`, `principal_41.csv`, `limma_32.csv`, or `limma_41.csv`, selecting the `ENSEMBL_ID` column. Check the identifier mappings returned by STRING before inspecting enrichment. The 9,825-gene offline background is not used for this website comparison.

Repeating the website procedure against a changed database version may produce different values.

## Random gene-set comparisons

The following values were recorded in the original STRING comparison for Principal Response (G1+G2+G4) genes and 20 random gene sets per temperature. The manuscript describes the random sets as approximately degree-matched using Human Protein Atlas interaction degree. The values below are transcribed from the author's comparison tables; file numbers retain their original labels. The settings above apply to this comparison.

| Gene set | 32°C PPI enrichment p-value | 41°C PPI enrichment p-value |
|---|---|---|
| Principal Response (G1+G2+G4) | 5.86e-05 | 1.3e-08 |
| Random file 1 | 0.757 | 0.0164 |
| Random file 2 | 0.0274 | 0.012 |
| Random file 3 | 0.0771 | 0.339 |
| Random file 4 | 0.0358 | 0.175 |
| Random file 5 | 0.0133 | 0.0537 |
| Random file 6 | 0.205 | 0.0325 |
| Random file 7 | 0.00848 | 0.842 |
| Random file 8 | 0.0342 | 0.0776 |
| Random file 9 | 0.283 | 0.22 |
| Random file 10 | 0.0381 | 0.134 |
| Random file 11 | 0.409 | 0.459 |
| Random file 12 | 0.135 | 0.324 |
| Random file 13 | 0.644 | 0.187 |
| Random file 14 | 0.000852 | 0.0011 |
| Random file 15 | 0.0184 | 0.00456 |
| Random file 16 | 0.552 | 0.043 |
| Random file 17 | 1.14e-05 | 0.726 |
| Random file 18 | 0.2 | 0.219 |
| Random file 19 | 0.467 | 0.0117 |
| Random file 20 | 0.204 | 0.0747 |

Across the 20 random sets, STRING PPI enrichment p-values ranged from 1.14e-05 to 0.757 at 32°C and from 0.0011 to 0.842 at 41°C. Nine of 20 random sets at 32°C and seven of 20 at 41°C had p < 0.05. One random set at 32°C had a p-value below that of the Principal Response set; none did at 41°C.

These are comparisons across different random gene sets, not repeated queries of the same input. Matching on Human Protein Atlas degree does not ensure matching on STRING degree. The original random gene lists and STRING identifier mappings would be needed to rerun these queries.

The [CSV](../reference/string_random_controls.csv) contains the same 42 p-values. This descriptive comparison is separate from the framework's 1,000-set network randomizations and does not use a t-test on p-values.
