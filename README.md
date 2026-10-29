# Condition-Informed PPI Network Evaluation

**A Condition-Informed Framework for Evaluating Temperature-Responsive Gene Sets Using Protein Interaction Data**.

This algorithm evaluates selected differentially expressed gene (DEG) sets by constructing protein interaction networks in which shared intermediaries are restricted to the condition-specific differentially expressed genes. It compares these networks with networks allowing unrestricted intermediaries and with approximately degree-matched random gene sets.

## Quick start

Use Python **3.12**. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run.py --verify --plots
```

On Windows, activate the environment with `.venv\Scripts\activate`.

The run uses bundled, frozen inputs and makes **no network requests**. Outputs go to `results/`; source data and reference results are never overwritten. To choose another output directory:

```bash
python run.py --output my_results --verify --plots
```

`--verify` compares regenerated numerical results and randomization draws with the saved reference results. Graph semantics can also be checked independently:

```bash
python -m unittest discover -s tests -v
```

## Where to find things

| Folder/file | Contents |
|---|---|
| `run.py` | One entry point for the reported computational analyses |
| `src/` | Network construction, sampling, statistical tests, and plots |
| `data/` | Frozen interaction cache, gene sets, DEG bridge lists, and housekeeping IDs |
| `reference/` | Final numerical results and draws used for verification |
| `docs/analysis_guide.md` | Analysis-to-manuscript mapping and output definitions |
| `docs/data_sources.md` | Input provenance and identifier conventions |
| `docs/string_comparison.md` | Reported STRING results and original website settings |
| `archive/degree_shift/` | Figure S1 degree-shift analysis, runner, and figure panels |
| `tests/` | Independent checks of graph construction and LCC calculations |

## Included analyses

- Principal Response (G1+G2+G4), Support Group (G3), and both Principal-Response-sized and G3-sized limma subsets at 32°C and 41°C.
- Second-order-only and combined first- and second-order networks.
- Paired unrestricted/restricted comparisons using 1,000 identical random sets per gene set. ΔLCC is unrestricted LCC size minus DEG-restricted LCC size.
- Upper-tail empirical LCC-size and LCC-edge-count enrichment tests, plus the reported lower-tail comparisons.
- Iteration sensitivity at 100, 500, 1,000, and 1,500 nested draws from the same sequence.
- One-sided Fisher tests for G3 housekeeping allocation and background enrichment.
- G3 median-degree tests using 1,000 size-matched random sets, without degree matching, before and after removing housekeeping genes.
- Condition-specific descriptive degree statistics and regenerated degree-distribution plots.

The frozen cache contains **9,825 candidate genes**. Networks use 20 observed-degree quantile bins for approximate matching. All sets at a temperature use the same full DEG bridge list. Bridges establish edges between selected nodes; they are not added as extra graph nodes. First-order links are undirected if either endpoint records the other. Parallel edges and self-edges are not added. Activating/inhibitory effects and expression direction are not assigned.

## Findings represented by these files

At 32°C, DEG-restricted LCC sizes were larger than the random expectation for Principal Response, G3, and Principal-Response-sized limma sets. G3 also had enriched LCC edge counts. Principal Response genes preferentially retained connectivity under DEG restriction; the G3 and limma retention tests were not significant. Corresponding restricted-network and retention comparisons at 41°C were not significant. G3-sized limma subsets showed no significant restricted LCC-size, restricted edge-count, or connectivity-retention enrichment at either temperature; their unrestricted 32°C LCC edge counts were enriched. Housekeeping genes were disproportionately allocated to G3 at both temperatures; background enrichment was significant only at 41°C.


## Citation and licensing

Bradley Frishman, Jorge Gonzalez, and Valery Forbes. *A Condition-Informed Framework for Evaluating Temperature-Responsive Gene Sets Using Protein Interaction Data*. Working manuscript.

Code is distributed under the existing [MIT license](LICENSE). Third-party database content remains subject to its source terms; see [data sources](docs/data_sources.md).
