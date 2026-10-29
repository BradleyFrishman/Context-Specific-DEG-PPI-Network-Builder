# Degree-shifted reference scenarios (Figure S1)

## Analysis

[degree_shift_original.py](degree_shift_original.py) computes the four degree-shifted reference scenarios in Figure S1 using DEG-restricted second-order networks and 20 random gene sets per scenario.

The supplied configuration uses the 32°C Principal Response (G1+G2+G4) list, **20 iterations per scenario**, and 8 worker threads. For 41°C, select the corresponding Principal Response interaction file and full condition-specific DEG list in the user parameters.

| Scenario | Implemented change to degree-bin counts |
|---|---|
| Conservative Up | Move the entire count in the lowest occupied degree bin to degree + 1 |
| Aggressive Up | Move every bin's count to degree + 1 |
| Conservative Down | Move the entire count in the lowest occupied degree bin to degree − 1, floored at zero |
| Aggressive Down | Move every bin's count to degree − 1, floored at zero |

Here, bins are exact integer counts of recorded interaction partners, so a one-bin shift is one interaction. These are not the 20 quantile bins used by the main pipeline. The conservative implementation selects the lowest bin with a positive gene count; it does not explicitly exclude degree zero. The shifts alter target bin counts, rather than modifying gene interaction lists.

Each scenario samples background genes without replacement within a draw. If a target degree cannot be filled, sampling falls back to successively lower available degrees. If no lower degree remains, the script reports unsampled genes. An edge joins two selected genes when they share at least one recorded interaction partner in the full condition-specific DEG list. The network is second-order only and DEG-restricted; intermediary genes are not added as graph nodes. Only endpoints of qualifying edges enter the graph. The script prints the observed LCC and the mean randomized LCC for each scenario.

## Run both temperatures

From the repository root, install the dependencies and run:

```bash
python -m pip install -r requirements.txt
python archive/degree_shift/run_degree_shift.py
```

The runner prepares the Principal Response interaction tables and headerless DEG lists from the bundled data, then executes the original script's algorithm with 20 iterations and 8 workers. Only the three input paths are substituted in memory. Results are written to `results/degree_shift/comparison.csv`, with a console log for each temperature. Use `--output PATH` to select another output directory.

The original script can also be run directly after setting `observed_file`, `background_file`, and `long_list_file`. The first two require `ENSEMBL_ID` and semicolon-separated `KNOWN_INTERACTIONS` columns; the DEG list is a headerless CSV with identifiers in its first column.

## Figure S1 and rerun comparison

`32C.png` and `41C.png` are Figure S1 panels A and B. Their numerical values are in `reported_means.csv`; `run.py --plots` uses this table to render the figure.

A fresh 20-draw run is recorded in `rerun_comparison.csv`. It reproduces the observed LCC sizes exactly (133 at 32°C and 140 at 41°C). Randomized means vary between runs because the original script samples without a fixed seed. The comparison table gives both the Figure S1 values and the fresh estimates; it does not replace the Figure S1 values.

The degree-shift scenarios use 20 draws per scenario. The main enrichment tests use 1,000 draws; the separate iteration-sensitivity analysis evaluates 100, 500, 1,000, and 1,500 draws.
