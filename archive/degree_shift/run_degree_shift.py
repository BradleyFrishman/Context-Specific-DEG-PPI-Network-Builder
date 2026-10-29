"""Run Figure S1 degree-shift scenarios with the bundled interaction data."""
import argparse
import ast
import contextlib
import io
from pathlib import Path
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "degree_shift")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source_path = Path(__file__).with_name("degree_shift_original.py")
    background = pd.read_csv(ROOT / "data/background_interactions.csv.gz", dtype=str).fillna("")
    reported = pd.read_csv(Path(__file__).with_name("reported_means.csv"))
    labels = {
        "Conservative Up": "Conservative upshift",
        "Aggressive Up": "Aggressive upshift",
        "Conservative Down": "Conservative downshift",
        "Aggressive Down": "Aggressive downshift",
    }
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        background_path = temp / "background.csv"
        background.to_csv(background_path, index=False)
        for temperature in (32, 41):
            ids = pd.read_csv(ROOT / f"data/gene_sets/principal_{temperature}.csv")["ENSEMBL_ID"]
            observed = background.set_index("ENSEMBL_ID").reindex(ids).dropna().reset_index()
            observed = observed[observed["KNOWN_INTERACTIONS"] != ""]
            observed_path = temp / f"observed_{temperature}.csv"
            observed.to_csv(observed_path, index=False)
            bridge_path = temp / f"bridges_{temperature}.csv"
            pd.read_csv(ROOT / f"data/deg_bridges/{temperature}.csv").to_csv(
                bridge_path, header=False, index=False
            )
            # Replace only input paths in memory; preserve the original algorithm/settings.
            replacements = {
                "observed_file": str(observed_path),
                "background_file": str(background_path),
                "long_list_file": str(bridge_path),
            }
            tree = ast.parse(source_path.read_text(encoding="utf-8"))
            for node in tree.body:
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    target = node.targets[0]
                    if isinstance(target, ast.Name) and target.id in replacements:
                        node.value = ast.Constant(replacements[target.id])
            ast.fix_missing_locations(tree)
            namespace = {"__name__": "__degree_shift__"}
            log = io.StringIO()
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                exec(compile(tree, str(source_path), "exec"), namespace)
            (args.output / f"{temperature}C.log").write_text(log.getvalue(), encoding="utf-8")
            values = {"Observed": namespace["obs_lcc"]}
            values.update({labels[label]: stats["mean_lcc"] for label, stats in namespace["results"].items()})
            for scenario, value in values.items():
                original = reported.loc[
                    (reported["temperature"] == temperature) & (reported["scenario"] == scenario),
                    "lcc_size",
                ].item()
                rows.append({
                    "temperature": temperature, "scenario": scenario,
                    "rerun_lcc_size": value, "figure_s1_lcc_size": original,
                    "difference": value - original,
                    "iterations": namespace["num_iterations"] if scenario != "Observed" else 0,
                })
    output = pd.DataFrame(rows)
    output.to_csv(args.output / "comparison.csv", index=False, float_format="%.2f")
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
