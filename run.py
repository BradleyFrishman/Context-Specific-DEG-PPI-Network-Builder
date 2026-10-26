#!/usr/bin/env python3
"""Reproduce the manuscript's computational analyses from frozen inputs."""
import argparse
import json
import platform
from pathlib import Path
import numpy
import scipy
from src.common import ROOT, verify_inputs
from src.analyses import network_analysis, housekeeping_and_degree, verify_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results', help='Output directory (default: results/)')
    parser.add_argument('--verify', action='store_true', help='Compare regenerated numerical outputs with reference files')
    parser.add_argument('--plots', action='store_true', help='Also regenerate summary and degree plots')
    args = parser.parse_args()
    out = args.output.resolve()
    if out == ROOT or ROOT / 'data' == out or ROOT / 'reference' == out or ROOT / 'src' == out:
        parser.error('Choose a separate output directory')
    for protected in ('data', 'reference', 'src', 'docs', 'archive', 'tests'):
        if (ROOT / protected) in out.parents:
            parser.error('Output must not be inside an input/source directory')
    out.mkdir(parents=True, exist_ok=True)
    manifest = verify_inputs()
    network_analysis(out)
    housekeeping_and_degree(out)
    matches = verify_results(out) if args.verify else None
    if args.plots:
        from src.plots import make_plots
        make_plots(out)
    (out / 'run_metadata.json').write_text(json.dumps(dict(
        python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
        network_seed=20261011, network_iterations=1000,
        sensitivity_checkpoints=[100, 500, 1000, 1500], degree_iterations=1000,
        reference_values_verified=matches, inputs=manifest), indent=2) + '\n')
    print(f'Finished: {out}', flush=True)


if __name__ == '__main__':
    main()
