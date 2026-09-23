from __future__ import annotations

import argparse

from .two_cell_calibration import write_two_cell_external_fidelity_benchmark_summary_artifact


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a two-cell external-fidelity MuJoCo benchmark run.")
    parser.add_argument("--run", required=True, help="two_cell_external_fidelity_matrix_mjcf_run.json")
    parser.add_argument("--comparison", required=True, help="two_cell_fidelity_matrix_measurement_comparison.json")
    parser.add_argument("--calibration", default="", help="two_cell_fidelity_matrix_parameter_calibration.json")
    parser.add_argument("--out", default="outputs/two_cell_external_fidelity_benchmark_summary")
    args = parser.parse_args()

    files = write_two_cell_external_fidelity_benchmark_summary_artifact(
        args.run,
        args.comparison,
        args.calibration or None,
        args.out,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
