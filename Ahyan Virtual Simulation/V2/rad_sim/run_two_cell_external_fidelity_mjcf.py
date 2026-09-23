from __future__ import annotations

import argparse
import json
from pathlib import Path

from .two_cell_calibration import (
    write_two_cell_external_fidelity_benchmark_summary_artifact,
    write_two_cell_fidelity_matrix_measurement_comparison_artifacts,
)
from .two_cell_bench import (
    two_cell_external_fidelity_mjcf_measurements_csv,
    two_cell_external_fidelity_mjcf_run_report,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the dense two-cell external-fidelity MJCF matrix with MuJoCo.")
    parser.add_argument(
        "--index",
        default="outputs/two_cell_bench_packet/two_cell_external_fidelity_matrix_mjcf_index.json",
        help="Path to two_cell_external_fidelity_matrix_mjcf_index.json.",
    )
    parser.add_argument("--out", default="outputs/two_cell_external_fidelity_mjcf_run")
    parser.add_argument("--steps", type=int, default=120)
    parser.add_argument("--timestep", type=float, default=0.001)
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument(
        "--case-id",
        action="append",
        default=[],
        help="Specific caseId to run. May be repeated.",
    )
    parser.add_argument(
        "--case-ids",
        default="",
        help="Comma-separated caseIds to run.",
    )
    parser.add_argument(
        "--assume-mujoco-unavailable",
        action="store_true",
        help="Write the not-run audit artifact without importing MuJoCo.",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Compare emitted measurement rows against the reduced fidelity matrix.",
    )
    parser.add_argument(
        "--compare-out",
        default="",
        help="Optional comparison output directory. Defaults to <out>_comparison.",
    )
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)
    availability = {"mujoco": False} if args.assume_mujoco_unavailable else None
    case_ids = [case_id.strip() for case_id in args.case_ids.split(",") if case_id.strip()]
    case_ids.extend(args.case_id)
    report = two_cell_external_fidelity_mjcf_run_report(
        args.index,
        output_dir=output_dir,
        steps=args.steps,
        timestep=args.timestep,
        engine_availability=availability,
        max_cases=args.max_cases,
        case_ids=case_ids,
    )
    target = output_dir / "two_cell_external_fidelity_matrix_mjcf_run.json"
    target.write_text(json.dumps(report, indent=2), encoding="utf-8")
    measurements = output_dir / "two_cell_external_fidelity_matrix_mjcf_measurements.csv"
    measurements.write_text(two_cell_external_fidelity_mjcf_measurements_csv(report), encoding="utf-8")
    print(f"run: {target}")
    print(f"measurements: {measurements}")
    if args.compare:
        comparison_dir = Path(args.compare_out) if args.compare_out else Path(f"{output_dir}_comparison")
        comparison_files = write_two_cell_fidelity_matrix_measurement_comparison_artifacts(
            measurements,
            comparison_dir,
            tolerance=args.tolerance,
        )
        for label, path in comparison_files.items():
            print(f"{label}: {path}")
        benchmark_files = write_two_cell_external_fidelity_benchmark_summary_artifact(
            target,
            comparison_files["comparison"],
            comparison_files["calibration"],
            comparison_dir,
        )
        for label, path in benchmark_files.items():
            print(f"{label}: {path}")


if __name__ == "__main__":
    main()
