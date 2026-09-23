from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cell_geometry import (
    RADHardwareProfile,
    config_with_hardware_profile,
    hardware_profile_from_json,
)
from .experiments import (
    CalibrationExperimentProtocol,
    build_calibration_experiment_protocol,
    calibration_bench_execution_validation,
    calibration_experiment_measurements_from_json,
    calibration_experiment_protocol_from_json,
    export_calibration_bench_execution_validation_csv,
    export_calibration_model_profile_holdout_validation_csv,
)
from .models import LatticeConfig


def _readme_text(filenames: dict[str, str]) -> str:
    return "\n".join(
        [
            "# RAD Calibration Bench Execution Validation",
            "",
            "This folder contains reports generated after ingesting filled fit",
            "and holdout calibration measurement files. It is a validation",
            "bookkeeping artifact, not a proof of physical accuracy.",
            "",
            "## Files",
            "",
            f"- `{filenames['report']}`: complete execution-validation report.",
            f"- `{filenames['summary_csv']}`: one-row validation summary.",
            f"- `{filenames['holdout_validation']}`: detailed fit/holdout residual validation.",
            f"- `{filenames['holdout_csv']}`: two-row fit/holdout residual table.",
            f"- `{filenames['model_profile']}`: model profile used for holdout replay.",
            f"- `{filenames['fit_comparison_report']}`: fit-file comparison and parameter estimates.",
            "",
            "## Interpretation",
            "",
            "- `residualValidationPass` means the bounded simulator profile improved both fit and holdout residual bookkeeping.",
            "- `independentValidationPass` additionally requires complete split provenance metadata.",
            "- Physical contact, friction, stiffness, gravity, and material-law accuracy remain unvalidated until measured separately.",
            "",
        ]
    )


def write_calibration_bench_execution_validation_artifacts(
    out_dir: str | Path,
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    fit_results_json: str,
    holdout_results_json: str,
    *,
    profile: dict[str, object] | None = None,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> dict[str, str]:
    """Write executed calibration validation artifacts for filled bench files."""

    fit_measurements = calibration_experiment_measurements_from_json(fit_results_json)
    holdout_measurements = calibration_experiment_measurements_from_json(
        holdout_results_json
    )
    report = calibration_bench_execution_validation(
        config,
        protocol,
        fit_measurements,
        holdout_measurements,
        profile,
        tolerance=tolerance,
        min_improvement=min_improvement,
    )
    holdout_validation = report["holdoutValidation"]
    filenames = {
        "report": "calibration_bench_execution_validation.json",
        "summary_csv": "calibration_bench_execution_validation_summary.csv",
        "holdout_validation": "calibration_model_profile_holdout_validation.json",
        "holdout_csv": "calibration_model_profile_holdout_validation.csv",
        "model_profile": "calibration_model_profile.json",
        "fit_comparison_report": "calibration_fit_comparison_report.json",
        "readme": "README.md",
    }
    payloads = {
        "report": json.dumps(report, indent=2),
        "summary_csv": export_calibration_bench_execution_validation_csv(report),
        "holdout_validation": json.dumps(holdout_validation, indent=2),
        "holdout_csv": export_calibration_model_profile_holdout_validation_csv(
            holdout_validation
        ),
        "model_profile": json.dumps(report["modelProfile"], indent=2),
        "fit_comparison_report": json.dumps(report["fitComparisonReport"], indent=2),
        "readme": _readme_text(filenames),
    }
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}
    for key, filename in filenames.items():
        target = out_path / filename
        target.write_text(payloads[key], encoding="utf-8")
        written[key] = str(target)
    return written


def _parse_cell(value: str) -> tuple[int, int]:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("cell coordinates must use row,col format")
    try:
        row, col = (int(parts[0]), int(parts[1]))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("cell coordinates must be integers") from exc
    return (row, col)


def _config_from_args(args: argparse.Namespace) -> LatticeConfig:
    return LatticeConfig(
        rows=args.rows,
        cols=args.cols,
        backlash=args.backlash,
        cell_size=args.cell_size,
        initial_alpha=args.initial_alpha,
        coupling_gain=args.coupling_gain,
        z_coupling_gain=args.z_coupling_gain,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare filled RAD calibration fit/holdout bench files."
    )
    parser.add_argument("--fit", required=True, help="filled fit results JSON")
    parser.add_argument("--holdout", required=True, help="filled holdout results JSON")
    parser.add_argument("--out", required=True, help="output folder for validation reports")
    parser.add_argument(
        "--protocol",
        help="calibration_experiment_protocol.json from the bench packet; defaults to a generated protocol",
    )
    parser.add_argument("--profile", help="optional calibration model profile JSON")
    parser.add_argument(
        "--hardware-profile",
        help="optional rad-sim.hardware-profile.v1 JSON file with measured dimensions",
    )
    parser.add_argument("--rows", type=int, default=3, help="lattice row count")
    parser.add_argument("--cols", type=int, default=3, help="lattice column count")
    parser.add_argument("--backlash", type=float, default=0.1, help="dead-zone backlash")
    parser.add_argument("--cell-size", type=float, default=1.0, help="normalized cell size")
    parser.add_argument(
        "--initial-alpha",
        type=float,
        default=1.0,
        help="initial normalized alpha/dilation value",
    )
    parser.add_argument(
        "--coupling-gain",
        type=float,
        default=0.55,
        help="alpha propagation coupling gain",
    )
    parser.add_argument(
        "--z-coupling-gain",
        type=float,
        default=0.32,
        help="vertical residual propagation gain",
    )
    parser.add_argument("--pin-radius", type=float, default=0.18, help="normalized pin radius")
    parser.add_argument("--hole-radius", type=float, default=0.225, help="normalized hole radius")
    parser.add_argument(
        "--center-cell",
        type=_parse_cell,
        help="generated protocol center cell as row,col; ignored when --protocol is supplied",
    )
    parser.add_argument("--alpha-step", type=float, default=-0.25, help="generated protocol alpha step")
    parser.add_argument("--z-step", type=float, default=0.30, help="generated protocol z step")
    parser.add_argument("--repeat-count", type=int, default=3, help="generated protocol repeat count")
    parser.add_argument(
        "--no-lock-control",
        action="store_true",
        help="omit generated locked-cell control step",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="measurement comparison tolerance",
    )
    parser.add_argument(
        "--min-improvement",
        type=float,
        default=0.0,
        help="minimum residual-score improvement required for profile pass",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    hardware_profile: RADHardwareProfile | None = (
        hardware_profile_from_json(Path(args.hardware_profile).read_text(encoding="utf-8"))
        if args.hardware_profile
        else None
    )
    config = _config_from_args(args)
    effective_config = (
        config_with_hardware_profile(config, hardware_profile)
        if hardware_profile is not None
        else config
    )
    protocol = (
        calibration_experiment_protocol_from_json(
            Path(args.protocol).read_text(encoding="utf-8")
        )
        if args.protocol
        else build_calibration_experiment_protocol(
            effective_config,
            center_cell=args.center_cell,
            alpha_step=args.alpha_step,
            z_step=args.z_step,
            hardware_profile=hardware_profile,
            repeat_count=args.repeat_count,
            include_lock_control=not args.no_lock_control,
        )
    )
    profile = json.loads(Path(args.profile).read_text(encoding="utf-8")) if args.profile else None
    written = write_calibration_bench_execution_validation_artifacts(
        args.out,
        effective_config,
        protocol,
        Path(args.fit).read_text(encoding="utf-8"),
        Path(args.holdout).read_text(encoding="utf-8"),
        profile=profile,
        tolerance=args.tolerance,
        min_improvement=args.min_improvement,
    )
    print("Wrote RAD calibration execution validation artifacts:")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
