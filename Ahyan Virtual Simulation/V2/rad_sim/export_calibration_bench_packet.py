from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cell_geometry import (
    RADHardwareProfile,
    config_with_hardware_profile,
    export_hardware_profile_json,
    hardware_profile_from_config,
    hardware_profile_from_json,
)
from .experiments import (
    CalibrationExperimentProtocol,
    build_calibration_experiment_protocol,
    calibration_bench_packet,
)
from .models import LatticeConfig


def _readme_text(filenames: dict[str, str]) -> str:
    return "\n".join(
        [
            "# RAD Calibration Bench Packet",
            "",
            "This folder contains simulator-derived planning artifacts for RAD",
            "calibration experiments. Treat it as a bench handoff, not as",
            "calibrated hardware validation.",
            "",
            "## Files",
            "",
            f"- `{filenames['packet']}`: complete JSON packet with notebook, protocol, templates, and validation guidance.",
            f"- `{filenames['notebook']}`: bench notebook with instruments, phases, dataset plan, and scenario rows.",
            f"- `{filenames['notebook_csv']}`: one-row-per-scenario CSV table for a lab notebook or spreadsheet.",
            f"- `{filenames['protocol']}`: machine-readable calibration protocol used by simulator comparison code.",
            f"- `{filenames['fit_template']}`: blank fit measurement rows. Fill this before fitting a profile.",
            f"- `{filenames['holdout_template']}`: blank holdout rows tied to the frozen profile ID.",
            f"- `{filenames['hardware_profile']}`: hardware dimensions used to derive the effective normalized model.",
            "",
            "## Bench Notes",
            "",
            "- Keep the fit and holdout raw files separate.",
            "- Do not inspect or tune against the holdout rows before freezing the profile.",
            "- Set `profileFrozenAt` after the fit profile is generated and before holdout collection.",
            "- Replace null measured fields with bench data; do not copy simulator predictions into measured rows.",
            "- Passing residual checks supports calibration bookkeeping only, not contact, friction, stiffness, or material-law truth.",
            "",
        ]
    )


def write_calibration_bench_packet_artifacts(
    out_dir: str | Path,
    config: LatticeConfig | None = None,
    protocol: CalibrationExperimentProtocol | None = None,
    *,
    center_cell: tuple[int, int] | None = None,
    alpha_step: float = -0.25,
    z_step: float = 0.30,
    repeat_count: int = 3,
    include_lock_control: bool = True,
    fit_dataset_id: str = "fit-run-001",
    holdout_dataset_id: str = "holdout-run-001",
    profile_id: str = "calibration-profile-v1",
    profile_frozen_at: str | None = None,
    hardware_profile: RADHardwareProfile | None = None,
) -> dict[str, str]:
    """Write calibration bench packet artifacts and return their file paths."""

    resolved_config = LatticeConfig(rows=3, cols=3) if config is None else config
    effective_config = (
        config_with_hardware_profile(resolved_config, hardware_profile)
        if hardware_profile is not None
        else resolved_config
    )
    resolved_protocol = (
        build_calibration_experiment_protocol(
            effective_config,
            center_cell=center_cell,
            alpha_step=alpha_step,
            z_step=z_step,
            hardware_profile=hardware_profile,
            repeat_count=repeat_count,
            include_lock_control=include_lock_control,
        )
        if protocol is None
        else protocol
    )
    packet = calibration_bench_packet(
        effective_config,
        resolved_protocol,
        fit_dataset_id=fit_dataset_id,
        holdout_dataset_id=holdout_dataset_id,
        profile_id=profile_id,
        profile_frozen_at=profile_frozen_at,
    )
    profile = (
        hardware_profile
        if hardware_profile is not None
        else hardware_profile_from_config(effective_config)
    )
    filenames = dict(packet["filenames"])  # type: ignore[arg-type]
    filenames["hardware_profile"] = "hardware_profile.json"
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    payloads: dict[str, str] = {
        "packet": json.dumps(packet, indent=2),
        "notebook": json.dumps(packet["benchNotebook"], indent=2),
        "notebook_csv": str(packet["benchNotebookCsv"]),
        "protocol": json.dumps(packet["experimentProtocol"], indent=2),
        "fit_template": json.dumps(packet["fitResultsTemplate"], indent=2),
        "holdout_template": json.dumps(packet["holdoutResultsTemplate"], indent=2),
        "hardware_profile": export_hardware_profile_json(profile),
        "readme": _readme_text(filenames),
    }
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
        description="Export a RAD calibration bench packet folder."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="output folder for calibration bench JSON, CSV, and README artifacts",
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
        help="protocol center cell as row,col; defaults to grid center",
    )
    parser.add_argument("--alpha-step", type=float, default=-0.25, help="alpha contraction step")
    parser.add_argument("--z-step", type=float, default=0.30, help="vertical actuation step")
    parser.add_argument(
        "--repeat-count",
        type=int,
        default=3,
        help="number of repeated measurements requested in the protocol",
    )
    parser.add_argument(
        "--no-lock-control",
        action="store_true",
        help="omit the locked-cell control protocol step",
    )
    parser.add_argument("--fit-dataset-id", default="fit-run-001", help="fit dataset ID")
    parser.add_argument(
        "--holdout-dataset-id",
        default="holdout-run-001",
        help="holdout dataset ID",
    )
    parser.add_argument(
        "--profile-id",
        default="calibration-profile-v1",
        help="profile ID that the holdout template should reference",
    )
    parser.add_argument(
        "--profile-frozen-at",
        help="profile freeze timestamp to embed in the holdout template",
    )
    parser.add_argument(
        "--hardware-profile",
        help="optional rad-sim.hardware-profile.v1 JSON file with measured dimensions",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    hardware_profile = (
        hardware_profile_from_json(Path(args.hardware_profile).read_text(encoding="utf-8"))
        if args.hardware_profile
        else None
    )
    written = write_calibration_bench_packet_artifacts(
        args.out,
        config=_config_from_args(args),
        center_cell=args.center_cell,
        alpha_step=args.alpha_step,
        z_step=args.z_step,
        repeat_count=args.repeat_count,
        include_lock_control=not args.no_lock_control,
        fit_dataset_id=args.fit_dataset_id,
        holdout_dataset_id=args.holdout_dataset_id,
        profile_id=args.profile_id,
        profile_frozen_at=args.profile_frozen_at,
        hardware_profile=hardware_profile,
    )
    print("Wrote RAD calibration bench packet artifacts:")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
