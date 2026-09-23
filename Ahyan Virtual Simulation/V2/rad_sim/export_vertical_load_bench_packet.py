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
from .models import LatticeConfig, LatticeState, LoadCase
from .operators import (
    export_vertical_load_physical_preview_report_csv,
    vertical_load_bench_packet,
)
from .unit_scale import physical_unit_scale_metadata


def _default_fixed_cells(rows: int) -> tuple[tuple[int, int], ...]:
    return ((max(0, rows // 2), 0),)


def _readme_text(filenames: dict[str, str]) -> str:
    return "\n".join(
        [
            "# RAD Vertical-Load Bench Packet",
            "",
            "This folder contains simulator-derived planning artifacts for the",
            "vertical-load RAD bench experiments. Treat the files as a measurement",
            "handoff, not as calibrated hardware validation.",
            "",
            "## Files",
            "",
            f"- `{filenames['packet']}`: complete packet with preview, protocol, template, and comparison guidance.",
            f"- `{filenames['preview_json']}`: spring-hinge physical-preview scenarios.",
            f"- `{filenames['preview_csv']}`: bench-planning table for the same scenarios.",
            f"- `{filenames['measurement_template']}`: blank measurement rows to fill after each bench run.",
            f"- `{filenames['experiment_protocol']}`: step-by-step collection protocol.",
            f"- `{filenames['unit_scale_metadata']}`: normalized-to-millimeter unit scale metadata with Lean theorem links.",
            f"- `{filenames['hardware_profile']}`: hardware dimensions used to derive the effective normalized model.",
            "",
            "## Bench Notes",
            "",
            "- Do not copy predictedHeight into measuredHeight.",
            "- Record fixture compliance, visible slip, friction, and backlash changes.",
            "- Keep simulator predictions and measured rows separate until comparison.",
            "- Re-run this export whenever grid, clearance, stiffness, or load settings change.",
            "",
        ]
    )


def write_vertical_load_bench_packet_artifacts(
    out_dir: str | Path,
    config: LatticeConfig | None = None,
    state: LatticeState | None = None,
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    repeat_count: int = 3,
    include_fields: bool = False,
    hardware_profile: RADHardwareProfile | None = None,
) -> dict[str, str]:
    """Write vertical-load bench packet artifacts and return their file paths."""

    resolved_config = LatticeConfig(rows=1, cols=4) if config is None else config
    effective_config = (
        config_with_hardware_profile(resolved_config, hardware_profile)
        if hardware_profile is not None
        else resolved_config
    )
    resolved_state = (
        LatticeState.uniform(effective_config)
        if state is None
        else state.normalized(effective_config)
    )
    resolved_load = (
        LoadCase(fixed_cells=_default_fixed_cells(effective_config.rows))
        if load_case is None
        else load_case
    )
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    packet = vertical_load_bench_packet(
        resolved_config,
        state=resolved_state,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=resolved_load,
        repeat_count=repeat_count,
        include_fields=include_fields,
        hardware_profile=hardware_profile,
    )
    preview = packet["previewReport"]
    measurement_template = packet["measurementTemplate"]
    experiment_protocol = packet["experimentProtocol"]
    unit_scale_metadata = packet.get(
        "unitScaleMetadata",
        physical_unit_scale_metadata(
            resolved_config,
            hardware_profile=hardware_profile,
        ),
    )
    profile = (
        hardware_profile
        if hardware_profile is not None
        else hardware_profile_from_config(effective_config)
    )

    filenames = {
        "packet": "vertical_load_bench_packet.json",
        "preview_json": "vertical_load_physical_preview_report.json",
        "preview_csv": "vertical_load_physical_preview_report.csv",
        "measurement_template": "vertical_load_energy_measurement_template.json",
        "experiment_protocol": "vertical_load_energy_experiment_protocol.json",
        "unit_scale_metadata": "physical_unit_scale_metadata.json",
        "hardware_profile": "hardware_profile.json",
        "readme": "README.md",
    }
    payloads: dict[str, str] = {
        "packet": json.dumps(packet, indent=2),
        "preview_json": json.dumps(preview, indent=2),
        "preview_csv": export_vertical_load_physical_preview_report_csv(preview),
        "measurement_template": json.dumps(measurement_template, indent=2),
        "experiment_protocol": json.dumps(experiment_protocol, indent=2),
        "unit_scale_metadata": json.dumps(unit_scale_metadata, indent=2),
        "hardware_profile": export_hardware_profile_json(profile),
        "readme": _readme_text(filenames),
    }

    written: dict[str, str] = {}
    for key, filename in filenames.items():
        target = out_path / filename
        target.write_text(payloads[key], encoding="utf-8")
        written[key] = str(target)
    return written


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


def _load_case_from_args(args: argparse.Namespace) -> LoadCase:
    fixed_cells: tuple[tuple[int, int], ...]
    if args.fixed_cell is None:
        fixed_cells = _default_fixed_cells(args.rows)
    else:
        fixed_cells = tuple(_parse_cell(value) for value in args.fixed_cell)
    return LoadCase(
        fixed_cells=fixed_cells,
        axial_stiffness=args.axial_stiffness,
        hinge_stiffness=args.hinge_stiffness,
        lock_stiffness=args.lock_stiffness,
        maxiter=args.maxiter,
    )


def _parse_cell(value: str) -> tuple[int, int]:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("cell coordinates must use row,col format")
    try:
        row, col = (int(parts[0]), int(parts[1]))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("cell coordinates must be integers") from exc
    return (row, col)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export a RAD vertical-load bench packet folder."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="output folder for JSON, CSV, and README bench artifacts",
    )
    parser.add_argument("--rows", type=int, default=1, help="lattice row count")
    parser.add_argument("--cols", type=int, default=4, help="lattice column count")
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
        "--fixed-cell",
        action="append",
        help="fixed cell as row,col; may be repeated; defaults to middle-left cell",
    )
    parser.add_argument(
        "--axial-stiffness",
        type=float,
        default=25.0,
        help="spring-hinge axial stiffness",
    )
    parser.add_argument(
        "--hinge-stiffness",
        type=float,
        default=1.0,
        help="spring-hinge rotational stiffness",
    )
    parser.add_argument(
        "--lock-stiffness",
        type=float,
        default=100.0,
        help="fixed-cell penalty stiffness",
    )
    parser.add_argument("--maxiter", type=int, default=500, help="SciPy optimizer iteration cap")
    parser.add_argument(
        "--contact-stiffness",
        type=float,
        default=1.0,
        help="uncalibrated contact proxy stiffness",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="measurement comparison tolerance",
    )
    parser.add_argument(
        "--repeat-count",
        type=int,
        default=3,
        help="number of repeated measurements requested in the protocol",
    )
    parser.add_argument(
        "--include-fields",
        action="store_true",
        help="include field arrays in JSON previews",
    )
    parser.add_argument(
        "--hardware-profile",
        help="optional rad-sim.hardware-profile.v1 JSON file with measured dimensions",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    config = _config_from_args(args)
    load_case = _load_case_from_args(args)
    hardware_profile = (
        hardware_profile_from_json(Path(args.hardware_profile).read_text(encoding="utf-8"))
        if args.hardware_profile
        else None
    )
    written = write_vertical_load_bench_packet_artifacts(
        args.out,
        config=config,
        load_case=load_case,
        tolerance=args.tolerance,
        contact_stiffness=args.contact_stiffness,
        repeat_count=args.repeat_count,
        include_fields=args.include_fields,
        hardware_profile=hardware_profile,
    )
    print("Wrote RAD vertical-load bench packet artifacts:")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
