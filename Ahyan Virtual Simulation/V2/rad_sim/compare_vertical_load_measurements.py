from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

from .cell_geometry import (
    RADHardwareProfile,
    config_with_hardware_profile,
    export_hardware_profile_json,
    hardware_profile_from_config,
    hardware_profile_from_json,
)
from .models import LatticeConfig, LoadCase
from .operators import (
    compare_vertical_load_energy_measurement_results,
    export_vertical_load_energy_comparison_report_json,
    vertical_load_energy_measurement_results_from_json,
)
from .unit_scale import physical_unit_scale_metadata


def _grid_payload(results: dict[str, object]) -> dict[str, object]:
    grid = results.get("grid", {})
    return grid if isinstance(grid, dict) else {}


def _optional_float(value: float | None, fallback: object, default: float) -> float:
    if value is not None:
        return float(value)
    if fallback is None:
        return float(default)
    return float(fallback)


def config_from_vertical_load_measurement_results(
    results: dict[str, object],
    rows: int | None = None,
    cols: int | None = None,
    backlash: float = 0.1,
    cell_size: float = 1.0,
    initial_alpha: float = 1.0,
    coupling_gain: float = 0.55,
    z_coupling_gain: float | None = None,
    pin_radius: float | None = None,
    hole_radius: float | None = None,
) -> LatticeConfig:
    """Infer a dimensionless config from filled vertical-load measurement JSON."""

    grid = _grid_payload(results)
    inferred_rows = int(rows if rows is not None else grid.get("rows", 1))
    inferred_cols = int(cols if cols is not None else grid.get("cols", 1))
    inferred_z_gain = _optional_float(
        z_coupling_gain,
        grid.get("zCouplingGain"),
        0.32,
    )
    clearance = grid.get("pinHoleClearance", None)
    if pin_radius is None and hole_radius is None and clearance is not None:
        inferred_pin = 0.0
        inferred_hole = float(clearance)
    else:
        inferred_pin = 0.18 if pin_radius is None else float(pin_radius)
        inferred_hole = (
            inferred_pin + float(clearance)
            if hole_radius is None and clearance is not None
            else 0.225
            if hole_radius is None
            else float(hole_radius)
        )
    return LatticeConfig(
        rows=inferred_rows,
        cols=inferred_cols,
        backlash=backlash,
        cell_size=cell_size,
        initial_alpha=initial_alpha,
        coupling_gain=coupling_gain,
        z_coupling_gain=inferred_z_gain,
        pin_radius=inferred_pin,
        hole_radius=inferred_hole,
    )


def _comparison_summary_csv(report: dict[str, object]) -> str:
    output = io.StringIO()
    fieldnames = [
        "scenario",
        "passes_tolerance",
        "all_load_cells_measured",
        "missing_measurements",
        "required_cell_count",
        "measured_cell_count",
        "max_abs_work_or_contact_error",
    ]
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    for scenario in report.get("scenarios", []):
        if not isinstance(scenario, dict):
            continue
        validation = scenario.get("validation", {})
        summary = (
            validation.get("summary", {})
            if isinstance(validation, dict)
            else {}
        )
        if not isinstance(summary, dict):
            summary = {}
        writer.writerow(
            {
                "scenario": scenario.get("name", ""),
                "passes_tolerance": summary.get("passesTolerance", False),
                "all_load_cells_measured": summary.get("allLoadCellsMeasured", False),
                "missing_measurements": summary.get("missingMeasurementCount", 0),
                "required_cell_count": summary.get("requiredCellCount", 0),
                "measured_cell_count": summary.get("measuredCellCount", 0),
                "max_abs_work_or_contact_error": summary.get(
                    "maxAbsWorkOrContactError",
                    0.0,
                ),
            }
        )
    return output.getvalue()


def _readme_text(filenames: dict[str, str]) -> str:
    return "\n".join(
        [
            "# RAD Vertical-Load Measurement Comparison",
            "",
            "This folder contains the simulator comparison for a filled",
            "vertical-load measurement JSON file.",
            "",
            "## Files",
            "",
            f"- `{filenames['report']}`: full claim-labeled comparison report.",
            f"- `{filenames['summary_csv']}`: one-row-per-scenario summary table.",
            f"- `{filenames['unit_scale_metadata']}`: normalized-to-millimeter unit scale metadata with Lean theorem links.",
            f"- `{filenames['hardware_profile']}`: hardware dimensions used to derive the effective normalized model.",
            "",
            "## Interpretation",
            "",
            "- A passing row means the filled measurements match the current normalized simulator equations within tolerance.",
            "- It does not prove calibrated gravity, friction, fixture compliance, or rigid-body contact.",
            "- Keep the original filled measurement JSON with this folder for traceability.",
            "",
        ]
    )


def write_vertical_load_energy_comparison_artifacts(
    input_path: str | Path,
    out_dir: str | Path,
    config: LatticeConfig | None = None,
    load_case: LoadCase | None = None,
    tolerance: float | None = None,
    contact_stiffness: float | None = None,
    hardware_profile: RADHardwareProfile | None = None,
) -> dict[str, str]:
    """Read filled vertical-load measurements and write comparison artifacts."""

    source = Path(input_path)
    results = vertical_load_energy_measurement_results_from_json(
        source.read_text(encoding="utf-8")
    )
    base_config = (
        config
        if config is not None
        else config_from_vertical_load_measurement_results(results)
    )
    resolved_config = (
        config_with_hardware_profile(base_config, hardware_profile)
        if hardware_profile is not None
        else base_config
    )
    report = compare_vertical_load_energy_measurement_results(
        resolved_config,
        results,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=load_case,
    )
    unit_scale_metadata = physical_unit_scale_metadata(
        base_config,
        hardware_profile=hardware_profile,
    )
    report["unitScaleMetadata"] = unit_scale_metadata
    profile = (
        hardware_profile
        if hardware_profile is not None
        else hardware_profile_from_config(resolved_config)
    )
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    filenames = {
        "report": "vertical_load_energy_comparison_report.json",
        "summary_csv": "vertical_load_energy_comparison_summary.csv",
        "unit_scale_metadata": "physical_unit_scale_metadata.json",
        "hardware_profile": "hardware_profile.json",
        "readme": "README.md",
    }
    payloads = {
        "report": export_vertical_load_energy_comparison_report_json(report),
        "summary_csv": _comparison_summary_csv(report),
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


def _parse_cell(value: str) -> tuple[int, int]:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("cell coordinates must use row,col format")
    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("cell coordinates must be integers") from exc


def _load_case_from_args(args: argparse.Namespace) -> LoadCase:
    fixed_cells = (
        tuple(_parse_cell(value) for value in args.fixed_cell)
        if args.fixed_cell
        else ((0, 0),)
    )
    return LoadCase(
        fixed_cells=fixed_cells,
        axial_stiffness=args.axial_stiffness,
        hinge_stiffness=args.hinge_stiffness,
        lock_stiffness=args.lock_stiffness,
        maxiter=args.maxiter,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare filled RAD vertical-load bench measurements."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="filled vertical_load_energy_measurement_template JSON file",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="output folder for comparison report and summary CSV",
    )
    parser.add_argument("--rows", type=int, default=None, help="override row count")
    parser.add_argument("--cols", type=int, default=None, help="override column count")
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
        default=None,
        help="override vertical residual propagation gain",
    )
    parser.add_argument(
        "--pin-radius",
        type=float,
        default=None,
        help="override normalized pin radius",
    )
    parser.add_argument(
        "--hole-radius",
        type=float,
        default=None,
        help="override normalized hole radius",
    )
    parser.add_argument(
        "--fixed-cell",
        action="append",
        help="fixed cell as row,col; may be repeated; defaults to 0,0",
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
        default=None,
        help="override contact proxy stiffness from the measurement file",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="override tolerance from the measurement file",
    )
    parser.add_argument(
        "--hardware-profile",
        help="optional rad-sim.hardware-profile.v1 JSON file with measured dimensions",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    results = vertical_load_energy_measurement_results_from_json(
        Path(args.input).read_text(encoding="utf-8")
    )
    config = config_from_vertical_load_measurement_results(
        results,
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
    hardware_profile = (
        hardware_profile_from_json(Path(args.hardware_profile).read_text(encoding="utf-8"))
        if args.hardware_profile
        else None
    )
    written = write_vertical_load_energy_comparison_artifacts(
        args.input,
        args.out,
        config=config,
        load_case=_load_case_from_args(args),
        tolerance=args.tolerance,
        contact_stiffness=args.contact_stiffness,
        hardware_profile=hardware_profile,
    )
    print("Wrote RAD vertical-load measurement comparison artifacts:")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
