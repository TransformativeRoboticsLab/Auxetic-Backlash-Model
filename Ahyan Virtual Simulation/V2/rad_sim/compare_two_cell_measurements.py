from __future__ import annotations

import argparse
from pathlib import Path

from .models import LatticeConfig
from .two_cell_calibration import (
    export_two_cell_measurement_comparison_json,
    export_two_cell_parameter_calibration_json,
    two_cell_measurements_from_csv,
    two_cell_measurements_from_json,
)


def load_two_cell_measurements(path: str | Path):
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if source.suffix.lower() == ".json":
        return two_cell_measurements_from_json(text)
    return two_cell_measurements_from_csv(text)


def write_two_cell_measurement_comparison_artifacts(
    input_path: str | Path,
    out: str | Path,
    *,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    measurements = load_two_cell_measurements(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    files = {
        "comparison": output_dir / "two_cell_measurement_comparison.json",
        "calibration": output_dir / "two_cell_parameter_calibration.json",
    }
    files["comparison"].write_text(
        export_two_cell_measurement_comparison_json(measurements, config),
        encoding="utf-8",
    )
    files["calibration"].write_text(
        export_two_cell_parameter_calibration_json(measurements, config),
        encoding="utf-8",
    )
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare measured two-cell RAD data against the reduced proxy.")
    parser.add_argument("--input", required=True, help="Filled two-cell measurement CSV or JSON.")
    parser.add_argument("--out", default="outputs/two_cell_measurement_comparison")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    args = parser.parse_args()
    files = write_two_cell_measurement_comparison_artifacts(
        args.input,
        args.out,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
