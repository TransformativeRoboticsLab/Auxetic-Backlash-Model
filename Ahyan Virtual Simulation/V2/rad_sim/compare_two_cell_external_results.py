from __future__ import annotations

import argparse
from pathlib import Path

from .models import LatticeConfig
from .two_cell_bench import (
    compare_two_cell_external_results,
    export_two_cell_external_comparison_json,
    two_cell_external_results_from_csv,
    two_cell_external_results_from_json,
)


def load_two_cell_external_results(path: str | Path) -> list[dict[str, object]]:
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if source.suffix.lower() == ".json":
        return two_cell_external_results_from_json(text)
    return two_cell_external_results_from_csv(text)


def write_two_cell_external_comparison_artifacts(
    input_path: str | Path,
    out: str | Path,
    *,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
    tolerance: float = 1e-6,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = load_two_cell_external_results(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    files = {"comparison": output_dir / "two_cell_external_comparison.json"}
    files["comparison"].write_text(
        export_two_cell_external_comparison_json(results, config, tolerance=tolerance),
        encoding="utf-8",
    )
    return files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare external two-cell RAD results against the quasistatic solver."
    )
    parser.add_argument("--input", required=True, help="Filled two-cell external result CSV or JSON.")
    parser.add_argument("--out", default="outputs/two_cell_external_comparison")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    files = write_two_cell_external_comparison_artifacts(
        args.input,
        args.out,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
        tolerance=args.tolerance,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
