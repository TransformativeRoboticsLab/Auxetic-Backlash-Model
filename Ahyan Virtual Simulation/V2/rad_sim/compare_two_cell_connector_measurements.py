from __future__ import annotations

import argparse

from .two_cell_calibration import (
    write_two_cell_connector_measurement_comparison_artifacts,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare measured two-cell RAD connector data against the A360 CAD-derived proxy."
    )
    parser.add_argument("--input", required=True, help="Filled connector measurement CSV or JSON.")
    parser.add_argument("--out", default="outputs/two_cell_connector_measurement_comparison")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    files = write_two_cell_connector_measurement_comparison_artifacts(
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
