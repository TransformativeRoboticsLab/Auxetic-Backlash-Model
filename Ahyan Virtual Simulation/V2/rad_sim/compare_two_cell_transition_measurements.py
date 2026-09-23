from __future__ import annotations

import argparse

from .two_cell_calibration import (
    write_two_cell_radius_backlash_transition_comparison_artifacts,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare filled two-cell RAD radius/backlash transition measurements "
            "against the reduced transition-bracket report."
        )
    )
    parser.add_argument("--input", required=True, help="Filled transition measurement CSV or JSON.")
    parser.add_argument("--out", default="outputs/two_cell_radius_backlash_transition_comparison")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    parser.add_argument("--hole-steps", type=int, default=9)
    parser.add_argument("--backlash-steps", type=int, default=9)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    files = write_two_cell_radius_backlash_transition_comparison_artifacts(
        args.input,
        args.out,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
        hole_steps=args.hole_steps,
        backlash_steps=args.backlash_steps,
        tolerance=args.tolerance,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
