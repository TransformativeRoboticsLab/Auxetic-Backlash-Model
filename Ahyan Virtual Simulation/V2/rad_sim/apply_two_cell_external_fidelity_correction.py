from __future__ import annotations

import argparse

from .two_cell_calibration import write_two_cell_external_fidelity_correction_application_artifacts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply a two-cell external-fidelity correction profile to the dense prediction table."
    )
    parser.add_argument("--profile", required=True, help="two_cell_external_fidelity_correction_profile.json")
    parser.add_argument("--measurements", default="", help="Optional filled measurement CSV or JSON for evaluation.")
    parser.add_argument("--out", default="outputs/two_cell_external_fidelity_correction_application")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    args = parser.parse_args()

    files = write_two_cell_external_fidelity_correction_application_artifacts(
        args.profile,
        args.out,
        measurements_path=args.measurements or None,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
