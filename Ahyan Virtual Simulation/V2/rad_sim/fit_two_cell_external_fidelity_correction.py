from __future__ import annotations

import argparse

from .two_cell_calibration import write_two_cell_external_fidelity_correction_profile_artifact


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fit a train/holdout correction profile from two-cell external-fidelity measurements."
    )
    parser.add_argument("--input", required=True, help="Filled fidelity-matrix measurement CSV or JSON.")
    parser.add_argument("--out", default="outputs/two_cell_external_fidelity_correction_profile")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    parser.add_argument("--holdout-stride", type=int, default=5)
    parser.add_argument(
        "--holdout-case-id",
        action="append",
        default=[],
        help="Specific caseId to reserve for holdout. May be repeated.",
    )
    args = parser.parse_args()

    files = write_two_cell_external_fidelity_correction_profile_artifact(
        args.input,
        args.out,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
        holdout_stride=args.holdout_stride,
        holdout_case_ids=args.holdout_case_id,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
