from __future__ import annotations

import json
from fractions import Fraction

from .cell_geometry import (
    PAPER_RAD_REFERENCE,
    PaperRADReference,
    RADHardwareProfile,
    calibrate_paper_rad_config,
    calibration_readiness,
    config_with_hardware_profile,
    hardware_profile_from_config,
)
from .models import LatticeConfig


UNIT_SCALE_LEAN_THEOREMS: tuple[str, ...] = (
    "measurementUnitScaleNat_zero",
    "measurementUnitScaleNat_preserves_denominator",
    "measurementUnitScaleInt_zero",
    "measurementUnitScaleInt_preserves_denominator",
    "measurementUnitScaleResidualNat_zero_when_equal",
    "measurementUnitScaleResidualNat_preserves_denominator",
    "measurementUnitScaleResidualInt_zero_when_equal",
    "measurementUnitScaleResidualInt_preserves_denominator",
)
HARDWARE_PROFILE_LEAN_THEOREMS: tuple[str, ...] = (
    "hardwareProfileCoverageComplete_true_when_equal",
    "hardwareProfileCoverageMissing_zero_when_complete",
    "hardwareProfileCoverageMissing_preserves_total",
)


def _scaled_value(value: float, max_denominator: int) -> dict[str, object]:
    fraction = Fraction(float(value)).limit_denominator(max_denominator)
    return {
        "value": float(value),
        "numerator": int(fraction.numerator),
        "denominator": int(fraction.denominator),
    }


def _scaled_mm_quantity(
    value_mm: float,
    max_denominator: int,
) -> dict[str, object]:
    payload = _scaled_value(value_mm, max_denominator)
    payload["unit"] = "mm"
    return payload


def physical_unit_scale_metadata(
    config: LatticeConfig,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
    hardware_profile: RADHardwareProfile | None = None,
    max_denominator: int = 1_000_000,
) -> dict[str, object]:
    """Return claim-labeled normalized-to-physical unit-scale metadata."""

    if max_denominator <= 0:
        raise ValueError("max_denominator must be positive")
    effective_reference = (
        hardware_profile.to_reference() if hardware_profile is not None else reference
    )
    effective_config = (
        config_with_hardware_profile(config, hardware_profile)
        if hardware_profile is not None
        else config
    )
    calibration = calibrate_paper_rad_config(effective_config, effective_reference)
    profile = (
        hardware_profile
        if hardware_profile is not None
        else hardware_profile_from_config(effective_config, effective_reference)
    )
    readiness = calibration_readiness(profile)
    mm_per_model_unit = _scaled_value(
        calibration.mm_per_model_unit,
        max_denominator,
    )
    model_units_per_mm = _scaled_value(
        1.0 / calibration.mm_per_model_unit,
        max_denominator,
    )
    return {
        "schema": "rad-sim.physical-unit-scale-metadata.v1",
        "method": (
            "normalized simulator length quantities are exported with explicit "
            "millimeter scale metadata and Lean theorem links for zero and "
            "denominator-preservation invariants"
        ),
        "claimLabels": {
            "unitScaleInvariants": "Lean-proven theorem scaffold",
            "paperReferenceScale": "paper-supported assumption",
            "hardwareProfileValues": "experimentally unvalidated physical assumption",
            "simulatorConfiguration": "simulator-derived empirical law",
        },
        "formalizationTarget": {
            "id": "measurement_unit_scale_invariants",
            "leanTheorems": list(UNIT_SCALE_LEAN_THEOREMS),
            "hardwareProfileLeanTheorems": list(HARDWARE_PROFILE_LEAN_THEOREMS),
            "claimLimit": (
                "The theorems prove finite unit-scale metadata invariants only; "
                "they do not prove that the selected physical scale is calibrated."
            ),
        },
        "grid": {
            "rows": effective_config.rows,
            "cols": effective_config.cols,
            "cellSize": effective_config.cell_size,
            "backlash": effective_config.backlash,
            "pinRadius": effective_config.pin_radius,
            "holeRadius": effective_config.hole_radius,
            "pinHoleClearance": effective_config.pin_hole_clearance,
            "zCouplingGain": effective_config.z_coupling_gain,
        },
        "profileApplication": {
            "applied": hardware_profile is not None,
            "inputBacklash": config.backlash,
            "effectiveBacklash": effective_config.backlash,
            "inputPinRadius": config.pin_radius,
            "effectivePinRadius": effective_config.pin_radius,
            "inputHoleRadius": config.hole_radius,
            "effectiveHoleRadius": effective_config.hole_radius,
        },
        "scale": {
            "modelUnitToMillimeter": {
                **mm_per_model_unit,
                "unit": "mm/model-unit",
            },
            "millimeterToModelUnit": {
                **model_units_per_mm,
                "unit": "model-unit/mm",
            },
            "maxDenominator": int(max_denominator),
        },
        "paperReference": {
            "sideLengthMm": reference.side_length_mm,
            "normalizedBacklash": reference.normalized_backlash,
            "referenceBacklashMm": reference.reference_backlash_mm,
            "fabricationHoleToleranceMm": reference.fabrication_hole_tolerance_mm,
        },
        "configuredPhysicalQuantities": {
            "sideLengthMm": _scaled_mm_quantity(
                calibration.side_length_mm,
                max_denominator,
            ),
            "configuredBacklashMm": _scaled_mm_quantity(
                calibration.configured_backlash_mm,
                max_denominator,
            ),
            "pinRadiusMm": _scaled_mm_quantity(
                calibration.pin_radius_mm,
                max_denominator,
            ),
            "holeRadiusMm": _scaled_mm_quantity(
                calibration.hole_radius_mm,
                max_denominator,
            ),
            "pinHoleClearanceMm": _scaled_mm_quantity(
                calibration.pin_hole_clearance_mm,
                max_denominator,
            ),
            "fabricationHoleToleranceMm": _scaled_mm_quantity(
                calibration.fabrication_hole_tolerance_mm,
                max_denominator,
            ),
        },
        "hardwareProfile": profile.to_dict(),
        "zeroAndResidualExamples": {
            "zeroModelLengthToMm": {
                "modelValue": 0.0,
                "scaledNumerator": 0,
                "denominator": mm_per_model_unit["denominator"],
                "theorem": "measurementUnitScaleNat_zero",
            },
            "zeroSignedCommandToPhysicalUnit": {
                "modelValue": 0,
                "scaledNumerator": 0,
                "denominator": mm_per_model_unit["denominator"],
                "theorem": "measurementUnitScaleInt_zero",
            },
            "equalUnsignedResidualToPhysicalUnit": {
                "residualNumerator": 0,
                "denominator": mm_per_model_unit["denominator"],
                "theorem": "measurementUnitScaleResidualNat_zero_when_equal",
            },
            "equalSignedResidualToPhysicalUnit": {
                "residualNumerator": 0,
                "denominator": mm_per_model_unit["denominator"],
                "theorem": "measurementUnitScaleResidualInt_zero_when_equal",
            },
        },
        "calibrationReadiness": {
            "profileName": readiness.profile_name,
            "level": readiness.level,
            "measuredFields": list(readiness.measured_fields),
            "missingFields": list(readiness.missing_fields),
            "solverGaps": list(readiness.solver_gaps),
            "coverageRatio": readiness.coverage_ratio,
            "summary": readiness.summary,
        },
        "artifactLinks": {
            "benchPacketFunction": "write_vertical_load_bench_packet_artifacts",
            "comparisonFunction": "write_vertical_load_energy_comparison_artifacts",
            "calibrationFunction": "calibrate_paper_rad_config",
        },
        "limitations": [
            "This is unit-scale metadata, not calibrated physical validation.",
            "Millimeter values are derived from the selected paper/reference scale unless a measured profile is introduced.",
            "Force, stiffness, contact, friction, and gravity scales remain uncalibrated.",
        ],
    }


def export_physical_unit_scale_metadata_json(
    config: LatticeConfig,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
    hardware_profile: RADHardwareProfile | None = None,
    max_denominator: int = 1_000_000,
) -> str:
    return json.dumps(
        physical_unit_scale_metadata(
            config,
            reference=reference,
            hardware_profile=hardware_profile,
            max_denominator=max_denominator,
        ),
        indent=2,
    )
