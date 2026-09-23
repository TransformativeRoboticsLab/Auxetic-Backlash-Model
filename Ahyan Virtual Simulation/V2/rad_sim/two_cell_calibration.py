from __future__ import annotations

import csv
import io
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .models import LatticeConfig
from .two_cell_bench import (
    TWO_CELL_PACKET_SCHEMA,
    TwoCellBenchControls,
    simulate_two_cell_bench,
    two_cell_connector_measurement_template,
    two_cell_fidelity_matrix_measurement_template,
    two_cell_physical_test_packet,
    two_cell_radius_backlash_transition_report,
)


TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA = "rad-sim.two-cell-measurement-comparison.v1"
TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA = "rad-sim.two-cell-connector-measurement-comparison.v1"
TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA = (
    "rad-sim.two-cell-fidelity-matrix-measurement-comparison.v1"
)
TWO_CELL_PARAMETER_CALIBRATION_SCHEMA = "rad-sim.two-cell-parameter-calibration.v1"
TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA = (
    "rad-sim.two-cell-fidelity-matrix-parameter-calibration.v1"
)
TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA = (
    "rad-sim.two-cell-radius-backlash-transition-comparison.v1"
)
TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA = (
    "rad-sim.two-cell-radius-backlash-transition-rerun.v1"
)
TWO_CELL_EXTERNAL_FIDELITY_BENCHMARK_SCHEMA = "rad-sim.two-cell-external-fidelity-benchmark-summary.v1"
TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA = "rad-sim.two-cell-external-fidelity-correction-profile.v1"
TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA = (
    "rad-sim.two-cell-external-fidelity-correction-application.v1"
)


@dataclass(frozen=True)
class TwoCellMeasurement:
    case_id: str
    left_x: float | None = None
    left_y: float | None = None
    left_z: float | None = None
    right_x: float | None = None
    right_y: float | None = None
    right_z: float | None = None
    right_alpha: float | None = None
    right_theta: float | None = None
    alpha_command: float | None = None
    z_command: float | None = None
    pin_radius: float | None = None
    hole_radius: float | None = None
    contact_mode_observed: str = ""
    notes: str = ""

    @property
    def observed_values(self) -> dict[str, float]:
        values = {
            "leftX": self.left_x,
            "leftY": self.left_y,
            "leftZ": self.left_z,
            "rightX": self.right_x,
            "rightY": self.right_y,
            "rightZ": self.right_z,
            "rightAlpha": self.right_alpha,
            "rightTheta": self.right_theta,
        }
        return {key: value for key, value in values.items() if value is not None}


@dataclass(frozen=True)
class TwoCellConnectorMeasurement:
    case_id: str
    connector: str
    left_x_mm: float | None = None
    left_y_mm: float | None = None
    left_z_mm: float | None = None
    right_x_mm: float | None = None
    right_y_mm: float | None = None
    right_z_mm: float | None = None
    lateral_slip_mm: float | None = None
    vertical_slip_mm: float | None = None
    total_slip_mm: float | None = None
    contact_mode_observed: str = ""
    lock_held_observed: str = ""
    alpha_command: float | None = None
    z_command: float | None = None
    backlash: float | None = None
    pin_radius: float | None = None
    hole_radius: float | None = None
    measurement_source: str = ""
    notes: str = ""

    @property
    def observed_marker_values(self) -> dict[str, float]:
        values = {
            "leftXmm": self.left_x_mm,
            "leftYmm": self.left_y_mm,
            "leftZmm": self.left_z_mm,
            "rightXmm": self.right_x_mm,
            "rightYmm": self.right_y_mm,
            "rightZmm": self.right_z_mm,
        }
        return {key: value for key, value in values.items() if value is not None}

    @property
    def observed_slip_values(self) -> dict[str, float]:
        values = {
            "lateralSlipMm": self.lateral_slip_mm,
            "verticalSlipMm": self.vertical_slip_mm,
            "totalSlipMm": self.total_slip_mm,
        }
        return {key: value for key, value in values.items() if value is not None}


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    numeric = float(value)
    if not math.isfinite(numeric):
        return None
    return numeric


def _get(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    return None


def two_cell_measurements_from_rows(rows: Iterable[dict[str, Any]]) -> list[TwoCellMeasurement]:
    measurements: list[TwoCellMeasurement] = []
    for row in rows:
        case_id = str(_get(row, "caseId", "case_id") or "").strip()
        if not case_id:
            continue
        measurements.append(
            TwoCellMeasurement(
                case_id=case_id,
                left_x=_optional_float(_get(row, "leftX", "left_x")),
                left_y=_optional_float(_get(row, "leftY", "left_y")),
                left_z=_optional_float(_get(row, "leftZ", "left_z")),
                right_x=_optional_float(_get(row, "rightX", "right_x")),
                right_y=_optional_float(_get(row, "rightY", "right_y")),
                right_z=_optional_float(_get(row, "rightZ", "right_z")),
                right_alpha=_optional_float(_get(row, "rightAlpha", "right_alpha")),
                right_theta=_optional_float(_get(row, "rightTheta", "right_theta")),
                alpha_command=_optional_float(_get(row, "alphaCommand", "alpha_command")),
                z_command=_optional_float(_get(row, "zCommand", "z_command")),
                pin_radius=_optional_float(_get(row, "pinRadius", "pin_radius")),
                hole_radius=_optional_float(_get(row, "holeRadius", "hole_radius")),
                contact_mode_observed=str(_get(row, "contactModeObserved", "contact_mode_observed") or ""),
                notes=str(_get(row, "notes") or ""),
            )
        )
    return measurements


def two_cell_measurements_from_csv(text: str) -> list[TwoCellMeasurement]:
    return two_cell_measurements_from_rows(csv.DictReader(io.StringIO(text)))


def two_cell_measurements_from_json(text: str) -> list[TwoCellMeasurement]:
    payload = json.loads(text)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get("measurements"), list):
        rows = payload["measurements"]
    elif isinstance(payload, dict) and payload.get("schema") == TWO_CELL_PACKET_SCHEMA:
        rows = [
            {
                "caseId": case["caseId"],
                "leftX": case["bench"]["cells"][0]["center"]["x"],
                "leftY": case["bench"]["cells"][0]["center"]["y"],
                "leftZ": case["bench"]["cells"][0]["center"]["z"],
                "rightX": case["bench"]["cells"][1]["center"]["x"],
                "rightY": case["bench"]["cells"][1]["center"]["y"],
                "rightZ": case["bench"]["cells"][1]["center"]["z"],
                "rightAlpha": case["bench"]["cells"][1]["alpha"],
                "rightTheta": case["bench"]["cells"][1]["theta"],
                "contactModeObserved": case["bench"]["connector"]["contactMode"],
            }
            for case in payload.get("cases", [])
        ]
    else:
        rows = []
    return two_cell_measurements_from_rows(rows)


def two_cell_connector_measurements_from_rows(
    rows: Iterable[dict[str, Any]],
) -> list[TwoCellConnectorMeasurement]:
    measurements: list[TwoCellConnectorMeasurement] = []
    for row in rows:
        case_id = str(_get(row, "caseId", "case_id") or "").strip()
        connector = str(_get(row, "connector", "connectorName", "connector_name") or "").strip()
        if not case_id or not connector:
            continue
        measurements.append(
            TwoCellConnectorMeasurement(
                case_id=case_id,
                connector=connector,
                left_x_mm=_optional_float(_get(row, "observedLeftXmm", "leftXmm", "left_x_mm")),
                left_y_mm=_optional_float(_get(row, "observedLeftYmm", "leftYmm", "left_y_mm")),
                left_z_mm=_optional_float(_get(row, "observedLeftZmm", "leftZmm", "left_z_mm")),
                right_x_mm=_optional_float(_get(row, "observedRightXmm", "rightXmm", "right_x_mm")),
                right_y_mm=_optional_float(_get(row, "observedRightYmm", "rightYmm", "right_y_mm")),
                right_z_mm=_optional_float(_get(row, "observedRightZmm", "rightZmm", "right_z_mm")),
                lateral_slip_mm=_optional_float(
                    _get(row, "observedLateralSlipMm", "lateralSlipMm", "lateral_slip_mm")
                ),
                vertical_slip_mm=_optional_float(
                    _get(row, "observedVerticalSlipMm", "verticalSlipMm", "vertical_slip_mm")
                ),
                total_slip_mm=_optional_float(_get(row, "observedTotalSlipMm", "totalSlipMm", "total_slip_mm")),
                contact_mode_observed=str(_get(row, "observedContactMode", "contactModeObserved") or ""),
                lock_held_observed=str(_get(row, "lockHeldObserved", "lock_held_observed") or ""),
                alpha_command=_optional_float(_get(row, "alphaCommand", "alpha_command")),
                z_command=_optional_float(_get(row, "zCommand", "z_command")),
                backlash=_optional_float(_get(row, "backlash")),
                pin_radius=_optional_float(_get(row, "pinRadius", "pin_radius")),
                hole_radius=_optional_float(_get(row, "holeRadius", "hole_radius")),
                measurement_source=str(_get(row, "measurementSource", "measurement_source") or ""),
                notes=str(_get(row, "notes") or ""),
            )
        )
    return measurements


def two_cell_connector_measurements_from_csv(text: str) -> list[TwoCellConnectorMeasurement]:
    return two_cell_connector_measurements_from_rows(csv.DictReader(io.StringIO(text)))


def two_cell_connector_measurements_from_json(text: str) -> list[TwoCellConnectorMeasurement]:
    payload = json.loads(text)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get("rows"), list):
        rows = payload["rows"]
    elif isinstance(payload, dict) and isinstance(payload.get("measurements"), list):
        rows = payload["measurements"]
    else:
        rows = []
    return two_cell_connector_measurements_from_rows(rows)


def default_two_cell_case_controls(base: TwoCellBenchControls | None = None) -> dict[str, TwoCellBenchControls]:
    base = base or TwoCellBenchControls()
    return {
        "free_contract_lift": base,
        "free_expand_pushdown": TwoCellBenchControls(
            alpha_command=0.35,
            z_command=-0.35,
            hole_sweep_max=base.hole_sweep_max,
            hole_sweep_steps=base.hole_sweep_steps,
            left_position_locked=base.left_position_locked,
        ),
        "right_state_locked": TwoCellBenchControls(
            alpha_command=base.alpha_command,
            z_command=base.z_command,
            hole_sweep_max=base.hole_sweep_max,
            hole_sweep_steps=base.hole_sweep_steps,
            left_position_locked=base.left_position_locked,
            right_locked=True,
        ),
        "right_position_locked": TwoCellBenchControls(
            alpha_command=base.alpha_command,
            z_command=base.z_command,
            hole_sweep_max=base.hole_sweep_max,
            hole_sweep_steps=base.hole_sweep_steps,
            left_position_locked=base.left_position_locked,
            right_position_locked=True,
        ),
        "left_free_neighbor_residual": TwoCellBenchControls(
            alpha_command=base.alpha_command,
            z_command=base.z_command,
            hole_sweep_max=base.hole_sweep_max,
            hole_sweep_steps=base.hole_sweep_steps,
            left_position_locked=False,
        ),
    }


def _controls_for_measurement(
    measurement: TwoCellMeasurement,
    controls_by_case: dict[str, TwoCellBenchControls],
) -> TwoCellBenchControls:
    controls = controls_by_case.get(measurement.case_id, TwoCellBenchControls())
    return TwoCellBenchControls(
        alpha_command=controls.alpha_command if measurement.alpha_command is None else measurement.alpha_command,
        z_command=controls.z_command if measurement.z_command is None else measurement.z_command,
        hole_sweep_max=controls.hole_sweep_max,
        hole_sweep_steps=controls.hole_sweep_steps,
        left_position_locked=controls.left_position_locked,
        right_locked=controls.right_locked,
        right_position_locked=controls.right_position_locked,
    )


def _predicted_values(bench: dict[str, Any]) -> dict[str, float]:
    left = bench["cells"][0]
    right = bench["cells"][1]
    return {
        "leftX": left["center"]["x"],
        "leftY": left["center"]["y"],
        "leftZ": left["center"]["z"],
        "rightX": right["center"]["x"],
        "rightY": right["center"]["y"],
        "rightZ": right["center"]["z"],
        "rightAlpha": right["alpha"],
        "rightTheta": right["theta"],
    }


def compare_two_cell_measurements(
    measurements: Iterable[TwoCellMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls_by_case = default_two_cell_case_controls(base_controls)
    rows = []
    squared = 0.0
    count = 0
    max_abs = 0.0
    contact_matches = 0
    contact_count = 0
    missing_cases = []
    for measurement in measurements:
        controls = _controls_for_measurement(measurement, controls_by_case)
        hole_radius = measurement.hole_radius if measurement.hole_radius is not None else config.hole_radius
        pin_radius = measurement.pin_radius if measurement.pin_radius is not None else config.pin_radius
        bench = simulate_two_cell_bench(
            config,
            controls,
            hole_radius=hole_radius,
            pin_radius=pin_radius,
        )
        predicted = _predicted_values(bench)
        residuals = {}
        for key, observed in measurement.observed_values.items():
            residual = observed - predicted[key]
            residuals[key] = residual
            squared += residual * residual
            count += 1
            max_abs = max(max_abs, abs(residual))
        observed_contact = measurement.contact_mode_observed.strip()
        if observed_contact:
            contact_count += 1
            if observed_contact == bench["connector"]["contactMode"]:
                contact_matches += 1
        if not measurement.observed_values:
            missing_cases.append(measurement.case_id)
        rows.append(
            {
                "caseId": measurement.case_id,
                "predicted": predicted,
                "observed": measurement.observed_values,
                "residual": residuals,
                "contactModePredicted": bench["connector"]["contactMode"],
                "contactModeObserved": observed_contact,
                "availableFields": sorted(measurement.observed_values),
            }
        )
    rms = math.sqrt(squared / count) if count else None
    return {
        "schema": TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA,
        "model": "cad-derived-two-cell-rigid-contact-proxy",
        "sampleCount": len(rows),
        "observedScalarCount": count,
        "rmsError": rms,
        "maxAbsError": max_abs if count else None,
        "contactModeAccuracy": contact_matches / contact_count if contact_count else None,
        "missingObservationCases": missing_cases,
        "rows": rows,
        "acceptance": {
            "readyForParameterFit": count >= 8,
            "requiresMoreData": count < 8,
            "requiresSegmentedPhysics": True,
        },
    }


_TRUE_LOCK_VALUES = {"1", "true", "t", "yes", "y", "held", "pass", "passed", "locked", "lock-held"}
_FALSE_LOCK_VALUES = {
    "0",
    "false",
    "f",
    "no",
    "n",
    "fell",
    "fail",
    "failed",
    "released",
    "not-held",
    "unheld",
    "lock-fell",
}


def _optional_bool(value: str) -> bool | None:
    normalized = value.strip().lower().replace(" ", "-").replace("_", "-")
    if not normalized:
        return None
    if normalized in _TRUE_LOCK_VALUES:
        return True
    if normalized in _FALSE_LOCK_VALUES:
        return False
    return None


def _connector_prediction_map(
    config: LatticeConfig,
    base_controls: TwoCellBenchControls | None,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, Any]]:
    template = two_cell_connector_measurement_template(config, base_controls)
    prediction_map = {
        (str(row["caseId"]), str(row["connector"])): row
        for row in template["rows"]
    }
    return prediction_map, template


def _implied_slip(measurement: TwoCellConnectorMeasurement) -> dict[str, float]:
    marker = measurement.observed_marker_values
    implied: dict[str, float] = {}
    if "leftXmm" in marker and "rightXmm" in marker:
        implied["slipXmm"] = marker["rightXmm"] - marker["leftXmm"]
    if "leftYmm" in marker and "rightYmm" in marker:
        implied["slipYmm"] = marker["rightYmm"] - marker["leftYmm"]
    if "leftZmm" in marker and "rightZmm" in marker:
        implied["slipZmm"] = marker["rightZmm"] - marker["leftZmm"]
        implied["verticalSlipMm"] = abs(implied["slipZmm"])
    if "slipXmm" in implied and "slipYmm" in implied:
        implied["lateralSlipMm"] = math.hypot(implied["slipXmm"], implied["slipYmm"])
    if "lateralSlipMm" in implied and "verticalSlipMm" in implied:
        implied["totalSlipMm"] = math.hypot(implied["lateralSlipMm"], implied["verticalSlipMm"])
    return implied


def _lock_expected(row: dict[str, Any]) -> bool:
    lock_mode = str(row.get("lockMode", "")).strip().lower()
    return lock_mode in {"right_state_locked", "right_position_locked"}


def compare_two_cell_connector_measurements(
    measurements: Iterable[TwoCellConnectorMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    measurements = list(measurements)
    prediction_map, template = _connector_prediction_map(config, base_controls)
    rows = []
    marker_squared = 0.0
    marker_count = 0
    marker_max_abs = 0.0
    slip_squared = 0.0
    slip_count = 0
    slip_max_abs = 0.0
    contact_count = 0
    contact_matches = 0
    lock_count = 0
    lock_matches = 0
    lock_expected_count = 0
    missing_predictions = []
    missing_observation_rows = []
    missing_marker_rows = []
    missing_slip_rows = []
    contact_mismatch_rows = []
    lock_mismatch_rows = []
    tolerance_failure_rows = []

    marker_prediction_keys = {
        "leftXmm": "predictedLeftXmm",
        "leftYmm": "predictedLeftYmm",
        "leftZmm": "predictedLeftZmm",
        "rightXmm": "predictedRightXmm",
        "rightYmm": "predictedRightYmm",
        "rightZmm": "predictedRightZmm",
    }
    slip_prediction_keys = {
        "lateralSlipMm": "predictedLateralSlipMm",
        "verticalSlipMm": "predictedVerticalSlipMm",
        "totalSlipMm": "predictedTotalSlipMm",
    }

    for measurement in measurements:
        row_key = (measurement.case_id, measurement.connector)
        predicted = prediction_map.get(row_key)
        row_label = f"{measurement.case_id}:{measurement.connector}"
        if predicted is None:
            missing_predictions.append(row_label)
            rows.append(
                {
                    "caseId": measurement.case_id,
                    "connector": measurement.connector,
                    "status": "missing-prediction",
                    "observed": {
                        **measurement.observed_marker_values,
                        **measurement.observed_slip_values,
                    },
                }
            )
            continue

        marker_residuals = {}
        marker_row_squared = 0.0
        marker_row_count = 0
        marker_row_max_abs = 0.0
        for observed_key, predicted_key in marker_prediction_keys.items():
            observed = measurement.observed_marker_values.get(observed_key)
            if observed is None:
                continue
            residual = observed - float(predicted[predicted_key])
            marker_residuals[observed_key] = residual
            marker_squared += residual * residual
            marker_count += 1
            marker_max_abs = max(marker_max_abs, abs(residual))
            marker_row_squared += residual * residual
            marker_row_count += 1
            marker_row_max_abs = max(marker_row_max_abs, abs(residual))

        implied_slip = _implied_slip(measurement)
        observed_slip = dict(implied_slip)
        observed_slip.update(measurement.observed_slip_values)
        slip_residuals = {}
        slip_row_squared = 0.0
        slip_row_count = 0
        slip_row_max_abs = 0.0
        for observed_key, predicted_key in slip_prediction_keys.items():
            observed = observed_slip.get(observed_key)
            if observed is None:
                continue
            residual = observed - float(predicted[predicted_key])
            slip_residuals[observed_key] = residual
            slip_squared += residual * residual
            slip_count += 1
            slip_max_abs = max(slip_max_abs, abs(residual))
            slip_row_squared += residual * residual
            slip_row_count += 1
            slip_row_max_abs = max(slip_row_max_abs, abs(residual))

        observed_contact = measurement.contact_mode_observed.strip()
        contact_match = None
        if observed_contact:
            contact_count += 1
            contact_match = observed_contact.lower() == str(predicted["predictedContactMode"]).lower()
            if contact_match:
                contact_matches += 1
            else:
                contact_mismatch_rows.append(row_label)

        expected_lock = _lock_expected(predicted)
        if expected_lock:
            lock_expected_count += 1
        observed_lock = _optional_bool(measurement.lock_held_observed)
        lock_match = None
        if expected_lock and observed_lock is not None:
            lock_count += 1
            lock_match = observed_lock is True
            if lock_match:
                lock_matches += 1
            else:
                lock_mismatch_rows.append(row_label)

        available_fields = sorted(
            [
                *measurement.observed_marker_values,
                *measurement.observed_slip_values,
                *(["observedContactMode"] if observed_contact else []),
                *(["lockHeldObserved"] if measurement.lock_held_observed.strip() else []),
            ]
        )
        if not available_fields:
            missing_observation_rows.append(row_label)
        if not measurement.observed_marker_values:
            missing_marker_rows.append(row_label)
        if not observed_slip:
            missing_slip_rows.append(row_label)

        row_max_abs = max(marker_row_max_abs, slip_row_max_abs)
        if row_max_abs > tolerance:
            tolerance_failure_rows.append(row_label)

        rows.append(
            {
                "caseId": measurement.case_id,
                "connector": measurement.connector,
                "status": "matched",
                "lockMode": predicted["lockMode"],
                "predictedContactMode": predicted["predictedContactMode"],
                "observedContactMode": observed_contact,
                "contactModeMatches": contact_match,
                "lockExpected": expected_lock,
                "lockHeldObserved": observed_lock,
                "lockHeldMatches": lock_match,
                "predictedSlipMm": {
                    "lateral": predicted["predictedLateralSlipMm"],
                    "vertical": predicted["predictedVerticalSlipMm"],
                    "total": predicted["predictedTotalSlipMm"],
                },
                "observedSlipMm": {
                    key: value
                    for key, value in observed_slip.items()
                    if key in {"lateralSlipMm", "verticalSlipMm", "totalSlipMm"}
                },
                "impliedSlipMm": implied_slip,
                "markerResidualMm": marker_residuals,
                "slipResidualMm": slip_residuals,
                "markerRmsErrorMm": math.sqrt(marker_row_squared / marker_row_count)
                if marker_row_count
                else None,
                "slipRmsErrorMm": math.sqrt(slip_row_squared / slip_row_count) if slip_row_count else None,
                "maxAbsResidualMm": row_max_abs if marker_row_count or slip_row_count else None,
                "availableFields": available_fields,
                "measurementSource": measurement.measurement_source,
                "notes": measurement.notes,
            }
        )

    matched_count = len([row for row in rows if row.get("status") == "matched"])
    observed_scalar_count = marker_count + slip_count
    contact_accuracy = contact_matches / contact_count if contact_count else None
    lock_accuracy = lock_matches / lock_count if lock_count else None
    passes_numeric_tolerance = observed_scalar_count > 0 and not tolerance_failure_rows
    passes_contact = contact_accuracy in (None, 1.0)
    passes_locks = lock_accuracy in (None, 1.0)
    passes_tolerance = (
        matched_count > 0
        and not missing_predictions
        and passes_numeric_tolerance
        and passes_contact
        and passes_locks
    )
    missing_evidence = []
    if not measurements:
        missing_evidence.append("filledConnectorMeasurementRows")
    if missing_observation_rows:
        missing_evidence.append("observedConnectorRows")
    if marker_count == 0:
        missing_evidence.append("observedConnectorMarkerPositions")
    if slip_count == 0:
        missing_evidence.append("observedConnectorSlip")
    if contact_count == 0:
        missing_evidence.append("observedContactModes")
    if lock_expected_count and lock_count == 0:
        missing_evidence.append("observedLockHeldFlags")
    if missing_predictions:
        missing_evidence.append("matchedCaseConnectorRows")
    if observed_scalar_count and tolerance_failure_rows:
        missing_evidence.append("connectorToleranceAgreement")
    if contact_mismatch_rows:
        missing_evidence.append("contactModeAgreement")
    if lock_mismatch_rows:
        missing_evidence.append("lockHeldAgreement")

    return {
        "schema": TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA,
        "model": "a360-cad-derived-two-cell-connector-proxy",
        "cadReference": template["cadReference"],
        "toleranceMm": tolerance,
        "sampleCount": len(measurements),
        "templateRowCount": len(template["rows"]),
        "matchedRowCount": matched_count,
        "observedConnectorRowCount": len(measurements) - len(missing_observation_rows),
        "observedScalarCount": observed_scalar_count,
        "markerResidualCount": marker_count,
        "slipResidualCount": slip_count,
        "markerRmsErrorMm": math.sqrt(marker_squared / marker_count) if marker_count else None,
        "markerMaxAbsErrorMm": marker_max_abs if marker_count else None,
        "slipRmsErrorMm": math.sqrt(slip_squared / slip_count) if slip_count else None,
        "slipMaxAbsErrorMm": slip_max_abs if slip_count else None,
        "contactModeAccuracy": contact_accuracy,
        "lockHeldAccuracy": lock_accuracy,
        "missingPredictionRows": missing_predictions,
        "missingObservationRows": missing_observation_rows,
        "missingMarkerRows": missing_marker_rows,
        "missingSlipRows": missing_slip_rows,
        "contactMismatchRows": contact_mismatch_rows,
        "lockMismatchRows": lock_mismatch_rows,
        "toleranceFailureRows": tolerance_failure_rows,
        "missingEvidence": sorted(set(missing_evidence)),
        "rows": rows,
        "acceptance": {
            "readyForConnectorCalibration": observed_scalar_count >= 6 and matched_count > 0,
            "coversFullConnectorSuite": matched_count == len(template["rows"]) and not missing_observation_rows,
            "passesTolerance": passes_tolerance,
            "requiresMoreData": bool(missing_evidence),
            "requiresSegmentedPhysics": True,
            "physicalAccuracyValidated": False,
        },
    }


def two_cell_fidelity_matrix_measurements_from_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    measurements: list[dict[str, Any]] = []
    for row in rows:
        case_id = str(_get(row, "caseId", "case_id") or "").strip()
        connector = str(_get(row, "connector", "connectorName", "connector_name") or "").strip()
        if not case_id or not connector:
            continue
        measurements.append(dict(row))
    return measurements


def two_cell_fidelity_matrix_measurements_from_csv(text: str) -> list[dict[str, Any]]:
    return two_cell_fidelity_matrix_measurements_from_rows(csv.DictReader(io.StringIO(text)))


def two_cell_fidelity_matrix_measurements_from_json(text: str) -> list[dict[str, Any]]:
    payload = json.loads(text)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get("rows"), list):
        rows = payload["rows"]
    elif isinstance(payload, dict) and isinstance(payload.get("measurements"), list):
        rows = payload["measurements"]
    else:
        rows = []
    return two_cell_fidelity_matrix_measurements_from_rows(rows)


def _fidelity_prediction_map(
    config: LatticeConfig,
    base_controls: TwoCellBenchControls | None,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, Any]]:
    template = two_cell_fidelity_matrix_measurement_template(config, base_controls)
    prediction_map = {
        (str(row["caseId"]), str(row["connector"])): row
        for row in template["rows"]
    }
    return prediction_map, template


def compare_two_cell_fidelity_matrix_measurements(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    measurements = two_cell_fidelity_matrix_measurements_from_rows(measurements)
    prediction_map, template = _fidelity_prediction_map(config, base_controls)
    numeric_pairs = {
        "observedLeftCellX": "predictedLeftCellX",
        "observedLeftCellY": "predictedLeftCellY",
        "observedLeftCellZ": "predictedLeftCellZ",
        "observedRightCellX": "predictedRightCellX",
        "observedRightCellY": "predictedRightCellY",
        "observedRightCellZ": "predictedRightCellZ",
        "observedLeftAlpha": "predictedLeftAlpha",
        "observedRightAlpha": "predictedRightAlpha",
        "observedLeftTheta": "predictedLeftTheta",
        "observedRightTheta": "predictedRightTheta",
        "observedLeftXmm": "predictedLeftXmm",
        "observedLeftYmm": "predictedLeftYmm",
        "observedLeftZmm": "predictedLeftZmm",
        "observedRightXmm": "predictedRightXmm",
        "observedRightYmm": "predictedRightYmm",
        "observedRightZmm": "predictedRightZmm",
        "observedLateralSlipMm": "predictedLateralSlipMm",
        "observedVerticalSlipMm": "predictedVerticalSlipMm",
        "observedTotalSlipMm": "predictedTotalSlipMm",
    }
    cell_fields = {key for key in numeric_pairs if "Cell" in key or key.endswith("Alpha") or key.endswith("Theta")}
    connector_fields = set(numeric_pairs) - cell_fields
    rows = []
    squared = 0.0
    count = 0
    max_abs = 0.0
    cell_squared = 0.0
    cell_count = 0
    connector_squared = 0.0
    connector_count = 0
    contact_count = 0
    contact_matches = 0
    lock_count = 0
    lock_matches = 0
    missing_predictions = []
    missing_observation_rows = []
    tolerance_failures = []
    contact_mismatches = []
    lock_mismatches = []
    for measurement in measurements:
        case_id = str(_get(measurement, "caseId", "case_id") or "").strip()
        connector = str(_get(measurement, "connector", "connectorName", "connector_name") or "").strip()
        label = f"{case_id}:{connector}"
        predicted = prediction_map.get((case_id, connector))
        if predicted is None:
            missing_predictions.append(label)
            rows.append({"caseId": case_id, "connector": connector, "status": "missing-prediction"})
            continue

        residuals = {}
        row_max_abs = 0.0
        row_count = 0
        for observed_key, predicted_key in numeric_pairs.items():
            observed = _optional_float(_get(measurement, observed_key))
            if observed is None:
                continue
            residual = observed - float(predicted[predicted_key])
            residuals[observed_key] = residual
            squared += residual * residual
            count += 1
            row_count += 1
            max_abs = max(max_abs, abs(residual))
            row_max_abs = max(row_max_abs, abs(residual))
            if observed_key in cell_fields:
                cell_squared += residual * residual
                cell_count += 1
            if observed_key in connector_fields:
                connector_squared += residual * residual
                connector_count += 1
        observed_contact = str(_get(measurement, "observedContactMode", "contactModeObserved") or "").strip()
        contact_match = None
        if observed_contact:
            contact_count += 1
            contact_match = observed_contact.lower() == str(predicted["predictedContactMode"]).lower()
            if contact_match:
                contact_matches += 1
            else:
                contact_mismatches.append(label)
        expected_lock = _lock_expected(predicted)
        observed_lock = _optional_bool(str(_get(measurement, "lockHeldObserved", "lock_held_observed") or ""))
        lock_match = None
        if expected_lock and observed_lock is not None:
            lock_count += 1
            lock_match = observed_lock is True
            if lock_match:
                lock_matches += 1
            else:
                lock_mismatches.append(label)
        if not residuals and not observed_contact and observed_lock is None:
            missing_observation_rows.append(label)
        if row_max_abs > tolerance:
            tolerance_failures.append(label)
        rows.append(
            {
                "caseId": case_id,
                "matrixIndex": predicted["matrixIndex"],
                "connector": connector,
                "status": "matched",
                "actuationCase": predicted["actuationCase"],
                "lockMode": predicted["lockMode"],
                "predictedContactMode": predicted["predictedContactMode"],
                "observedContactMode": observed_contact,
                "contactModeMatches": contact_match,
                "lockExpected": expected_lock,
                "lockHeldObserved": observed_lock,
                "lockHeldMatches": lock_match,
                "residuals": residuals,
                "maxAbsResidual": row_max_abs if row_count else None,
                "observedScalarCount": row_count,
                "measurementSource": str(_get(measurement, "measurementSource", "measurement_source") or ""),
                "notes": str(_get(measurement, "notes") or ""),
            }
        )
    matched_count = sum(1 for row in rows if row.get("status") == "matched")
    contact_accuracy = contact_matches / contact_count if contact_count else None
    lock_accuracy = lock_matches / lock_count if lock_count else None
    missing_evidence = []
    if not measurements:
        missing_evidence.append("filledFidelityMatrixMeasurementRows")
    if missing_observation_rows:
        missing_evidence.append("observedMatrixRows")
    if count == 0:
        missing_evidence.append("observedNumericFields")
    if cell_count == 0:
        missing_evidence.append("observedCellCoordinates")
    if connector_count == 0:
        missing_evidence.append("observedConnectorCoordinatesOrSlip")
    if contact_count == 0:
        missing_evidence.append("observedContactModes")
    if lock_count == 0:
        missing_evidence.append("observedLockHeldFlags")
    if missing_predictions:
        missing_evidence.append("matchedMatrixCaseConnectorRows")
    if tolerance_failures:
        missing_evidence.append("fidelityMatrixToleranceAgreement")
    if contact_mismatches:
        missing_evidence.append("contactModeAgreement")
    if lock_mismatches:
        missing_evidence.append("lockHeldAgreement")
    passes_tolerance = (
        matched_count > 0
        and count > 0
        and not missing_predictions
        and not tolerance_failures
        and not contact_mismatches
        and not lock_mismatches
    )
    return {
        "schema": TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA,
        "model": "two-cell-fidelity-matrix-measured-vs-reduced-proxy",
        "cadReference": template["cadReference"],
        "tolerance": tolerance,
        "sampleCount": len(measurements),
        "templateRowCount": len(template["rows"]),
        "matchedRowCount": matched_count,
        "observedScalarCount": count,
        "observedCellScalarCount": cell_count,
        "observedConnectorScalarCount": connector_count,
        "rmsError": math.sqrt(squared / count) if count else None,
        "maxAbsError": max_abs if count else None,
        "cellRmsError": math.sqrt(cell_squared / cell_count) if cell_count else None,
        "connectorRmsError": math.sqrt(connector_squared / connector_count) if connector_count else None,
        "contactModeAccuracy": contact_accuracy,
        "lockHeldAccuracy": lock_accuracy,
        "missingPredictionRows": missing_predictions,
        "missingObservationRows": missing_observation_rows,
        "toleranceFailureRows": tolerance_failures,
        "contactMismatchRows": contact_mismatches,
        "lockMismatchRows": lock_mismatches,
        "missingEvidence": sorted(set(missing_evidence)),
        "rows": rows,
        "acceptance": {
            "readyForMatrixCalibration": count >= 24 and matched_count > 0,
            "coversFullFidelityMatrix": matched_count == len(template["rows"]) and not missing_observation_rows,
            "passesTolerance": passes_tolerance,
            "requiresMoreData": bool(missing_evidence),
            "requiresSegmentedPhysics": True,
            "physicalAccuracyValidated": False,
        },
    }


def _fit_scale(samples: list[tuple[float, float]]) -> dict[str, Any]:
    usable = [(predicted, observed) for predicted, observed in samples if abs(predicted) > 1e-12]
    if not usable:
        return {
            "estimate": None,
            "sampleCount": 0,
            "status": "insufficient-excitation",
            "meanResidual": None,
            "rmsResidual": None,
        }
    numerator = sum(predicted * observed for predicted, observed in usable)
    denominator = sum(predicted * predicted for predicted, _ in usable)
    estimate = numerator / denominator if denominator else None
    residuals = [observed - predicted for predicted, observed in usable]
    squared = sum(value * value for value in residuals)
    return {
        "estimate": estimate,
        "sampleCount": len(usable),
        "status": "estimated" if estimate is not None else "insufficient-excitation",
        "meanResidual": sum(residuals) / len(residuals),
        "rmsResidual": math.sqrt(squared / len(residuals)),
    }


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _bounded(value: float | None, lower: float, upper: float, fallback: float) -> float:
    if value is None or not math.isfinite(value):
        return fallback
    return max(lower, min(upper, value))


def _axis_residual_groups(
    measurements: list[dict[str, Any]],
    prediction_map: dict[tuple[str, str], dict[str, Any]],
    *,
    axis: str,
) -> list[dict[str, Any]]:
    buckets: dict[float, dict[str, Any]] = {}
    for measurement in measurements:
        case_id = str(_get(measurement, "caseId", "case_id") or "").strip()
        connector = str(_get(measurement, "connector", "connectorName", "connector_name") or "").strip()
        predicted = prediction_map.get((case_id, connector))
        if predicted is None:
            continue
        axis_value = _optional_float(_get(measurement, axis))
        if axis_value is None:
            axis_value = _optional_float(predicted.get(axis))
        if axis_value is None:
            continue
        bucket = buckets.setdefault(
            float(axis_value),
            {
                axis: float(axis_value),
                "count": 0,
                "rightZResiduals": [],
                "rightAlphaResiduals": [],
                "verticalSlipResiduals": [],
                "lateralSlipResiduals": [],
            },
        )
        for observed_key, predicted_key, residual_key in (
            ("observedRightCellZ", "predictedRightCellZ", "rightZResiduals"),
            ("observedRightAlpha", "predictedRightAlpha", "rightAlphaResiduals"),
            ("observedVerticalSlipMm", "predictedVerticalSlipMm", "verticalSlipResiduals"),
            ("observedLateralSlipMm", "predictedLateralSlipMm", "lateralSlipResiduals"),
        ):
            observed = _optional_float(_get(measurement, observed_key))
            if observed is None:
                continue
            bucket[residual_key].append(observed - float(predicted[predicted_key]))
            bucket["count"] += 1
    grouped = []
    for value in sorted(buckets):
        bucket = buckets[value]
        grouped.append(
            {
                axis: bucket[axis],
                "observedScalarCount": bucket["count"],
                "meanRightZResidual": _mean(bucket["rightZResiduals"]),
                "meanRightAlphaResidual": _mean(bucket["rightAlphaResiduals"]),
                "meanVerticalSlipResidualMm": _mean(bucket["verticalSlipResiduals"]),
                "meanLateralSlipResidualMm": _mean(bucket["lateralSlipResiduals"]),
            }
        )
    return grouped


def calibrate_two_cell_fidelity_matrix_parameters(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Estimate reduced-proxy calibration factors from the dense two-cell matrix."""

    config = config or LatticeConfig(rows=1, cols=2)
    measurements = two_cell_fidelity_matrix_measurements_from_rows(measurements)
    prediction_map, template = _fidelity_prediction_map(config, base_controls)
    comparison = compare_two_cell_fidelity_matrix_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
    )
    alpha_samples: list[tuple[float, float]] = []
    theta_samples: list[tuple[float, float]] = []
    z_samples: list[tuple[float, float]] = []
    lateral_slip_samples: list[tuple[float, float]] = []
    vertical_slip_samples: list[tuple[float, float]] = []
    total_slip_samples: list[tuple[float, float]] = []
    alpha_biases: list[float] = []
    z_biases: list[float] = []
    vertical_slip_biases: list[float] = []
    for measurement in measurements:
        case_id = str(_get(measurement, "caseId", "case_id") or "").strip()
        connector = str(_get(measurement, "connector", "connectorName", "connector_name") or "").strip()
        predicted = prediction_map.get((case_id, connector))
        if predicted is None:
            continue
        for observed_key, predicted_key in (
            ("observedLeftAlpha", "predictedLeftAlpha"),
            ("observedRightAlpha", "predictedRightAlpha"),
        ):
            observed = _optional_float(_get(measurement, observed_key))
            if observed is not None:
                predicted_delta = float(predicted[predicted_key]) - config.initial_alpha
                observed_delta = observed - config.initial_alpha
                alpha_samples.append((predicted_delta, observed_delta))
                alpha_biases.append(observed - float(predicted[predicted_key]))
        for observed_key, predicted_key in (
            ("observedLeftTheta", "predictedLeftTheta"),
            ("observedRightTheta", "predictedRightTheta"),
        ):
            observed = _optional_float(_get(measurement, observed_key))
            if observed is not None:
                theta_samples.append((float(predicted[predicted_key]), observed))
        for observed_key, predicted_key in (
            ("observedLeftCellZ", "predictedLeftCellZ"),
            ("observedRightCellZ", "predictedRightCellZ"),
        ):
            observed = _optional_float(_get(measurement, observed_key))
            if observed is not None:
                z_samples.append((float(predicted[predicted_key]), observed))
                z_biases.append(observed - float(predicted[predicted_key]))
        for observed_key, predicted_key, samples in (
            ("observedLateralSlipMm", "predictedLateralSlipMm", lateral_slip_samples),
            ("observedVerticalSlipMm", "predictedVerticalSlipMm", vertical_slip_samples),
            ("observedTotalSlipMm", "predictedTotalSlipMm", total_slip_samples),
        ):
            observed = _optional_float(_get(measurement, observed_key))
            if observed is not None:
                predicted_value = float(predicted[predicted_key])
                samples.append((predicted_value, observed))
                if observed_key == "observedVerticalSlipMm":
                    vertical_slip_biases.append(observed - predicted_value)
    alpha_fit = _fit_scale(alpha_samples)
    theta_fit = _fit_scale(theta_samples)
    z_fit = _fit_scale(z_samples)
    lateral_slip_fit = _fit_scale(lateral_slip_samples)
    vertical_slip_fit = _fit_scale(vertical_slip_samples)
    total_slip_fit = _fit_scale(total_slip_samples)
    alpha_scale = _bounded(alpha_fit["estimate"], 0.25, 4.0, 1.0)
    z_scale = _bounded(z_fit["estimate"], 0.25, 4.0, 1.0)
    vertical_slip_scale = _bounded(vertical_slip_fit["estimate"], 0.25, 4.0, 1.0)
    current_clearance = max(0.0, config.pin_hole_clearance)
    proposed_clearance = max(0.0, current_clearance * vertical_slip_scale)
    proposed = {
        "couplingGain": _bounded(config.coupling_gain * alpha_scale, 0.0, 1.0, config.coupling_gain),
        "zCouplingGain": _bounded(config.z_coupling_gain * z_scale, 0.0, 1.0, config.z_coupling_gain),
        "holeRadius": max(config.pin_radius, config.pin_radius + proposed_clearance),
        "pinHoleClearance": proposed_clearance,
        "applyAutomatically": False,
        "reason": "reduced-proxy estimate only; require holdout bench and segmented contact validation before applying to physical claims",
    }
    ready = bool(comparison["acceptance"]["readyForMatrixCalibration"])
    return {
        "schema": TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA,
        "model": "two-cell-fidelity-matrix-reduced-proxy-calibration-estimate",
        "cadReference": template["cadReference"],
        "comparisonSummary": {
            "sampleCount": comparison["sampleCount"],
            "matchedRowCount": comparison["matchedRowCount"],
            "observedScalarCount": comparison["observedScalarCount"],
            "rmsError": comparison["rmsError"],
            "maxAbsError": comparison["maxAbsError"],
            "passesTolerance": comparison["acceptance"]["passesTolerance"],
            "missingEvidence": comparison["missingEvidence"],
        },
        "estimates": {
            "alphaResponseScale": alpha_fit,
            "thetaScale": theta_fit,
            "zResponseScale": z_fit,
            "lateralSlipScale": lateral_slip_fit,
            "verticalSlipScale": vertical_slip_fit,
            "totalSlipScale": total_slip_fit,
            "meanAlphaBias": _mean(alpha_biases),
            "meanZBias": _mean(z_biases),
            "meanVerticalSlipBiasMm": _mean(vertical_slip_biases),
        },
        "axisDiagnostics": {
            "byHoleRadius": _axis_residual_groups(measurements, prediction_map, axis="holeRadius"),
            "byBacklash": _axis_residual_groups(measurements, prediction_map, axis="backlash"),
        },
        "proposedReducedProxyUpdates": proposed,
        "acceptance": {
            "readyForReducedProxyCalibration": ready,
            "requiresMoreData": not ready,
            "requiresHoldoutValidation": True,
            "requiresSegmentedPhysics": True,
            "physicalAccuracyValidated": False,
        },
        "claimLabels": {
            "calibration": "fits scale/bias terms for the reduced two-cell proxy only",
            "holeRadius": "vertical slip scale suggests an effective clearance update, not an exact CAD hole dimension",
            "physicalAccuracy": "not validated until segmented CAD/contact and independent bench measurements agree",
        },
    }


def two_cell_radius_backlash_transition_measurements_from_rows(
    rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    measurements: list[dict[str, Any]] = []
    for row in rows:
        case_id = str(_get(row, "caseId", "case_id") or "").strip()
        phase_case_id = str(_get(row, "phaseDiagramCaseId", "phase_diagram_case_id") or "").strip()
        if not case_id and not phase_case_id:
            continue
        measurements.append(dict(row))
    return measurements


def two_cell_radius_backlash_transition_measurements_from_csv(text: str) -> list[dict[str, Any]]:
    return two_cell_radius_backlash_transition_measurements_from_rows(csv.DictReader(io.StringIO(text)))


def two_cell_radius_backlash_transition_measurements_from_json(text: str) -> list[dict[str, Any]]:
    payload = json.loads(text)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get("rows"), list):
        rows = payload["rows"]
    elif isinstance(payload, dict) and isinstance(payload.get("measurements"), list):
        rows = payload["measurements"]
    elif isinstance(payload, dict) and isinstance(payload.get("measurementCases"), list):
        rows = payload["measurementCases"]
    else:
        rows = []
    return two_cell_radius_backlash_transition_measurements_from_rows(rows)


def _transition_prediction_maps(
    config: LatticeConfig,
    base_controls: TwoCellBenchControls | None,
    report_overrides: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, Any]]:
    report = two_cell_radius_backlash_transition_report(config, base_controls, **report_overrides)
    by_case_id = {str(row["caseId"]): row for row in report["measurementCases"]}
    by_phase_case_id = {str(row["phaseDiagramCaseId"]): row for row in report["measurementCases"]}
    brackets_by_phase_case_id: dict[str, list[dict[str, Any]]] = {}
    for bracket in report["transitionBrackets"]:
        for side_key, side in (("leftCaseId", "lower"), ("rightCaseId", "upper")):
            case_id = str(bracket.get(side_key, ""))
            if not case_id:
                continue
            brackets_by_phase_case_id.setdefault(case_id, []).append({**bracket, "transitionSide": side})
    return by_case_id, by_phase_case_id, brackets_by_phase_case_id, report


def _phase_label(value: Any) -> str:
    normalized = str(value or "").strip().lower()
    normalized = normalized.replace("+", " and ").replace("_", " ").replace("-", " ")
    normalized = " ".join(normalized.split())
    aliases = {
        "axial vertical contact": "axial-and-vertical-contact",
        "axial and vertical contact": "axial-and-vertical-contact",
        "axial contact": "axial-contact",
        "vertical contact": "vertical-contact",
        "position locked": "position-locked",
        "state locked": "state-locked",
        "state locked vertical free": "state-locked-vertical-free",
        "free play": "free-play",
        "gravity sag": "gravity-sag",
    }
    return aliases.get(normalized, normalized.replace(" ", "-"))


def _transition_shift_vote(
    predicted: dict[str, Any],
    observed_phase: str,
    brackets: list[dict[str, Any]],
) -> str:
    if not observed_phase or observed_phase == str(predicted.get("phase", "")):
        return "within-predicted-bracket"
    votes: list[str] = []
    for bracket in brackets:
        side = str(bracket.get("transitionSide", ""))
        if side == "lower" and observed_phase == str(bracket.get("toPhase", "")):
            votes.append(f"{bracket['axis']}-boundary-lower-than-predicted")
        elif side == "upper" and observed_phase == str(bracket.get("fromPhase", "")):
            votes.append(f"{bracket['axis']}-boundary-higher-than-predicted")
    return votes[0] if votes else "phase-mismatch-not-adjacent"


def compare_two_cell_radius_backlash_transition_measurements(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    measurements = two_cell_radius_backlash_transition_measurements_from_rows(measurements)
    by_case_id, by_phase_case_id, brackets_by_phase_case_id, report = _transition_prediction_maps(
        config,
        base_controls,
        dict(report_overrides),
    )
    numeric_pairs = {
        "observedRightZ": "rightZ",
        "observedRightAlpha": "rightAlpha",
        "observedConnectorMaxVerticalSlipMm": "connectorMaxVerticalSlipMm",
        "observedConnectorMaxVerticalExcessMm": "connectorMaxVerticalExcessMm",
        "observedConnectorTotalContactPenalty": "connectorTotalContactPenalty",
    }
    rows: list[dict[str, Any]] = []
    squared = 0.0
    count = 0
    max_abs = 0.0
    phase_count = 0
    phase_matches = 0
    lock_count = 0
    lock_matches = 0
    missing_predictions: list[str] = []
    missing_observation_rows: list[str] = []
    phase_mismatches: list[str] = []
    tolerance_failures: list[str] = []
    shift_votes: list[str] = []
    z_samples: list[tuple[float, float]] = []
    alpha_samples: list[tuple[float, float]] = []
    vertical_slip_samples: list[tuple[float, float]] = []
    vertical_excess_samples: list[tuple[float, float]] = []
    for measurement in measurements:
        case_id = str(_get(measurement, "caseId", "case_id") or "").strip()
        phase_case_id = str(_get(measurement, "phaseDiagramCaseId", "phase_diagram_case_id") or "").strip()
        label = phase_case_id or case_id
        predicted = by_phase_case_id.get(phase_case_id) if phase_case_id else None
        predicted = predicted or by_case_id.get(case_id)
        if predicted is None:
            missing_predictions.append(label)
            rows.append({"caseId": case_id, "phaseDiagramCaseId": phase_case_id, "status": "missing-prediction"})
            continue
        observed_phase = _phase_label(_get(measurement, "observedPhase", "observedContactPhase"))
        predicted_phase = str(predicted["phase"])
        phase_match = None
        shift_vote = "unobserved-phase"
        if observed_phase:
            phase_count += 1
            phase_match = observed_phase == predicted_phase
            if phase_match:
                phase_matches += 1
            else:
                phase_mismatches.append(label)
            shift_vote = _transition_shift_vote(
                predicted,
                observed_phase,
                brackets_by_phase_case_id.get(str(predicted.get("phaseDiagramCaseId", "")), []),
            )
            shift_votes.append(shift_vote)

        residuals: dict[str, float] = {}
        row_count = 0
        row_max_abs = 0.0
        for observed_key, predicted_key in numeric_pairs.items():
            observed = _optional_float(_get(measurement, observed_key))
            if observed is None:
                continue
            predicted_value = float(predicted[predicted_key])
            residual = observed - predicted_value
            residuals[observed_key] = residual
            squared += residual * residual
            count += 1
            row_count += 1
            max_abs = max(max_abs, abs(residual))
            row_max_abs = max(row_max_abs, abs(residual))
            if observed_key == "observedRightZ":
                z_samples.append((predicted_value, observed))
            elif observed_key == "observedRightAlpha":
                alpha_samples.append((predicted_value - config.initial_alpha, observed - config.initial_alpha))
            elif observed_key == "observedConnectorMaxVerticalSlipMm":
                vertical_slip_samples.append((predicted_value, observed))
            elif observed_key == "observedConnectorMaxVerticalExcessMm":
                vertical_excess_samples.append((predicted_value, observed))

        observed_lock = _optional_bool(str(_get(measurement, "observedLockHeld", "lockHeldObserved") or ""))
        lock_match = None
        if observed_lock is not None and _lock_expected(predicted):
            lock_count += 1
            lock_match = observed_lock is True and bool(predicted.get("lockPass", False))
            if lock_match:
                lock_matches += 1
        if not residuals and not observed_phase and observed_lock is None:
            missing_observation_rows.append(label)
        if row_count and row_max_abs > tolerance:
            tolerance_failures.append(label)
        rows.append(
            {
                "caseId": predicted["caseId"],
                "phaseDiagramCaseId": predicted["phaseDiagramCaseId"],
                "matrixIndex": predicted["matrixIndex"],
                "status": "matched",
                "category": predicted["category"],
                "predictedPhase": predicted_phase,
                "observedPhase": observed_phase,
                "phaseMatches": phase_match,
                "transitionShiftVote": shift_vote,
                "backlash": predicted["backlash"],
                "holeRadius": predicted["holeRadius"],
                "effectiveBacklash": predicted.get("effectiveBacklash", predicted["backlash"]),
                "effectiveHoleRadius": predicted.get("effectiveHoleRadius", predicted["holeRadius"]),
                "effectivePinHoleClearanceMm": predicted.get(
                    "effectivePinHoleClearanceMm",
                    predicted["pinHoleClearanceMm"],
                ),
                "predictedRightZ": predicted["rightZ"],
                "predictedRightAlpha": predicted["rightAlpha"],
                "predictedConnectorMaxVerticalSlipMm": predicted["connectorMaxVerticalSlipMm"],
                "predictedConnectorMaxVerticalExcessMm": predicted["connectorMaxVerticalExcessMm"],
                "observedScalarCount": row_count,
                "maxAbsResidual": row_max_abs if row_count else None,
                "residuals": residuals,
                "lockExpected": _lock_expected(predicted),
                "lockHeldObserved": observed_lock,
                "lockHeldMatches": lock_match,
                "notes": str(_get(measurement, "notes") or ""),
            }
        )
    phase_accuracy = phase_matches / phase_count if phase_count else None
    lock_accuracy = lock_matches / lock_count if lock_count else None
    rms_error = math.sqrt(squared / count) if count else None
    vote_counts = dict(sorted(Counter(shift_votes).items()))
    hole_lower_votes = sum(value for key, value in vote_counts.items() if key == "holeRadius-boundary-lower-than-predicted")
    hole_higher_votes = sum(value for key, value in vote_counts.items() if key == "holeRadius-boundary-higher-than-predicted")
    if hole_lower_votes and not hole_higher_votes:
        transition_direction = "effective-hole-transition-lower-than-reduced-model"
    elif hole_higher_votes and not hole_lower_votes:
        transition_direction = "effective-hole-transition-higher-than-reduced-model"
    elif hole_lower_votes or hole_higher_votes:
        transition_direction = "mixed-transition-shift"
    else:
        transition_direction = "not-enough-transition-phase-evidence"
    hole_bracket_widths = [
        float(item["bracketWidth"])
        for item in report["transitionBrackets"]
        if str(item.get("axis", "")) == "holeRadius"
    ]
    backlash_lower_votes = sum(
        value for key, value in vote_counts.items() if key == "backlash-boundary-lower-than-predicted"
    )
    backlash_higher_votes = sum(
        value for key, value in vote_counts.items() if key == "backlash-boundary-higher-than-predicted"
    )
    backlash_bracket_widths = [
        float(item["bracketWidth"])
        for item in report["transitionBrackets"]
        if str(item.get("axis", "")) == "backlash"
    ]
    hole_bias = 0.0
    if hole_lower_votes and not hole_higher_votes:
        hole_bias = 0.5 * _mean(hole_bracket_widths) if hole_bracket_widths else 0.0
    elif hole_higher_votes and not hole_lower_votes:
        hole_bias = -0.5 * _mean(hole_bracket_widths) if hole_bracket_widths else 0.0
    backlash_bias = 0.0
    if backlash_lower_votes and not backlash_higher_votes:
        backlash_bias = 0.5 * _mean(backlash_bracket_widths) if backlash_bracket_widths else 0.0
    elif backlash_higher_votes and not backlash_lower_votes:
        backlash_bias = -0.5 * _mean(backlash_bracket_widths) if backlash_bracket_widths else 0.0
    z_fit = _fit_scale(z_samples)
    alpha_fit = _fit_scale(alpha_samples)
    vertical_slip_fit = _fit_scale(vertical_slip_samples)
    vertical_excess_fit = _fit_scale(vertical_excess_samples)
    z_scale = _bounded(z_fit["estimate"], 0.25, 4.0, 1.0)
    alpha_scale = _bounded(alpha_fit["estimate"], 0.25, 4.0, 1.0)
    vertical_slip_scale = _bounded(vertical_slip_fit["estimate"], 0.25, 4.0, 1.0)
    current_clearance = max(0.0, config.pin_hole_clearance)
    proposed_clearance = max(0.0, current_clearance * vertical_slip_scale + hole_bias)
    proposed = {
        "couplingGain": _bounded(config.coupling_gain * alpha_scale, 0.0, 1.0, config.coupling_gain),
        "zCouplingGain": _bounded(config.z_coupling_gain * z_scale, 0.0, 1.0, config.z_coupling_gain),
        "backlash": max(0.0, config.backlash + backlash_bias),
        "holeRadius": max(config.pin_radius, config.pin_radius + proposed_clearance),
        "pinHoleClearance": proposed_clearance,
        "effectiveHoleRadiusBias": hole_bias,
        "effectiveBacklashBias": backlash_bias,
        "transitionBoundaryDirection": transition_direction,
        "applyAutomatically": False,
        "reason": "reduced-proxy transition estimate only; require holdout bench and segmented contact validation before physical claims",
    }
    missing_evidence = []
    if not measurements:
        missing_evidence.append("filledTransitionMeasurementRows")
    if missing_predictions:
        missing_evidence.append("matchedTransitionCaseRows")
    if missing_observation_rows:
        missing_evidence.append("observedTransitionRows")
    if phase_count == 0:
        missing_evidence.append("observedTransitionPhases")
    if count == 0:
        missing_evidence.append("observedTransitionNumericFields")
    if phase_mismatches:
        missing_evidence.append("transitionPhaseAgreement")
    if tolerance_failures:
        missing_evidence.append("transitionNumericTolerance")
    passes_tolerance = (
        bool(measurements)
        and not missing_predictions
        and not tolerance_failures
        and (phase_count == 0 or not phase_mismatches)
        and count > 0
    )
    return {
        "schema": TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA,
        "model": "two-cell-radius-backlash-transition-measured-vs-reduced-proxy-comparison",
        "cadReference": report["cadReference"],
        "transitionReportSummary": report["summary"],
        "rows": rows,
        "summary": {
            "providedRowCount": len(measurements),
            "matchedRowCount": sum(1 for row in rows if row.get("status") == "matched"),
            "missingPredictionCount": len(missing_predictions),
            "missingObservationRowCount": len(missing_observation_rows),
            "observedScalarCount": count,
            "rmsError": rms_error,
            "maxAbsError": max_abs if count else None,
            "phaseObservationCount": phase_count,
            "phaseAccuracy": phase_accuracy,
            "phaseMismatchCount": len(phase_mismatches),
            "lockObservationCount": lock_count,
            "lockHeldAccuracy": lock_accuracy,
            "transitionShiftVoteCounts": vote_counts,
            "transitionShiftDirection": transition_direction,
            "passesTolerance": passes_tolerance,
            "physicalAccuracyValidated": False,
            "missingEvidence": sorted(set(missing_evidence)),
        },
        "calibrationEstimate": {
            "rightZScale": z_fit,
            "rightAlphaResponseScale": alpha_fit,
            "verticalSlipScale": vertical_slip_fit,
            "verticalExcessScale": vertical_excess_fit,
            "transitionBoundaryDirection": transition_direction,
            "proposedReducedProxyUpdates": proposed,
            "applyAutomatically": False,
            "reason": "transition rows estimate reduced-proxy correction only; real geometry still needs segmented CAD/contact and holdout measurements",
        },
        "acceptance": {
            "readyForReducedProxyTransitionCalibration": count >= 4 and phase_count >= 2 and not missing_predictions,
            "passesTolerance": passes_tolerance,
            "requiresMoreData": bool(missing_evidence),
            "requiresSegmentedPhysics": True,
            "requiresHoldoutValidation": True,
            "physicalAccuracyValidated": False,
        },
        "claimBoundary": {
            "allowedClaim": "compares filled transition-adjacent bench rows against reduced-model boundary predictions",
            "blockedClaim": "exact real RAD transition law without segmented CAD, contact parameters, and independent physical holdout sweeps",
            "physicalAccuracyValidated": False,
        },
    }


def _transition_midpoint_value(bracket: dict[str, Any]) -> float | None:
    if str(bracket.get("axis", "")) == "holeRadius":
        return _finite_float(bracket.get("midHoleRadius"))
    if str(bracket.get("axis", "")) == "backlash":
        return _finite_float(bracket.get("midBacklash"))
    return None


def _transition_fixed_value(bracket: dict[str, Any]) -> float | None:
    if str(bracket.get("axis", "")) == "holeRadius":
        return _finite_float(bracket.get("fixedBacklash"))
    if str(bracket.get("axis", "")) == "backlash":
        return _finite_float(bracket.get("fixedHoleRadius"))
    return None


def _transition_boundary_movement_rows(
    baseline: list[dict[str, Any]],
    calibrated: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    movement_rows: list[dict[str, Any]] = []
    remaining = list(calibrated)
    for index, before in enumerate(baseline, start=1):
        before_axis = str(before.get("axis", ""))
        before_pair = (str(before.get("fromPhase", "")), str(before.get("toPhase", "")))
        before_midpoint = _transition_midpoint_value(before)
        before_fixed = _transition_fixed_value(before)
        best_index: int | None = None
        best_score: tuple[float, float, float] | None = None
        for candidate_index, candidate in enumerate(remaining):
            axis_penalty = 0.0 if str(candidate.get("axis", "")) == before_axis else 1.0
            pair_penalty = (
                0.0
                if (str(candidate.get("fromPhase", "")), str(candidate.get("toPhase", ""))) == before_pair
                else 1.0
            )
            candidate_fixed = _transition_fixed_value(candidate)
            candidate_midpoint = _transition_midpoint_value(candidate)
            fixed_delta = (
                abs(float(candidate_fixed) - float(before_fixed))
                if candidate_fixed is not None and before_fixed is not None
                else 1e6
            )
            midpoint_delta = (
                abs(float(candidate_midpoint) - float(before_midpoint))
                if candidate_midpoint is not None and before_midpoint is not None
                else 1e6
            )
            score = (axis_penalty + pair_penalty, fixed_delta, midpoint_delta)
            if best_score is None or score < best_score:
                best_score = score
                best_index = candidate_index
        after = remaining.pop(best_index) if best_index is not None and remaining else None
        after_midpoint = _transition_midpoint_value(after) if after else None
        movement_rows.append(
            {
                "index": index,
                "status": "matched-nearest" if after else "missing-after-calibration",
                "baselineTransitionId": before.get("transitionId", ""),
                "calibratedTransitionId": after.get("transitionId", "") if after else "",
                "axis": before_axis,
                "baselineFromPhase": before.get("fromPhase", ""),
                "baselineToPhase": before.get("toPhase", ""),
                "calibratedFromPhase": after.get("fromPhase", "") if after else "",
                "calibratedToPhase": after.get("toPhase", "") if after else "",
                "baselineMidpoint": before_midpoint,
                "calibratedMidpoint": after_midpoint,
                "midpointDelta": (
                    float(after_midpoint) - float(before_midpoint)
                    if after_midpoint is not None and before_midpoint is not None
                    else None
                ),
                "baselineFixedValue": before_fixed,
                "calibratedFixedValue": _transition_fixed_value(after) if after else None,
                "physicalAccuracyValidated": False,
            }
        )
    start = len(movement_rows) + 1
    for offset, after in enumerate(remaining):
        movement_rows.append(
            {
                "index": start + offset,
                "status": "new-after-calibration",
                "baselineTransitionId": "",
                "calibratedTransitionId": after.get("transitionId", ""),
                "axis": after.get("axis", ""),
                "baselineFromPhase": "",
                "baselineToPhase": "",
                "calibratedFromPhase": after.get("fromPhase", ""),
                "calibratedToPhase": after.get("toPhase", ""),
                "baselineMidpoint": None,
                "calibratedMidpoint": _transition_midpoint_value(after),
                "midpointDelta": None,
                "baselineFixedValue": None,
                "calibratedFixedValue": _transition_fixed_value(after),
                "physicalAccuracyValidated": False,
            }
        )
    return movement_rows


def two_cell_radius_backlash_transition_rerun_report(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> dict[str, Any]:
    """Rerun the reduced transition diagram using measured boundary-bias estimates."""

    config = config or LatticeConfig(rows=1, cols=2)
    measurements = two_cell_radius_backlash_transition_measurements_from_rows(measurements)
    baseline_comparison = compare_two_cell_radius_backlash_transition_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
        **report_overrides,
    )
    updates = baseline_comparison["calibrationEstimate"]["proposedReducedProxyUpdates"]
    hole_bias = float(updates.get("effectiveHoleRadiusBias", 0.0))
    backlash_bias = float(updates.get("effectiveBacklashBias", 0.0))
    calibrated_overrides = dict(report_overrides)
    calibrated_overrides["hole_radius_bias"] = hole_bias
    calibrated_overrides["backlash_bias"] = backlash_bias
    baseline_report = two_cell_radius_backlash_transition_report(config, base_controls, **report_overrides)
    calibrated_report = two_cell_radius_backlash_transition_report(
        config,
        base_controls,
        **calibrated_overrides,
    )
    calibrated_comparison = compare_two_cell_radius_backlash_transition_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
        **calibrated_overrides,
    )
    movement_rows = _transition_boundary_movement_rows(
        baseline_report["transitionBrackets"],
        calibrated_report["transitionBrackets"],
    )
    finite_deltas = [
        abs(float(row["midpointDelta"]))
        for row in movement_rows
        if row.get("midpointDelta") is not None
    ]
    return {
        "schema": TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA,
        "model": "two-cell-radius-backlash-transition-calibrated-rerun-report",
        "cadReference": baseline_report["cadReference"],
        "axisBiases": {
            "effectiveHoleRadiusBias": hole_bias,
            "effectiveBacklashBias": backlash_bias,
            "transitionBoundaryDirection": updates.get("transitionBoundaryDirection", ""),
        },
        "proposedReducedProxyUpdates": updates,
        "baseline": baseline_report,
        "calibrated": calibrated_report,
        "baselineComparison": baseline_comparison,
        "calibratedComparison": calibrated_comparison,
        "transitionMovement": movement_rows,
        "summary": {
            "status": "calibrated-transition-rerun-ready-needs-holdout-validation",
            "providedMeasurementRowCount": len(measurements),
            "baselineTransitionCount": baseline_report["summary"]["transitionBracketCount"],
            "calibratedTransitionCount": calibrated_report["summary"]["transitionBracketCount"],
            "transitionCountDelta": (
                calibrated_report["summary"]["transitionBracketCount"]
                - baseline_report["summary"]["transitionBracketCount"]
            ),
            "baselinePhaseAccuracy": baseline_comparison["summary"]["phaseAccuracy"],
            "calibratedPhaseAccuracy": calibrated_comparison["summary"]["phaseAccuracy"],
            "transitionMovementCount": len(movement_rows),
            "maxAbsMidpointDelta": max(finite_deltas) if finite_deltas else 0.0,
            "transitionShiftDirection": updates.get("transitionBoundaryDirection", ""),
            "physicalAccuracyValidated": False,
            "remainingEvidence": [
                "independent holdout transition measurements after applying the proxy correction",
                "segmented CAD bodies with pin-hole contact surfaces",
                "external rigid-body contact results using the same nominal geometry",
                "bench repeatability across at least two physical assemblies",
            ],
        },
        "claimBoundary": {
            "allowedClaim": "shows how measured transition votes move the reduced-model effective boundary",
            "blockedClaim": "validated physical transition law until holdout sweeps and segmented CAD/contact agree",
            "physicalAccuracyValidated": False,
        },
    }


def calibrate_two_cell_parameters(
    measurements: Iterable[TwoCellMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    backlash_candidates: Iterable[float] | None = None,
    hole_radius_candidates: Iterable[float] | None = None,
    coupling_gain_candidates: Iterable[float] | None = None,
    z_coupling_gain_candidates: Iterable[float] | None = None,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    measurements = list(measurements)
    backlash_candidates = list(backlash_candidates or [config.backlash * t for t in (0.5, 0.75, 1.0, 1.25, 1.5)])
    hole_radius_candidates = list(hole_radius_candidates or [config.hole_radius, config.pin_radius + 0.5 * config.pin_hole_clearance, config.pin_radius + 1.5 * config.pin_hole_clearance])
    coupling_gain_candidates = list(coupling_gain_candidates or [config.coupling_gain])
    z_coupling_gain_candidates = list(z_coupling_gain_candidates or [config.z_coupling_gain * t for t in (0.5, 0.75, 1.0, 1.25, 1.5)])
    candidates = []
    for backlash in backlash_candidates:
        for hole_radius in hole_radius_candidates:
            if hole_radius < config.pin_radius:
                continue
            for coupling_gain in coupling_gain_candidates:
                for z_coupling_gain in z_coupling_gain_candidates:
                    trial = LatticeConfig(
                        rows=1,
                        cols=2,
                        backlash=max(0.0, backlash),
                        cell_size=config.cell_size,
                        initial_alpha=config.initial_alpha,
                        coupling_gain=max(0.0, min(1.0, coupling_gain)),
                        z_coupling_gain=max(0.0, min(1.0, z_coupling_gain)),
                        pin_radius=config.pin_radius,
                        hole_radius=hole_radius,
                        alpha_min=config.alpha_min,
                        alpha_max=config.alpha_max,
                    )
                    comparison = compare_two_cell_measurements(measurements, trial, base_controls)
                    candidates.append(
                        {
                            "backlash": trial.backlash,
                            "holeRadius": trial.hole_radius,
                            "pinHoleClearance": trial.pin_hole_clearance,
                            "couplingGain": trial.coupling_gain,
                            "zCouplingGain": trial.z_coupling_gain,
                            "rmsError": comparison["rmsError"],
                            "maxAbsError": comparison["maxAbsError"],
                            "observedScalarCount": comparison["observedScalarCount"],
                        }
                    )
    valid = [candidate for candidate in candidates if candidate["rmsError"] is not None]
    best = min(valid, key=lambda item: (item["rmsError"], item["maxAbsError"])) if valid else None
    return {
        "schema": TWO_CELL_PARAMETER_CALIBRATION_SCHEMA,
        "candidateCount": len(candidates),
        "observedScalarCount": best["observedScalarCount"] if best else 0,
        "best": best,
        "topCandidates": sorted(valid, key=lambda item: (item["rmsError"], item["maxAbsError"]))[:10],
        "baseline": compare_two_cell_measurements(measurements, config, base_controls),
        "measurementTemplate": two_cell_physical_test_packet(config, base_controls)["measurementTemplateColumns"],
        "interpretation": {
            "calibrationScope": "reduced two-cell proxy only",
            "requiresSegmentedPhysics": True,
            "nextStep": "measure the packet cases, fit this proxy, then compare against an external rigid-body contact model",
        },
    }


def export_two_cell_measurement_comparison_json(
    measurements: Iterable[TwoCellMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(compare_two_cell_measurements(measurements, config, base_controls), indent=2)


def export_two_cell_connector_measurement_comparison_json(
    measurements: Iterable[TwoCellConnectorMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> str:
    return json.dumps(
        compare_two_cell_connector_measurements(
            measurements,
            config,
            base_controls,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_two_cell_connector_measurement_comparison_csv(
    measurements: Iterable[TwoCellConnectorMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> str:
    report = compare_two_cell_connector_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
    )
    fields = [
        "caseId",
        "connector",
        "status",
        "lockMode",
        "predictedContactMode",
        "observedContactMode",
        "contactModeMatches",
        "lockExpected",
        "lockHeldObserved",
        "lockHeldMatches",
        "markerRmsErrorMm",
        "slipRmsErrorMm",
        "maxAbsResidualMm",
        "availableFields",
        "measurementSource",
        "notes",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in report["rows"]:
        writer.writerow(
            {
                field: ";".join(row[field]) if field == "availableFields" and isinstance(row.get(field), list)
                else row.get(field, "")
                for field in fields
            }
        )
    return output.getvalue()


def export_two_cell_fidelity_matrix_measurement_comparison_json(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> str:
    return json.dumps(
        compare_two_cell_fidelity_matrix_measurements(
            measurements,
            config,
            base_controls,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_two_cell_fidelity_matrix_measurement_comparison_csv(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> str:
    report = compare_two_cell_fidelity_matrix_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
    )
    fields = [
        "caseId",
        "matrixIndex",
        "connector",
        "status",
        "actuationCase",
        "lockMode",
        "predictedContactMode",
        "observedContactMode",
        "contactModeMatches",
        "lockExpected",
        "lockHeldObserved",
        "lockHeldMatches",
        "observedScalarCount",
        "maxAbsResidual",
        "measurementSource",
        "notes",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in report["rows"]:
        writer.writerow({field: row.get(field, "") for field in fields})
    return output.getvalue()


def export_two_cell_fidelity_matrix_parameter_calibration_json(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> str:
    return json.dumps(
        calibrate_two_cell_fidelity_matrix_parameters(
            measurements,
            config,
            base_controls,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_two_cell_radius_backlash_transition_comparison_json(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> str:
    return json.dumps(
        compare_two_cell_radius_backlash_transition_measurements(
            measurements,
            config,
            base_controls,
            tolerance=tolerance,
            **report_overrides,
        ),
        indent=2,
    )


def export_two_cell_radius_backlash_transition_comparison_csv(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> str:
    report = compare_two_cell_radius_backlash_transition_measurements(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
        **report_overrides,
    )
    fields = [
        "caseId",
        "phaseDiagramCaseId",
        "matrixIndex",
        "status",
        "category",
        "predictedPhase",
        "observedPhase",
        "phaseMatches",
        "transitionShiftVote",
        "backlash",
        "holeRadius",
        "effectiveBacklash",
        "effectiveHoleRadius",
        "effectivePinHoleClearanceMm",
        "predictedRightZ",
        "predictedRightAlpha",
        "predictedConnectorMaxVerticalSlipMm",
        "predictedConnectorMaxVerticalExcessMm",
        "observedScalarCount",
        "maxAbsResidual",
        "lockExpected",
        "lockHeldObserved",
        "lockHeldMatches",
        "notes",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in report["rows"]:
        writer.writerow({field: row.get(field, "") for field in fields})
    return output.getvalue()


def export_two_cell_radius_backlash_transition_rerun_json(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> str:
    return json.dumps(
        two_cell_radius_backlash_transition_rerun_report(
            measurements,
            config,
            base_controls,
            tolerance=tolerance,
            **report_overrides,
        ),
        indent=2,
    )


def export_two_cell_radius_backlash_transition_rerun_csv(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
    **report_overrides: Any,
) -> str:
    report = two_cell_radius_backlash_transition_rerun_report(
        measurements,
        config,
        base_controls,
        tolerance=tolerance,
        **report_overrides,
    )
    fields = [
        "index",
        "status",
        "baselineTransitionId",
        "calibratedTransitionId",
        "axis",
        "baselineFromPhase",
        "baselineToPhase",
        "calibratedFromPhase",
        "calibratedToPhase",
        "baselineMidpoint",
        "calibratedMidpoint",
        "midpointDelta",
        "baselineFixedValue",
        "calibratedFixedValue",
        "physicalAccuracyValidated",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in report["transitionMovement"]:
        writer.writerow({field: row.get(field, "") for field in fields})
    return output.getvalue()


def load_two_cell_connector_measurements(path: str | Path) -> list[TwoCellConnectorMeasurement]:
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if source.suffix.lower() == ".json":
        return two_cell_connector_measurements_from_json(text)
    return two_cell_connector_measurements_from_csv(text)


def load_two_cell_fidelity_matrix_measurements(path: str | Path) -> list[dict[str, Any]]:
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if source.suffix.lower() == ".json":
        return two_cell_fidelity_matrix_measurements_from_json(text)
    return two_cell_fidelity_matrix_measurements_from_csv(text)


def load_two_cell_radius_backlash_transition_measurements(path: str | Path) -> list[dict[str, Any]]:
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    if source.suffix.lower() == ".json":
        return two_cell_radius_backlash_transition_measurements_from_json(text)
    return two_cell_radius_backlash_transition_measurements_from_csv(text)


def write_two_cell_connector_measurement_comparison_artifacts(
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
    measurements = load_two_cell_connector_measurements(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    files = {
        "comparison": output_dir / "two_cell_connector_measurement_comparison.json",
        "comparison_csv": output_dir / "two_cell_connector_measurement_comparison.csv",
    }
    files["comparison"].write_text(
        export_two_cell_connector_measurement_comparison_json(
            measurements,
            config,
            tolerance=tolerance,
        ),
        encoding="utf-8",
    )
    files["comparison_csv"].write_text(
        export_two_cell_connector_measurement_comparison_csv(
            measurements,
            config,
            tolerance=tolerance,
        ),
        encoding="utf-8",
    )
    return files


def write_two_cell_fidelity_matrix_measurement_comparison_artifacts(
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
    measurements = load_two_cell_fidelity_matrix_measurements(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    files = {
        "comparison": output_dir / "two_cell_fidelity_matrix_measurement_comparison.json",
        "comparison_csv": output_dir / "two_cell_fidelity_matrix_measurement_comparison.csv",
        "calibration": output_dir / "two_cell_fidelity_matrix_parameter_calibration.json",
    }
    files["comparison"].write_text(
        export_two_cell_fidelity_matrix_measurement_comparison_json(
            measurements,
            config,
            tolerance=tolerance,
        ),
        encoding="utf-8",
    )
    files["comparison_csv"].write_text(
        export_two_cell_fidelity_matrix_measurement_comparison_csv(
            measurements,
            config,
            tolerance=tolerance,
        ),
        encoding="utf-8",
    )
    files["calibration"].write_text(
        export_two_cell_fidelity_matrix_parameter_calibration_json(
            measurements,
            config,
            tolerance=tolerance,
        ),
        encoding="utf-8",
    )
    return files


def write_two_cell_radius_backlash_transition_comparison_artifacts(
    input_path: str | Path,
    out: str | Path,
    *,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
    tolerance: float = 1e-6,
    hole_steps: int = 9,
    backlash_steps: int = 9,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    measurements = load_two_cell_radius_backlash_transition_measurements(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    overrides = {"hole_steps": hole_steps, "backlash_steps": backlash_steps}
    files = {
        "comparison": output_dir / "two_cell_radius_backlash_transition_comparison.json",
        "comparison_csv": output_dir / "two_cell_radius_backlash_transition_comparison.csv",
        "rerun": output_dir / "two_cell_radius_backlash_transition_rerun.json",
        "rerun_csv": output_dir / "two_cell_radius_backlash_transition_rerun.csv",
    }
    files["comparison"].write_text(
        export_two_cell_radius_backlash_transition_comparison_json(
            measurements,
            config,
            tolerance=tolerance,
            **overrides,
        ),
        encoding="utf-8",
    )
    files["comparison_csv"].write_text(
        export_two_cell_radius_backlash_transition_comparison_csv(
            measurements,
            config,
            tolerance=tolerance,
            **overrides,
        ),
        encoding="utf-8",
    )
    files["rerun"].write_text(
        export_two_cell_radius_backlash_transition_rerun_json(
            measurements,
            config,
            tolerance=tolerance,
            **overrides,
        ),
        encoding="utf-8",
    )
    files["rerun_csv"].write_text(
        export_two_cell_radius_backlash_transition_rerun_csv(
            measurements,
            config,
            tolerance=tolerance,
            **overrides,
        ),
        encoding="utf-8",
    )
    return files


def _load_json_file(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _finite_float(value: Any) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _rms(values: list[float]) -> float | None:
    return math.sqrt(sum(value * value for value in values) / len(values)) if values else None


def _mode_counts(values: Iterable[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value)
        if not key:
            continue
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _group_external_rows(
    measurement_rows: list[dict[str, Any]],
    comparison_by_key: dict[tuple[str, str], dict[str, Any]],
    axis: str,
) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in measurement_rows:
        groups.setdefault(str(row.get(axis, "")), []).append(row)
    summaries = []
    for value, rows in sorted(groups.items(), key=lambda item: item[0]):
        comparison_rows = [
            comparison_by_key[(str(row.get("caseId", "")), str(row.get("connector", "")))]
            for row in rows
            if (str(row.get("caseId", "")), str(row.get("connector", ""))) in comparison_by_key
        ]
        lateral = [_finite_float(row.get("observedLateralSlipMm")) for row in rows]
        vertical = [_finite_float(row.get("observedVerticalSlipMm")) for row in rows]
        total = [_finite_float(row.get("observedTotalSlipMm")) for row in rows]
        lateral_values = [value for value in lateral if value is not None]
        vertical_values = [value for value in vertical if value is not None]
        total_values = [value for value in total if value is not None]
        residuals = [
            abs(float(residual))
            for comparison in comparison_rows
            for residual in dict(comparison.get("residuals", {})).values()
            if _finite_float(residual) is not None
        ]
        contact_flags = [
            bool(comparison.get("contactModeMatches"))
            for comparison in comparison_rows
            if comparison.get("contactModeMatches") is not None
        ]
        lock_flags = [
            bool(comparison.get("lockHeldMatches"))
            for comparison in comparison_rows
            if comparison.get("lockHeldMatches") is not None
        ]
        summaries.append(
            {
                "axis": axis,
                "value": value,
                "rowCount": len(rows),
                "caseCount": len({str(row.get("caseId", "")) for row in rows}),
                "meanObservedLateralSlipMm": _mean(lateral_values),
                "maxObservedLateralSlipMm": max(lateral_values, default=None),
                "meanObservedVerticalSlipMm": _mean(vertical_values),
                "maxObservedVerticalSlipMm": max(vertical_values, default=None),
                "meanObservedTotalSlipMm": _mean(total_values),
                "maxObservedTotalSlipMm": max(total_values, default=None),
                "residualRms": _rms(residuals),
                "maxAbsResidual": max(residuals, default=None),
                "contactModeAccuracy": _mean([1.0 if flag else 0.0 for flag in contact_flags]),
                "lockHeldAccuracy": _mean([1.0 if flag else 0.0 for flag in lock_flags]),
                "observedContactModes": _mode_counts(row.get("observedContactMode", "") for row in rows),
            }
        )
    return summaries


def two_cell_external_fidelity_benchmark_summary(
    run_report: dict[str, Any],
    comparison_report: dict[str, Any],
    calibration_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    calibration_report = calibration_report or {}
    measurement_rows = [
        dict(row)
        for row in run_report.get("results", {}).get("measurementRows", [])
        if isinstance(row, dict)
    ]
    comparison_rows = [
        dict(row)
        for row in comparison_report.get("rows", [])
        if isinstance(row, dict)
    ]
    comparison_by_key = {
        (str(row.get("caseId", "")), str(row.get("connector", ""))): row
        for row in comparison_rows
    }
    axis_summaries = {
        "byHoleRadius": _group_external_rows(measurement_rows, comparison_by_key, "holeRadius"),
        "byBacklash": _group_external_rows(measurement_rows, comparison_by_key, "backlash"),
        "byLockMode": _group_external_rows(measurement_rows, comparison_by_key, "lockMode"),
        "byActuationCase": _group_external_rows(measurement_rows, comparison_by_key, "actuationCase"),
    }
    all_axis_rows = [row for rows in axis_summaries.values() for row in rows]
    worst_axis_rows = sorted(
        [row for row in all_axis_rows if row.get("residualRms") is not None],
        key=lambda row: float(row["residualRms"]),
        reverse=True,
    )[:8]
    acceptance = dict(comparison_report.get("acceptance", {}))
    run_summary = dict(run_report.get("summary", {}))
    comparison_summary = {
        "sampleCount": comparison_report.get("sampleCount", 0),
        "templateRowCount": comparison_report.get("templateRowCount", 0),
        "matchedRowCount": comparison_report.get("matchedRowCount", 0),
        "observedScalarCount": comparison_report.get("observedScalarCount", 0),
        "rmsError": comparison_report.get("rmsError"),
        "maxAbsError": comparison_report.get("maxAbsError"),
        "cellRmsError": comparison_report.get("cellRmsError"),
        "connectorRmsError": comparison_report.get("connectorRmsError"),
        "contactModeAccuracy": comparison_report.get("contactModeAccuracy"),
        "lockHeldAccuracy": comparison_report.get("lockHeldAccuracy"),
        "passesTolerance": acceptance.get("passesTolerance", False),
        "coversFullFidelityMatrix": acceptance.get("coversFullFidelityMatrix", False),
        "missingEvidence": list(comparison_report.get("missingEvidence", [])),
    }
    remaining_evidence = sorted(
        set(
            [
                *list(run_summary.get("missingEvidence", [])),
                *list(comparison_summary.get("missingEvidence", [])),
                "segmentedCadContactModel",
                "benchCoordinateHoldout",
                "measuredFrictionAndStiffness",
            ]
        )
    )
    return {
        "schema": TWO_CELL_EXTERNAL_FIDELITY_BENCHMARK_SCHEMA,
        "model": "two-cell-external-mujoco-proxy-full-matrix-benchmark",
        "runSummary": run_summary,
        "comparisonSummary": comparison_summary,
        "calibrationSummary": calibration_report.get("comparisonSummary", {}),
        "calibrationEstimates": calibration_report.get("estimates", {}),
        "proposedReducedProxyUpdates": calibration_report.get("proposedReducedProxyUpdates", {}),
        "axisSummaries": axis_summaries,
        "worstAxisResiduals": worst_axis_rows,
        "observedContactModes": _mode_counts(row.get("observedContactMode", "") for row in measurement_rows),
        "summary": {
            "externalRunComplete": bool(run_summary.get("externalFidelityMjcfRunComplete")),
            "coversFullFidelityMatrix": bool(acceptance.get("coversFullFidelityMatrix")),
            "passesReducedModelTolerance": bool(acceptance.get("passesTolerance")),
            "readyForReducedProxyCalibration": bool(acceptance.get("readyForMatrixCalibration")),
            "physicalAccuracyValidated": False,
            "remainingEvidence": remaining_evidence,
            "remainingEvidenceCount": len(remaining_evidence),
        },
        "claimLabels": {
            "externalRun": "real MuJoCo execution of generated proxy MJCF cases",
            "calibration": "reduced-proxy calibration evidence, not a physical validation",
            "physicalAccuracy": "false until segmented CAD and bench holdout evidence pass",
        },
    }


def write_two_cell_external_fidelity_benchmark_summary_artifact(
    run_path: str | Path,
    comparison_path: str | Path,
    calibration_path: str | Path | None,
    out: str | Path,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    calibration = _load_json_file(calibration_path) if calibration_path is not None else {}
    summary = two_cell_external_fidelity_benchmark_summary(
        _load_json_file(run_path),
        _load_json_file(comparison_path),
        calibration,
    )
    files = {"benchmark": output_dir / "two_cell_external_fidelity_benchmark_summary.json"}
    files["benchmark"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return files


FIDELITY_CORRECTION_FIELD_PAIRS = {
    "observedLeftCellX": "predictedLeftCellX",
    "observedLeftCellY": "predictedLeftCellY",
    "observedLeftCellZ": "predictedLeftCellZ",
    "observedRightCellX": "predictedRightCellX",
    "observedRightCellY": "predictedRightCellY",
    "observedRightCellZ": "predictedRightCellZ",
    "observedLeftXmm": "predictedLeftXmm",
    "observedLeftYmm": "predictedLeftYmm",
    "observedLeftZmm": "predictedLeftZmm",
    "observedRightXmm": "predictedRightXmm",
    "observedRightYmm": "predictedRightYmm",
    "observedRightZmm": "predictedRightZmm",
    "observedLateralSlipMm": "predictedLateralSlipMm",
    "observedVerticalSlipMm": "predictedVerticalSlipMm",
    "observedTotalSlipMm": "predictedTotalSlipMm",
}


def _case_sort_key(case_id: str) -> tuple[int, str]:
    match = re.search(r"(\d+)$", str(case_id))
    return (int(match.group(1)) if match else 10**9, str(case_id))


def _split_case_ids_for_holdout(
    rows: list[dict[str, Any]],
    *,
    holdout_stride: int,
    holdout_case_ids: Iterable[str] | None,
) -> tuple[set[str], set[str]]:
    case_ids = sorted({str(row.get("caseId", "")) for row in rows if str(row.get("caseId", ""))}, key=_case_sort_key)
    explicit = {str(case_id).strip() for case_id in list(holdout_case_ids or []) if str(case_id).strip()}
    if explicit:
        holdout = explicit & set(case_ids)
    else:
        stride = max(2, int(holdout_stride))
        holdout = {case_id for index, case_id in enumerate(case_ids) if index % stride == stride - 1}
        if not holdout and case_ids:
            holdout = {case_ids[-1]}
    train = set(case_ids) - holdout
    if not train and holdout:
        moved = sorted(holdout, key=_case_sort_key)[0]
        holdout.remove(moved)
        train.add(moved)
    return train, holdout


def _fit_affine_correction(samples: list[tuple[float, float]]) -> dict[str, Any]:
    if not samples:
        return {
            "status": "missing-observations",
            "sampleCount": 0,
            "slope": 1.0,
            "intercept": 0.0,
            "baselineRms": None,
            "correctedRms": None,
            "improvement": None,
        }
    xs = [sample[0] for sample in samples]
    ys = [sample[1] for sample in samples]
    mean_x = _mean(xs) or 0.0
    mean_y = _mean(ys) or 0.0
    denominator = sum((x - mean_x) ** 2 for x in xs)
    if denominator > 1e-12:
        slope = sum((x - mean_x) * (y - mean_y) for x, y in samples) / denominator
        intercept = mean_y - slope * mean_x
        status = "affine-fit"
    else:
        slope = 1.0
        intercept = _mean([y - x for x, y in samples]) or 0.0
        status = "bias-only-fit"
    baseline_residuals = [y - x for x, y in samples]
    corrected_residuals = [y - (slope * x + intercept) for x, y in samples]
    baseline_rms = _rms(baseline_residuals)
    corrected_rms = _rms(corrected_residuals)
    improvement = None
    if baseline_rms is not None and baseline_rms > 1e-12 and corrected_rms is not None:
        improvement = 1.0 - corrected_rms / baseline_rms
    return {
        "status": status,
        "sampleCount": len(samples),
        "slope": slope,
        "intercept": intercept,
        "baselineRms": baseline_rms,
        "correctedRms": corrected_rms,
        "improvement": improvement,
    }


def _evaluate_corrections(
    rows: list[dict[str, Any]],
    prediction_map: dict[tuple[str, str], dict[str, Any]],
    correction_fields: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    baseline_residuals: list[float] = []
    corrected_residuals: list[float] = []
    cell_baseline: list[float] = []
    cell_corrected: list[float] = []
    connector_baseline: list[float] = []
    connector_corrected: list[float] = []
    contact_count = 0
    contact_matches = 0
    lock_count = 0
    lock_matches = 0
    for row in rows:
        case_id = str(row.get("caseId", ""))
        connector = str(row.get("connector", ""))
        predicted = prediction_map.get((case_id, connector))
        if predicted is None:
            continue
        for observed_key, predicted_key in FIDELITY_CORRECTION_FIELD_PAIRS.items():
            observed = _finite_float(_get(row, observed_key))
            if observed is None:
                continue
            predicted_value = _finite_float(predicted.get(predicted_key))
            if predicted_value is None:
                continue
            correction = correction_fields.get(observed_key, {})
            corrected_value = float(correction.get("slope", 1.0)) * predicted_value + float(
                correction.get("intercept", 0.0)
            )
            baseline = observed - predicted_value
            corrected = observed - corrected_value
            baseline_residuals.append(baseline)
            corrected_residuals.append(corrected)
            if "Cell" in observed_key:
                cell_baseline.append(baseline)
                cell_corrected.append(corrected)
            else:
                connector_baseline.append(baseline)
                connector_corrected.append(corrected)
        observed_contact = str(_get(row, "observedContactMode", "contactModeObserved") or "").strip()
        if observed_contact:
            contact_count += 1
            if observed_contact.lower() == str(predicted.get("predictedContactMode", "")).lower():
                contact_matches += 1
        expected_lock = _lock_expected(predicted)
        observed_lock = _optional_bool(str(_get(row, "lockHeldObserved", "lock_held_observed") or ""))
        if expected_lock and observed_lock is not None:
            lock_count += 1
            if observed_lock is True:
                lock_matches += 1
    baseline_rms = _rms(baseline_residuals)
    corrected_rms = _rms(corrected_residuals)
    return {
        "rowCount": len(rows),
        "caseCount": len({str(row.get("caseId", "")) for row in rows if str(row.get("caseId", ""))}),
        "observedScalarCount": len(baseline_residuals),
        "baselineRms": baseline_rms,
        "correctedRms": corrected_rms,
        "baselineMaxAbs": max((abs(value) for value in baseline_residuals), default=None),
        "correctedMaxAbs": max((abs(value) for value in corrected_residuals), default=None),
        "cellBaselineRms": _rms(cell_baseline),
        "cellCorrectedRms": _rms(cell_corrected),
        "connectorBaselineRms": _rms(connector_baseline),
        "connectorCorrectedRms": _rms(connector_corrected),
        "improvement": 1.0 - corrected_rms / baseline_rms
        if baseline_rms is not None and baseline_rms > 1e-12 and corrected_rms is not None
        else None,
        "contactModeAccuracy": contact_matches / contact_count if contact_count else None,
        "lockHeldAccuracy": lock_matches / lock_count if lock_count else None,
    }


def fit_two_cell_external_fidelity_correction_profile(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    *,
    holdout_stride: int = 5,
    holdout_case_ids: Iterable[str] | None = None,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    rows = two_cell_fidelity_matrix_measurements_from_rows(measurements)
    prediction_map, template = _fidelity_prediction_map(config, base_controls)
    train_case_ids, holdout_case_id_set = _split_case_ids_for_holdout(
        rows,
        holdout_stride=holdout_stride,
        holdout_case_ids=holdout_case_ids,
    )
    train_rows = [row for row in rows if str(row.get("caseId", "")) in train_case_ids]
    holdout_rows = [row for row in rows if str(row.get("caseId", "")) in holdout_case_id_set]

    correction_fields: dict[str, dict[str, Any]] = {}
    for observed_key, predicted_key in FIDELITY_CORRECTION_FIELD_PAIRS.items():
        samples: list[tuple[float, float]] = []
        for row in train_rows:
            predicted = prediction_map.get((str(row.get("caseId", "")), str(row.get("connector", ""))))
            if predicted is None:
                continue
            observed = _finite_float(_get(row, observed_key))
            predicted_value = _finite_float(predicted.get(predicted_key))
            if observed is not None and predicted_value is not None:
                samples.append((predicted_value, observed))
        correction_fields[observed_key] = {
            "predictedField": predicted_key,
            **_fit_affine_correction(samples),
        }

    train_eval = _evaluate_corrections(train_rows, prediction_map, correction_fields)
    holdout_eval = _evaluate_corrections(holdout_rows, prediction_map, correction_fields)
    holdout_improves = (
        holdout_eval["improvement"] is not None and float(holdout_eval["improvement"]) > 0.0
    )
    return {
        "schema": TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA,
        "model": "two-cell-external-mujoco-proxy-affine-correction-profile",
        "cadReference": template["cadReference"],
        "split": {
            "holdoutStride": holdout_stride,
            "trainCaseCount": len(train_case_ids),
            "holdoutCaseCount": len(holdout_case_id_set),
            "trainCaseIds": sorted(train_case_ids, key=_case_sort_key),
            "holdoutCaseIds": sorted(holdout_case_id_set, key=_case_sort_key),
        },
        "fieldCorrections": correction_fields,
        "train": train_eval,
        "holdout": holdout_eval,
        "summary": {
            "readyForProxyPreview": bool(holdout_improves and holdout_eval["observedScalarCount"] > 0),
            "holdoutImprovesRms": bool(holdout_improves),
            "holdoutImprovement": holdout_eval["improvement"],
            "coversFullFidelityMatrix": len(rows) == len(template["rows"]),
            "physicalAccuracyValidated": False,
            "remainingEvidence": [
                "segmentedCadContactModel",
                "benchCoordinateHoldout",
                "measuredFrictionAndStiffness",
            ],
        },
        "claimLabels": {
            "correction": "affine surrogate from reduced model predictions to MuJoCo proxy observations",
            "validation": "case-id holdout validation against generated external-engine proxy data",
            "physicalAccuracy": "false until independent bench and segmented CAD evidence pass",
        },
    }


def export_two_cell_external_fidelity_correction_profile_json(
    measurements: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
    **kwargs: Any,
) -> str:
    return json.dumps(
        fit_two_cell_external_fidelity_correction_profile(
            measurements,
            config,
            base_controls,
            **kwargs,
        ),
        indent=2,
    )


def write_two_cell_external_fidelity_correction_profile_artifact(
    input_path: str | Path,
    out: str | Path,
    *,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
    holdout_stride: int = 5,
    holdout_case_ids: Iterable[str] | None = None,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    measurements = load_two_cell_fidelity_matrix_measurements(input_path)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    files = {"profile": output_dir / "two_cell_external_fidelity_correction_profile.json"}
    files["profile"].write_text(
        export_two_cell_external_fidelity_correction_profile_json(
            measurements,
            config,
            holdout_stride=holdout_stride,
            holdout_case_ids=holdout_case_ids,
        ),
        encoding="utf-8",
    )
    return files


def _corrected_field_name(observed_key: str) -> str:
    return "corrected" + observed_key[len("observed") :] if observed_key.startswith("observed") else f"corrected{observed_key}"


def _residual_field_name(prefix: str, observed_key: str) -> str:
    suffix = observed_key[len("observed") :] if observed_key.startswith("observed") else observed_key
    return f"{prefix}Residual{suffix}"


def apply_two_cell_external_fidelity_correction_profile(
    profile: dict[str, Any],
    measurements: Iterable[dict[str, Any]] | None = None,
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    prediction_map, template = _fidelity_prediction_map(config, base_controls)
    measurement_rows = (
        two_cell_fidelity_matrix_measurements_from_rows(measurements)
        if measurements is not None
        else []
    )
    measurement_map = {
        (str(row.get("caseId", "")), str(row.get("connector", ""))): row
        for row in measurement_rows
    }
    corrections = {
        key: dict(value)
        for key, value in dict(profile.get("fieldCorrections", {})).items()
        if isinstance(value, dict)
    }
    applied_rows: list[dict[str, Any]] = []
    for template_row in template["rows"]:
        case_id = str(template_row["caseId"])
        connector = str(template_row["connector"])
        measurement = measurement_map.get((case_id, connector), {})
        row: dict[str, Any] = {
            "caseId": case_id,
            "matrixIndex": template_row["matrixIndex"],
            "connector": connector,
            "leftSite": template_row["leftSite"],
            "rightSite": template_row["rightSite"],
            "actuationCase": template_row["actuationCase"],
            "lockMode": template_row["lockMode"],
            "backlash": template_row["backlash"],
            "alphaCommand": template_row["alphaCommand"],
            "zCommand": template_row["zCommand"],
            "gravityForce": template_row["gravityForce"],
            "pinRadius": template_row["pinRadius"],
            "holeRadius": template_row["holeRadius"],
            "pinHoleClearanceMm": template_row["pinHoleClearanceMm"],
            "predictedContactMode": template_row["predictedContactMode"],
            "observedContactMode": measurement.get("observedContactMode", ""),
            "lockHeldObserved": measurement.get("lockHeldObserved", ""),
            "measurementSource": measurement.get("measurementSource", ""),
        }
        for observed_key, predicted_key in FIDELITY_CORRECTION_FIELD_PAIRS.items():
            predicted_value = _finite_float(template_row.get(predicted_key))
            correction = corrections.get(observed_key, {})
            corrected_key = _corrected_field_name(observed_key)
            row[predicted_key] = template_row[predicted_key]
            if predicted_value is None:
                row[corrected_key] = ""
                continue
            corrected_value = float(correction.get("slope", 1.0)) * predicted_value + float(
                correction.get("intercept", 0.0)
            )
            row[corrected_key] = corrected_value
            observed = _finite_float(measurement.get(observed_key))
            if observed is not None:
                row[observed_key] = observed
                row[_residual_field_name("baseline", observed_key)] = observed - predicted_value
                row[_residual_field_name("corrected", observed_key)] = observed - corrected_value
            else:
                row[observed_key] = ""
                row[_residual_field_name("baseline", observed_key)] = ""
                row[_residual_field_name("corrected", observed_key)] = ""
        applied_rows.append(row)

    evaluation = _evaluate_corrections(measurement_rows, prediction_map, corrections) if measurement_rows else {
        "rowCount": 0,
        "caseCount": 0,
        "observedScalarCount": 0,
        "baselineRms": None,
        "correctedRms": None,
        "improvement": None,
    }
    ready = bool(
        evaluation.get("improvement") is not None
        and float(evaluation["improvement"]) > 0.0
        and evaluation.get("observedScalarCount", 0) > 0
    )
    return {
        "schema": TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA,
        "model": "two-cell-external-fidelity-corrected-preview-table",
        "profileSchema": profile.get("schema", ""),
        "cadReference": template["cadReference"],
        "rows": applied_rows,
        "evaluation": evaluation,
        "summary": {
            "rowCount": len(applied_rows),
            "measurementRowCount": len(measurement_rows),
            "correctedFieldCount": len(corrections),
            "readyForProxyPreview": ready,
            "physicalAccuracyValidated": False,
            "remainingEvidence": [
                "segmentedCadContactModel",
                "benchCoordinateHoldout",
                "measuredFrictionAndStiffness",
            ],
        },
        "claimLabels": {
            "correctedRows": "MuJoCo-proxy corrected preview predictions",
            "use": "can seed simulator/UI previews but not hardware execution",
            "physicalAccuracy": "false until independent bench and segmented CAD evidence pass",
        },
    }


def export_two_cell_external_fidelity_correction_application_json(
    profile: dict[str, Any],
    measurements: Iterable[dict[str, Any]] | None = None,
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(
        apply_two_cell_external_fidelity_correction_profile(
            profile,
            measurements,
            config,
            base_controls,
        ),
        indent=2,
    )


def export_two_cell_external_fidelity_correction_application_csv(application: dict[str, Any]) -> str:
    rows = [dict(row) for row in application.get("rows", []) if isinstance(row, dict)]
    fieldnames = list(rows[0]) if rows else []
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def write_two_cell_external_fidelity_correction_application_artifacts(
    profile_path: str | Path,
    out: str | Path,
    *,
    measurements_path: str | Path | None = None,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    profile = _load_json_file(profile_path)
    measurements = load_two_cell_fidelity_matrix_measurements(measurements_path) if measurements_path is not None else None
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    application = apply_two_cell_external_fidelity_correction_profile(profile, measurements, config)
    files = {
        "application": output_dir / "two_cell_external_fidelity_correction_application.json",
        "application_csv": output_dir / "two_cell_external_fidelity_correction_application.csv",
    }
    files["application"].write_text(json.dumps(application, indent=2), encoding="utf-8")
    files["application_csv"].write_text(
        export_two_cell_external_fidelity_correction_application_csv(application),
        encoding="utf-8",
    )
    return files


def export_two_cell_parameter_calibration_json(
    measurements: Iterable[TwoCellMeasurement],
    config: LatticeConfig | None = None,
    base_controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(calibrate_two_cell_parameters(measurements, config, base_controls), indent=2)
