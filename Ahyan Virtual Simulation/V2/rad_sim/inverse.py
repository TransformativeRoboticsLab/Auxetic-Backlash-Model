from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable

import numpy as np

try:
    from scipy.optimize import lsq_linear
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only in minimal runtimes.
    if exc.name not in {"scipy", "scipy.optimize"}:
        raise
    _INVERSE_SCIPY_IMPORT_ERROR = exc

    def lsq_linear(*args, **kwargs):
        raise ImportError(
            "SciPy is required for inverse design. Install project dependencies "
            "before using solve_inverse_design."
        ) from _INVERSE_SCIPY_IMPORT_ERROR

from .experiments import (
    ResponseMatrix,
    SourceCommand,
    build_response_matrix,
    characterize_response,
)
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .operators import lattice_topology_diagnostic
try:
    from .spring_hinge import solve_spring_hinge_3d
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only in minimal runtimes.
    if exc.name not in {"scipy", "scipy.optimize"}:
        raise
    _SPRING_HINGE_IMPORT_ERROR = exc

    def solve_spring_hinge_3d(*args, **kwargs):
        raise ImportError(
            "SciPy is required for spring-hinge mechanics. Install project dependencies "
            "before using solve_spring_hinge_3d."
        ) from _SPRING_HINGE_IMPORT_ERROR


@dataclass(frozen=True)
class InverseDesignResult:
    matrix: ResponseMatrix
    coefficients: np.ndarray
    commands: tuple[SourceCommand, ...]
    state: LatticeState
    baseline: SimulationResult
    result: SimulationResult
    rms_error_before: float
    rms_error_after: float
    max_abs_error_before: float
    max_abs_error_after: float
    linearized_rms_error: float
    success: bool
    message: str
    target_alpha: np.ndarray | None = None
    target_height: np.ndarray | None = None
    alpha_residual: np.ndarray | None = None
    height_residual: np.ndarray | None = None
    reachable_alpha_mask: np.ndarray | None = None
    reachable_height_mask: np.ndarray | None = None
    underactuated_alpha_mask: np.ndarray | None = None
    underactuated_height_mask: np.ndarray | None = None
    positive_height_reachable_mask: np.ndarray | None = None
    negative_height_reachable_mask: np.ndarray | None = None
    positive_height_underactuated_mask: np.ndarray | None = None
    negative_height_underactuated_mask: np.ndarray | None = None
    topology_reachable_mask: np.ndarray | None = None
    topology_blocked_alpha_mask: np.ndarray | None = None
    topology_blocked_height_mask: np.ndarray | None = None
    max_command_multiplier: float = 0.0
    saturated_column_count: int = 0
    near_saturated_column_count: int = 0

    @property
    def active_actuator_count(self) -> int:
        return len(self.commands)

    @property
    def alpha_underactuated_cells(self) -> int:
        if self.underactuated_alpha_mask is None:
            return 0
        return int(np.count_nonzero(self.underactuated_alpha_mask))

    @property
    def height_underactuated_cells(self) -> int:
        if self.underactuated_height_mask is None:
            return 0
        return int(np.count_nonzero(self.underactuated_height_mask))

    @property
    def positive_height_reachable_cells(self) -> int:
        if self.positive_height_reachable_mask is None:
            return 0
        return int(np.count_nonzero(self.positive_height_reachable_mask))

    @property
    def negative_height_reachable_cells(self) -> int:
        if self.negative_height_reachable_mask is None:
            return 0
        return int(np.count_nonzero(self.negative_height_reachable_mask))

    @property
    def positive_height_underactuated_cells(self) -> int:
        if self.positive_height_underactuated_mask is None:
            return 0
        return int(np.count_nonzero(self.positive_height_underactuated_mask))

    @property
    def negative_height_underactuated_cells(self) -> int:
        if self.negative_height_underactuated_mask is None:
            return 0
        return int(np.count_nonzero(self.negative_height_underactuated_mask))

    @property
    def topology_blocked_alpha_cells(self) -> int:
        if self.topology_blocked_alpha_mask is None:
            return 0
        return int(np.count_nonzero(self.topology_blocked_alpha_mask))

    @property
    def topology_blocked_height_cells(self) -> int:
        if self.topology_blocked_height_mask is None:
            return 0
        return int(np.count_nonzero(self.topology_blocked_height_mask))

    @property
    def saturated_column_fraction(self) -> float:
        if self.coefficients.size == 0:
            return 0.0
        return float(self.saturated_column_count / self.coefficients.size)

    @property
    def near_saturated_column_fraction(self) -> float:
        if self.coefficients.size == 0:
            return 0.0
        return float(self.near_saturated_column_count / self.coefficients.size)

    @property
    def alpha_rms_residual(self) -> float:
        return 0.0 if self.alpha_residual is None else _rms(self.alpha_residual.reshape(-1))

    @property
    def height_rms_residual(self) -> float:
        return 0.0 if self.height_residual is None else _rms(self.height_residual.reshape(-1))

    @property
    def max_abs_alpha_residual(self) -> float:
        if self.alpha_residual is None or self.alpha_residual.size == 0:
            return 0.0
        return float(np.max(np.abs(self.alpha_residual)))

    @property
    def max_abs_height_residual(self) -> float:
        if self.height_residual is None or self.height_residual.size == 0:
            return 0.0
        return float(np.max(np.abs(self.height_residual)))

    def to_dict(
        self,
        *,
        tolerance: float = 1e-9,
        include_response_matrix: bool = True,
    ) -> dict[str, object]:
        return {
            "schema": "rad-sim.inverse-design-result.v1",
            "strategy": "bounded-linearized-response-fit",
            "grid": {"rows": self.matrix.cell_shape[0], "cols": self.matrix.cell_shape[1]},
            "success": self.success,
            "message": self.message,
            "commands": [
                {
                    "row": command.cell[0],
                    "col": command.cell[1],
                    "alpha": command.alpha,
                    "z": command.z,
                }
                for command in self.commands
            ],
            "coefficients": self.coefficients.tolist(),
            "metrics": {
                "activeActuatorCount": self.active_actuator_count,
                "rmsErrorBefore": self.rms_error_before,
                "rmsErrorAfter": self.rms_error_after,
                "maxAbsErrorBefore": self.max_abs_error_before,
                "maxAbsErrorAfter": self.max_abs_error_after,
                "linearizedRmsError": self.linearized_rms_error,
                "alphaRmsResidual": self.alpha_rms_residual,
                "heightRmsResidual": self.height_rms_residual,
                "maxAbsAlphaResidual": self.max_abs_alpha_residual,
                "maxAbsHeightResidual": self.max_abs_height_residual,
                "alphaUnderactuatedCells": self.alpha_underactuated_cells,
                "heightUnderactuatedCells": self.height_underactuated_cells,
                "positiveHeightReachableCells": self.positive_height_reachable_cells,
                "negativeHeightReachableCells": self.negative_height_reachable_cells,
                "positiveHeightUnderactuatedCells": self.positive_height_underactuated_cells,
                "negativeHeightUnderactuatedCells": self.negative_height_underactuated_cells,
                "topologyBlockedAlphaCells": self.topology_blocked_alpha_cells,
                "topologyBlockedHeightCells": self.topology_blocked_height_cells,
                "saturatedColumnCount": self.saturated_column_count,
                "nearSaturatedColumnCount": self.near_saturated_column_count,
                "saturatedColumnFraction": self.saturated_column_fraction,
                "nearSaturatedColumnFraction": self.near_saturated_column_fraction,
                "maxCommandMultiplier": self.max_command_multiplier,
            },
            "target": {
                "alpha": _array_or_none(self.target_alpha),
                "height": _array_or_none(self.target_height),
            },
            "residual": {
                "alpha": _array_or_none(self.alpha_residual),
                "height": _array_or_none(self.height_residual),
            },
            "reachability": {
                "alpha": _array_or_none(self.reachable_alpha_mask),
                "height": _array_or_none(self.reachable_height_mask),
                "underactuatedAlpha": _array_or_none(self.underactuated_alpha_mask),
                "underactuatedHeight": _array_or_none(self.underactuated_height_mask),
                "positiveHeight": _array_or_none(self.positive_height_reachable_mask),
                "negativeHeight": _array_or_none(self.negative_height_reachable_mask),
                "underactuatedPositiveHeight": _array_or_none(self.positive_height_underactuated_mask),
                "underactuatedNegativeHeight": _array_or_none(self.negative_height_underactuated_mask),
                "topological": _array_or_none(self.topology_reachable_mask),
                "topologyBlockedAlpha": _array_or_none(self.topology_blocked_alpha_mask),
                "topologyBlockedHeight": _array_or_none(self.topology_blocked_height_mask),
            },
            "responseMatrix": (
                self.matrix.to_dict(tolerance=tolerance)
                if include_response_matrix
                else {"schema": "rad-sim.response-matrix.v1", "diagnostics": self.matrix.to_dict(tolerance=tolerance)["diagnostics"]}
            ),
        }


@dataclass(frozen=True)
class InversePhysicalValidation:
    inverse: InverseDesignResult
    physical_baseline: SimulationResult
    physical_result: SimulationResult
    height_residual_before: np.ndarray | None
    height_residual_after: np.ndarray | None
    alpha_residual_after: np.ndarray | None
    height_model_error: np.ndarray
    center_model_error: np.ndarray
    physical_rms_height_error_before: float
    physical_rms_height_error_after: float
    physical_max_abs_height_error_before: float
    physical_max_abs_height_error_after: float
    height_rms_model_error: float
    center_rms_model_error: float
    max_abs_height_model_error: float
    max_abs_center_model_error: float
    physical_energy: float
    physical_success: bool

    @property
    def physical_height_error_improvement(self) -> float:
        return self.physical_rms_height_error_before - self.physical_rms_height_error_after

    @property
    def model_agreement_score(self) -> float:
        return 1.0 / (1.0 + self.center_rms_model_error)

    def to_dict(self) -> dict[str, object]:
        center_norm = np.linalg.norm(self.center_model_error, axis=2)
        return {
            "schema": "rad-sim.inverse-physical-validation.v1",
            "physicalSuccess": self.physical_success,
            "metrics": {
                "physicalRmsHeightErrorBefore": self.physical_rms_height_error_before,
                "physicalRmsHeightErrorAfter": self.physical_rms_height_error_after,
                "physicalHeightErrorImprovement": self.physical_height_error_improvement,
                "physicalMaxAbsHeightErrorBefore": self.physical_max_abs_height_error_before,
                "physicalMaxAbsHeightErrorAfter": self.physical_max_abs_height_error_after,
                "heightRmsModelError": self.height_rms_model_error,
                "centerRmsModelError": self.center_rms_model_error,
                "maxAbsHeightModelError": self.max_abs_height_model_error,
                "maxAbsCenterModelError": self.max_abs_center_model_error,
                "modelAgreementScore": self.model_agreement_score,
                "physicalEnergy": self.physical_energy,
            },
            "fields": {
                "heightResidualBefore": _array_or_none(self.height_residual_before),
                "heightResidualAfter": _array_or_none(self.height_residual_after),
                "alphaResidualAfter": _array_or_none(self.alpha_residual_after),
                "heightModelError": _array_or_none(self.height_model_error),
                "centerModelErrorNorm": center_norm.tolist(),
            },
        }


def _target_array(
    value: np.ndarray | Iterable[Iterable[float]] | None,
    shape: tuple[int, int],
    name: str,
) -> np.ndarray | None:
    if value is None:
        return None
    target = np.asarray(value, dtype=float)
    if target.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    return target


def _array_or_none(value: np.ndarray | None) -> list | None:
    return None if value is None else value.tolist()


def _baseline_state(
    config: LatticeConfig,
    locked_cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]] = (),
) -> LatticeState:
    state = LatticeState.uniform(config)
    removed = {(int(r), int(c)) for r, c in removed_cells}
    for r, c in removed:
        state.removed_mask[r, c] = True
    for r, c in locked_cells:
        if (int(r), int(c)) in removed:
            continue
        state.locked_mask[r, c] = True
    return state


def _unlocked_cells(
    config: LatticeConfig,
    actuator_cells: Iterable[tuple[int, int]] | None,
    locked_cells: tuple[tuple[int, int], ...],
    removed_cells: tuple[tuple[int, int], ...] = (),
) -> tuple[tuple[int, int], ...]:
    locked = set(locked_cells)
    removed = set(removed_cells)
    cells = (
        actuator_cells
        if actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    return tuple(
        (int(r), int(c))
        for r, c in cells
        if (int(r), int(c)) not in locked and (int(r), int(c)) not in removed
    )


def _topology_reachable_mask(
    config: LatticeConfig,
    actuator_cells: tuple[tuple[int, int], ...],
    removed_cells: tuple[tuple[int, int], ...],
) -> np.ndarray:
    state = _baseline_state(config, (), removed_cells)
    topology = lattice_topology_diagnostic(config, state)
    labels = topology["component_labels"]
    active_components = {
        int(labels[r, c])
        for r, c in actuator_cells
        if 0 <= r < config.rows and 0 <= c < config.cols and int(labels[r, c]) >= 0
    }
    return np.vectorize(lambda label: int(label) in active_components)(labels)


def _topology_blocked_mask(
    target: np.ndarray | None,
    baseline: np.ndarray,
    topology_reachable: np.ndarray,
    tolerance: float,
) -> np.ndarray | None:
    if target is None:
        return None
    requested_change = np.abs(target - baseline) > tolerance
    return requested_change & ~topology_reachable


def _stack_system(
    matrix: ResponseMatrix,
    baseline: SimulationResult,
    target_alpha: np.ndarray | None,
    target_height: np.ndarray | None,
    alpha_weight: float,
    height_weight: float,
) -> tuple[np.ndarray, np.ndarray]:
    a_blocks: list[np.ndarray] = []
    b_blocks: list[np.ndarray] = []
    if target_alpha is not None and alpha_weight > 0:
        a_blocks.append(matrix.alpha * alpha_weight)
        b_blocks.append((target_alpha - baseline.alpha).reshape(-1) * alpha_weight)
    if target_height is not None and height_weight > 0:
        a_blocks.append(matrix.height * height_weight)
        b_blocks.append(
            (target_height - baseline.metadata["height"]).reshape(-1) * height_weight
        )
    if not a_blocks:
        raise ValueError("at least one positive-weight target is required")
    return np.vstack(a_blocks), np.concatenate(b_blocks)


def _error_vector(
    result: SimulationResult,
    target_alpha: np.ndarray | None,
    target_height: np.ndarray | None,
    alpha_weight: float,
    height_weight: float,
) -> np.ndarray:
    blocks = []
    if target_alpha is not None and alpha_weight > 0:
        blocks.append((target_alpha - result.alpha).reshape(-1) * alpha_weight)
    if target_height is not None and height_weight > 0:
        blocks.append((target_height - result.metadata["height"]).reshape(-1) * height_weight)
    if not blocks:
        return np.zeros(0, dtype=float)
    return np.concatenate(blocks)


def _rms(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(values**2)))


def _max_abs(values: np.ndarray | None) -> float:
    if values is None or values.size == 0:
        return 0.0
    return float(np.max(np.abs(values)))


def _profile_update_map(profile: dict[str, object]) -> dict[str, dict[str, object]]:
    updates = profile.get("recommendedUpdates", [])
    if not isinstance(updates, list):
        return {}
    mapped: dict[str, dict[str, object]] = {}
    for update in updates:
        if not isinstance(update, dict):
            continue
        name = update.get("name")
        if isinstance(name, str):
            mapped[name] = update
    return mapped


def _profile_update_value(
    updates: dict[str, dict[str, object]],
    name: str,
) -> float | None:
    value = updates.get(name, {}).get("proposed")
    try:
        numeric = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return numeric if np.isfinite(numeric) else None


def _safe_profile_update_count(updates: dict[str, dict[str, object]]) -> int:
    return sum(1 for update in updates.values() if update.get("safeToApply") is True)


def _profile_scale(value: float | None, uncertainty: float | None, tolerance: float) -> float:
    candidates = [abs(float(tolerance)), 1e-12]
    if value is not None:
        candidates.append(abs(float(value)))
    if uncertainty is not None:
        candidates.append(abs(float(uncertainty)))
    return max(candidates)


def _target_mask(
    target: np.ndarray | None,
    baseline: np.ndarray,
    tolerance: float,
) -> np.ndarray:
    if target is None:
        return np.zeros_like(baseline, dtype=bool)
    return np.abs(target - baseline) > tolerance


def _band_failure_count(
    residual: np.ndarray | None,
    band: float,
    tolerance: float,
) -> int:
    if residual is None:
        return 0
    return int(np.count_nonzero(np.abs(residual) > band + tolerance))


def _csv_scalar(value: object) -> str:
    text = str(value)
    if any(mark in text for mark in (",", "\n", '"')):
        return '"' + text.replace('"', '""') + '"'
    return text


def _dict_section(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _report_float(value: object) -> float | None:
    try:
        numeric = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return numeric if np.isfinite(numeric) else None


def _report_int(value: object) -> int:
    numeric = _report_float(value)
    return int(numeric) if numeric is not None else 0


def _aggregate_commands(
    commands: tuple[SourceCommand, ...],
    coefficients: np.ndarray,
    tolerance: float,
) -> tuple[SourceCommand, ...]:
    by_cell: dict[tuple[int, int], list[float]] = {}
    for command, coefficient in zip(commands, coefficients, strict=True):
        alpha = command.alpha * float(coefficient)
        z = command.z * float(coefficient)
        if abs(alpha) <= tolerance and abs(z) <= tolerance:
            continue
        values = by_cell.setdefault(command.cell, [0.0, 0.0])
        values[0] += alpha
        values[1] += z
    return tuple(
        SourceCommand(cell=cell, alpha=values[0], z=values[1])
        for cell, values in sorted(by_cell.items())
        if abs(values[0]) > tolerance or abs(values[1]) > tolerance
    )


def _state_with_commands(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]] = (),
) -> LatticeState:
    removed = {(int(r), int(c)) for r, c in removed_cells}
    state = _baseline_state(config, locked_cells, removed)
    for command in commands:
        r, c = command.cell
        if (int(r), int(c)) in removed:
            continue
        state.actuator_grid[r, c] += command.alpha
        state.z_actuator_grid[r, c] += command.z
    return state


def _reachable_mask(matrix: ResponseMatrix, family: str, tolerance: float) -> np.ndarray:
    response = matrix.alpha if family == "alpha" else matrix.height
    if response.shape[1] == 0:
        flat = np.zeros(response.shape[0], dtype=bool)
    else:
        flat = np.any(np.abs(response) > tolerance, axis=1)
    return flat.reshape(matrix.cell_shape)


def _signed_height_reachable_masks(
    config: LatticeConfig,
    actuator_cells: tuple[tuple[int, int], ...],
    locked_cells: tuple[tuple[int, int], ...],
    removed_cells: tuple[tuple[int, int], ...],
    *,
    alpha_step: float,
    z_step: float,
    include_alpha: bool,
    include_z: bool,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray]:
    positive = np.zeros(config.rows * config.cols, dtype=bool)
    negative = np.zeros(config.rows * config.cols, dtype=bool)
    alpha_probe = abs(float(alpha_step))
    z_probe = abs(float(z_step))
    for cell in actuator_cells:
        commands: list[SourceCommand] = []
        if include_alpha and alpha_probe > tolerance:
            commands.extend(
                (
                    SourceCommand(cell=cell, alpha=alpha_probe),
                    SourceCommand(cell=cell, alpha=-alpha_probe),
                )
            )
        if include_z and z_probe > tolerance:
            commands.extend(
                (
                    SourceCommand(cell=cell, z=z_probe),
                    SourceCommand(cell=cell, z=-z_probe),
                )
            )
        for command in commands:
            response = characterize_response(
                config,
                (command,),
                locked_cells,
                tolerance,
                removed_cells=removed_cells,
            )
            delta = response.height_delta.reshape(-1)
            positive |= delta > tolerance
            negative |= delta < -tolerance
    return positive.reshape(config.rows, config.cols), negative.reshape(config.rows, config.cols)


def _underactuated_mask(
    target: np.ndarray | None,
    baseline: np.ndarray,
    reachable: np.ndarray,
    tolerance: float,
) -> np.ndarray | None:
    if target is None:
        return None
    requested_change = np.abs(target - baseline) > tolerance
    return requested_change & ~reachable


def _signed_height_underactuated_masks(
    target: np.ndarray | None,
    baseline: np.ndarray,
    positive_reachable: np.ndarray,
    negative_reachable: np.ndarray,
    tolerance: float,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    if target is None:
        return None, None
    residual = target - baseline
    positive_requested = residual > tolerance
    negative_requested = residual < -tolerance
    return positive_requested & ~positive_reachable, negative_requested & ~negative_reachable


def _saturation_counts(
    coefficients: np.ndarray,
    max_command_multiplier: float,
    tolerance: float,
) -> tuple[int, int]:
    if coefficients.size == 0:
        return 0, 0
    margin = max(tolerance, 1e-6)
    magnitudes = np.abs(coefficients)
    saturated = magnitudes >= max_command_multiplier - margin
    near_saturated = magnitudes >= 0.9 * max_command_multiplier
    return int(np.count_nonzero(saturated)), int(np.count_nonzero(near_saturated))


def solve_inverse_design(
    config: LatticeConfig,
    *,
    target_height: np.ndarray | Iterable[Iterable[float]] | None = None,
    target_alpha: np.ndarray | Iterable[Iterable[float]] | None = None,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    locked_cells: Iterable[tuple[int, int]] = (),
    removed_cells: Iterable[tuple[int, int]] = (),
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    alpha_weight: float = 0.35,
    height_weight: float = 1.0,
    regularization: float = 1e-4,
    max_command_multiplier: float = 4.0,
    tolerance: float = 1e-9,
) -> InverseDesignResult:
    """Fit actuator commands to target alpha/height fields with a linearized model.

    The columns are finite command responses from ``build_response_matrix``.
    Coefficients scale those command columns, so this is a local inverse-design
    approximation around the unactuated baseline, not a full nonlinear optimizer.
    """
    shape = (config.rows, config.cols)
    target_alpha_array = _target_array(target_alpha, shape, "target_alpha")
    target_height_array = _target_array(target_height, shape, "target_height")
    if target_alpha_array is None and target_height_array is None:
        raise ValueError("target_alpha or target_height is required")
    if regularization < 0:
        raise ValueError("regularization must be non-negative")
    if max_command_multiplier <= 0:
        raise ValueError("max_command_multiplier must be positive")

    locked = tuple((int(r), int(c)) for r, c in locked_cells)
    removed = tuple((int(r), int(c)) for r, c in removed_cells)
    candidates = _unlocked_cells(config, actuator_cells, locked, removed)
    matrix = build_response_matrix(
        config,
        actuator_cells=candidates,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=locked,
        removed_cells=removed,
        tolerance=tolerance,
    )
    baseline = simulate_kinematic(config, _baseline_state(config, locked, removed))
    a_matrix, b_vector = _stack_system(
        matrix,
        baseline,
        target_alpha_array,
        target_height_array,
        alpha_weight,
        height_weight,
    )
    if regularization > 0:
        damp = np.sqrt(regularization) * np.eye(a_matrix.shape[1])
        a_solve = np.vstack([a_matrix, damp])
        b_solve = np.concatenate([b_vector, np.zeros(a_matrix.shape[1], dtype=float)])
    else:
        a_solve = a_matrix
        b_solve = b_vector

    if a_matrix.shape[1] == 0:
        coefficients = np.zeros(0, dtype=float)
        success = False
        message = "no unlocked actuator candidates"
    else:
        fit = lsq_linear(
            a_solve,
            b_solve,
            bounds=(-max_command_multiplier, max_command_multiplier),
            lsmr_tol="auto",
        )
        coefficients = np.asarray(fit.x, dtype=float)
        success = bool(fit.success)
        message = str(fit.message)

    commands = _aggregate_commands(matrix.commands, coefficients, tolerance)
    state = _state_with_commands(config, commands, locked, removed)
    result = simulate_kinematic(config, state)
    before = _error_vector(
        baseline, target_alpha_array, target_height_array, alpha_weight, height_weight
    )
    after = _error_vector(
        result, target_alpha_array, target_height_array, alpha_weight, height_weight
    )
    linearized = a_matrix @ coefficients - b_vector if coefficients.size else -b_vector
    alpha_residual = (
        None if target_alpha_array is None else target_alpha_array - result.alpha
    )
    height_residual = (
        None if target_height_array is None else target_height_array - result.metadata["height"]
    )
    reachable_alpha = _reachable_mask(matrix, "alpha", tolerance)
    reachable_height = _reachable_mask(matrix, "height", tolerance)
    positive_height_reachable, negative_height_reachable = _signed_height_reachable_masks(
        config,
        candidates,
        locked,
        removed,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        tolerance=tolerance,
    )
    (
        positive_height_underactuated,
        negative_height_underactuated,
    ) = _signed_height_underactuated_masks(
        target_height_array,
        baseline.metadata["height"],
        positive_height_reachable,
        negative_height_reachable,
        tolerance,
    )
    saturated_count, near_saturated_count = _saturation_counts(
        coefficients,
        max_command_multiplier,
        tolerance,
    )
    topology_reachable = _topology_reachable_mask(config, candidates, removed)
    return InverseDesignResult(
        matrix=matrix,
        coefficients=coefficients,
        commands=commands,
        state=state,
        baseline=baseline,
        result=result,
        rms_error_before=_rms(before),
        rms_error_after=_rms(after),
        max_abs_error_before=float(np.max(np.abs(before))) if before.size else 0.0,
        max_abs_error_after=float(np.max(np.abs(after))) if after.size else 0.0,
        linearized_rms_error=_rms(linearized),
        success=success,
        message=message,
        target_alpha=target_alpha_array,
        target_height=target_height_array,
        alpha_residual=alpha_residual,
        height_residual=height_residual,
        reachable_alpha_mask=reachable_alpha,
        reachable_height_mask=reachable_height,
        underactuated_alpha_mask=_underactuated_mask(
            target_alpha_array,
            baseline.alpha,
            reachable_alpha,
            tolerance,
        ),
        underactuated_height_mask=_underactuated_mask(
            target_height_array,
            baseline.metadata["height"],
            reachable_height,
            tolerance,
        ),
        positive_height_reachable_mask=positive_height_reachable,
        negative_height_reachable_mask=negative_height_reachable,
        positive_height_underactuated_mask=positive_height_underactuated,
        negative_height_underactuated_mask=negative_height_underactuated,
        topology_reachable_mask=topology_reachable,
        topology_blocked_alpha_mask=_topology_blocked_mask(
            target_alpha_array,
            baseline.alpha,
            topology_reachable,
            tolerance,
        ),
        topology_blocked_height_mask=_topology_blocked_mask(
            target_height_array,
            baseline.metadata["height"],
            topology_reachable,
            tolerance,
        ),
        max_command_multiplier=max_command_multiplier,
        saturated_column_count=saturated_count,
        near_saturated_column_count=near_saturated_count,
    )


def validate_inverse_design_physical(
    inverse: InverseDesignResult,
    load_case: LoadCase | None = None,
) -> InversePhysicalValidation:
    """Check an inverse design with the 3D spring-hinge relaxation model.

    This is an explicit validation pass for the kinematic inverse proposal. It
    does not re-optimize commands; it measures how the same command state behaves
    under the spring-hinge preview and reports target residuals plus disagreement
    from the kinematic solution.
    """

    config = inverse.result.config
    physical_baseline = solve_spring_hinge_3d(config, inverse.baseline.state, load_case)
    physical_result = solve_spring_hinge_3d(config, inverse.state, load_case)
    physical_height = physical_result.metadata["height"]
    baseline_height = physical_baseline.metadata["height"]
    kinematic_height = inverse.result.metadata["height"]
    height_residual_before = (
        None if inverse.target_height is None else inverse.target_height - baseline_height
    )
    height_residual_after = (
        None if inverse.target_height is None else inverse.target_height - physical_height
    )
    alpha_residual_after = (
        None if inverse.target_alpha is None else inverse.target_alpha - physical_result.alpha
    )
    height_model_error = physical_height - kinematic_height
    center_model_error = physical_result.deformed_centers_3d - inverse.result.deformed_centers_3d
    center_norm = np.linalg.norm(center_model_error, axis=2)
    return InversePhysicalValidation(
        inverse=inverse,
        physical_baseline=physical_baseline,
        physical_result=physical_result,
        height_residual_before=height_residual_before,
        height_residual_after=height_residual_after,
        alpha_residual_after=alpha_residual_after,
        height_model_error=height_model_error,
        center_model_error=center_model_error,
        physical_rms_height_error_before=(
            0.0 if height_residual_before is None else _rms(height_residual_before.reshape(-1))
        ),
        physical_rms_height_error_after=(
            0.0 if height_residual_after is None else _rms(height_residual_after.reshape(-1))
        ),
        physical_max_abs_height_error_before=_max_abs(height_residual_before),
        physical_max_abs_height_error_after=_max_abs(height_residual_after),
        height_rms_model_error=_rms(height_model_error.reshape(-1)),
        center_rms_model_error=_rms(center_norm.reshape(-1)),
        max_abs_height_model_error=_max_abs(height_model_error),
        max_abs_center_model_error=_max_abs(center_norm),
        physical_energy=float(physical_result.metadata.get("energy", np.nan)),
        physical_success=bool(
            physical_baseline.metadata.get("success", False)
            and physical_result.metadata.get("success", False)
        ),
    )


def inverse_design_report(
    inverse: InverseDesignResult,
    validation: InversePhysicalValidation | None = None,
    *,
    tolerance: float = 1e-9,
    include_response_matrix: bool = True,
) -> dict[str, object]:
    return {
        "schema": "rad-sim.inverse-design-report.v1",
        "inverse": inverse.to_dict(
            tolerance=tolerance,
            include_response_matrix=include_response_matrix,
        ),
        "physicalValidation": None if validation is None else validation.to_dict(),
        "assumptions": {
            "inverseModel": "linearized finite response columns around the unactuated baseline",
            "physicalValidation": "spring-hinge pass validates the proposed commands but does not re-optimize them",
        },
    }


def export_inverse_design_report_json(
    inverse: InverseDesignResult,
    validation: InversePhysicalValidation | None = None,
    *,
    tolerance: float = 1e-9,
    include_response_matrix: bool = True,
) -> str:
    return json.dumps(
        inverse_design_report(
            inverse,
            validation,
            tolerance=tolerance,
            include_response_matrix=include_response_matrix,
        ),
        indent=2,
    )


def reachable_equilibrium_profile_inverse_report(
    inverse: InverseDesignResult,
    empirical_profile: dict[str, object],
    validation: InversePhysicalValidation | None = None,
    *,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Score an inverse solve against a reachable-equilibrium empirical profile.

    The empirical profile is calibration metadata. This report deliberately does
    not mutate ``LatticeConfig`` or re-run optimization; it records how the
    existing inverse residuals compare with bounded profile scales.
    """

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    updates = _profile_update_map(empirical_profile)
    alpha_scale = _profile_update_value(updates, "alphaResponseScale")
    height_scale = _profile_update_value(updates, "heightResponseScale")
    topology_tolerance = _profile_update_value(updates, "topologyLeakageTolerance")
    group_sequence_tolerance = _profile_update_value(updates, "groupSequenceTolerance")
    uncertainty_budget = _profile_update_value(updates, "uncertaintyBudget")
    alpha_band = _profile_scale(alpha_scale, uncertainty_budget, tolerance)
    height_band = _profile_scale(height_scale, uncertainty_budget, tolerance)
    if topology_tolerance is not None:
        height_band = max(height_band, abs(float(topology_tolerance)))

    baseline_alpha = inverse.baseline.alpha
    baseline_height = inverse.baseline.metadata["height"]
    alpha_target_mask = _target_mask(inverse.target_alpha, baseline_alpha, tolerance)
    height_target_mask = _target_mask(inverse.target_height, baseline_height, tolerance)
    target_mask = alpha_target_mask | height_target_mask
    target_cell_count = int(np.count_nonzero(target_mask))
    alpha_rms = inverse.alpha_rms_residual
    height_rms = inverse.height_rms_residual
    alpha_max = inverse.max_abs_alpha_residual
    height_max = inverse.max_abs_height_residual
    normalized_alpha_rms = alpha_rms / alpha_band
    normalized_height_rms = height_rms / height_band
    normalized_alpha_max = alpha_max / alpha_band
    normalized_height_max = height_max / height_band
    alpha_failures = _band_failure_count(inverse.alpha_residual, alpha_band, tolerance)
    height_failures = _band_failure_count(inverse.height_residual, height_band, tolerance)

    profile_summary = empirical_profile.get("summary", {})
    if not isinstance(profile_summary, dict):
        profile_summary = {}
    empirical_ready = bool(profile_summary.get("empiricalProfileReady"))
    safe_count = _safe_profile_update_count(updates)
    finite_parameters = {
        "alphaResponseScale": alpha_scale,
        "heightResponseScale": height_scale,
        "topologyLeakageTolerance": topology_tolerance,
        "groupSequenceTolerance": group_sequence_tolerance,
        "uncertaintyBudget": uncertainty_budget,
    }
    missing_evidence: list[str] = []
    if not empirical_ready:
        missing_evidence.append("empiricalProfileReady")
    if safe_count <= 0:
        missing_evidence.append("safeBoundedProposals")
    for name, value in finite_parameters.items():
        if value is None:
            missing_evidence.append(name)
    if inverse.target_alpha is None and inverse.target_height is None:
        missing_evidence.append("inverseTarget")
    if target_cell_count <= 0:
        missing_evidence.append("targetCells")
    if not inverse.success:
        missing_evidence.append("inverseSolveSuccess")
    profile_ready = len(missing_evidence) == 0
    score_terms = [
        normalized_alpha_rms if inverse.target_alpha is not None else None,
        normalized_height_rms if inverse.target_height is not None else None,
        float(alpha_failures + height_failures),
    ]
    finite_score_terms = [float(value) for value in score_terms if value is not None]
    weighted_score = (
        float(np.sqrt(np.mean(np.square(finite_score_terms))))
        if finite_score_terms
        else 0.0
    )

    return {
        "schema": "rad-sim.reachable-equilibrium-profile-inverse.v1",
        "method": (
            "profile-aware diagnostic scoring of an already solved inverse plan; "
            "the empirical profile is read-only calibration metadata"
        ),
        "profile": {
            "schema": empirical_profile.get("schema", ""),
            "empiricalProfileReady": empirical_ready,
            "safeBoundedProposalCount": safe_count,
            "proposalNames": sorted(updates),
        },
        "parameters": {
            **finite_parameters,
            "alphaResidualBand": alpha_band,
            "heightResidualBand": height_band,
            "tolerance": float(tolerance),
        },
        "target": {
            "targetCellCount": target_cell_count,
            "alphaTargetCellCount": int(np.count_nonzero(alpha_target_mask)),
            "heightTargetCellCount": int(np.count_nonzero(height_target_mask)),
        },
        "inverse": {
            "success": bool(inverse.success),
            "message": inverse.message,
            "activeActuatorCount": inverse.active_actuator_count,
            "rmsErrorBefore": inverse.rms_error_before,
            "rmsErrorAfter": inverse.rms_error_after,
            "maxAbsErrorBefore": inverse.max_abs_error_before,
            "maxAbsErrorAfter": inverse.max_abs_error_after,
            "saturatedColumnCount": inverse.saturated_column_count,
            "nearSaturatedColumnCount": inverse.near_saturated_column_count,
            "alphaUnderactuatedCells": inverse.alpha_underactuated_cells,
            "heightUnderactuatedCells": inverse.height_underactuated_cells,
            "topologyBlockedAlphaCells": inverse.topology_blocked_alpha_cells,
            "topologyBlockedHeightCells": inverse.topology_blocked_height_cells,
        },
        "residuals": {
            "alphaRmsResidual": alpha_rms,
            "heightRmsResidual": height_rms,
            "maxAbsAlphaResidual": alpha_max,
            "maxAbsHeightResidual": height_max,
            "normalizedAlphaRmsResidual": normalized_alpha_rms,
            "normalizedHeightRmsResidual": normalized_height_rms,
            "normalizedMaxAbsAlphaResidual": normalized_alpha_max,
            "normalizedMaxAbsHeightResidual": normalized_height_max,
            "alphaProfileBandFailureCount": alpha_failures,
            "heightProfileBandFailureCount": height_failures,
        },
        "physicalValidation": (
            None
            if validation is None
            else {
                "schema": "rad-sim.inverse-physical-validation.v1",
                "physicalSuccess": validation.physical_success,
                "physicalRmsHeightErrorAfter": validation.physical_rms_height_error_after,
                "physicalMaxAbsHeightErrorAfter": validation.physical_max_abs_height_error_after,
                "heightRmsModelError": validation.height_rms_model_error,
                "centerRmsModelError": validation.center_rms_model_error,
            }
        ),
        "summary": {
            "status": (
                "reachable-equilibrium-profile-inverse-ready"
                if profile_ready
                else "needs-profile-aware-inverse-review"
            ),
            "profileInverseReady": profile_ready,
            "profileWeightedResidualScore": weighted_score,
            "profileBandFailureCount": alpha_failures + height_failures,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_profile_inverse_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumProfileInverseNat",
            "leanPredicate": "reachableEquilibriumProfileInverseReadyNat",
            "schema": "rad-sim.reachable-equilibrium-profile-inverse.v1",
        },
        "claimLabels": {
            "profile": "bench-measured empirical law",
            "inverse": "linearized finite-response inverse plan",
            "diagnostic": "read-only calibration-aware residual certificate",
            "physicalAccuracy": "not a physical-law proof and not an optimizer mutation",
        },
        "limitations": [
            "The report scores residuals against empirical bands but does not prove reachability.",
            "Profile scales are fitted metadata and should be replaced by independent hardware measurements.",
            "Spring-hinge or rigid-body contact validation remains a separate check.",
        ],
    }


def export_reachable_equilibrium_profile_inverse_json(
    inverse: InverseDesignResult,
    empirical_profile: dict[str, object],
    validation: InversePhysicalValidation | None = None,
    *,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        reachable_equilibrium_profile_inverse_report(
            inverse,
            empirical_profile,
            validation,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_reachable_equilibrium_profile_inverse_csv(
    report: dict[str, object],
) -> str:
    summary = report.get("summary", {})
    profile = report.get("profile", {})
    residuals = report.get("residuals", {})
    target = report.get("target", {})
    inverse = report.get("inverse", {})
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(profile, dict):
        profile = {}
    if not isinstance(residuals, dict):
        residuals = {}
    if not isinstance(target, dict):
        target = {}
    if not isinstance(inverse, dict):
        inverse = {}
    rows = [
        [
            "schema",
            "profile_inverse_ready",
            "empirical_profile_ready",
            "target_cells",
            "active_actuators",
            "profile_weighted_residual_score",
            "alpha_rms_residual",
            "height_rms_residual",
            "alpha_band_failures",
            "height_band_failures",
            "missing_evidence",
        ],
        [
            report.get("schema", ""),
            summary.get("profileInverseReady", ""),
            profile.get("empiricalProfileReady", ""),
            target.get("targetCellCount", ""),
            inverse.get("activeActuatorCount", ""),
            summary.get("profileWeightedResidualScore", ""),
            residuals.get("alphaRmsResidual", ""),
            residuals.get("heightRmsResidual", ""),
            residuals.get("alphaProfileBandFailureCount", ""),
            residuals.get("heightProfileBandFailureCount", ""),
            ";".join(str(item) for item in summary.get("missingEvidence", []))
            if isinstance(summary.get("missingEvidence"), list)
            else "",
        ],
    ]
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)


def reachable_equilibrium_profile_inverse_acceptance_report(
    profile_inverse_report: dict[str, object],
    *,
    max_weighted_residual_score: float = 1.0,
    max_band_failures: int = 0,
    max_active_actuators: int = 24,
    allow_underactuated: bool = False,
    require_physical_validation: bool = False,
    score_scale: int = 1000,
) -> dict[str, object]:
    """Gate a profile-aware inverse diagnostic before it is accepted for preview."""

    if max_weighted_residual_score < 0:
        raise ValueError("max_weighted_residual_score must be non-negative")
    if max_band_failures < 0:
        raise ValueError("max_band_failures must be non-negative")
    if max_active_actuators < 0:
        raise ValueError("max_active_actuators must be non-negative")
    if score_scale <= 0:
        raise ValueError("score_scale must be positive")

    summary = _dict_section(profile_inverse_report.get("summary"))
    inverse = _dict_section(profile_inverse_report.get("inverse"))
    residuals = _dict_section(profile_inverse_report.get("residuals"))
    physical = profile_inverse_report.get("physicalValidation")
    physical_section = _dict_section(physical)
    score = _report_float(summary.get("profileWeightedResidualScore"))
    band_failures = _report_int(summary.get("profileBandFailureCount"))
    active_actuators = _report_int(inverse.get("activeActuatorCount"))
    underactuated_cells = sum(
        _report_int(inverse.get(name))
        for name in (
            "alphaUnderactuatedCells",
            "heightUnderactuatedCells",
            "topologyBlockedAlphaCells",
            "topologyBlockedHeightCells",
        )
    )
    profile_inverse_ready = bool(summary.get("profileInverseReady"))
    inverse_success = bool(inverse.get("success"))
    physical_pass = (
        bool(physical_section.get("physicalSuccess"))
        if isinstance(physical, dict)
        else not require_physical_validation
    )

    missing_evidence: list[str] = []
    if profile_inverse_report.get("schema") != "rad-sim.reachable-equilibrium-profile-inverse.v1":
        missing_evidence.append("profileInverseReportSchema")
    if not profile_inverse_ready:
        missing_evidence.append("profileInverseReady")
    if not inverse_success:
        missing_evidence.append("inverseSolveSuccess")
    if score is None:
        missing_evidence.append("profileWeightedResidualScore")
    if require_physical_validation and not isinstance(physical, dict):
        missing_evidence.append("physicalValidation")

    failed_criteria: list[str] = []
    if score is not None and score > max_weighted_residual_score:
        failed_criteria.append("profileWeightedResidualScore")
    if band_failures > max_band_failures:
        failed_criteria.append("profileBandFailureCount")
    if active_actuators > max_active_actuators:
        failed_criteria.append("activeActuatorBudget")
    if underactuated_cells > 0 and not allow_underactuated:
        failed_criteria.append("underactuatedOrTopologyBlockedCells")
    if require_physical_validation and not physical_pass:
        failed_criteria.append("physicalValidation")

    if missing_evidence:
        decision = "reject-missing-evidence"
    elif failed_criteria:
        decision = "review-required"
    else:
        decision = "accept-for-preview"
    accepted = decision == "accept-for-preview"
    score_scaled = int(np.ceil((score or 0.0) * score_scale))
    limit_scaled = int(np.floor(max_weighted_residual_score * score_scale))
    return {
        "schema": "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
        "sourceReportSchema": profile_inverse_report.get("schema", ""),
        "method": (
            "thresholded acceptance gate for a read-only profile-aware inverse "
            "residual certificate"
        ),
        "criteria": {
            "maxWeightedResidualScore": float(max_weighted_residual_score),
            "maxBandFailures": int(max_band_failures),
            "maxActiveActuators": int(max_active_actuators),
            "allowUnderactuated": bool(allow_underactuated),
            "requirePhysicalValidation": bool(require_physical_validation),
            "scoreScale": int(score_scale),
        },
        "metrics": {
            "profileInverseReady": profile_inverse_ready,
            "inverseSolveSuccess": inverse_success,
            "profileWeightedResidualScore": score,
            "profileWeightedResidualScoreScaled": score_scaled,
            "maxWeightedResidualScoreScaled": limit_scaled,
            "profileBandFailureCount": band_failures,
            "activeActuatorCount": active_actuators,
            "underactuatedOrTopologyBlockedCellCount": underactuated_cells,
            "physicalValidationPass": bool(physical_pass),
            "alphaProfileBandFailureCount": _report_int(
                residuals.get("alphaProfileBandFailureCount")
            ),
            "heightProfileBandFailureCount": _report_int(
                residuals.get("heightProfileBandFailureCount")
            ),
        },
        "decision": {
            "decision": decision,
            "acceptedForPreview": accepted,
            "reviewRequired": not accepted,
            "failedCriteria": failed_criteria,
        },
        "summary": {
            "status": (
                "profile-aware-inverse-accepted"
                if accepted
                else "profile-aware-inverse-not-accepted"
            ),
            "profileInverseAcceptanceReady": accepted,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_profile_inverse_acceptance_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumProfileInverseAcceptanceNat",
            "leanPredicate": "reachableEquilibriumProfileInverseAcceptanceReadyNat",
            "schema": "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
        },
        "claimLabels": {
            "gate": "finite thresholded acceptance predicate",
            "profile": "read-only empirical metadata",
            "inverse": "accepted only for preview, not automatic hardware execution",
            "physicalAccuracy": "not a physical-law or convergence proof",
        },
        "limitations": [
            "Acceptance means the residual certificate fits declared thresholds; it does not prove nonlinear reachability.",
            "Underactuated and topology-blocked targets default to review unless explicitly allowed.",
            "Physical validation remains optional unless required by the caller.",
        ],
    }


def export_reachable_equilibrium_profile_inverse_acceptance_json(
    profile_inverse_report: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_profile_inverse_acceptance_report(
            profile_inverse_report,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_profile_inverse_acceptance_csv(
    report: dict[str, object],
) -> str:
    summary = _dict_section(report.get("summary"))
    decision = _dict_section(report.get("decision"))
    metrics = _dict_section(report.get("metrics"))
    rows = [
        [
            "schema",
            "acceptance_ready",
            "decision",
            "accepted_for_preview",
            "profile_inverse_ready",
            "residual_score",
            "band_failures",
            "active_actuators",
            "underactuated_or_topology_blocked_cells",
            "failed_criteria",
            "missing_evidence",
        ],
        [
            report.get("schema", ""),
            summary.get("profileInverseAcceptanceReady", ""),
            decision.get("decision", ""),
            decision.get("acceptedForPreview", ""),
            metrics.get("profileInverseReady", ""),
            metrics.get("profileWeightedResidualScore", ""),
            metrics.get("profileBandFailureCount", ""),
            metrics.get("activeActuatorCount", ""),
            metrics.get("underactuatedOrTopologyBlockedCellCount", ""),
            ";".join(str(item) for item in decision.get("failedCriteria", []))
            if isinstance(decision.get("failedCriteria"), list)
            else "",
            ";".join(str(item) for item in summary.get("missingEvidence", []))
            if isinstance(summary.get("missingEvidence"), list)
            else "",
        ],
    ]
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)


def _command_preview_record(index: int, command: SourceCommand) -> dict[str, object]:
    row, col = command.cell
    return {
        "index": int(index),
        "row": int(row),
        "col": int(col),
        "alpha": float(command.alpha),
        "z": float(command.z),
        "commandType": "profile-inverse-actuator-command",
    }


def reachable_equilibrium_profile_inverse_preview_packet(
    inverse: InverseDesignResult,
    profile_inverse_report: dict[str, object],
    acceptance_report: dict[str, object],
    *,
    packet_id: str | None = None,
    notes: str = "",
) -> dict[str, object]:
    """Package an accepted profile-aware inverse plan for preview/lab handoff."""

    acceptance_summary = _dict_section(acceptance_report.get("summary"))
    acceptance_decision = _dict_section(acceptance_report.get("decision"))
    profile_summary = _dict_section(profile_inverse_report.get("summary"))
    target = _dict_section(profile_inverse_report.get("target"))
    residuals = _dict_section(profile_inverse_report.get("residuals"))
    commands = [
        _command_preview_record(index, command)
        for index, command in enumerate(inverse.commands)
    ]
    preview_events = [
        {
            "index": command["index"],
            "type": "set-actuator-command",
            "row": command["row"],
            "col": command["col"],
            "alpha": command["alpha"],
            "z": command["z"],
            "source": "accepted-profile-aware-inverse",
        }
        for command in commands
    ]
    missing_evidence: list[str] = []
    if profile_inverse_report.get("schema") != "rad-sim.reachable-equilibrium-profile-inverse.v1":
        missing_evidence.append("profileInverseReportSchema")
    if acceptance_report.get("schema") != "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1":
        missing_evidence.append("acceptanceReportSchema")
    if not profile_summary.get("profileInverseReady"):
        missing_evidence.append("profileInverseReady")
    if not acceptance_summary.get("profileInverseAcceptanceReady"):
        missing_evidence.append("profileInverseAcceptanceReady")
    if acceptance_decision.get("acceptedForPreview") is not True:
        missing_evidence.append("acceptedForPreview")
    if not commands:
        missing_evidence.append("commandRecords")
    if _report_int(target.get("targetCellCount")) <= 0:
        missing_evidence.append("targetRecords")
    if (
        _report_float(residuals.get("heightRmsResidual")) is None
        and _report_float(residuals.get("alphaRmsResidual")) is None
    ):
        missing_evidence.append("residualRecords")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
        "packetId": packet_id or "profile-inverse-preview",
        "sourceReportSchemas": {
            "profileInverse": profile_inverse_report.get("schema", ""),
            "acceptance": acceptance_report.get("schema", ""),
            "inverse": "rad-sim.inverse-design-result.v1",
        },
        "method": (
            "read-only packet of accepted profile-aware inverse commands for "
            "visual preview, bench review, or supplemental artifacts"
        ),
        "acceptance": {
            "decision": acceptance_decision.get("decision", ""),
            "acceptedForPreview": bool(acceptance_decision.get("acceptedForPreview")),
            "failedCriteria": acceptance_decision.get("failedCriteria", []),
        },
        "target": {
            "targetCellCount": target.get("targetCellCount", 0),
            "alphaTargetCellCount": target.get("alphaTargetCellCount", 0),
            "heightTargetCellCount": target.get("heightTargetCellCount", 0),
        },
        "residuals": {
            "profileWeightedResidualScore": profile_summary.get(
                "profileWeightedResidualScore"
            ),
            "profileBandFailureCount": profile_summary.get("profileBandFailureCount"),
            "alphaRmsResidual": residuals.get("alphaRmsResidual"),
            "heightRmsResidual": residuals.get("heightRmsResidual"),
        },
        "commands": commands,
        "previewEvents": preview_events,
        "reviewProtocol": {
            "mode": "preview-only",
            "operator": "accepted-profile-aware-inverse",
            "requiresHumanReviewBeforeHardware": True,
            "notes": notes,
            "steps": [
                "Load the packet in the browser or notebook.",
                "Inspect target residual, unreachable masks, and actuator budget.",
                "Run spring-hinge or rigid-body contact validation before physical actuation.",
                "Record any manual overrides as a separate event sequence.",
            ],
        },
        "summary": {
            "status": (
                "profile-aware-inverse-preview-packet-ready"
                if ready
                else "profile-aware-inverse-preview-packet-not-ready"
            ),
            "profileInversePreviewPacketReady": ready,
            "commandCount": len(commands),
            "eventCount": len(preview_events),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_profile_inverse_preview_packet_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumProfileInversePreviewPacketNat",
            "leanPredicate": "reachableEquilibriumProfileInversePreviewPacketReadyNat",
            "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
        },
        "claimLabels": {
            "packet": "finite preview handoff artifact",
            "commands": "accepted inverse command proposal",
            "hardware": "requires separate human and physical validation",
            "physicalAccuracy": "not a hardware execution proof",
        },
        "limitations": [
            "The packet is read-only and does not apply commands.",
            "Accepted preview packets still require physical validation before hardware use.",
            "The commands come from the existing inverse solve and are not re-optimized here.",
        ],
    }


def export_reachable_equilibrium_profile_inverse_preview_packet_json(
    inverse: InverseDesignResult,
    profile_inverse_report: dict[str, object],
    acceptance_report: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_profile_inverse_preview_packet(
            inverse,
            profile_inverse_report,
            acceptance_report,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_profile_inverse_preview_packet_csv(
    packet: dict[str, object],
) -> str:
    summary = _dict_section(packet.get("summary"))
    target = _dict_section(packet.get("target"))
    acceptance = _dict_section(packet.get("acceptance"))
    residuals = _dict_section(packet.get("residuals"))
    rows = [
        [
            "schema",
            "packet_ready",
            "decision",
            "target_cells",
            "command_count",
            "event_count",
            "profile_weighted_residual_score",
            "profile_band_failures",
            "missing_evidence",
        ],
        [
            packet.get("schema", ""),
            summary.get("profileInversePreviewPacketReady", ""),
            acceptance.get("decision", ""),
            target.get("targetCellCount", ""),
            summary.get("commandCount", ""),
            summary.get("eventCount", ""),
            residuals.get("profileWeightedResidualScore", ""),
            residuals.get("profileBandFailureCount", ""),
            ";".join(str(item) for item in summary.get("missingEvidence", []))
            if isinstance(summary.get("missingEvidence"), list)
            else "",
        ],
    ]
    for command in packet.get("commands", []):
        if not isinstance(command, dict):
            continue
        rows.append(
            [
                "command",
                command.get("index", ""),
                command.get("row", ""),
                command.get("col", ""),
                command.get("alpha", ""),
                command.get("z", ""),
                command.get("commandType", ""),
                "",
                "",
            ]
        )
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)


def _packet_command_sources(
    config: LatticeConfig,
    packet: dict[str, object],
    *,
    locked_cells: Iterable[tuple[int, int]] = (),
    removed_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> tuple[tuple[SourceCommand, ...], list[dict[str, object]]]:
    locked = {(int(r), int(c)) for r, c in locked_cells}
    removed = {(int(r), int(c)) for r, c in removed_cells}
    commands: list[SourceCommand] = []
    invalid: list[dict[str, object]] = []
    for index, raw in enumerate(packet.get("commands", [])):
        if not isinstance(raw, dict):
            invalid.append({"index": index, "reason": "command record is not an object"})
            continue
        row_value = raw.get("row")
        col_value = raw.get("col")
        alpha = _report_float(raw.get("alpha"))
        z = _report_float(raw.get("z"))
        try:
            row = int(row_value)  # type: ignore[arg-type]
            col = int(col_value)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            invalid.append({"index": index, "reason": "row/col must be integers"})
            continue
        if not (0 <= row < config.rows and 0 <= col < config.cols):
            invalid.append({"index": index, "row": row, "col": col, "reason": "cell outside grid"})
            continue
        if (row, col) in removed:
            invalid.append({"index": index, "row": row, "col": col, "reason": "cell is removed"})
            continue
        if (row, col) in locked:
            invalid.append({"index": index, "row": row, "col": col, "reason": "cell is locked"})
            continue
        if alpha is None and z is None:
            invalid.append({"index": index, "row": row, "col": col, "reason": "missing finite alpha/z command"})
            continue
        alpha_value = 0.0 if alpha is None else alpha
        z_value = 0.0 if z is None else z
        if abs(alpha_value) <= tolerance and abs(z_value) <= tolerance:
            continue
        commands.append(SourceCommand((row, col), alpha=alpha_value, z=z_value))
    return tuple(commands), invalid


def _state_from_packet_commands(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    *,
    locked_cells: Iterable[tuple[int, int]] = (),
    removed_cells: Iterable[tuple[int, int]] = (),
) -> LatticeState:
    state = LatticeState.uniform(config)
    for r, c in removed_cells:
        if 0 <= int(r) < config.rows and 0 <= int(c) < config.cols:
            state.removed_mask[int(r), int(c)] = True
    for r, c in locked_cells:
        if 0 <= int(r) < config.rows and 0 <= int(c) < config.cols:
            state.locked_mask[int(r), int(c)] = True
    for command in commands:
        row, col = command.cell
        state.actuator_grid[row, col] += command.alpha
        state.z_actuator_grid[row, col] += command.z
    return state


def reachable_equilibrium_profile_inverse_preview_replay_report(
    config: LatticeConfig,
    packet: dict[str, object],
    *,
    target_alpha: np.ndarray | Iterable[Iterable[float]] | None = None,
    target_height: np.ndarray | Iterable[Iterable[float]] | None = None,
    locked_cells: Iterable[tuple[int, int]] = (),
    removed_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
    residual_agreement_tolerance: float = 1e-9,
) -> dict[str, object]:
    """Replay a preview packet deterministically through the kinematic simulator."""

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if residual_agreement_tolerance < 0:
        raise ValueError("residual_agreement_tolerance must be non-negative")
    target_alpha_array = _target_array(
        target_alpha,
        (config.rows, config.cols),
        "target_alpha",
    )
    target_height_array = _target_array(
        target_height,
        (config.rows, config.cols),
        "target_height",
    )
    packet_summary = _dict_section(packet.get("summary"))
    packet_residuals = _dict_section(packet.get("residuals"))
    commands, invalid_commands = _packet_command_sources(
        config,
        packet,
        locked_cells=locked_cells,
        removed_cells=removed_cells,
        tolerance=tolerance,
    )
    state = _state_from_packet_commands(
        config,
        commands,
        locked_cells=locked_cells,
        removed_cells=removed_cells,
    )
    result = simulate_kinematic(config, state)
    alpha_residual = (
        None if target_alpha_array is None else target_alpha_array - result.alpha
    )
    height_residual = (
        None
        if target_height_array is None
        else target_height_array - result.metadata["height"]
    )
    alpha_rms = 0.0 if alpha_residual is None else _rms(alpha_residual.reshape(-1))
    height_rms = 0.0 if height_residual is None else _rms(height_residual.reshape(-1))
    packet_alpha_rms = _report_float(packet_residuals.get("alphaRmsResidual"))
    packet_height_rms = _report_float(packet_residuals.get("heightRmsResidual"))
    residual_diffs = []
    if packet_alpha_rms is not None and alpha_residual is not None:
        residual_diffs.append(abs(alpha_rms - packet_alpha_rms))
    if packet_height_rms is not None and height_residual is not None:
        residual_diffs.append(abs(height_rms - packet_height_rms))
    max_residual_disagreement = max(residual_diffs, default=0.0)
    comparable_residuals = bool(residual_diffs)
    residual_agreement_pass = (
        comparable_residuals
        and max_residual_disagreement <= max(tolerance, residual_agreement_tolerance)
    )
    event_count = len(packet.get("previewEvents", [])) if isinstance(packet.get("previewEvents"), list) else 0
    expected_command_count = _report_int(packet_summary.get("commandCount"))
    command_event_coverage_pass = event_count >= len(commands)
    target_cell_count = int(
        np.count_nonzero(np.abs(target_alpha_array - config.initial_alpha) > tolerance)
    ) if target_alpha_array is not None else 0
    if target_height_array is not None:
        target_cell_count += int(np.count_nonzero(np.abs(target_height_array) > tolerance))
    if target_cell_count <= 0:
        packet_target = _dict_section(packet.get("target"))
        target_cell_count = _report_int(packet_target.get("targetCellCount"))

    missing_evidence: list[str] = []
    if packet.get("schema") != "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1":
        missing_evidence.append("packetSchema")
    if not packet_summary.get("profileInversePreviewPacketReady"):
        missing_evidence.append("packetReady")
    if not commands:
        missing_evidence.append("replayedCommands")
    if invalid_commands:
        missing_evidence.append("validCommandRecords")
    if expected_command_count and expected_command_count != len(commands):
        missing_evidence.append("commandCountAgreement")
    if not command_event_coverage_pass:
        missing_evidence.append("eventCoverage")
    if target_alpha_array is None and target_height_array is None:
        missing_evidence.append("targetArrays")
    if target_cell_count <= 0:
        missing_evidence.append("targetRecords")
    if not comparable_residuals:
        missing_evidence.append("residualComparison")
    if comparable_residuals and not residual_agreement_pass:
        missing_evidence.append("residualAgreement")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
        "sourcePacketSchema": packet.get("schema", ""),
        "method": (
            "deterministic kinematic replay of a read-only accepted inverse "
            "preview packet"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
        },
        "commands": {
            "expectedCommandCount": expected_command_count,
            "replayedCommandCount": len(commands),
            "invalidCommandCount": len(invalid_commands),
            "invalidCommands": invalid_commands,
            "eventCount": event_count,
            "commandEventCoveragePass": command_event_coverage_pass,
        },
        "simulation": {
            "model": result.metadata.get("model", "kinematic"),
            "meanAlpha": result.metadata.get("mean_alpha"),
            "maxAbsHeight": _max_abs(result.metadata["height"]),
            "finiteDieOffCells": int(np.count_nonzero(np.isfinite(result.metadata["die_off"]))),
        },
        "residuals": {
            "targetCellCount": target_cell_count,
            "alphaRmsResidual": alpha_rms if target_alpha_array is not None else None,
            "heightRmsResidual": height_rms if target_height_array is not None else None,
            "packetAlphaRmsResidual": packet_alpha_rms,
            "packetHeightRmsResidual": packet_height_rms,
            "maxResidualDisagreement": max_residual_disagreement,
            "residualAgreementPass": residual_agreement_pass,
        },
        "summary": {
            "status": (
                "profile-aware-inverse-preview-replay-ready"
                if ready
                else "profile-aware-inverse-preview-replay-not-ready"
            ),
            "profileInversePreviewReplayReady": ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_profile_inverse_preview_replay_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumProfileInversePreviewReplayNat",
            "leanPredicate": "reachableEquilibriumProfileInversePreviewReplayReadyNat",
            "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
        },
        "claimLabels": {
            "replay": "deterministic simulator replay certificate",
            "packet": "read-only accepted inverse command packet",
            "physicalAccuracy": "kinematic replay only; not hardware execution or contact validation",
        },
        "limitations": [
            "Replay uses the current kinematic simulator and supplied config.",
            "Residual agreement is only checked when target arrays are supplied.",
            "Passing replay does not imply spring-hinge or rigid-body contact validity.",
        ],
    }


def export_reachable_equilibrium_profile_inverse_preview_replay_json(
    config: LatticeConfig,
    packet: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_profile_inverse_preview_replay_report(
            config,
            packet,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_profile_inverse_preview_replay_csv(
    report: dict[str, object],
) -> str:
    summary = _dict_section(report.get("summary"))
    commands = _dict_section(report.get("commands"))
    residuals = _dict_section(report.get("residuals"))
    rows = [
        [
            "schema",
            "replay_ready",
            "replayed_commands",
            "invalid_commands",
            "event_count",
            "target_cells",
            "height_rms_residual",
            "max_residual_disagreement",
            "residual_agreement_pass",
            "missing_evidence",
        ],
        [
            report.get("schema", ""),
            summary.get("profileInversePreviewReplayReady", ""),
            commands.get("replayedCommandCount", ""),
            commands.get("invalidCommandCount", ""),
            commands.get("eventCount", ""),
            residuals.get("targetCellCount", ""),
            residuals.get("heightRmsResidual", ""),
            residuals.get("maxResidualDisagreement", ""),
            residuals.get("residualAgreementPass", ""),
            ";".join(str(item) for item in summary.get("missingEvidence", []))
            if isinstance(summary.get("missingEvidence"), list)
            else "",
        ],
    ]
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)


def reachable_equilibrium_profile_inverse_preview_physical_report(
    config: LatticeConfig,
    packet: dict[str, object],
    *,
    target_alpha: np.ndarray | Iterable[Iterable[float]] | None = None,
    target_height: np.ndarray | Iterable[Iterable[float]] | None = None,
    locked_cells: Iterable[tuple[int, int]] = (),
    removed_cells: Iterable[tuple[int, int]] = (),
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
    residual_agreement_tolerance: float = 1e-9,
    max_height_model_error: float | None = None,
    max_center_model_error: float | None = None,
) -> dict[str, object]:
    """Run an accepted preview packet through the 3D spring-hinge preview."""

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if residual_agreement_tolerance < 0:
        raise ValueError("residual_agreement_tolerance must be non-negative")
    if max_height_model_error is not None and max_height_model_error < 0:
        raise ValueError("max_height_model_error must be non-negative")
    if max_center_model_error is not None and max_center_model_error < 0:
        raise ValueError("max_center_model_error must be non-negative")

    target_alpha_array = _target_array(
        target_alpha,
        (config.rows, config.cols),
        "target_alpha",
    )
    target_height_array = _target_array(
        target_height,
        (config.rows, config.cols),
        "target_height",
    )
    replay = reachable_equilibrium_profile_inverse_preview_replay_report(
        config,
        packet,
        target_alpha=target_alpha_array,
        target_height=target_height_array,
        locked_cells=locked_cells,
        removed_cells=removed_cells,
        tolerance=tolerance,
        residual_agreement_tolerance=residual_agreement_tolerance,
    )
    commands, invalid_commands = _packet_command_sources(
        config,
        packet,
        locked_cells=locked_cells,
        removed_cells=removed_cells,
        tolerance=tolerance,
    )
    state = _state_from_packet_commands(
        config,
        commands,
        locked_cells=locked_cells,
        removed_cells=removed_cells,
    )
    baseline_state = _baseline_state(config, locked_cells, removed_cells)
    kinematic = simulate_kinematic(config, state)
    physical_baseline = solve_spring_hinge_3d(config, baseline_state, load_case)
    physical_result = solve_spring_hinge_3d(config, state, load_case)
    physical_height = physical_result.metadata["height"]
    kinematic_height = kinematic.metadata["height"]
    height_model_error = physical_height - kinematic_height
    center_model_error = physical_result.deformed_centers_3d - kinematic.deformed_centers_3d
    center_norm = np.linalg.norm(center_model_error, axis=2)
    alpha_residual = (
        None if target_alpha_array is None else target_alpha_array - physical_result.alpha
    )
    height_residual = (
        None if target_height_array is None else target_height_array - physical_height
    )
    target_cell_count = int(
        np.count_nonzero(
            np.abs(target_alpha_array - config.initial_alpha) > tolerance
        )
    ) if target_alpha_array is not None else 0
    if target_height_array is not None:
        target_cell_count += int(
            np.count_nonzero(np.abs(target_height_array) > tolerance)
        )
    if target_cell_count <= 0:
        target_cell_count = _report_int(
            _dict_section(packet.get("target")).get("targetCellCount")
        )

    height_rms_model_error = _rms(height_model_error.reshape(-1))
    center_rms_model_error = _rms(center_norm.reshape(-1))
    max_abs_height_model_error = _max_abs(height_model_error)
    max_abs_center_model_error = _max_abs(center_norm)
    physical_energy = _report_float(physical_result.metadata.get("energy"))
    finite_model_comparison = all(
        np.isfinite(value)
        for value in (
            height_rms_model_error,
            center_rms_model_error,
            max_abs_height_model_error,
            max_abs_center_model_error,
        )
    )
    model_error_pass = True
    if max_height_model_error is not None:
        model_error_pass = model_error_pass and (
            max_abs_height_model_error <= max_height_model_error + tolerance
        )
    if max_center_model_error is not None:
        model_error_pass = model_error_pass and (
            max_abs_center_model_error <= max_center_model_error + tolerance
        )
    physical_success = bool(
        physical_baseline.metadata.get("success", False)
        and physical_result.metadata.get("success", False)
    )
    replay_ready = bool(
        _dict_section(replay.get("summary")).get("profileInversePreviewReplayReady")
    )
    missing_evidence: list[str] = []
    if not replay_ready:
        missing_evidence.append("previewReplayReady")
    if not commands:
        missing_evidence.append("commandRecords")
    if invalid_commands:
        missing_evidence.append("validCommandRecords")
    if not physical_success:
        missing_evidence.append("physicalSolverSuccess")
    if target_alpha_array is None and target_height_array is None:
        missing_evidence.append("targetArrays")
    if target_cell_count <= 0:
        missing_evidence.append("targetRecords")
    if not finite_model_comparison:
        missing_evidence.append("modelComparisonRecords")
    if physical_energy is None:
        missing_evidence.append("energyRecords")
    if not model_error_pass:
        missing_evidence.append("modelErrorThreshold")
    ready = len(missing_evidence) == 0
    load = load_case or LoadCase()
    return {
        "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
        "sourcePacketSchema": packet.get("schema", ""),
        "sourceReplaySchema": replay.get("schema", ""),
        "method": (
            "read-only spring-hinge 3D physical preview of an accepted inverse "
            "command packet"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
        },
        "commands": {
            "replayedCommandCount": len(commands),
            "invalidCommandCount": len(invalid_commands),
        },
        "loadCase": {
            "fixedCells": [list(cell) for cell in load.fixed_cells],
            "prescribedDisplacementCount": len(load.prescribed_displacements),
            "externalForceCount": len(load.external_forces),
            "axialStiffness": load.axial_stiffness,
            "hingeStiffness": load.hinge_stiffness,
            "lockStiffness": load.lock_stiffness,
            "maxiter": load.maxiter,
        },
        "physical": {
            "model": physical_result.metadata.get("model", "spring_hinge_3d"),
            "baselineSuccess": bool(physical_baseline.metadata.get("success", False)),
            "resultSuccess": bool(physical_result.metadata.get("success", False)),
            "physicalSuccess": physical_success,
            "iterations": int(physical_result.metadata.get("iterations", 0)),
            "physicalEnergy": physical_energy,
            "storedEnergy": _report_float(physical_result.metadata.get("stored_energy")),
            "springEnergy": _report_float(physical_result.metadata.get("spring_energy")),
            "hingeEnergy": _report_float(physical_result.metadata.get("hinge_energy")),
            "lockPenaltyEnergy": _report_float(
                physical_result.metadata.get("lock_penalty_energy")
            ),
            "springEdges": int(physical_result.metadata.get("spring_edges", 0)),
            "hingeTriples": int(physical_result.metadata.get("hinge_triples", 0)),
            "removedCells": int(physical_result.metadata.get("removed_cells", 0)),
        },
        "comparison": {
            "heightRmsModelError": height_rms_model_error,
            "centerRmsModelError": center_rms_model_error,
            "maxAbsHeightModelError": max_abs_height_model_error,
            "maxAbsCenterModelError": max_abs_center_model_error,
            "modelAgreementScore": 1.0 / (1.0 + center_rms_model_error),
            "finiteModelComparison": finite_model_comparison,
            "modelErrorThresholdPass": model_error_pass,
            "maxHeightModelErrorLimit": max_height_model_error,
            "maxCenterModelErrorLimit": max_center_model_error,
        },
        "residuals": {
            "targetCellCount": target_cell_count,
            "physicalAlphaRmsResidual": (
                None if alpha_residual is None else _rms(alpha_residual.reshape(-1))
            ),
            "physicalHeightRmsResidual": (
                None if height_residual is None else _rms(height_residual.reshape(-1))
            ),
            "physicalMaxAbsHeightResidual": _max_abs(height_residual),
        },
        "summary": {
            "status": (
                "profile-aware-inverse-preview-physical-ready"
                if ready
                else "profile-aware-inverse-preview-physical-not-ready"
            ),
            "profileInversePreviewPhysicalReady": ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_profile_inverse_preview_physical_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumProfileInversePreviewPhysicalNat",
            "leanPredicate": "reachableEquilibriumProfileInversePreviewPhysicalReadyNat",
            "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
        },
        "claimLabels": {
            "gate": "finite physical-preview evidence predicate",
            "physics": "simulator-derived spring-hinge preview",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
            "hardware": "not hardware execution authorization",
        },
        "limitations": [
            "The spring-hinge preview is normalized and uncalibrated.",
            "Passing this gate does not prove rigid-body contact, friction, gravity, or material accuracy.",
            "The report replays a command packet and does not re-optimize inverse design commands.",
        ],
    }


def export_reachable_equilibrium_profile_inverse_preview_physical_json(
    config: LatticeConfig,
    packet: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_profile_inverse_preview_physical_report(
            config,
            packet,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_profile_inverse_preview_physical_csv(
    report: dict[str, object],
) -> str:
    summary = _dict_section(report.get("summary"))
    commands = _dict_section(report.get("commands"))
    physical = _dict_section(report.get("physical"))
    comparison = _dict_section(report.get("comparison"))
    residuals = _dict_section(report.get("residuals"))
    rows = [
        [
            "schema",
            "physical_ready",
            "physical_success",
            "replayed_commands",
            "target_cells",
            "physical_height_rms_residual",
            "height_rms_model_error",
            "center_rms_model_error",
            "physical_energy",
            "missing_evidence",
        ],
        [
            report.get("schema", ""),
            summary.get("profileInversePreviewPhysicalReady", ""),
            physical.get("physicalSuccess", ""),
            commands.get("replayedCommandCount", ""),
            residuals.get("targetCellCount", ""),
            residuals.get("physicalHeightRmsResidual", ""),
            comparison.get("heightRmsModelError", ""),
            comparison.get("centerRmsModelError", ""),
            physical.get("physicalEnergy", ""),
            ";".join(str(item) for item in summary.get("missingEvidence", []))
            if isinstance(summary.get("missingEvidence"), list)
            else "",
        ],
    ]
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)
