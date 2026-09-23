from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.optimize import lsq_linear

from .experiments import (
    ResponseMatrix,
    SourceCommand,
    build_response_matrix,
    characterize_response,
)
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .spring_hinge import solve_spring_hinge_3d


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
) -> LatticeState:
    state = LatticeState.uniform(config)
    for r, c in locked_cells:
        state.locked_mask[r, c] = True
    return state


def _unlocked_cells(
    config: LatticeConfig,
    actuator_cells: Iterable[tuple[int, int]] | None,
    locked_cells: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    locked = set(locked_cells)
    cells = (
        actuator_cells
        if actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    return tuple((int(r), int(c)) for r, c in cells if (int(r), int(c)) not in locked)


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
) -> LatticeState:
    state = _baseline_state(config, locked_cells)
    for command in commands:
        r, c = command.cell
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
            response = characterize_response(config, (command,), locked_cells, tolerance)
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
    candidates = _unlocked_cells(config, actuator_cells, locked)
    matrix = build_response_matrix(
        config,
        actuator_cells=candidates,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=locked,
        tolerance=tolerance,
    )
    baseline = simulate_kinematic(config, _baseline_state(config, locked))
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
    state = _state_with_commands(config, commands, locked)
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
