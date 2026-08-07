from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.optimize import lsq_linear

from .experiments import ResponseMatrix, SourceCommand, build_response_matrix
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, SimulationResult


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

    @property
    def active_actuator_count(self) -> int:
        return len(self.commands)


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
    )
