from __future__ import annotations

import numpy as np

try:
    from scipy.optimize import minimize
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only in minimal runtimes.
    if exc.name not in {"scipy", "scipy.optimize"}:
        raise
    _SPRING_HINGE_IMPORT_ERROR = exc

    def minimize(*args, **kwargs):
        raise ImportError(
            "SciPy is required for spring-hinge mechanics. Install project dependencies "
            "before using solve_spring_hinge or solve_spring_hinge_3d."
        ) from _SPRING_HINGE_IMPORT_ERROR

from .geometry import cell_corners, grid_centers
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult


def _active_cell(active_mask: np.ndarray | None, r: int, c: int) -> bool:
    return active_mask is None or bool(active_mask[r, c])


def _edges(
    rows: int, cols: int, active_mask: np.ndarray | None = None
) -> list[tuple[int, int]]:
    idx = lambda r, c: r * cols + c
    out: list[tuple[int, int]] = []
    for r in range(rows):
        for c in range(cols):
            if not _active_cell(active_mask, r, c):
                continue
            if r + 1 < rows and _active_cell(active_mask, r + 1, c):
                out.append((idx(r, c), idx(r + 1, c)))
            if c + 1 < cols and _active_cell(active_mask, r, c + 1):
                out.append((idx(r, c), idx(r, c + 1)))
    return out


def _hinge_triples(
    rows: int, cols: int, active_mask: np.ndarray | None = None
) -> list[tuple[int, int, int]]:
    idx = lambda r, c: r * cols + c
    triples: list[tuple[int, int, int]] = []
    for r in range(rows):
        for c in range(1, cols - 1):
            if not (
                _active_cell(active_mask, r, c - 1)
                and _active_cell(active_mask, r, c)
                and _active_cell(active_mask, r, c + 1)
            ):
                continue
            triples.append((idx(r, c - 1), idx(r, c), idx(r, c + 1)))
    for c in range(cols):
        for r in range(1, rows - 1):
            if not (
                _active_cell(active_mask, r - 1, c)
                and _active_cell(active_mask, r, c)
                and _active_cell(active_mask, r + 1, c)
            ):
                continue
            triples.append((idx(r - 1, c), idx(r, c), idx(r + 1, c)))
    return triples


def spring_energy(points: np.ndarray, edges: list[tuple[int, int]], rest: np.ndarray, k: float) -> float:
    energy = 0.0
    for n, (i, j) in enumerate(edges):
        length = np.linalg.norm(points[i] - points[j])
        energy += 0.5 * k * (length - rest[n]) ** 2
    return float(energy)


def hinge_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u = a - b
    v = c - b
    denom = np.linalg.norm(u) * np.linalg.norm(v)
    if denom < 1e-12:
        return 0.0
    cosang = np.clip(np.dot(u, v) / denom, -1.0, 1.0)
    return float(np.arccos(cosang))


def hinge_energy(
    points: np.ndarray,
    triples: list[tuple[int, int, int]],
    rest_angles: np.ndarray,
    k_theta: float,
) -> float:
    energy = 0.0
    for n, (i, j, k) in enumerate(triples):
        angle = hinge_angle(points[i], points[j], points[k])
        energy += 0.5 * k_theta * (angle - rest_angles[n]) ** 2
    return float(energy)


def _lock_penalty_energy(
    points: np.ndarray,
    target: np.ndarray,
    penalty_mask: np.ndarray,
    k_lock: float,
) -> float:
    energy = 0.0
    for idx, enabled in enumerate(penalty_mask):
        if enabled:
            energy += 0.5 * k_lock * float(np.sum((points[idx] - target[idx]) ** 2))
    return float(energy)


def _external_potential_energy(
    points: np.ndarray,
    initial: np.ndarray,
    force: dict[int, np.ndarray],
) -> float:
    energy = 0.0
    for idx, f in force.items():
        energy -= float(np.dot(f, points[idx] - initial[idx]))
    return float(energy)


def _energy_breakdown(
    points: np.ndarray,
    initial: np.ndarray,
    target: np.ndarray,
    edges: list[tuple[int, int]],
    rest_lengths: np.ndarray,
    triples: list[tuple[int, int, int]],
    rest_angles: np.ndarray,
    penalty_mask: np.ndarray,
    force: dict[int, np.ndarray],
    load_case: LoadCase,
) -> dict[str, float]:
    spring = spring_energy(points, edges, rest_lengths, load_case.axial_stiffness)
    hinge = hinge_energy(points, triples, rest_angles, load_case.hinge_stiffness)
    lock = _lock_penalty_energy(points, target, penalty_mask, load_case.lock_stiffness)
    external = _external_potential_energy(points, initial, force)
    stored = spring + hinge + lock
    return {
        "spring_energy": float(spring),
        "hinge_energy": float(hinge),
        "lock_penalty_energy": float(lock),
        "stored_energy": float(stored),
        "external_potential_energy": float(external),
        "objective_energy": float(stored + external),
    }


def _vector(value: tuple[float, ...] | np.ndarray, dimensions: int) -> np.ndarray:
    vector = np.asarray(value, dtype=float)
    if vector.ndim != 1 or vector.size > dimensions:
        raise ValueError(f"load vector must have at most {dimensions} components")
    out = np.zeros(dimensions, dtype=float)
    out[: vector.size] = vector
    return out


def solve_spring_hinge(
    config: LatticeConfig,
    state: LatticeState,
    load_case: LoadCase | None = None,
) -> SimulationResult:
    load_case = load_case or LoadCase()
    kin = simulate_kinematic(config, state)
    rows, cols = config.rows, config.cols
    normalized = state.normalized(config)
    active_mask = ~normalized.removed_mask
    initial = kin.original_centers.reshape((-1, 2))
    target = kin.deformed_centers.reshape((-1, 2))
    edges = _edges(rows, cols, active_mask)
    rest_lengths = np.array([np.linalg.norm(initial[i] - initial[j]) for i, j in edges])
    triples = _hinge_triples(rows, cols, active_mask)
    rest_angles = np.array(
        [hinge_angle(initial[i], initial[j], initial[k]) for i, j, k in triples]
    )

    position_fixed = set(
        map(int, np.flatnonzero((normalized.position_locked_mask & active_mask).reshape(-1)))
    )
    fixed: set[int] = {r * cols + c for r, c in load_case.fixed_cells} | position_fixed
    prescribed = {
        r * cols + c: np.asarray(delta, dtype=float)
        for (r, c), delta in load_case.prescribed_displacements.items()
    }
    force = {
        r * cols + c: np.asarray(value, dtype=float)
        for (r, c), value in load_case.external_forces.items()
    }
    lock_mask = normalized.locked_mask.reshape(-1)
    actuator_mask = np.abs(normalized.actuator_grid.reshape(-1)) > 1e-12
    penalty_mask = lock_mask | actuator_mask

    def unpack(q: np.ndarray) -> np.ndarray:
        points = q.reshape((-1, 2)).copy()
        for idx, delta in prescribed.items():
            points[idx] = initial[idx] + delta
        for idx in fixed:
            points[idx] = initial[idx]
        return points

    def objective(q: np.ndarray) -> float:
        points = unpack(q)
        return _energy_breakdown(
            points,
            initial,
            target,
            edges,
            rest_lengths,
            triples,
            rest_angles,
            penalty_mask,
            force,
            load_case,
        )["objective_energy"]

    q0 = target.reshape(-1)
    result = minimize(
        objective,
        q0,
        method="L-BFGS-B",
        options={"maxiter": load_case.maxiter, "ftol": 1e-10},
    )
    centers = unpack(result.x).reshape((rows, cols, 2))
    final_points = centers.reshape((-1, 2))
    energy = _energy_breakdown(
        final_points,
        initial,
        target,
        edges,
        rest_lengths,
        triples,
        rest_angles,
        penalty_mask,
        force,
        load_case,
    )
    corners = cell_corners(config, centers, kin.alpha)
    height = kin.deformed_centers_3d[..., 2]
    centers_3d = np.dstack([centers[..., 0], centers[..., 1], height])
    corners_3d = np.zeros((rows, cols, 5, 3), dtype=float)
    corners_3d[..., :2] = corners
    corners_3d[..., 2] = height[..., None]
    return SimulationResult(
        config=config,
        state=kin.state,
        alpha=kin.alpha,
        theta_degrees=kin.theta_degrees,
        original_centers=kin.original_centers,
        deformed_centers=centers,
        original_corners=kin.original_corners,
        deformed_corners=corners,
        original_centers_3d=kin.original_centers_3d,
        deformed_centers_3d=centers_3d,
        original_corners_3d=kin.original_corners_3d,
        deformed_corners_3d=corners_3d,
        complex_original=kin.complex_original,
        complex_deformed=centers[..., 0] + 1j * centers[..., 1],
        metadata={
            "model": "spring_hinge",
            "success": bool(result.success),
            "message": str(result.message),
            "energy": float(result.fun),
            "objective_energy": energy["objective_energy"],
            "stored_energy": energy["stored_energy"],
            "spring_energy": energy["spring_energy"],
            "hinge_energy": energy["hinge_energy"],
            "lock_penalty_energy": energy["lock_penalty_energy"],
            "external_potential_energy": energy["external_potential_energy"],
            "iterations": int(result.nit),
            "position_locked_cells": len(position_fixed),
        },
    )


def solve_spring_hinge_3d(
    config: LatticeConfig,
    state: LatticeState,
    load_case: LoadCase | None = None,
) -> SimulationResult:
    load_case = load_case or LoadCase()
    kin = simulate_kinematic(config, state)
    rows, cols = config.rows, config.cols
    normalized = state.normalized(config)
    active_mask = ~normalized.removed_mask
    initial = kin.original_centers_3d.reshape((-1, 3))
    target = kin.deformed_centers_3d.reshape((-1, 3))
    edges = _edges(rows, cols, active_mask)
    rest_lengths = np.array([np.linalg.norm(initial[i] - initial[j]) for i, j in edges])
    triples = _hinge_triples(rows, cols, active_mask)
    rest_angles = np.array([hinge_angle(initial[i], initial[j], initial[k]) for i, j, k in triples])

    position_fixed = set(
        map(int, np.flatnonzero((normalized.position_locked_mask & active_mask).reshape(-1)))
    )
    fixed: set[int] = {r * cols + c for r, c in load_case.fixed_cells} | position_fixed
    prescribed = {
        r * cols + c: _vector(delta, 3)
        for (r, c), delta in load_case.prescribed_displacements.items()
    }
    force = {
        r * cols + c: _vector(value, 3)
        for (r, c), value in load_case.external_forces.items()
    }
    lock_mask = normalized.locked_mask.reshape(-1)
    alpha_actuator_mask = np.abs(normalized.actuator_grid.reshape(-1)) > 1e-12
    z_actuator_mask = np.abs(normalized.z_actuator_grid.reshape(-1)) > 1e-12
    penalty_mask = lock_mask | alpha_actuator_mask | z_actuator_mask

    def unpack(q: np.ndarray) -> np.ndarray:
        points = q.reshape((-1, 3)).copy()
        for idx, delta in prescribed.items():
            points[idx] = initial[idx] + delta
        for idx in fixed:
            points[idx] = initial[idx]
        return points

    def objective(q: np.ndarray) -> float:
        points = unpack(q)
        return _energy_breakdown(
            points,
            initial,
            target,
            edges,
            rest_lengths,
            triples,
            rest_angles,
            penalty_mask,
            force,
            load_case,
        )["objective_energy"]

    q0 = target.reshape(-1)
    result = minimize(
        objective,
        q0,
        method="L-BFGS-B",
        options={"maxiter": load_case.maxiter, "ftol": 1e-10},
    )
    points = unpack(result.x).reshape((rows, cols, 3))
    final_points = points.reshape((-1, 3))
    energy = _energy_breakdown(
        final_points,
        initial,
        target,
        edges,
        rest_lengths,
        triples,
        rest_angles,
        penalty_mask,
        force,
        load_case,
    )
    centers = points[..., :2]
    height = points[..., 2]
    corners = cell_corners(config, centers, kin.alpha)
    corners_3d = np.zeros((rows, cols, 5, 3), dtype=float)
    corners_3d[..., :2] = corners
    corners_3d[..., 2] = height[..., None]
    target_error = points.reshape((-1, 3)) - target
    return SimulationResult(
        config=config,
        state=kin.state,
        alpha=kin.alpha,
        theta_degrees=kin.theta_degrees,
        original_centers=kin.original_centers,
        deformed_centers=centers,
        original_corners=kin.original_corners,
        deformed_corners=corners,
        original_centers_3d=kin.original_centers_3d,
        deformed_centers_3d=points,
        original_corners_3d=kin.original_corners_3d,
        deformed_corners_3d=corners_3d,
        complex_original=kin.complex_original,
        complex_deformed=centers[..., 0] + 1j * centers[..., 1],
        metadata={
            "model": "spring_hinge_3d",
            "success": bool(result.success),
            "message": str(result.message),
            "energy": float(result.fun),
            "objective_energy": energy["objective_energy"],
            "stored_energy": energy["stored_energy"],
            "spring_energy": energy["spring_energy"],
            "hinge_energy": energy["hinge_energy"],
            "lock_penalty_energy": energy["lock_penalty_energy"],
            "external_potential_energy": energy["external_potential_energy"],
            "iterations": int(result.nit),
            "height": height,
            "kinematic_height": kin.metadata["height"],
            "target_rms_error": float(np.sqrt(np.mean(target_error**2))),
            "spring_edges": len(edges),
            "hinge_triples": len(triples),
            "removed_cells": int(np.count_nonzero(normalized.removed_mask)),
            "position_locked_cells": len(position_fixed),
        },
    )
