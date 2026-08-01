from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .geometry import cell_corners, grid_centers
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult


def _edges(rows: int, cols: int) -> list[tuple[int, int]]:
    idx = lambda r, c: r * cols + c
    out: list[tuple[int, int]] = []
    for r in range(rows):
        for c in range(cols):
            if r + 1 < rows:
                out.append((idx(r, c), idx(r + 1, c)))
            if c + 1 < cols:
                out.append((idx(r, c), idx(r, c + 1)))
    return out


def _hinge_triples(rows: int, cols: int) -> list[tuple[int, int, int]]:
    idx = lambda r, c: r * cols + c
    triples: list[tuple[int, int, int]] = []
    for r in range(rows):
        for c in range(1, cols - 1):
            triples.append((idx(r, c - 1), idx(r, c), idx(r, c + 1)))
    for c in range(cols):
        for r in range(1, rows - 1):
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


def solve_spring_hinge(
    config: LatticeConfig,
    state: LatticeState,
    load_case: LoadCase | None = None,
) -> SimulationResult:
    load_case = load_case or LoadCase()
    kin = simulate_kinematic(config, state)
    rows, cols = config.rows, config.cols
    initial = kin.original_centers.reshape((-1, 2))
    target = kin.deformed_centers.reshape((-1, 2))
    edges = _edges(rows, cols)
    rest_lengths = np.array([np.linalg.norm(initial[i] - initial[j]) for i, j in edges])
    triples = _hinge_triples(rows, cols)
    rest_angles = np.array([hinge_angle(initial[i], initial[j], initial[k]) for i, j, k in triples])

    fixed: set[int] = {r * cols + c for r, c in load_case.fixed_cells}
    prescribed = {
        r * cols + c: np.asarray(delta, dtype=float)
        for (r, c), delta in load_case.prescribed_displacements.items()
    }
    force = {
        r * cols + c: np.asarray(value, dtype=float)
        for (r, c), value in load_case.external_forces.items()
    }
    lock_mask = state.normalized(config).locked_mask.reshape(-1)
    actuator_mask = np.abs(state.normalized(config).actuator_grid.reshape(-1)) > 1e-12
    penalty_mask = lock_mask | actuator_mask

    def unpack(q: np.ndarray) -> np.ndarray:
        points = q.reshape((-1, 2)).copy()
        for idx in fixed:
            points[idx] = initial[idx]
        for idx, delta in prescribed.items():
            points[idx] = initial[idx] + delta
        return points

    def objective(q: np.ndarray) -> float:
        points = unpack(q)
        value = spring_energy(points, edges, rest_lengths, load_case.axial_stiffness)
        value += hinge_energy(points, triples, rest_angles, load_case.hinge_stiffness)
        for idx, enabled in enumerate(penalty_mask):
            if enabled:
                value += 0.5 * load_case.lock_stiffness * np.sum((points[idx] - target[idx]) ** 2)
        for idx, f in force.items():
            value -= float(np.dot(f, points[idx] - initial[idx]))
        return float(value)

    q0 = target.reshape(-1)
    result = minimize(
        objective,
        q0,
        method="L-BFGS-B",
        options={"maxiter": load_case.maxiter, "ftol": 1e-10},
    )
    centers = unpack(result.x).reshape((rows, cols, 2))
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
            "iterations": int(result.nit),
        },
    )
