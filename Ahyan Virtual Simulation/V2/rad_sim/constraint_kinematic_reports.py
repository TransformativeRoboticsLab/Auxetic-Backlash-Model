from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .constraint_kinematic import (
    ConstraintKinematicSettings,
    active_neighbor_edges,
    solve_constraint_kinematic,
)
from .models import LatticeConfig, LatticeState, SimulationResult


CONSTRAINT_REPORT_SCHEMA = "rad-sim.constraint-kinematic-scenario-report.v1"


def _result_summary(name: str, result: SimulationResult) -> dict[str, Any]:
    height = result.deformed_centers_3d[..., 2]
    metadata = result.metadata
    return {
        "name": name,
        "rows": result.config.rows,
        "cols": result.config.cols,
        "cellSize": result.config.cell_size,
        "backlash": result.config.backlash,
        "success": bool(metadata.get("success", False)),
        "iterations": int(metadata.get("iterations", 0)),
        "edgeCount": int(metadata.get("edge_count", 0)),
        "thetaDeadZoneDegrees": float(metadata.get("theta_dead_zone_degrees", 0.0)),
        "maxEdgeError": float(metadata.get("max_edge_error", 0.0)),
        "meanEdgeError": float(metadata.get("mean_edge_error", 0.0)),
        "maxBacklashResidual": float(metadata.get("max_backlash_residual", 0.0)),
        "meanBacklashResidual": float(metadata.get("mean_backlash_residual", 0.0)),
        "boundaryContactCount": int(metadata.get("boundary_contact_count", 0)),
        "maxHeight": float(np.max(height)) if height.size else 0.0,
        "minHeight": float(np.min(height)) if height.size else 0.0,
        "maxAbsHeight": float(np.max(np.abs(height))) if height.size else 0.0,
        "meanAlpha": float(np.mean(result.alpha)),
        "minTheta": float(np.min(result.theta_degrees)),
        "maxTheta": float(np.max(result.theta_degrees)),
    }


def _actuated_state(config: LatticeConfig, command: tuple[int, int, float]) -> LatticeState:
    state = LatticeState.uniform(config)
    r, c, alpha_command = command
    state.actuator_grid[r, c] = alpha_command
    return state


def _fixed_end_state(config: LatticeConfig, command: tuple[int, int, float]) -> LatticeState:
    state = _actuated_state(config, command)
    state.position_locked_mask[0, 0] = True
    state.position_locked_mask[config.rows - 1, config.cols - 1] = True
    state.position_locked_mask[0, config.cols - 1] = True
    state.position_locked_mask[config.rows - 1, 0] = True
    return state


def constraint_solver_scenario_reports() -> dict[str, Any]:
    """Run compact benchmark scenarios for the paper-grounded constraint solver."""

    scenarios: list[tuple[str, LatticeConfig, LatticeState, dict[str, Any]]] = []
    scenarios.append(
        (
            "1x2 free pair",
            LatticeConfig(rows=1, cols=2, backlash=0.08),
            _actuated_state(LatticeConfig(rows=1, cols=2, backlash=0.08), (0, 0, -0.22)),
            {},
        )
    )
    scenarios.append(
        (
            "1x15 fixed-end strand",
            LatticeConfig(rows=1, cols=15, backlash=0.08),
            _fixed_end_state(LatticeConfig(rows=1, cols=15, backlash=0.08), (0, 7, -0.22)),
            {},
        )
    )
    scenarios.append(
        (
            "15x1 fixed-end strand",
            LatticeConfig(rows=15, cols=1, backlash=0.08),
            _fixed_end_state(LatticeConfig(rows=15, cols=1, backlash=0.08), (7, 0, -0.22)),
            {},
        )
    )
    config_3x3 = LatticeConfig(rows=3, cols=3, backlash=0.08)
    scenarios.append(("3x3 center actuation", config_3x3, _actuated_state(config_3x3, (1, 1, -0.2)), {}))
    config_15x15 = LatticeConfig(rows=15, cols=15, backlash=0.08)
    scenarios.append(("15x15 center actuation", config_15x15, _actuated_state(config_15x15, (7, 7, -0.16)), {}))

    summaries = []
    report_settings = ConstraintKinematicSettings(max_nfev=80)
    for name, config, state, boundary in scenarios:
        normalized = state.normalized(config)
        result = solve_constraint_kinematic(config, normalized, boundary=boundary, settings=report_settings)
        summary = _result_summary(name, result)
        summary["activeEdges"] = len(active_neighbor_edges(config, normalized))
        summaries.append(summary)

    return {
        "schema": CONSTRAINT_REPORT_SCHEMA,
        "model": "constraint_kinematic",
        "governingEquations": {
            "normalizedBacklash": "b_norm = b / L",
            "angleDeadZone": "Delta phi = asin(b / L)",
            "deadZone": "f(x) = max(0, x - b) + min(x + b, 0)",
            "alphaTheta": "theta[degrees] = 70 alpha - 60",
            "pairDistance": "D_ij = L(cos(theta_i) + cos(theta_j))",
        },
        "scenarios": summaries,
    }


def export_constraint_solver_scenario_reports_json(report: dict[str, Any] | None = None) -> str:
    return json.dumps(report or constraint_solver_scenario_reports(), indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export RAD constraint-kinematic scenario diagnostics.")
    parser.add_argument("--output", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    payload = export_constraint_solver_scenario_reports_json()
    if args.output is None:
        print(payload)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
