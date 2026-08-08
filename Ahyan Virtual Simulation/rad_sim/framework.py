from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .experiments import (
    OperatorInteractionGraph,
    PhysicalResponseComparison,
    ResponseCharacterization,
    ResponseDecayProfile,
    ResponseMatrix,
    SourceCommand,
    build_response_matrix,
    characterize_cluster,
    characterize_pairwise_interactions,
    characterize_response,
    compare_physical_cluster,
    response_decay_profile,
)
from .models import LatticeConfig, LatticeState, LoadCase
from .operators import (
    ProgrammableDiscontinuityEvent,
    SequenceOrderDiagnostic,
    compare_sequence_order,
)


@dataclass(frozen=True)
class ProgrammableDiscontinuityDiagnostic:
    """Executable summary of locality, reachability, and operator composition.

    The dead-zone law and alpha/theta kinematics are paper-supported. Treating
    superposition residuals as evidence of non-additive operator interaction is
    a modeling diagnostic introduced here for the programmable-discontinuity
    framework.
    """

    commands: tuple[SourceCommand, ...]
    locked_cells: tuple[tuple[int, int], ...]
    combined: ResponseCharacterization
    response_matrix: ResponseMatrix
    alpha_superposition_error: float
    height_superposition_error: float
    tolerance: float
    event_sequence: tuple[ProgrammableDiscontinuityEvent, ...] = ()
    sequence_order: SequenceOrderDiagnostic | None = None
    decay_profile: ResponseDecayProfile | None = None
    pairwise_interactions: OperatorInteractionGraph | None = None
    physical_validation: PhysicalResponseComparison | None = None

    @property
    def active_operator_count(self) -> int:
        return sum(
            1
            for command in self.commands
            if abs(command.alpha) > self.tolerance or abs(command.z) > self.tolerance
        )

    @property
    def alpha_locality_radius(self) -> int:
        return self.combined.effective_alpha_die_off

    @property
    def z_locality_radius(self) -> int:
        return self.combined.effective_z_die_off

    @property
    def alpha_decay_ratio(self) -> float:
        return 0.0 if self.decay_profile is None else self.decay_profile.alpha_ratio

    @property
    def z_decay_ratio(self) -> float:
        return 0.0 if self.decay_profile is None else self.decay_profile.z_ratio

    @property
    def alpha_decay_length(self) -> float:
        return 0.0 if self.decay_profile is None else self.decay_profile.alpha_length

    @property
    def z_decay_length(self) -> float:
        return 0.0 if self.decay_profile is None else self.decay_profile.z_length

    @property
    def reachable_alpha_cells(self) -> int:
        return self.response_matrix.reachable_alpha_cells(self.tolerance)

    @property
    def reachable_height_cells(self) -> int:
        return self.response_matrix.reachable_height_cells(self.tolerance)

    @property
    def alpha_rank(self) -> int:
        return self.response_matrix.alpha_rank

    @property
    def height_rank(self) -> int:
        return self.response_matrix.height_rank

    @property
    def nonadditive(self) -> bool:
        return (
            self.alpha_superposition_error > self.tolerance
            or self.height_superposition_error > self.tolerance
        )

    @property
    def order_sensitive(self) -> bool:
        return bool(self.sequence_order and self.sequence_order.order_sensitive)

    @property
    def noncommuting_adjacent_pairs(self) -> int:
        return 0 if self.sequence_order is None else self.sequence_order.noncommuting_adjacent_pairs

    @property
    def max_order_error(self) -> float:
        return 0.0 if self.sequence_order is None else self.sequence_order.max_order_error

    @property
    def nonadditive_pair_count(self) -> int:
        if self.pairwise_interactions is None:
            return 0
        return self.pairwise_interactions.nonadditive_pair_count

    @property
    def max_pairwise_interaction_error(self) -> float:
        if self.pairwise_interactions is None:
            return 0.0
        return self.pairwise_interactions.max_interaction_error

    @property
    def pairwise_interaction_hotspot_map(self) -> np.ndarray:
        if self.pairwise_interactions is None:
            return np.zeros_like(self.combined.alpha_delta)
        return self.pairwise_interactions.interaction_hotspot_map

    @property
    def max_pairwise_hotspot_error(self) -> float:
        if self.pairwise_interactions is None:
            return 0.0
        return self.pairwise_interactions.max_hotspot_error

    @property
    def pairwise_interaction_degree_map(self) -> np.ndarray:
        if self.pairwise_interactions is None:
            return np.zeros_like(self.combined.alpha_delta, dtype=int)
        return self.pairwise_interactions.interaction_degree_map

    @property
    def max_pairwise_interaction_degree(self) -> int:
        if self.pairwise_interactions is None:
            return 0
        return self.pairwise_interactions.max_interaction_degree

    @property
    def pairwise_interaction_density(self) -> float:
        if self.pairwise_interactions is None:
            return 0.0
        return self.pairwise_interactions.interaction_density

    @property
    def pairwise_interactions_truncated(self) -> bool:
        return bool(self.pairwise_interactions and self.pairwise_interactions.truncated)

    @property
    def total_cells(self) -> int:
        return self.combined.alpha_delta.size

    @property
    def alpha_underactuated_cells(self) -> int:
        return max(0, self.total_cells - self.reachable_alpha_cells)

    @property
    def height_underactuated_cells(self) -> int:
        return max(0, self.total_cells - self.reachable_height_cells)

    def to_dict(
        self,
        config: LatticeConfig | None = None,
        *,
        include_response_matrix: bool = True,
        include_fields: bool = True,
    ) -> dict[str, object]:
        return programmable_discontinuity_report(
            self,
            config=config,
            include_response_matrix=include_response_matrix,
            include_fields=include_fields,
        )


def _unique_command_cells(commands: tuple[SourceCommand, ...]) -> tuple[tuple[int, int], ...]:
    cells: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for command in commands:
        cell = (int(command.cell[0]), int(command.cell[1]))
        if cell in seen:
            continue
        seen.add(cell)
        cells.append(cell)
    return tuple(cells)


def _superposition_errors(
    config: LatticeConfig,
    commands: tuple[SourceCommand, ...],
    combined: ResponseCharacterization,
    locked_cells: tuple[tuple[int, int], ...],
    tolerance: float,
) -> tuple[float, float]:
    if len(commands) <= 1:
        return 0.0, 0.0
    alpha_expected = np.zeros_like(combined.alpha_delta)
    height_expected = np.zeros_like(combined.height_delta)
    for command in commands:
        single = characterize_response(
            config,
            (command,),
            locked_cells=locked_cells,
            tolerance=tolerance,
        )
        alpha_expected += single.alpha_delta
        height_expected += single.height_delta
    return (
        float(np.max(np.abs(combined.alpha_delta - alpha_expected))),
        float(np.max(np.abs(combined.height_delta - height_expected))),
    )


def _cell_dict(cell: tuple[int, int]) -> dict[str, int]:
    return {"row": int(cell[0]), "col": int(cell[1])}


def _command_dict(command: SourceCommand, index: int | None = None) -> dict[str, object]:
    payload: dict[str, object] = {
        "row": int(command.cell[0]),
        "col": int(command.cell[1]),
        "alpha": float(command.alpha),
        "z": float(command.z),
    }
    if index is not None:
        payload["index"] = index
    return payload


def _event_dict(
    event: ProgrammableDiscontinuityEvent, index: int | None = None
) -> dict[str, object]:
    payload: dict[str, object] = {
        "kind": event.kind,
        "cell": None if event.cell is None else _cell_dict(event.cell),
        "alpha": float(event.alpha),
        "z": float(event.z),
    }
    if index is not None:
        payload["index"] = index
    return payload


def _config_dict(config: LatticeConfig | None) -> dict[str, object] | None:
    if config is None:
        return None
    return {
        "rows": config.rows,
        "cols": config.cols,
        "cellSize": config.cell_size,
        "initialAlpha": config.initial_alpha,
        "backlash": config.backlash,
        "couplingGain": config.coupling_gain,
        "zCouplingGain": config.z_coupling_gain,
        "pinRadius": config.pin_radius,
        "holeRadius": config.hole_radius,
        "pinHoleClearance": config.pin_hole_clearance,
        "maxCouplingSteps": config.max_coupling_steps,
        "alphaMin": config.alpha_min,
        "alphaMax": config.alpha_max,
    }


def _response_dict(
    response: ResponseCharacterization,
    *,
    include_fields: bool,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "commands": [
            _command_dict(command, index) for index, command in enumerate(response.commands)
        ],
        "alphaReach": response.alpha_reach,
        "zReach": response.z_reach,
        "effectiveAlphaDieOff": response.effective_alpha_die_off,
        "effectiveZDieOff": response.effective_z_die_off,
        "maxAbsAlphaDelta": response.max_abs_alpha_delta,
        "maxAbsHeightDelta": response.max_abs_height_delta,
    }
    if include_fields:
        payload["fields"] = {
            "alphaDelta": response.alpha_delta.tolist(),
            "heightDelta": response.height_delta.tolist(),
            "actuatorInfluence": response.actuator_influence.tolist(),
            "zResidual": response.z_residual.tolist(),
            "alphaDieOff": response.alpha_die_off.tolist(),
            "zDieOff": response.z_die_off.tolist(),
        }
    return payload


def _decay_dict(decay: ResponseDecayProfile | None) -> dict[str, object] | None:
    if decay is None:
        return None
    return {
        "model": decay.model,
        "alphaRatio": decay.alpha_ratio,
        "zRatio": decay.z_ratio,
        "alphaLength": decay.alpha_length,
        "zLength": decay.z_length,
        "alphaShells": decay.alpha_shells,
        "zShells": decay.z_shells,
        "alphaReach": decay.alpha_reach,
        "zReach": decay.z_reach,
        "alphaFirst": decay.alpha_first,
        "zFirst": decay.z_first,
        "alphaLast": decay.alpha_last,
        "zLast": decay.z_last,
    }


def _state_distance_dict(distance) -> dict[str, object]:
    return {
        "modeChanges": distance.mode_changes,
        "commandAlphaError": distance.command_alpha_error,
        "commandZError": distance.command_z_error,
        "commandError": distance.command_error,
        "alphaGridError": distance.alpha_grid_error,
        "finalAlphaError": distance.final_alpha_error,
        "finalHeightError": distance.final_height_error,
        "finalError": distance.final_error,
    }


def _sequence_order_dict(
    sequence_order: SequenceOrderDiagnostic | None,
) -> dict[str, object] | None:
    if sequence_order is None:
        return None
    return {
        "eventCount": sequence_order.event_count,
        "adjacentPairCount": sequence_order.adjacent_pair_count,
        "noncommutingAdjacentPairs": sequence_order.noncommuting_adjacent_pairs,
        "orderSensitive": sequence_order.order_sensitive,
        "maxAdjacentAlphaError": sequence_order.max_adjacent_alpha_error,
        "maxAdjacentHeightError": sequence_order.max_adjacent_height_error,
        "maxAdjacentCommandError": sequence_order.max_adjacent_command_error,
        "maxAdjacentAlphaGridError": sequence_order.max_adjacent_alpha_grid_error,
        "maxOrderError": sequence_order.max_order_error,
        "reverse": _state_distance_dict(sequence_order.reverse),
        "adjacent": [
            {
                "index": adjacent.index,
                "firstKind": adjacent.first_kind,
                "secondKind": adjacent.second_kind,
                "sensitive": adjacent.sensitive,
                "distance": _state_distance_dict(adjacent.distance),
            }
            for adjacent in sequence_order.adjacent
        ],
    }


def _law_candidate(
    identifier: str,
    operator_class: str,
    property_name: str,
    statement: str,
    supported: bool,
    evidence: dict[str, object],
) -> dict[str, object]:
    return {
        "id": identifier,
        "operatorClass": operator_class,
        "property": property_name,
        "statement": statement,
        "supportedByDiagnostic": bool(supported),
        "status": "simulator-diagnostic",
        "evidence": evidence,
    }


def _framework_operator_law_candidates(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
) -> dict[str, object]:
    max_superposition_error = max(
        diagnostic.alpha_superposition_error,
        diagnostic.height_superposition_error,
    )
    total_cells = diagnostic.total_cells
    locality_supported = (
        diagnostic.alpha_locality_radius < total_cells
        and diagnostic.z_locality_radius < total_cells
    )
    rank_limited = (
        diagnostic.alpha_underactuated_cells > 0
        or diagnostic.height_underactuated_cells > 0
    )
    laws = [
        _law_candidate(
            "bounded_locality",
            "actuation and clearance operators",
            "locality",
            (
                "The measured response is bounded by finite alpha and height "
                "die-off radii at the report tolerance."
            ),
            locality_supported,
            {
                "alphaLocalityRadius": diagnostic.alpha_locality_radius,
                "zLocalityRadius": diagnostic.z_locality_radius,
                "alphaDecayRatio": diagnostic.alpha_decay_ratio,
                "zDecayRatio": diagnostic.z_decay_ratio,
                "tolerance": diagnostic.tolerance,
            },
        ),
        _law_candidate(
            "rank_limited_reachability",
            "finite actuator set",
            "reachable set",
            (
                "The selected actuator operators span a finite response subspace; "
                "unreached cells mark underactuated regions for the current command basis."
            ),
            rank_limited,
            {
                "alphaRank": diagnostic.alpha_rank,
                "heightRank": diagnostic.height_rank,
                "reachableAlphaCells": diagnostic.reachable_alpha_cells,
                "reachableHeightCells": diagnostic.reachable_height_cells,
                "alphaUnderactuatedCells": diagnostic.alpha_underactuated_cells,
                "heightUnderactuatedCells": diagnostic.height_underactuated_cells,
                "totalCells": total_cells,
            },
        ),
        _law_candidate(
            "composition_nonadditivity",
            "actuation composition",
            "nonadditivity",
            (
                "Operator composition is non-additive when the combined response "
                "exceeds the sum of individual responses by more than tolerance."
            ),
            diagnostic.nonadditive,
            {
                "alphaSuperpositionError": diagnostic.alpha_superposition_error,
                "heightSuperpositionError": diagnostic.height_superposition_error,
                "maxSuperpositionError": max_superposition_error,
                "tolerance": diagnostic.tolerance,
            },
        ),
        _law_candidate(
            "event_order_noncommutativity",
            "lock and actuation sequence",
            "noncommutativity",
            (
                "Lock, release, and actuation events are noncommutative when "
                "reversal or adjacent swaps change the final state by more than tolerance."
            ),
            diagnostic.order_sensitive,
            {
                "eventCount": len(diagnostic.event_sequence),
                "noncommutingAdjacentPairs": diagnostic.noncommuting_adjacent_pairs,
                "maxOrderError": diagnostic.max_order_error,
                "tolerance": diagnostic.tolerance,
            },
        ),
    ]
    return {
        "schema": "rad-sim.framework-law-candidates.v1",
        "method": "thresholded diagnostic predicates over locality, reachability, composition, and event-order metrics",
        "laws": laws,
    }


def _formalization_target(
    identifier: str,
    statement: str,
    source: str,
    status: str,
    ready_for_lean: bool,
    dependencies: tuple[str, ...],
    evidence: dict[str, object],
) -> dict[str, object]:
    return {
        "id": identifier,
        "statement": statement,
        "source": source,
        "status": status,
        "readyForLean": bool(ready_for_lean),
        "dependencies": list(dependencies),
        "evidence": evidence,
    }


def _lean_tooling_status() -> dict[str, object]:
    lean_path = shutil.which("lean")
    lake_path = shutil.which("lake")
    available = bool(lean_path and lake_path)
    return {
        "engine": "Lean",
        "leanPath": lean_path,
        "lakePath": lake_path,
        "available": available,
        "status": "available" if available else "not-found-on-path",
    }


def _framework_formalization_targets(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
) -> dict[str, object]:
    tooling = _lean_tooling_status()
    lean_ready = bool(tooling["available"])
    basic_status = "ready-for-lean" if lean_ready else "pending-lean-tooling"
    witness_status = (
        basic_status if diagnostic.order_sensitive else "pending-numeric-witness"
    )
    targets = [
        _formalization_target(
            "dead_zone_zero_inside_backlash",
            (
                "For b >= 0, f_b(x)=max(0,x-b)+min(x+b,0) equals 0 "
                "whenever -b <= x <= b."
            ),
            "paper-supported",
            basic_status,
            lean_ready,
            ("real max/min lemmas", "nonnegative backlash premise"),
            {"formula": "f(x)=max(0,x-b)+min(x+b,0)"},
        ),
        _formalization_target(
            "dead_zone_piecewise_linear_outside_gap",
            (
                "For b >= 0, f_b(x)=x-b when x >= b and f_b(x)=x+b "
                "when x <= -b."
            ),
            "paper-supported",
            basic_status,
            lean_ready,
            ("real max/min lemmas", "case split on backlash thresholds"),
            {"formula": "f(x)=max(0,x-b)+min(x+b,0)"},
        ),
        _formalization_target(
            "lock_projection_idempotent",
            "Applying the same lock projection twice is equivalent to applying it once.",
            "simulator-operator",
            basic_status,
            lean_ready,
            ("finite grid state model", "lock projection definition"),
            {"lockedCellCount": len(diagnostic.locked_cells)},
        ),
        _formalization_target(
            "finite_response_rank_bound",
            "The rank of a finite response matrix is bounded by its command-column count.",
            "linear-algebra-diagnostic",
            basic_status,
            lean_ready,
            ("finite matrix rank theorem", "response matrix column count"),
            {
                "alphaRank": diagnostic.alpha_rank,
                "heightRank": diagnostic.height_rank,
                "commandCount": len(diagnostic.response_matrix.commands),
            },
        ),
        _formalization_target(
            "noncommutativity_witness_from_order_error",
            (
                "If sequence-order distance is greater than tolerance, the "
                "corresponding event compositions are not equal."
            ),
            "simulator-diagnostic",
            witness_status,
            lean_ready and diagnostic.order_sensitive,
            ("state distance definition", "event composition semantics"),
            {
                "orderSensitive": diagnostic.order_sensitive,
                "maxOrderError": diagnostic.max_order_error,
                "tolerance": diagnostic.tolerance,
            },
        ),
        _formalization_target(
            "bounded_locality_witness",
            (
                "If all response magnitudes outside a reported die-off radius "
                "are below tolerance, the diagnostic has a finite locality witness."
            ),
            "simulator-diagnostic",
            "requires-calibrated-premise",
            False,
            ("normed response field", "thresholded locality definition"),
            {
                "alphaLocalityRadius": diagnostic.alpha_locality_radius,
                "zLocalityRadius": diagnostic.z_locality_radius,
                "tolerance": diagnostic.tolerance,
            },
        ),
    ]
    return {
        "schema": "rad-sim.formalization-targets.v1",
        "method": (
            "candidate theorem manifest; no Lean proof is emitted until Lean/Lake "
            "are available and the premises are first-principles enough to formalize"
        ),
        "tooling": tooling,
        "targets": targets,
    }


def _pairwise_interactions_dict(
    pairwise: OperatorInteractionGraph | None,
    *,
    include_fields: bool,
) -> dict[str, object] | None:
    if pairwise is None:
        return None
    payload: dict[str, object] = {
        "commandCount": len(pairwise.commands),
        "totalPairCount": pairwise.total_pair_count,
        "evaluatedPairCount": pairwise.evaluated_pair_count,
        "nonadditivePairCount": pairwise.nonadditive_pair_count,
        "maxAlphaError": pairwise.max_alpha_error,
        "maxHeightError": pairwise.max_height_error,
        "maxInteractionError": pairwise.max_interaction_error,
        "maxHotspotError": pairwise.max_hotspot_error,
        "maxInteractionDegree": pairwise.max_interaction_degree,
        "interactionDensity": pairwise.interaction_density,
        "truncated": pairwise.truncated,
        "tolerance": pairwise.tolerance,
        "interactions": [
            {
                "firstIndex": interaction.first_index,
                "secondIndex": interaction.second_index,
                "first": _command_dict(interaction.first),
                "second": _command_dict(interaction.second),
                "manhattanDistance": interaction.manhattan_distance,
                "alphaSuperpositionError": interaction.alpha_superposition_error,
                "heightSuperpositionError": interaction.height_superposition_error,
                "maxError": interaction.max_error,
                "nonadditive": interaction.nonadditive,
            }
            for interaction in pairwise.interactions
        ],
    }
    if include_fields:
        payload["fields"] = {
            "alphaErrorMatrix": pairwise.alpha_error_matrix.tolist(),
            "heightErrorMatrix": pairwise.height_error_matrix.tolist(),
            "interactionHotspotMap": pairwise.interaction_hotspot_map.tolist(),
            "interactionDegreeMap": pairwise.interaction_degree_map.tolist(),
        }
    return payload


def _response_matrix_dict(
    matrix: ResponseMatrix,
    tolerance: float,
    include_response_matrix: bool,
) -> dict[str, object]:
    payload = matrix.to_dict(tolerance=tolerance)
    if include_response_matrix:
        return payload
    return {
        "schema": payload["schema"],
        "grid": payload["grid"],
        "commands": payload["commands"],
        "diagnostics": payload["diagnostics"],
    }


def _physical_validation_dict(
    comparison: PhysicalResponseComparison | None,
    *,
    include_fields: bool,
) -> dict[str, object] | None:
    if comparison is None:
        return None
    payload: dict[str, object] = {
        "schema": "rad-sim.physical-response-comparison.v1",
        "model": comparison.physical_result.metadata.get("model", "spring_hinge_3d"),
        "commands": [
            _command_dict(command, index) for index, command in enumerate(comparison.commands)
        ],
        "lockedCells": [_cell_dict(cell) for cell in comparison.locked_cells],
        "physicalSuccess": comparison.physical_success,
        "physicalEnergy": comparison.physical_energy,
        "iterations": int(comparison.physical_result.metadata.get("iterations", 0)),
        "springEdges": int(comparison.physical_result.metadata.get("spring_edges", 0)),
        "hingeTriples": int(comparison.physical_result.metadata.get("hinge_triples", 0)),
        "alphaRmsModelError": comparison.alpha_rms_error,
        "heightRmsModelError": comparison.height_rms_error,
        "centerRmsModelError": comparison.center_rms_error,
        "maxAbsHeightModelError": comparison.max_abs_height_error,
        "maxAbsCenterModelError": comparison.max_abs_center_error,
        "maxAbsPhysicalHeightDelta": float(np.max(np.abs(comparison.physical_height_delta))),
        "maxAbsPhysicalAlphaDelta": float(np.max(np.abs(comparison.physical_alpha_delta))),
    }
    if include_fields:
        height_model_error = comparison.physical_height_delta - comparison.kinematic.height_delta
        alpha_model_error = comparison.physical_alpha_delta - comparison.kinematic.alpha_delta
        center_norm = np.linalg.norm(comparison.physical_center_delta, axis=2)
        payload["fields"] = {
            "physicalAlphaDelta": comparison.physical_alpha_delta.tolist(),
            "physicalHeightDelta": comparison.physical_height_delta.tolist(),
            "physicalCenterDelta": comparison.physical_center_delta.tolist(),
            "physicalCenterDeltaNorm": center_norm.tolist(),
            "alphaModelError": alpha_model_error.tolist(),
            "heightModelError": height_model_error.tolist(),
        }
    return payload


def diagnose_programmable_discontinuity(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    *,
    locked_cells: Iterable[tuple[int, int]] = (),
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    pairwise_max_pairs: int | None = 64,
    include_physical: bool = False,
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> ProgrammableDiscontinuityDiagnostic:
    """Measure a command set as a programmable-discontinuity operator.

    Locality is measured by dead-zone die-off, reachability by response-matrix
    support/rank, and composition by comparing the combined response to the sum
    of isolated command responses. Pairwise interactions are bounded by
    ``pairwise_max_pairs`` so large command sets remain usable.
    """

    command_tuple = tuple(commands)
    locked_tuple = tuple((int(r), int(c)) for r, c in locked_cells)
    event_tuple = tuple(event_sequence or ())
    cells = tuple(actuator_cells) if actuator_cells is not None else _unique_command_cells(command_tuple)
    combined = characterize_cluster(
        config,
        command_tuple,
        locked_cells=locked_tuple,
        tolerance=tolerance,
    )
    response_matrix = build_response_matrix(
        config,
        actuator_cells=cells,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=locked_tuple,
        tolerance=tolerance,
    )
    alpha_error, height_error = _superposition_errors(
        config,
        command_tuple,
        combined,
        locked_tuple,
        tolerance,
    )
    event_base_state = LatticeState.uniform(config)
    for r, c in locked_tuple:
        event_base_state.locked_mask[r, c] = True
    sequence_order = (
        compare_sequence_order(config, event_base_state, event_tuple, tolerance=tolerance)
        if event_tuple
        else None
    )
    decay_profile = response_decay_profile(combined, tolerance=tolerance)
    pairwise_interactions = characterize_pairwise_interactions(
        config,
        command_tuple,
        locked_cells=locked_tuple,
        tolerance=tolerance,
        max_pairs=pairwise_max_pairs,
    )
    physical_validation = (
        compare_physical_cluster(
            config,
            command_tuple,
            locked_cells=locked_tuple,
            load_case=load_case,
            tolerance=tolerance,
        )
        if include_physical
        else None
    )
    return ProgrammableDiscontinuityDiagnostic(
        commands=command_tuple,
        locked_cells=locked_tuple,
        combined=combined,
        response_matrix=response_matrix,
        alpha_superposition_error=alpha_error,
        height_superposition_error=height_error,
        tolerance=tolerance,
        event_sequence=event_tuple,
        sequence_order=sequence_order,
        decay_profile=decay_profile,
        pairwise_interactions=pairwise_interactions,
        physical_validation=physical_validation,
    )


def programmable_discontinuity_report(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
    config: LatticeConfig | None = None,
    *,
    include_response_matrix: bool = True,
    include_fields: bool = True,
) -> dict[str, object]:
    """Serialize an operator diagnostic as a versioned research artifact.

    The report is intentionally descriptive: it records the paper-supported
    equations used by the simulator and labels the locality, rank,
    superposition, and order-sensitivity measurements as diagnostics introduced
    for the programmable-discontinuity framework.
    """

    return {
        "schema": "rad-sim.programmable-discontinuity-report.v1",
        "grid": {
            "rows": diagnostic.response_matrix.cell_shape[0],
            "cols": diagnostic.response_matrix.cell_shape[1],
            "totalCells": diagnostic.total_cells,
        },
        "config": _config_dict(config),
        "operators": {
            "commands": [
                _command_dict(command, index)
                for index, command in enumerate(diagnostic.commands)
            ],
            "lockedCells": [_cell_dict(cell) for cell in diagnostic.locked_cells],
            "eventSequence": [
                _event_dict(event, index)
                for index, event in enumerate(diagnostic.event_sequence)
            ],
            "activeOperatorCount": diagnostic.active_operator_count,
        },
        "paperSupportedAssumptions": [
            {
                "name": "backlash dead-zone activation",
                "formula": "f(x)=max(0,x-b)+min(x+b,0)",
                "implementation": "rad_sim.coupling.backlash_activation",
            },
            {
                "name": "normalized backlash",
                "formula": "b_norm=b/L",
                "implementation": "LatticeConfig.backlash is dimensionless in v1",
            },
            {
                "name": "rotating-square angle/dilation relation",
                "formula": "theta_degrees=70*alpha-60",
                "implementation": "rad_sim.coupling.alpha_to_theta",
            },
        ],
        "simulatorDiagnostics": [
            {
                "name": "response matrix reachability",
                "interpretation": "finite command columns approximate local reachable alpha/height directions",
            },
            {
                "name": "superposition residual",
                "interpretation": "nonzero residual marks non-additive operator composition caused by thresholds, locks, saturation, or coupling",
            },
            {
                "name": "shellwise locality fit",
                "interpretation": "log-linear shell maxima estimate die-off but are not a constitutive law",
            },
            {
                "name": "event-order sensitivity",
                "interpretation": "reversal and adjacent-swap differences test noncommutativity of lock and actuation operators",
            },
        ],
        "operatorLawCandidates": _framework_operator_law_candidates(diagnostic),
        "formalizationTargets": _framework_formalization_targets(diagnostic),
        "locality": {
            "alphaLocalityRadius": diagnostic.alpha_locality_radius,
            "zLocalityRadius": diagnostic.z_locality_radius,
            "alphaDecayRatio": diagnostic.alpha_decay_ratio,
            "zDecayRatio": diagnostic.z_decay_ratio,
            "alphaDecayLength": diagnostic.alpha_decay_length,
            "zDecayLength": diagnostic.z_decay_length,
            "decayProfile": _decay_dict(diagnostic.decay_profile),
        },
        "reachability": {
            "reachableAlphaCells": diagnostic.reachable_alpha_cells,
            "reachableHeightCells": diagnostic.reachable_height_cells,
            "alphaUnderactuatedCells": diagnostic.alpha_underactuated_cells,
            "heightUnderactuatedCells": diagnostic.height_underactuated_cells,
            "alphaRank": diagnostic.alpha_rank,
            "heightRank": diagnostic.height_rank,
        },
        "composition": {
            "nonadditive": diagnostic.nonadditive,
            "alphaSuperpositionError": diagnostic.alpha_superposition_error,
            "heightSuperpositionError": diagnostic.height_superposition_error,
            "orderSensitive": diagnostic.order_sensitive,
            "noncommutingAdjacentPairs": diagnostic.noncommuting_adjacent_pairs,
            "maxOrderError": diagnostic.max_order_error,
            "nonadditivePairCount": diagnostic.nonadditive_pair_count,
            "maxPairwiseInteractionError": diagnostic.max_pairwise_interaction_error,
            "maxPairwiseHotspotError": diagnostic.max_pairwise_hotspot_error,
            "maxPairwiseInteractionDegree": diagnostic.max_pairwise_interaction_degree,
            "pairwiseInteractionDensity": diagnostic.pairwise_interaction_density,
            "pairwiseInteractionsTruncated": diagnostic.pairwise_interactions_truncated,
        },
        "combinedResponse": _response_dict(
            diagnostic.combined,
            include_fields=include_fields,
        ),
        "responseMatrix": _response_matrix_dict(
            diagnostic.response_matrix,
            diagnostic.tolerance,
            include_response_matrix,
        ),
        "pairwiseInteractions": _pairwise_interactions_dict(
            diagnostic.pairwise_interactions,
            include_fields=include_fields,
        ),
        "physicalValidation": _physical_validation_dict(
            diagnostic.physical_validation,
            include_fields=include_fields,
        ),
        "sequenceOrder": _sequence_order_dict(diagnostic.sequence_order),
        "tolerance": diagnostic.tolerance,
    }


def export_programmable_discontinuity_report_json(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
    config: LatticeConfig | None = None,
    *,
    include_response_matrix: bool = True,
    include_fields: bool = True,
) -> str:
    return json.dumps(
        programmable_discontinuity_report(
            diagnostic,
            config=config,
            include_response_matrix=include_response_matrix,
            include_fields=include_fields,
        ),
        indent=2,
    )
