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
        "positiveZReachCells": response.positive_z_reach,
        "negativeZReachCells": response.negative_z_reach,
        "maxPositiveHeightDelta": response.max_positive_height_delta,
        "maxNegativeHeightDelta": response.max_negative_height_delta,
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
            "group_operator_support_decomposition",
            (
                "A list-supported group operator acts only on the declared "
                "support, empty group support is neutral, and disjoint group "
                "supports commute in the abstract mode model."
            ),
            "simulator-operator",
            "lean-proved-discrete",
            True,
            (
                "LocalModeOperator.groupOperator",
                "list membership",
                "disjoint support premise",
            ),
            {
                "leanTheorems": [
                    "LocalModeOperator.groupOperator_support_inside",
                    "LocalModeOperator.groupOperator_support_outside",
                    "LocalModeOperator.groupOperator_empty_apply",
                    "LocalModeOperator.groupOperator_append_support_left",
                    "LocalModeOperator.groupOperator_append_support_right",
                    "LocalModeOperator.groupOperators_with_disjoint_lists_commute",
                ],
                "activeOperatorCount": diagnostic.active_operator_count,
            },
        ),
        _formalization_target(
            "removed_cell_clears_group_supported_constraints",
            (
                "If a constraint touches a removed cell, cell-removal "
                "constraint deletion makes that constraint inactive even after "
                "a group operator is applied."
            ),
            "simulator-operator",
            "lean-proved-discrete",
            True,
            (
                "CellGraph.removeCellConstraints",
                "CellConstraintMap.touches",
                "LocalModeOperator.groupOperator",
            ),
            {
                "leanTheorems": [
                    "CellGraph.removed_cell_constraints_clear_group_operator",
                    "CellGraph.removed_cell_constraints_not_active_after_group_operator",
                ],
                "activeOperatorCount": diagnostic.active_operator_count,
            },
        ),
        _formalization_target(
            "vertical_clearance_gates_residual_contact",
            (
                "A vertical command whose magnitude is inside pin-hole clearance "
                "has zero discrete residual transmission, and the clearance-excess "
                "contact penalty is nonnegative."
            ),
            "contact-mechanics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.clearanceExcessNat",
                "Mechanics.verticalResidualStepNat",
                "Mechanics.contactPenaltyFromClearanceNat",
            ),
            {
                "leanTheorems": [
                    "clearanceExcessNat_zero_inside",
                    "verticalResidualStepNat_zero_inside",
                    "contactPenaltyFromClearanceNat_nonnegative",
                ],
                "formula": "0.5*k*max(|z|-clearance,0)^2 in simulator diagnostics",
                "claimLimit": "discrete theorem target plus uncalibrated normalized simulator penalty",
            },
        ),
        _formalization_target(
            "fixed_cell_load_work_proxy_zero",
            (
                "The finite load-work magnitude proxy is nonnegative, and a "
                "fixed cell with zero displacement contributes zero load work."
            ),
            "load-mechanics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.loadWorkMagnitudeNat",
                "fixed-cell zero-displacement premise",
            ),
            {
                "leanTheorems": [
                    "loadWorkMagnitudeNat_nonnegative",
                    "loadWorkMagnitudeNat_zero_fixed",
                ],
                "claimLimit": "finite magnitude theorem plus uncalibrated simulator load-work proxy",
            },
        ),
        _formalization_target(
            "signed_vertical_load_work_residual_int",
            (
                "The signed finite load-work scaffold records the simulator "
                "sign convention: zero load or zero displacement gives zero "
                "signed work, and signed measured-versus-simulated residuals "
                "vanish when the signed quantities are equal."
            ),
            "signed-load-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.signedLoadWorkInt",
                "Mechanics.signedWorkResidualInt",
                "vertical-load comparison report schema",
            ),
            {
                "leanTheorems": [
                    "signedLoadWorkInt_zero_displacement",
                    "signedLoadWorkInt_zero_load",
                    "absoluteErrorInt_self",
                    "signedWorkResidualInt_zero_when_equal",
                    "signedEnergyResidualTripleInt_nonnegative",
                    "signedEnergyResidualTripleInt_zero_when_equal",
                ],
                "pythonFunction": "write_vertical_load_energy_comparison_artifacts",
                "cliModule": "rad_sim.compare_vertical_load_measurements",
                "reportSchema": "rad-sim.vertical-load-energy-comparison-report.v1",
                "claimLimit": "signed discrete validation scaffold; external load sign and hardware work remain calibration assumptions",
            },
        ),
        _formalization_target(
            "integer_scaled_mechanics_scaffold",
            (
                "Integer-scaled mechanics quantities carry an explicit "
                "denominator while preserving proved numerator facts: scaled "
                "spring/contact energies are nonnegative, zero displacement or "
                "penetration gives zero numerator, and scaled signed residual "
                "triples vanish when measured and simulated signed quantities "
                "agree."
            ),
            "integer-scaled-mechanics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ScaledNatQuantity",
                "Mechanics.ScaledIntQuantity",
                "future rational/real mechanics targets",
            ),
            {
                "leanTheorems": [
                    "scaledSpringEnergyNat_numerator_nonnegative",
                    "scaledSpringEnergyNat_zero_displacement",
                    "scaledSpringEnergyNat_preserves_denominator",
                    "scaledContactPenaltyNat_numerator_nonnegative",
                    "scaledContactPenaltyNat_zero_penetration",
                    "scaledSignedLoadWorkInt_zero_displacement",
                    "scaledSignedLoadWorkInt_zero_load",
                    "scaledSignedLoadWorkInt_preserves_denominator",
                    "scaledSignedWorkResidualInt_numerator_nonnegative",
                    "scaledSignedWorkResidualInt_zero_when_equal",
                    "scaledSignedEnergyResidualTripleInt_numerator_nonnegative",
                    "scaledSignedEnergyResidualTripleInt_zero_when_equal",
                ],
                "claimLimit": "integer-scaled numerator/denominator scaffold; not yet a field-valued rational or real mechanics proof",
                "nextFormalStep": "replace denominator-carrying records with Rat/Real semantics once algebra tooling is available",
            },
        ),
        _formalization_target(
            "measurement_unit_scale_invariants",
            (
                "A finite measurement-unit scale maps normalized simulator "
                "values into denominator-carrying physical-unit records while "
                "preserving zero commands, zero equality residuals, and "
                "denominator metadata for both nonnegative and signed "
                "quantities."
            ),
            "calibration-unit-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.MeasurementUnitScale",
                "calibrate_paper_rad_config",
                "bench artifact unit metadata",
            ),
            {
                "leanTheorems": [
                    "measurementUnitScaleNat_zero",
                    "measurementUnitScaleNat_preserves_denominator",
                    "measurementUnitScaleInt_zero",
                    "measurementUnitScaleInt_preserves_denominator",
                    "measurementUnitScaleResidualNat_zero_when_equal",
                    "measurementUnitScaleResidualNat_preserves_denominator",
                    "measurementUnitScaleResidualInt_zero_when_equal",
                    "measurementUnitScaleResidualInt_preserves_denominator",
                    "hardwareProfileCoverageComplete_true_when_equal",
                    "hardwareProfileCoverageMissing_zero_when_complete",
                    "hardwareProfileCoverageMissing_preserves_total",
                ],
                "unitScaleFunction": "physical_unit_scale_metadata",
                "unitScaleSchema": "rad-sim.physical-unit-scale-metadata.v1",
                "hardwareProfileSchema": "rad-sim.hardware-profile.v1",
                "hardwareProfileFunctions": [
                    "RADHardwareProfile.to_dict",
                    "hardware_profile_from_json",
                    "export_hardware_profile_json",
                ],
                "pythonFunction": "calibrate_paper_rad_config",
                "benchPacketFunction": "write_vertical_load_bench_packet_artifacts",
                "comparisonFunction": "write_vertical_load_energy_comparison_artifacts",
                "claimLimit": "finite unit-scale metadata invariant; not a calibrated physical-units proof",
                "nextFormalStep": "connect MeasurementUnitScale to Rat/Real unit conversion and measured hardware profiles",
            },
        ),
        _formalization_target(
            "mechanics_energy_certificate_nonnegative_proxy",
            (
                "The mechanical energy certificate separates a nonnegative "
                "stored/proxy total from signed external-work terms in the "
                "spring-hinge physical preview."
            ),
            "variational-mechanics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.mechanicalStoredEnergyNat",
                "spring-hinge solver energy breakdown",
                "mechanics energy certificate schema",
            ),
            {
                "leanTheorems": [
                    "mechanicalStoredEnergyNat_nonnegative",
                    "mechanicalStoredEnergyNat_zero_components",
                    "springEnergyNat_nonnegative",
                    "hingeEnergyNat_nonnegative",
                    "contactPenaltyNat_nonnegative",
                    "loadWorkMagnitudeNat_nonnegative",
                ],
                "pythonFunction": "mechanics_energy_certificate_to_dict",
                "schema": "rad-sim.mechanics-energy-certificate.v1",
                "claimLimit": "nonnegative proxy scaffold; signed work and contact calibration remain physical assumptions",
            },
        ),
        _formalization_target(
            "calibration_parameter_estimate_residual_bookkeeping",
            (
                "A calibration parameter-estimate report records finite "
                "sample counts, raw residuals, fitted residuals, and "
                "tolerances so zero fitted residuals with at least one sample "
                "satisfy the finite pass predicate."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationFitResidualNat",
                "calibration_experiment_comparison_report",
                "rad-sim.calibration-parameter-estimates.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationFitResidualPassNat_zero",
                    "calibrationFitResidualNat_preserves_sample_count",
                ],
                "pythonFunction": "calibration_experiment_comparison_report",
                "browserFunction": "calibrationParameterEstimates",
                "reportSchema": "rad-sim.calibration-comparison-report.v1",
                "parameterEstimateSchema": "rad-sim.calibration-parameter-estimates.v1",
                "claimLimit": "finite residual bookkeeping only; fitted parameters remain empirical until separately validated",
            },
        ),
        _formalization_target(
            "calibration_model_profile_safe_update_bounds",
            (
                "A calibration-derived model-profile update is marked safe "
                "only when it has at least one finite sample and its proposed "
                "bounded simulator parameter lies between the declared lower "
                "and upper bounds."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationModelProfileUpdateNat",
                "calibration_model_profile_from_report",
                "rad-sim.calibration-model-profile.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationModelProfileUpdateSafeNat_intro",
                    "calibrationModelProfileUpdateSafeNat_bounds",
                ],
                "pythonFunction": "calibration_model_profile_from_report",
                "pythonApplyFunction": "apply_calibration_model_profile",
                "browserFunction": "calibrationModelProfile",
                "browserApplyFunction": "applyCalibrationModelProfile",
                "profileSchema": "rad-sim.calibration-model-profile.v1",
                "applicationSchema": "rad-sim.calibration-model-profile-application.v1",
                "claimLimit": "finite safety gate for bounded simulator updates; not a proof that fitted parameters are physical laws",
            },
        ),
        _formalization_target(
            "calibration_model_profile_selection_predicate",
            (
                "A calibration model-profile candidate is selectable only "
                "when it has an applied update, non-increased missing "
                "observations, and a non-increased finite residual score."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationModelProfileCandidateNat",
                "calibration_model_profile_residual_comparison",
                "rad-sim.calibration-model-profile-selection.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationModelProfileCandidateImprovesNat_intro",
                    "calibrationModelProfileCandidateImprovesNat_applied_updates",
                    "calibrationModelProfileCandidateImprovesNat_residual_nonincrease",
                    "calibrationModelProfileCandidateImprovesNat_missing_nonincrease",
                    "calibrationModelProfileCandidateScoreNat_le_before",
                ],
                "pythonFunction": "select_calibration_model_profile",
                "browserFunction": "selectCalibrationModelProfile",
                "comparisonSchema": "rad-sim.calibration-model-profile-residual-comparison.v1",
                "selectionSchema": "rad-sim.calibration-model-profile-selection.v1",
                "claimLimit": "finite residual-selection rule only; independent validation remains required before claiming physical calibration",
            },
        ),
        _formalization_target(
            "calibration_model_profile_holdout_predicate",
            (
                "A train/holdout calibration profile validation passes only "
                "when a profile has an applied update and both fit and "
                "holdout residual scores and missing-observation counts do "
                "not increase."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationModelProfileHoldoutNat",
                "calibration_model_profile_holdout_validation",
                "rad-sim.calibration-model-profile-holdout-validation.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationModelProfileHoldoutPassNat_intro",
                    "calibrationModelProfileHoldoutPassNat_applied_updates",
                    "calibrationModelProfileHoldoutPassNat_fit_residual_nonincrease",
                    "calibrationModelProfileHoldoutPassNat_holdout_residual_nonincrease",
                    "calibrationModelProfileHoldoutPassNat_missing_nonincrease",
                    "calibrationModelProfileHoldoutScoreNat_le_before",
                ],
                "pythonFunction": "calibration_model_profile_holdout_validation",
                "csvExportFunction": "export_calibration_model_profile_holdout_validation_csv",
                "browserFunction": "calibrationModelProfileHoldoutValidation",
                "schema": "rad-sim.calibration-model-profile-holdout-validation.v1",
                "claimLimit": "finite train/holdout residual gate; still not a proof of contact, friction, material, or actuator physics",
            },
        ),
        _formalization_target(
            "calibration_train_holdout_split_metadata",
            (
                "A proper finite calibration split has positive fit samples, "
                "positive holdout samples, zero overlap, and a frozen profile "
                "before holdout evaluation."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationTrainHoldoutSplitNat",
                "rad-sim.calibration-train-holdout-split.v1",
                "bench data provenance metadata",
            ),
            {
                "leanTheorems": [
                    "calibrationTrainHoldoutSplitReadyNat_intro",
                    "calibrationTrainHoldoutSplitReadyNat_has_fit_samples",
                    "calibrationTrainHoldoutSplitReadyNat_has_holdout_samples",
                    "calibrationTrainHoldoutSplitReadyNat_zero_overlap",
                    "calibrationTrainHoldoutSplitReadyNat_profile_frozen",
                ],
                "pythonFunction": "calibration_model_profile_holdout_validation",
                "csvExportFunction": "export_calibration_model_profile_holdout_validation_csv",
                "browserFunction": "calibrationModelProfileHoldoutValidation",
                "schema": "rad-sim.calibration-train-holdout-split.v1",
                "claimLimit": "finite split metadata predicate; actual file independence and frozen-profile timing remain bench-protocol obligations",
            },
        ),
        _formalization_target(
            "calibration_train_holdout_file_provenance",
            (
                "A documented finite calibration holdout provenance record "
                "has known fit and holdout dataset IDs, distinct dataset and "
                "source-file IDs, marked fit/holdout roles, a frozen profile, "
                "and a holdout profile ID matching the frozen profile."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationTrainHoldoutProvenanceNat",
                "rad-sim.calibration-dataset-provenance.v1",
                "rad-sim.calibration-train-holdout-split.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationTrainHoldoutProvenanceReadyNat_intro",
                    "calibrationTrainHoldoutProvenanceReadyNat_distinct_files",
                    "calibrationTrainHoldoutProvenanceReadyNat_profile_frozen",
                    "calibrationTrainHoldoutProvenanceReadyNat_profile_matches",
                ],
                "pythonFunction": "calibration_model_profile_holdout_validation",
                "browserFunction": "calibrationModelProfileHoldoutValidation",
                "datasetProvenanceSchema": "rad-sim.calibration-dataset-provenance.v1",
                "splitSchema": "rad-sim.calibration-train-holdout-split.v1",
                "claimLimit": "finite provenance predicate; true file provenance and collection timing remain external laboratory obligations",
            },
        ),
        _formalization_target(
            "calibration_bench_protocol_coverage",
            (
                "A calibration bench notebook is coverage-ready only when it "
                "contains at least one scenario, at least one instrument, at "
                "least one output artifact, at least one measurement column, "
                "and both fit and holdout dataset roles."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationBenchProtocolCoverageNat",
                "calibration_bench_notebook",
                "rad-sim.calibration-bench-notebook.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationBenchProtocolCoverageReadyNat_intro",
                    "calibrationBenchProtocolCoverageReadyNat_has_scenarios",
                    "calibrationBenchProtocolCoverageReadyNat_has_outputs",
                    "calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles",
                ],
                "pythonFunction": "calibration_bench_notebook",
                "csvExportFunction": "export_calibration_bench_notebook_csv",
                "browserFunction": "calibrationBenchNotebook",
                "schema": "rad-sim.calibration-bench-notebook.v1",
                "claimLimit": "finite protocol-coverage predicate; actual independent collection, calibration, and physical-law validity remain laboratory obligations",
            },
        ),
        _formalization_target(
            "calibration_bench_packet_completeness",
            (
                "A calibration bench packet is complete only when it contains "
                "notebook artifacts, protocol artifacts, fit and holdout "
                "templates, validation instructions, and a sufficiently "
                "populated artifact manifest."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationBenchPacketCompletenessNat",
                "calibration_bench_packet",
                "rad-sim.calibration-bench-packet.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationBenchPacketCompleteNat_intro",
                    "calibrationBenchPacketCompleteNat_has_notebook",
                    "calibrationBenchPacketCompleteNat_has_fit_template",
                    "calibrationBenchPacketCompleteNat_has_holdout_template",
                    "calibrationBenchPacketCompleteNat_has_manifest",
                ],
                "pythonFunction": "calibration_bench_packet",
                "writerFunction": "write_calibration_bench_packet_artifacts",
                "browserFunction": "calibrationBenchPacket",
                "schema": "rad-sim.calibration-bench-packet.v1",
                "claimLimit": "finite packet-completeness predicate; file independence, bench execution, and physical calibration remain external evidence",
            },
        ),
        _formalization_target(
            "calibration_bench_executed_validation_gate",
            (
                "A filled calibration bench execution is validation-ready only "
                "when fit and holdout samples are present, at least one bounded "
                "model-profile update was applied, residual and independent "
                "validation flags are positive, and no required provenance "
                "evidence is missing."
            ),
            "bench-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.CalibrationBenchExecutedValidationNat",
                "calibration_bench_execution_validation",
                "rad-sim.calibration-bench-execution-validation.v1",
            ),
            {
                "leanTheorems": [
                    "calibrationBenchExecutedValidationReadyNat_intro",
                    "calibrationBenchExecutedValidationReadyNat_has_fit_samples",
                    "calibrationBenchExecutedValidationReadyNat_has_holdout_samples",
                    "calibrationBenchExecutedValidationReadyNat_has_applied_updates",
                    "calibrationBenchExecutedValidationReadyNat_has_independent_validation",
                    "calibrationBenchExecutedValidationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "calibration_bench_execution_validation",
                "writerFunction": "write_calibration_bench_execution_validation_artifacts",
                "browserFunction": "calibrationBenchExecutionValidation",
                "csvExportFunction": "export_calibration_bench_execution_validation_csv",
                "schema": "rad-sim.calibration-bench-execution-validation.v1",
                "claimLimit": "finite execution-evidence predicate; lab procedure truth, measurement accuracy, and physical-law validity remain external evidence",
            },
        ),
        _formalization_target(
            "measured_vertical_load_work_validation_zero_residual",
            (
                "A measured-versus-simulated vertical load-work validation "
                "has zero finite residual when the measured finite work "
                "quantity equals the simulated quantity."
            ),
            "bench-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.measuredWorkResidualNat",
                "vertical-load energy validation schema",
                "load-cell height measurements",
            ),
            {
                "leanTheorems": [
                    "absoluteErrorNat_self",
                    "measuredWorkResidualNat_zero_when_equal",
                    "measuredEnergyResidualTripleNat_zero_when_equal",
                    "loadWorkMagnitudeNat_nonnegative",
                ],
                "pythonFunction": "validate_vertical_load_energy_measurements",
                "templateFunction": "vertical_load_energy_measurement_template",
                "protocolFunction": "vertical_load_energy_experiment_protocol",
                "packetFunction": "vertical_load_bench_packet",
                "compareFunction": "compare_vertical_load_energy_measurement_results",
                "schema": "rad-sim.vertical-load-energy-validation.v1",
                "resultsSchema": "rad-sim.vertical-load-energy-measurement-results.v1",
                "protocolSchema": "rad-sim.vertical-load-energy-experiment-protocol.v1",
                "packetSchema": "rad-sim.vertical-load-bench-packet.v1",
                "reportSchema": "rad-sim.vertical-load-energy-comparison-report.v1",
                "claimLimit": "bench-data comparison scaffold; real contact/load calibration remains experimental",
            },
        ),
        _formalization_target(
            "vertical_load_comparison_pass_predicate",
            (
                "A filled vertical-load scenario with no missing measurements "
                "and zero signed-work, magnitude, and contact-proxy errors "
                "passes the finite comparison predicate for any tolerance; a "
                "scenario with a missing measurement fails that predicate."
            ),
            "bench-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.verticalLoadScenarioPassNat",
                "Mechanics.measuredEnergyResidualTripleNat",
                "vertical-load comparison report schema",
            ),
            {
                "leanTheorems": [
                    "measuredEnergyResidualTripleNat_nonnegative",
                    "measuredEnergyResidualTripleNat_zero_when_equal",
                    "verticalLoadScenarioPassNat_zero_errors",
                    "verticalLoadScenarioPassNat_false_when_missing",
                ],
                "pythonFunction": "write_vertical_load_energy_comparison_artifacts",
                "cliModule": "rad_sim.compare_vertical_load_measurements",
                "summaryCsv": "vertical_load_energy_comparison_summary.csv",
                "reportSchema": "rad-sim.vertical-load-energy-comparison-report.v1",
                "claimLimit": "finite pass/fail predicate; physical accuracy still requires bench calibration",
            },
        ),
        _formalization_target(
            "physical_validation_readiness_gate",
            (
                "A physical validation evidence packet is readiness-complete "
                "only when executed calibration is present, vertical-load "
                "comparison scenarios are present and passing, load/contact "
                "proxy terms are present, clearance is configured, and no "
                "required evidence is missing."
            ),
            "bench-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.PhysicalValidationReadinessNat",
                "physical_validation_readiness_report",
                "rad-sim.physical-validation-readiness.v1",
            ),
            {
                "leanTheorems": [
                    "physicalValidationReadyNat_intro",
                    "physicalValidationReadyNat_has_calibration",
                    "physicalValidationReadyNat_has_vertical_load",
                    "physicalValidationReadyNat_has_load_proxy",
                    "physicalValidationReadyNat_has_contact_proxy",
                    "physicalValidationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "physical_validation_readiness_report",
                "browserFunction": "physicalValidationReadinessReport",
                "csvExportFunction": "export_physical_validation_readiness_csv",
                "schema": "rad-sim.physical-validation-readiness.v1",
                "claimLimit": "finite evidence-readiness predicate; contact, friction, stiffness, gravity, and material-law accuracy remain external physical validation",
            },
        ),
        _formalization_target(
            "contact_state_abstraction_gate",
            (
                "A contact-state abstraction is complete only when active "
                "bodies have matching pin, hole, and clearance records, "
                "contact-state and penalty records cover the active bodies, "
                "clearance is configured, and no abstraction evidence is "
                "missing."
            ),
            "contact-state-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ContactStateAbstractionNat",
                "contact_state_abstraction_report",
                "rad-sim.contact-state-abstraction.v1",
            ),
            {
                "leanTheorems": [
                    "contactStateAbstractionReadyNat_intro",
                    "contactStateAbstractionReadyNat_has_bodies",
                    "contactStateAbstractionReadyNat_pins_match_bodies",
                    "contactStateAbstractionReadyNat_holes_match_bodies",
                    "contactStateAbstractionReadyNat_has_contact_records",
                    "contactStateAbstractionReadyNat_has_penalty_terms",
                    "contactStateAbstractionReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "contact_state_abstraction_report",
                "browserFunction": "contactStateAbstractionReport",
                "csvExportFunction": "export_contact_state_abstraction_csv",
                "schema": "rad-sim.contact-state-abstraction.v1",
                "claimLimit": "finite contact-state bookkeeping predicate; rigid-body contact, friction, and stiffness remain uncalibrated",
            },
        ),
        _formalization_target(
            "contact_graph_consistency_gate",
            (
                "A contact graph is consistency-complete only when removed-cell "
                "incident edges are deleted, removed cells have no active "
                "contact, active bodies have contact records, group-support "
                "cells are represented, and no consistency evidence is missing."
            ),
            "contact-state-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ContactGraphConsistencyNat",
                "contact_graph_consistency_report",
                "rad-sim.contact-graph-consistency.v1",
            ),
            {
                "leanTheorems": [
                    "contactGraphConsistentNat_intro",
                    "contactGraphConsistentNat_has_active_bodies",
                    "contactGraphConsistentNat_contact_records_cover_bodies",
                    "contactGraphConsistentNat_removed_edges_clear",
                    "contactGraphConsistentNat_removed_contacts_clear",
                    "contactGraphConsistentNat_support_records_cover_support",
                    "contactGraphConsistentNat_zero_missing_evidence",
                ],
                "pythonFunction": "contact_graph_consistency_report",
                "browserFunction": "contactGraphConsistencyReport",
                "csvExportFunction": "export_contact_graph_consistency_csv",
                "schema": "rad-sim.contact-graph-consistency.v1",
                "claimLimit": "finite graph/contact/support bookkeeping predicate; physical contact remains uncalibrated",
            },
        ),
        _formalization_target(
            "physical_realization_map_gate",
            (
                "A physical realization map is readiness-complete only when at "
                "least one abstract operator is present, every abstract "
                "operator has a realized simulator state effect, support "
                "record, claim label, contact-graph evidence, and no missing "
                "realization evidence."
            ),
            "physical-realization-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.PhysicalRealizationMapNat",
                "physical_realization_map_report",
                "rad-sim.physical-realization-map.v1",
            ),
            {
                "leanTheorems": [
                    "physicalRealizationMapReadyNat_intro",
                    "physicalRealizationMapReadyNat_has_operators",
                    "physicalRealizationMapReadyNat_realizes_all",
                    "physicalRealizationMapReadyNat_has_support_records",
                    "physicalRealizationMapReadyNat_has_state_effects",
                    "physicalRealizationMapReadyNat_has_claim_labels",
                    "physicalRealizationMapReadyNat_has_contact_graph",
                    "physicalRealizationMapReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "physical_realization_map_report",
                "browserFunction": "physicalRealizationMapReport",
                "csvExportFunction": "export_physical_realization_map_csv",
                "schema": "rad-sim.physical-realization-map.v1",
                "claimLimit": "finite abstract-to-simulator realization predicate; hardware realization remains experimentally unvalidated",
            },
        ),
        _formalization_target(
            "external_physics_engine_audit_gate",
            (
                "An external physics engine audit is readiness-complete only "
                "when an independent rigid-body/contact tool candidate is "
                "available, required features and validation scenarios are "
                "recorded, contact-model evidence exists, and no audit "
                "evidence is missing."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalPhysicsEngineAuditNat",
                "external_physics_engine_audit_report",
                "rad-sim.external-physics-engine-audit.v1",
            ),
            {
                "leanTheorems": [
                    "externalPhysicsEngineAuditReadyNat_intro",
                    "externalPhysicsEngineAuditReadyNat_has_engines",
                    "externalPhysicsEngineAuditReadyNat_has_available",
                    "externalPhysicsEngineAuditReadyNat_has_features",
                    "externalPhysicsEngineAuditReadyNat_has_scenarios",
                    "externalPhysicsEngineAuditReadyNat_has_contact_model",
                    "externalPhysicsEngineAuditReadyNat_has_independent_tool",
                    "externalPhysicsEngineAuditReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "external_physics_engine_audit_report",
                "browserFunction": "externalPhysicsEngineAuditReport",
                "csvExportFunction": "export_external_physics_engine_audit_csv",
                "schema": "rad-sim.external-physics-engine-audit.v1",
                "claimLimit": "finite external-engine readiness predicate; external solver execution and hardware agreement remain unvalidated",
            },
        ),
        _formalization_target(
            "mujoco_model_export_gate",
            (
                "A MuJoCo model export is readiness-complete only when a "
                "coarse external-engine model, active body records, fixed-body "
                "records, gravity records, MJCF XML bytes, and zero missing "
                "export evidence are present."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalPhysicsModelExportNat",
                "mujoco_model_export_report",
                "rad-sim.mujoco-model-export.v1",
            ),
            {
                "leanTheorems": [
                    "externalPhysicsModelExportReadyNat_intro",
                    "externalPhysicsModelExportReadyNat_has_model",
                    "externalPhysicsModelExportReadyNat_has_bodies",
                    "externalPhysicsModelExportReadyNat_has_fixed_bodies",
                    "externalPhysicsModelExportReadyNat_has_gravity",
                    "externalPhysicsModelExportReadyNat_has_xml",
                    "externalPhysicsModelExportReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "mujoco_model_export_report",
                "browserFunction": "mujocoModelExportReport",
                "xmlExportFunction": "export_mujoco_model_xml",
                "schema": "rad-sim.mujoco-model-export.v1",
                "claimLimit": "finite MJCF export-completeness predicate; exported geometry is a coarse proxy, not validated CAD/contact mechanics",
            },
        ),
        _formalization_target(
            "mujoco_pin_hole_contact_geometry_gate",
            (
                "A MuJoCo pin-hole contact-geometry inventory is "
                "readiness-complete only when pin, hole, clearance, contact "
                "pair, active contact-pair, MJCF geometry, and zero missing "
                "geometry evidence records are present."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalContactGeometryNat",
                "mujoco_pin_hole_contact_geometry_report",
                "rad-sim.mujoco-pin-hole-contact-geometry.v1",
            ),
            {
                "leanTheorems": [
                    "externalContactGeometryReadyNat_intro",
                    "externalContactGeometryReadyNat_has_pins",
                    "externalContactGeometryReadyNat_holes_cover_pins",
                    "externalContactGeometryReadyNat_clearance_covers_pins",
                    "externalContactGeometryReadyNat_pairs_cover_pins",
                    "externalContactGeometryReadyNat_active_pairs_cover_pins",
                    "externalContactGeometryReadyNat_has_xml",
                    "externalContactGeometryReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "mujoco_pin_hole_contact_geometry_report",
                "browserFunction": "mujocoPinHoleContactGeometryReport",
                "csvExportFunction": "export_mujoco_pin_hole_contact_geometry_csv",
                "schema": "rad-sim.mujoco-pin-hole-contact-geometry.v1",
                "claimLimit": "finite pin-hole contact-geometry inventory; proxy geometry, friction, compliance, and hardware contact remain unvalidated",
            },
        ),
        _formalization_target(
            "mujoco_contact_parameter_profile_gate",
            (
                "A MuJoCo contact-parameter profile is readiness-complete "
                "only when contact pairs have parameter, friction, solver, "
                "stiffness, damping, XML-attribute, and zero missing evidence "
                "records."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalContactParameterProfileNat",
                "mujoco_contact_parameter_report",
                "rad-sim.mujoco-contact-parameter-profile.v1",
            ),
            {
                "leanTheorems": [
                    "externalContactParameterProfileReadyNat_intro",
                    "externalContactParameterProfileReadyNat_has_pairs",
                    "externalContactParameterProfileReadyNat_parameters_cover_pairs",
                    "externalContactParameterProfileReadyNat_friction_covers_pairs",
                    "externalContactParameterProfileReadyNat_solver_covers_pairs",
                    "externalContactParameterProfileReadyNat_stiffness_covers_pairs",
                    "externalContactParameterProfileReadyNat_damping_covers_pairs",
                    "externalContactParameterProfileReadyNat_has_xml_attributes",
                    "externalContactParameterProfileReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "mujoco_contact_parameter_report",
                "browserFunction": "mujocoContactParameterReport",
                "csvExportFunction": "export_mujoco_contact_parameter_csv",
                "schema": "rad-sim.mujoco-contact-parameter-profile.v1",
                "claimLimit": "finite contact-parameter completeness predicate; values remain uncalibrated unless backed by bench measurements",
            },
        ),
        _formalization_target(
            "contact_parameter_calibration_packet_completeness",
            (
                "A contact-parameter calibration packet is readiness-complete "
                "only when contact-pair records, parameter records, measurement "
                "columns, fit and holdout template rows, artifact-manifest "
                "entries, attached profile evidence, and zero missing evidence "
                "are present."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ContactParameterCalibrationPacketNat",
                "contact_parameter_calibration_packet",
                "rad-sim.contact-parameter-calibration-packet.v1",
            ),
            {
                "leanTheorems": [
                    "contactParameterCalibrationPacketCompleteNat_intro",
                    "contactParameterCalibrationPacketCompleteNat_has_pairs",
                    "contactParameterCalibrationPacketCompleteNat_parameters_cover_pairs",
                    "contactParameterCalibrationPacketCompleteNat_has_measurement_columns",
                    "contactParameterCalibrationPacketCompleteNat_fit_rows_cover_pairs",
                    "contactParameterCalibrationPacketCompleteNat_holdout_rows_cover_pairs",
                    "contactParameterCalibrationPacketCompleteNat_has_manifest",
                    "contactParameterCalibrationPacketCompleteNat_has_profile_evidence",
                    "contactParameterCalibrationPacketCompleteNat_zero_missing_evidence",
                ],
                "pythonFunction": "contact_parameter_calibration_packet",
                "browserFunction": "contactParameterCalibrationPacket",
                "csvExportFunction": "export_contact_parameter_calibration_packet_csv",
                "schema": "rad-sim.contact-parameter-calibration-packet.v1",
                "claimLimit": "finite contact-calibration packet completeness predicate; fitted contact physics remains external bench evidence",
            },
        ),
        _formalization_target(
            "contact_parameter_bench_validation_gate",
            (
                "A contact-parameter bench validation is readiness-complete "
                "only when a ready packet, fit rows, holdout rows, completed "
                "measurements, residual checks, fit pass, holdout pass, "
                "independent holdout pass, and zero missing evidence are present."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ContactParameterBenchValidationNat",
                "compare_contact_parameter_calibration_results",
                "rad-sim.contact-parameter-bench-validation.v1",
            ),
            {
                "leanTheorems": [
                    "contactParameterBenchValidationReadyNat_intro",
                    "contactParameterBenchValidationReadyNat_has_packet",
                    "contactParameterBenchValidationReadyNat_has_fit_rows",
                    "contactParameterBenchValidationReadyNat_has_holdout_rows",
                    "contactParameterBenchValidationReadyNat_measurements_cover_rows",
                    "contactParameterBenchValidationReadyNat_has_residuals",
                    "contactParameterBenchValidationReadyNat_has_fit_pass",
                    "contactParameterBenchValidationReadyNat_has_holdout_pass",
                    "contactParameterBenchValidationReadyNat_has_independent_holdout",
                    "contactParameterBenchValidationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "compare_contact_parameter_calibration_results",
                "browserFunction": "compareContactParameterCalibrationResults",
                "csvExportFunction": "export_contact_parameter_bench_validation_csv",
                "schema": "rad-sim.contact-parameter-bench-validation.v1",
                "resultsSchema": "rad-sim.contact-parameter-calibration-results.v1",
                "claimLimit": "finite contact-parameter bench comparison predicate; physical contact laws remain empirical and require independent measurements",
            },
        ),
        _formalization_target(
            "contact_parameter_interval_calibration_gate",
            (
                "A contact-parameter interval calibration is readiness-complete "
                "only when a passing bench validation, parameter intervals, "
                "accepted interval records, uncertainty records, holdout agreement "
                "records, simulator parameters inside bounds, and zero missing "
                "evidence are present."
            ),
            "bench-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ContactParameterIntervalCalibrationNat",
                "contact_parameter_interval_calibration_report",
                "rad-sim.contact-parameter-interval-calibration.v1",
            ),
            {
                "leanTheorems": [
                    "contactParameterIntervalCalibrationReadyNat_intro",
                    "contactParameterIntervalCalibrationReadyNat_has_validation",
                    "contactParameterIntervalCalibrationReadyNat_has_intervals",
                    "contactParameterIntervalCalibrationReadyNat_accepts_all_intervals",
                    "contactParameterIntervalCalibrationReadyNat_has_uncertainty",
                    "contactParameterIntervalCalibrationReadyNat_has_holdout_agreement",
                    "contactParameterIntervalCalibrationReadyNat_parameters_inside_bounds",
                    "contactParameterIntervalCalibrationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "contact_parameter_interval_calibration_report",
                "browserFunction": "contactParameterIntervalCalibrationReport",
                "csvExportFunction": "export_contact_parameter_interval_calibration_csv",
                "schema": "rad-sim.contact-parameter-interval-calibration.v1",
                "upstreamSchema": "rad-sim.contact-parameter-bench-validation.v1",
                "claimLimit": "finite empirical interval-calibration predicate; interval acceptance is not a constitutive contact proof",
            },
        ),
        _formalization_target(
            "mujoco_external_run_gate",
            (
                "A MuJoCo external run is readiness-complete only when the "
                "engine is available, a ready model export is attached, body "
                "records are present, result records cover body records, "
                "requested steps are positive and completed, and no run "
                "evidence is missing."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalPhysicsRunNat",
                "mujoco_external_run_report",
                "rad-sim.mujoco-external-run.v1",
            ),
            {
                "leanTheorems": [
                    "externalPhysicsRunReadyNat_intro",
                    "externalPhysicsRunReadyNat_has_engine",
                    "externalPhysicsRunReadyNat_has_export",
                    "externalPhysicsRunReadyNat_has_bodies",
                    "externalPhysicsRunReadyNat_results_cover_bodies",
                    "externalPhysicsRunReadyNat_has_steps",
                    "externalPhysicsRunReadyNat_steps_completed",
                    "externalPhysicsRunReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "mujoco_external_run_report",
                "browserFunction": "mujocoExternalRunReport",
                "schema": "rad-sim.mujoco-external-run.v1",
                "claimLimit": "finite external-run completeness predicate; external solver output still needs simulator and bench comparison",
            },
        ),
        _formalization_target(
            "mujoco_external_comparison_gate",
            (
                "A MuJoCo external comparison is readiness-complete only when "
                "a complete external run is attached, simulator records are "
                "present, external and matched records cover simulator "
                "records, a tolerance record exists, the comparison is within "
                "tolerance, and no comparison evidence is missing."
            ),
            "external-physics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ExternalPhysicsComparisonNat",
                "mujoco_external_comparison_report",
                "rad-sim.mujoco-external-comparison.v1",
            ),
            {
                "leanTheorems": [
                    "externalPhysicsComparisonReadyNat_intro",
                    "externalPhysicsComparisonReadyNat_has_run",
                    "externalPhysicsComparisonReadyNat_has_simulator_records",
                    "externalPhysicsComparisonReadyNat_external_covers_simulator",
                    "externalPhysicsComparisonReadyNat_matched_covers_simulator",
                    "externalPhysicsComparisonReadyNat_has_tolerance",
                    "externalPhysicsComparisonReadyNat_within_tolerance",
                    "externalPhysicsComparisonReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "mujoco_external_comparison_report",
                "browserFunction": "mujocoExternalComparisonReport",
                "csvExportFunction": "export_mujoco_external_comparison_csv",
                "schema": "rad-sim.mujoco-external-comparison.v1",
                "claimLimit": "finite simulator-vs-external comparison predicate; bench validation remains external evidence",
            },
        ),
        _formalization_target(
            "equilibrium_relation_gate",
            (
                "An equilibrium relation is readiness-complete only when "
                "physical realization evidence is ready, solver success and "
                "contact-state evidence are present, energy terms are covered "
                "by nonnegative terms, the residual is within tolerance, and "
                "no equilibrium evidence is missing."
            ),
            "variational-mechanics-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.EquilibriumRelationNat",
                "equilibrium_relation_report",
                "rad-sim.equilibrium-relation.v1",
            ),
            {
                "leanTheorems": [
                    "equilibriumRelationReadyNat_intro",
                    "equilibriumRelationReadyNat_has_realization",
                    "equilibriumRelationReadyNat_has_solver",
                    "equilibriumRelationReadyNat_has_contact",
                    "equilibriumRelationReadyNat_has_energy_terms",
                    "equilibriumRelationReadyNat_energy_terms_nonnegative",
                    "equilibriumRelationReadyNat_residual_within_tolerance",
                    "equilibriumRelationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "equilibrium_relation_report",
                "browserFunction": "equilibriumRelationReport",
                "csvExportFunction": "export_equilibrium_relation_csv",
                "schema": "rad-sim.equilibrium-relation.v1",
                "claimLimit": "finite equilibrium evidence predicate; continuous minimization and physical material/contact accuracy remain unproved",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_controllability_gate",
            (
                "A reachable-equilibrium controllability report is readiness-"
                "complete only when equilibrium evidence is ready, an actuator "
                "basis and response columns exist, target cells and reachable "
                "responses are recorded, topology evidence is present, optional "
                "full-target reachability policy is satisfied, and no evidence "
                "is missing."
            ),
            "reachability-controllability-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumControllabilityNat",
                "reachable_equilibrium_controllability_report",
                "rad-sim.reachable-equilibrium-controllability.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumControllabilityReadyNat_intro",
                    "reachableEquilibriumControllabilityReadyNat_has_equilibrium",
                    "reachableEquilibriumControllabilityReadyNat_has_actuator_basis",
                    "reachableEquilibriumControllabilityReadyNat_has_response_columns",
                    "reachableEquilibriumControllabilityReadyNat_has_target_cells",
                    "reachableEquilibriumControllabilityReadyNat_has_reachable_responses",
                    "reachableEquilibriumControllabilityReadyNat_has_topology",
                    "reachableEquilibriumControllabilityReadyNat_target_policy",
                    "reachableEquilibriumControllabilityReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_controllability_report",
                "browserFunction": "reachableEquilibriumControllabilityReport",
                "csvExportFunction": "export_reachable_equilibrium_controllability_csv",
                "schema": "rad-sim.reachable-equilibrium-controllability.v1",
                "claimLimit": "finite response-matrix reachability predicate; nonlinear controllability and hardware feasibility remain unproved",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_bench_protocol_gate",
            (
                "A reachable-equilibrium bench protocol is readiness-complete "
                "only when reachable-equilibrium evidence is attached, protocol "
                "steps and actuator-column trials exist, target observation "
                "cells and measurement columns are recorded, topology-blocked "
                "target policy is represented, pass/fail criteria are present, "
                "and no protocol evidence is missing."
            ),
            "physical-experiment-protocol-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumBenchProtocolNat",
                "reachable_equilibrium_bench_protocol",
                "rad-sim.reachable-equilibrium-bench-protocol.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumBenchProtocolReadyNat_intro",
                    "reachableEquilibriumBenchProtocolReadyNat_has_reachability",
                    "reachableEquilibriumBenchProtocolReadyNat_has_steps",
                    "reachableEquilibriumBenchProtocolReadyNat_has_actuator_trials",
                    "reachableEquilibriumBenchProtocolReadyNat_has_targets",
                    "reachableEquilibriumBenchProtocolReadyNat_has_measurement_columns",
                    "reachableEquilibriumBenchProtocolReadyNat_has_topology_policy",
                    "reachableEquilibriumBenchProtocolReadyNat_has_pass_fail",
                    "reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_bench_protocol",
                "browserFunction": "reachableEquilibriumBenchProtocol",
                "csvExportFunction": "export_reachable_equilibrium_bench_protocol_csv",
                "schema": "rad-sim.reachable-equilibrium-bench-protocol.v1",
                "claimLimit": "bench protocol coverage predicate; physical measurements and nonlinear controllability remain unproved",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_bench_validation_gate",
            (
                "A reachable-equilibrium bench validation is readiness-complete "
                "only when the protocol is ready, result rows and target "
                "measurements exist, completed measurements cover target "
                "measurements, reachability, topology-leakage, and group-"
                "sequence checks are present, the comparison passes, and no "
                "validation evidence is missing."
            ),
            "physical-experiment-validation-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumBenchValidationNat",
                "compare_reachable_equilibrium_bench_results",
                "rad-sim.reachable-equilibrium-bench-comparison.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumBenchValidationReadyNat_intro",
                    "reachableEquilibriumBenchValidationReadyNat_has_protocol",
                    "reachableEquilibriumBenchValidationReadyNat_has_rows",
                    "reachableEquilibriumBenchValidationReadyNat_has_targets",
                    "reachableEquilibriumBenchValidationReadyNat_measurements_complete",
                    "reachableEquilibriumBenchValidationReadyNat_has_reachability_checks",
                    "reachableEquilibriumBenchValidationReadyNat_has_topology_checks",
                    "reachableEquilibriumBenchValidationReadyNat_has_group_checks",
                    "reachableEquilibriumBenchValidationReadyNat_has_pass",
                    "reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "compare_reachable_equilibrium_bench_results",
                "templateFunction": "reachable_equilibrium_bench_results_template",
                "browserFunction": "compareReachableEquilibriumBenchResults",
                "csvExportFunction": "export_reachable_equilibrium_bench_comparison_csv",
                "schema": "rad-sim.reachable-equilibrium-bench-comparison.v1",
                "resultsSchema": "rad-sim.reachable-equilibrium-bench-results.v1",
                "claimLimit": "finite filled-result validation predicate; physical laws and nonlinear controllability remain unproved",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_amplitude_calibration_gate",
            (
                "A reachable-equilibrium amplitude calibration is readiness-"
                "complete only when bench validation is ready, repeated-trial "
                "groups, amplitude estimates, residual-field cells, uncertainty "
                "bands, topology-leakage bands, group-sequence residuals, and "
                "zero missing evidence are present."
            ),
            "physical-amplitude-calibration-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumAmplitudeCalibrationNat",
                "reachable_equilibrium_amplitude_calibration_report",
                "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_intro",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_validation",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_repeated_trials",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_amplitudes",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_residual_field",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_uncertainty",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_topology_bands",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_has_group_residuals",
                    "reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_amplitude_calibration_report",
                "browserFunction": "reachableEquilibriumAmplitudeCalibrationReport",
                "csvExportFunction": "export_reachable_equilibrium_amplitude_calibration_csv",
                "schema": "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
                "claimLimit": "finite repeated-trial amplitude artifact; statistical model, physical laws, and nonlinear controllability remain unproved",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_empirical_profile_gate",
            (
                "A reachable-equilibrium empirical profile is readiness-complete "
                "only when amplitude calibration evidence is ready, bounded and "
                "safe update proposals exist, tolerance and uncertainty records "
                "are present, a holdout hook is represented, and no profile "
                "evidence is missing."
            ),
            "physical-empirical-profile-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumEmpiricalProfileNat",
                "reachable_equilibrium_empirical_profile_from_amplitude",
                "rad-sim.reachable-equilibrium-empirical-profile.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumEmpiricalProfileReadyNat_intro",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_amplitude",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_bounded_proposals",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_safe_proposals",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_tolerances",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_uncertainty",
                    "reachableEquilibriumEmpiricalProfileReadyNat_has_holdout_hooks",
                    "reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_empirical_profile_from_amplitude",
                "browserFunction": "reachableEquilibriumEmpiricalProfileFromAmplitude",
                "csvExportFunction": "export_reachable_equilibrium_empirical_profile_csv",
                "schema": "rad-sim.reachable-equilibrium-empirical-profile.v1",
                "claimLimit": "bounded empirical profile proposal artifact; no v1 LatticeConfig mutation or physical-law proof",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_profile_inverse_gate",
            (
                "A profile-aware inverse diagnostic is readiness-complete only "
                "when the empirical profile is ready, safe proposals and target "
                "cells exist, the inverse solve succeeded, residual and score "
                "records are present, the profile is used read-only, and no "
                "evidence is missing."
            ),
            "profile-aware-inverse-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumProfileInverseNat",
                "reachable_equilibrium_profile_inverse_report",
                "rad-sim.reachable-equilibrium-profile-inverse.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumProfileInverseReadyNat_intro",
                    "reachableEquilibriumProfileInverseReadyNat_has_profile",
                    "reachableEquilibriumProfileInverseReadyNat_has_safe_proposals",
                    "reachableEquilibriumProfileInverseReadyNat_has_targets",
                    "reachableEquilibriumProfileInverseReadyNat_has_solve",
                    "reachableEquilibriumProfileInverseReadyNat_has_residuals",
                    "reachableEquilibriumProfileInverseReadyNat_has_scores",
                    "reachableEquilibriumProfileInverseReadyNat_read_only_profile_use",
                    "reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_profile_inverse_report",
                "browserFunction": "reachableEquilibriumProfileInverseReport",
                "csvExportFunction": "export_reachable_equilibrium_profile_inverse_csv",
                "schema": "rad-sim.reachable-equilibrium-profile-inverse.v1",
                "claimLimit": "profile-aware inverse residual certificate; no optimizer mutation, physical-law proof, or fabrication calibration",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_profile_inverse_acceptance_gate",
            (
                "A profile-aware inverse acceptance decision is ready only when "
                "the profile-inverse report is ready, score records exist, the "
                "scaled residual score, band-failure count, and actuator count "
                "fit declared limits, the plan is accepted for preview, and no "
                "evidence is missing."
            ),
            "profile-aware-inverse-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumProfileInverseAcceptanceNat",
                "reachable_equilibrium_profile_inverse_acceptance_report",
                "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_intro",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_has_profile_inverse",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_has_score_records",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_score_within_limit",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_band_failures_within_limit",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_actuators_within_limit",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_accepted",
                    "reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_profile_inverse_acceptance_report",
                "browserFunction": "reachableEquilibriumProfileInverseAcceptanceReport",
                "csvExportFunction": "export_reachable_equilibrium_profile_inverse_acceptance_csv",
                "schema": "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
                "claimLimit": "thresholded preview acceptance gate; not a hardware execution authorization or nonlinear controllability proof",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_profile_inverse_preview_packet_gate",
            (
                "An accepted profile-aware inverse preview packet is ready only "
                "when acceptance evidence is ready, command records exist, "
                "preview event records cover commands, target and residual "
                "records exist, the packet is marked read-only, and no evidence "
                "is missing."
            ),
            "profile-aware-inverse-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumProfileInversePreviewPacketNat",
                "reachable_equilibrium_profile_inverse_preview_packet",
                "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_intro",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_acceptance",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_commands",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_events_cover_commands",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_targets",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_residuals",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_read_only",
                    "reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_profile_inverse_preview_packet",
                "browserFunction": "reachableEquilibriumProfileInversePreviewPacket",
                "csvExportFunction": "export_reachable_equilibrium_profile_inverse_preview_packet_csv",
                "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
                "claimLimit": "read-only preview/lab handoff artifact; not hardware execution authorization",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_profile_inverse_preview_replay_gate",
            (
                "A profile-aware inverse preview replay is ready only when the "
                "packet is ready, command records replay exactly, preview events "
                "cover the commands, simulation and target-residual records "
                "exist, residuals agree with packet metadata, replay is read-"
                "only, and no evidence is missing."
            ),
            "profile-aware-inverse-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumProfileInversePreviewReplayNat",
                "reachable_equilibrium_profile_inverse_preview_replay_report",
                "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_intro",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_packet",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_commands",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_replays_all_commands",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_event_coverage",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_simulation",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_target_residuals",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_residual_agreement",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_read_only",
                    "reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_profile_inverse_preview_replay_report",
                "browserFunction": "reachableEquilibriumProfileInversePreviewReplayReport",
                "csvExportFunction": "export_reachable_equilibrium_profile_inverse_preview_replay_csv",
                "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
                "claimLimit": "deterministic kinematic replay certificate; not physical contact validation or hardware execution",
            },
        ),
        _formalization_target(
            "reachable_equilibrium_profile_inverse_preview_physical_gate",
            (
                "A profile-aware inverse preview physical check is ready only "
                "when the replay certificate is ready, the physical-preview "
                "solver succeeds, command and target-residual records exist, "
                "model-comparison and energy records exist, the check is read-"
                "only, and no evidence is missing."
            ),
            "profile-aware-inverse-scaffold",
            "lean-proved-discrete",
            True,
            (
                "Mechanics.ReachableEquilibriumProfileInversePreviewPhysicalNat",
                "reachable_equilibrium_profile_inverse_preview_physical_report",
                "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
            ),
            {
                "leanTheorems": [
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_intro",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_replay",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_solver",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_commands",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_targets",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_model_comparison",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_energy",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_read_only",
                    "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence",
                ],
                "pythonFunction": "reachable_equilibrium_profile_inverse_preview_physical_report",
                "browserFunction": "reachableEquilibriumProfileInversePreviewPhysicalReport",
                "csvExportFunction": "export_reachable_equilibrium_profile_inverse_preview_physical_csv",
                "schema": "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
                "claimLimit": "normalized spring-hinge/browser spring-preview check; not a calibrated contact, gravity, or hardware-execution proof",
            },
        ),
        _formalization_target(
            "spring_hinge_removed_topology_load_comparison",
            (
                "A removed-cell vertical residual/load case can be compared "
                "against the 3D spring-hinge physical-preview solver using the "
                "same deleted spring/hinge topology and fixed-cell load support."
            ),
            "spring-hinge-diagnostic",
            "simulator-diagnostic",
            False,
            (
                "compare_vertical_residual_spring_hinge_3d",
                "solve_spring_hinge_3d",
                "LoadCase.fixed_cells",
            ),
            {
                "pythonFunction": "compare_vertical_residual_spring_hinge_3d",
                "reportFunction": "build_vertical_load_physical_preview_report",
                "csvExportFunction": "export_vertical_load_physical_preview_report_csv",
                "comparisonSchema": "rad-sim.vertical-removal-physical-comparison.v1",
                "reportSchema": "rad-sim.vertical-load-physical-preview-report.v1",
                "solver": "solve_spring_hinge_3d",
                "claimLimit": "numerical physical preview; not a calibrated rigid-body contact proof",
            },
        ),
        _formalization_target(
            "browser_spring_preview_removed_edge_deletion",
            (
                "The browser spring-preview relaxation treats removed cells as "
                "deleted graph nodes: removed cells do not average with active "
                "neighbors, and incident link-strain entries are marked as "
                "absent spring edges."
            ),
            "spring-preview-diagnostic",
            "simulator-diagnostic",
            False,
            (
                "CellGraph.removeCell",
                "web/physics.js simulatePhysicalRelaxation",
                "link-strain absent-edge counts",
            ),
            {
                "leanTheorems": [
                    "CellGraph.removed_cell_not_present",
                    "CellGraph.removal_deletes_one_step_from_removed",
                    "CellGraph.removal_deletes_one_step_to_removed",
                ],
                "browserFunction": "RAD.simulatePhysicalRelaxation",
                "browserMetrics": [
                    "physicalActiveSpringEdges",
                    "physicalSkippedSpringEdges",
                ],
                "validation": "tests/validate_web_modules.js",
                "claimLimit": "browser relaxation topology consistency only; not calibrated rigid-body contact",
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


def formalization_target_manifest(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
) -> dict[str, object]:
    """Return the standalone theorem-target manifest for proof work."""

    return _framework_formalization_targets(diagnostic)


def export_formalization_target_manifest_json(
    diagnostic: ProgrammableDiscontinuityDiagnostic,
) -> str:
    return json.dumps(formalization_target_manifest(diagnostic), indent=2)


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
