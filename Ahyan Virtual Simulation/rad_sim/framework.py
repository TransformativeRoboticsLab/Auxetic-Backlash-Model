from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .experiments import (
    OperatorInteractionGraph,
    ResponseCharacterization,
    ResponseDecayProfile,
    ResponseMatrix,
    SourceCommand,
    build_response_matrix,
    characterize_cluster,
    characterize_pairwise_interactions,
    characterize_response,
    response_decay_profile,
)
from .models import LatticeConfig, LatticeState
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
    )
