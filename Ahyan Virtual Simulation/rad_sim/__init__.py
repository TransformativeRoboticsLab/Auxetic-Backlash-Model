"""RAD backlash lattice simulator."""

from .cell_geometry import (
    PAPER_RAD_REFERENCE,
    PaperRADCellGeometry,
    PaperRADReference,
    RADJointGeometry,
    build_paper_rad_cell_geometry,
)
from .coupling import alpha_to_theta, backlash_activation, theta_to_alpha
from .experiments import (
    PairCharacterization,
    ResponseCharacterization,
    ResponseMatrix,
    SourceCommand,
    build_response_matrix,
    characterize_cluster,
    characterize_pair,
    characterize_response,
    characterize_single_cell,
)
from .kinematic import simulate_kinematic
from .inverse import InverseDesignResult, solve_inverse_design
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .operators import (
    DeadZonePropagationOperator,
    PropagationResult,
    alpha_backlash_operator,
    evaluate_programmable_operators,
    lock_projection,
    vertical_clearance_operator,
)
from .spring_hinge import solve_spring_hinge, solve_spring_hinge_3d
from .visualization import interactive_app, plot_lattice

__all__ = [
    "LatticeConfig",
    "LatticeState",
    "LoadCase",
    "SimulationResult",
    "PaperRADReference",
    "RADJointGeometry",
    "PaperRADCellGeometry",
    "PAPER_RAD_REFERENCE",
    "PropagationResult",
    "DeadZonePropagationOperator",
    "SourceCommand",
    "ResponseCharacterization",
    "PairCharacterization",
    "ResponseMatrix",
    "InverseDesignResult",
    "build_paper_rad_cell_geometry",
    "build_response_matrix",
    "solve_inverse_design",
    "alpha_to_theta",
    "theta_to_alpha",
    "backlash_activation",
    "characterize_response",
    "characterize_single_cell",
    "characterize_pair",
    "characterize_cluster",
    "alpha_backlash_operator",
    "vertical_clearance_operator",
    "lock_projection",
    "evaluate_programmable_operators",
    "simulate_kinematic",
    "solve_spring_hinge",
    "solve_spring_hinge_3d",
    "plot_lattice",
    "interactive_app",
]
