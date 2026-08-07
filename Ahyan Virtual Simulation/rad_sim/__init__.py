"""RAD backlash lattice simulator."""

from .coupling import alpha_to_theta, backlash_activation, theta_to_alpha
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .spring_hinge import solve_spring_hinge
from .visualization import interactive_app, plot_lattice

__all__ = [
    "LatticeConfig",
    "LatticeState",
    "LoadCase",
    "SimulationResult",
    "alpha_to_theta",
    "theta_to_alpha",
    "backlash_activation",
    "simulate_kinematic",
    "solve_spring_hinge",
    "plot_lattice",
    "interactive_app",
]
