import argparse

from .models import LatticeConfig, LatticeState, LoadCase
from .spring_hinge import solve_spring_hinge
from .visualization import plot_lattice


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the RAD spring-hinge quasistatic demo.")
    parser.add_argument("--no-show", action="store_true", help="build the figure without opening a window")
    args = parser.parse_args()

    config = LatticeConfig(rows=5, cols=5, backlash=0.1, coupling_gain=0.5)
    state = LatticeState.uniform(config)
    state.actuator_grid[2, 2] = -0.32
    state.actuator_grid[1, 3] = 0.18
    load = LoadCase(fixed_cells=((0, 0), (0, 4)), lock_stiffness=150.0)
    result = solve_spring_hinge(config, state, load)
    print("model:", result.metadata["model"])
    print("success:", result.metadata["success"])
    print("energy:", round(result.metadata["energy"], 6))
    print("iterations:", result.metadata["iterations"])
    plot_lattice(result, show=not args.no_show)


if __name__ == "__main__":
    main()
