import argparse

import numpy as np

from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState
from .visualization import plot_lattice


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the RAD kinematic backlash demo.")
    parser.add_argument("--no-show", action="store_true", help="build the figure without opening a window")
    args = parser.parse_args()

    config = LatticeConfig(rows=7, cols=7, backlash=0.1, coupling_gain=0.55)
    state = LatticeState.uniform(config)
    state.actuator_grid[3, 3] = -0.35
    state.actuator_grid[1, 5] = 0.22
    state.locked_mask[0, :] = True
    state.alpha_grid[0, :] = np.linspace(0.9, 1.1, config.cols)
    result = simulate_kinematic(config, state)
    print("model:", result.metadata["model"])
    print("mean alpha:", round(result.metadata["mean_alpha"], 4))
    plot_lattice(result, show=not args.no_show)


if __name__ == "__main__":
    main()
