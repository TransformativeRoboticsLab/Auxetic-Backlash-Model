from __future__ import annotations

import numpy as np

from .models import LatticeConfig, LatticeState, SimulationResult


def _draw_cells_3d(ax, corners: np.ndarray, color: str, alpha: float, linewidth: float) -> None:
    rows, cols = corners.shape[:2]
    for r in range(rows):
        for c in range(cols):
            ax.plot(
                corners[r, c, :, 0],
                corners[r, c, :, 1],
                corners[r, c, :, 2],
                color=color,
                alpha=alpha,
                lw=linewidth,
            )


def _set_equal_3d(ax, points: np.ndarray) -> None:
    flat = points.reshape((-1, 3))
    mins = flat.min(axis=0)
    maxs = flat.max(axis=0)
    centers = (mins + maxs) / 2.0
    radius = max(float(np.max(maxs - mins)) / 2.0, 0.5)
    ax.set_xlim(centers[0] - radius, centers[0] + radius)
    ax.set_ylim(centers[1] - radius, centers[1] + radius)
    ax.set_zlim(centers[2] - radius * 0.45, centers[2] + radius * 0.45)


def plot_lattice(result: SimulationResult, show: bool = True):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(14, 9), constrained_layout=True)
    ax0 = fig.add_subplot(2, 2, 1, projection="3d")
    ax1 = fig.add_subplot(2, 2, 2, projection="3d")
    ax2 = fig.add_subplot(2, 2, 3, projection="3d")
    ax3 = fig.add_subplot(2, 2, 4, projection="3d")

    _draw_cells_3d(ax0, result.original_corners_3d, "#8a8a8a", 0.28, 1.0)
    _draw_cells_3d(ax0, result.deformed_corners_3d, "#0f6b8f", 0.95, 1.8)
    ax0.scatter(
        result.deformed_centers_3d[..., 0],
        result.deformed_centers_3d[..., 1],
        result.deformed_centers_3d[..., 2],
        s=18,
        c="#202020",
    )
    ax0.set_title("3D Backlash Lattice")
    _set_equal_3d(ax0, np.concatenate([result.original_centers_3d, result.deformed_centers_3d]))

    surface = ax1.plot_surface(
        result.deformed_centers_3d[..., 0],
        result.deformed_centers_3d[..., 1],
        result.deformed_centers_3d[..., 2],
        facecolors=plt.cm.viridis((result.alpha - result.alpha.min()) / (np.ptp(result.alpha) + 1e-12)),
        rstride=1,
        cstride=1,
        linewidth=0.4,
        edgecolor="#333333",
        alpha=0.9,
    )
    surface.set_array(result.alpha.ravel())
    surface.set_cmap("viridis")
    ax1.set_title("Actuated Surface Height")
    fig.colorbar(surface, ax=ax1, fraction=0.046, label="alpha")

    z0 = result.complex_original.reshape(-1)
    z1 = result.complex_deformed.reshape(-1)
    map_height = np.abs(z1 - z0)
    ax2.scatter(z0.real, z0.imag, np.zeros_like(map_height), s=22, label="z", color="#8a8a8a")
    ax2.scatter(z1.real, z1.imag, map_height, s=28, label="w=f(z)", color="#ba3a3a")
    for a, b, h in zip(z0, z1, map_height, strict=False):
        ax2.plot([a.real, b.real], [a.imag, b.imag], [0.0, h], color="#444444", alpha=0.25, lw=0.8)
    ax2.set_title("3D Complex-Plane Displacement")
    ax2.legend()

    rows = np.arange(result.config.rows)
    cols = np.arange(result.config.cols)
    cc, rr = np.meshgrid(cols, rows)
    ax3.plot_surface(cc, rr, result.theta_degrees, cmap="plasma", edgecolor="#333333", linewidth=0.4)
    ax3.set_title("Cell Rotation Angle")
    ax3.set_xlabel("col")
    ax3.set_ylabel("row")
    ax3.set_zlabel("theta deg")

    for ax in (ax0, ax1, ax2, ax3):
        ax.grid(True, alpha=0.2)
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle(f"RAD {result.metadata.get('model', 'model')} simulator")
    if show:
        plt.show()
    return fig, (ax0, ax1, ax2, ax3)


def interactive_app():
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except Exception:
        from .kinematic import simulate_kinematic

        config = LatticeConfig()
        state = LatticeState.uniform(config)
        state.actuator_grid[config.rows // 2, config.cols // 2] = -0.25
        return plot_lattice(simulate_kinematic(config, state))

    from .kinematic import simulate_kinematic

    rows = widgets.IntSlider(value=5, min=2, max=12, description="rows")
    cols = widgets.IntSlider(value=5, min=2, max=12, description="cols")
    backlash = widgets.FloatSlider(value=0.1, min=0.0, max=0.5, step=0.01, description="backlash")
    act_r = widgets.IntSlider(value=2, min=0, max=11, description="act row")
    act_c = widgets.IntSlider(value=2, min=0, max=11, description="act col")
    command = widgets.FloatSlider(value=-0.25, min=-0.7, max=0.7, step=0.01, description="command")
    lock = widgets.Checkbox(value=False, description="lock actuated cell")
    output = widgets.Output()

    def redraw(*_):
        act_r.max = rows.value - 1
        act_c.max = cols.value - 1
        r = min(act_r.value, rows.value - 1)
        c = min(act_c.value, cols.value - 1)
        config = LatticeConfig(rows=rows.value, cols=cols.value, backlash=backlash.value)
        state = LatticeState.uniform(config)
        state.actuator_grid[r, c] = command.value
        state.locked_mask[r, c] = lock.value
        with output:
            output.clear_output(wait=True)
            plot_lattice(simulate_kinematic(config, state), show=True)

    for control in (rows, cols, backlash, act_r, act_c, command, lock):
        control.observe(redraw, names="value")

    display(widgets.VBox([widgets.HBox([rows, cols, backlash]), widgets.HBox([act_r, act_c, command, lock]), output]))
    redraw()
