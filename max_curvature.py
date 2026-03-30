"""
Recreate: Max Curvature — Variable Cell Size & Backlash

Style matches the Die-off Distance figure:
  - contourf fill (no scatter dots)
  - red dotted lines for the original data curves
  - gray "no-data zone" where curvature = 0 / outside hull
  - viridis colormap
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata


# ──────────────────────────────────────────────────────────────
# 1.  Curvature model  — swap your real formula in here
# ──────────────────────────────────────────────────────────────

def compute_max_curvature(
    normalized_backlash: np.ndarray,
    cell_size_mm: np.ndarray,
) -> np.ndarray:
    """κ [1/mm] = 4 · β_norm / L   (calibrated to match colorbar range)"""
    β = np.asarray(normalized_backlash, dtype=float)
    L = np.asarray(cell_size_mm,        dtype=float)
    return 4.0 * β / L


# ──────────────────────────────────────────────────────────────
# 2.  Sample the sparse design curves
# ──────────────────────────────────────────────────────────────

def sample_design_curves(
    n_curves:    int   = 13,
    abs_bl_min:  float = 0.20,
    abs_bl_max:  float = 2.10,
    cs_min:      float = 10.0,
    cs_max:      float = 50.0,
    n_per_curve: int   = 200,
) -> tuple:
    """
    Returns (bl_norm, cell_size, kappa) flat arrays for every
    sample on every constant-abs-backlash hyperbola.
    """
    abs_bl_vals = np.linspace(abs_bl_min, abs_bl_max, n_curves)
    cs_vals     = np.linspace(cs_min, cs_max, n_per_curve)

    bl_parts, cs_parts, k_parts = [], [], []
    for abs_bl in abs_bl_vals:
        bl = abs_bl / cs_vals
        κ  = compute_max_curvature(bl, cs_vals)
        bl_parts.append(bl)
        cs_parts.append(cs_vals.copy())
        k_parts.append(κ)

    return (
        np.concatenate(bl_parts),
        np.concatenate(cs_parts),
        np.concatenate(k_parts),
    )


# ──────────────────────────────────────────────────────────────
# 3.  Interpolate onto a regular grid
# ──────────────────────────────────────────────────────────────

def build_grid(
    bl_known:    np.ndarray,
    cs_known:    np.ndarray,
    kappa_known: np.ndarray,
    abs_bl_min:  float = 0.20,
    abs_bl_max:  float = 2.10,
    n_grid:      int   = 500,
    method:      str   = "cubic",
) -> tuple:
    """
    Returns (BL, CS, KAPPA) meshgrids.
    Points outside the two bounding hyperbolas are forced to NaN
    so nothing is rendered outside the actual data region.
    """
    bl_grid = np.linspace(bl_known.min(), bl_known.max(), n_grid)
    cs_grid = np.linspace(cs_known.min(), cs_known.max(), n_grid)
    BL, CS  = np.meshgrid(bl_grid, cs_grid)

    KAPPA = griddata(
        points=(bl_known, cs_known),
        values=kappa_known,
        xi=(BL, CS),
        method=method,
        rescale=True,
    )

    # Mask anything outside the left (min) and right (max) bounding curves
    left_bound  = abs_bl_min / CS   # β_norm at the innermost curve
    right_bound = abs_bl_max / CS   # β_norm at the outermost curve
    outside = (BL < left_bound) | (BL > right_bound)
    KAPPA[outside] = np.nan

    return BL, CS, KAPPA


# ──────────────────────────────────────────────────────────────
# 4.  Main plot
# ──────────────────────────────────────────────────────────────

def plot_max_curvature(
    n_curves:      int            = 13,
    n_per_curve:   int            = 200,
    n_grid:        int            = 500,
    n_levels:      int            = 40,
    abs_bl_min:    float          = 0.20,
    abs_bl_max:    float          = 2.10,
    cs_min:        float          = 10.0,
    cs_max:        float          = 50.0,
    interp_method: str            = "cubic",
    cmap:          str            = "viridis",
    figsize:       tuple          = (7.5, 6.5),
    save_path:     Optional[str]  = None,
) -> tuple:
    """
    Full pipeline: sample → interpolate → contourf + red curve overlay.
    Only the region strictly between the bounding hyperbolas is rendered.
    """

    # ── Data ────────────────────────────────────────────────────
    bl_s, cs_s, κ_s = sample_design_curves(
        n_curves=n_curves, n_per_curve=n_per_curve,
        abs_bl_min=abs_bl_min, abs_bl_max=abs_bl_max,
        cs_min=cs_min, cs_max=cs_max,
    )
    BL, CS, KAPPA = build_grid(
        bl_s, cs_s, κ_s,
        abs_bl_min=abs_bl_min, abs_bl_max=abs_bl_max,
        n_grid=n_grid, method=interp_method,
    )

    vmin   = np.nanmin(KAPPA)
    vmax   = np.nanmax(KAPPA)
    levels = np.linspace(vmin, vmax, n_levels)

    # ── Figure ───────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=figsize)

    # Filled contours — NaN cells are simply transparent
    cf = ax.contourf(BL, CS, KAPPA,
                     levels=levels,
                     cmap=cmap,
                     extend="neither")

    # Red dotted lines — one per original design curve
    abs_bl_vals = np.linspace(abs_bl_min, abs_bl_max, n_curves)
    cs_line     = np.linspace(cs_min, cs_max, n_per_curve)
    for abs_bl in abs_bl_vals:
        bl_line = abs_bl / cs_line
        ax.plot(bl_line, cs_line,
                color="red", linestyle="dotted",
                linewidth=1.0, alpha=0.85)

    # ── Colorbar ─────────────────────────────────────────────────
    cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Max Curvature [1/mm]", rotation=270,
                   labelpad=22, fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))

    # ── Axes ─────────────────────────────────────────────────────
    ax.set_xlabel("Normalized Backlash [n.d.]", fontsize=14)
    ax.set_ylabel("Cell Size [mm]", fontsize=14)
    ax.set_title("Max Curvature — Variable Cell Size & Backlash",
                 fontsize=15, fontweight="bold")
    ax.tick_params(axis="both", labelsize=12)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5))

    ax.set_xlim(0, bl_s.max() * 1.01)
    ax.set_ylim(cs_min, cs_max)

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"Saved → {save_path}")

    return fig, ax


# ──────────────────────────────────────────────────────────────
# 5.  Entry point
# ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    fig, ax = plot_max_curvature(
        n_curves=13,
        n_per_curve=200,
        n_grid=500,
        n_levels=40,
        interp_method="cubic",
        cmap="viridis",
        save_path="max_curvature.png",
    )
    plt.show()