# revised_relax_function_3D.py

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ---------------- Core Helpers ---------------- #

def _closest_point_on_segment(a, b, p):
    """Return closest point to p on segment ab in 3D."""
    ab = b - a
    denom = np.dot(ab, ab) + 1e-12
    t = np.dot(p - a, ab) / denom
    t = np.clip(t, 0.0, 1.0)
    return a + t * ab


def closest_point_on_polyline(poly, p):
    """Return closest point on polyline to point p in 3D."""
    best, best_d2 = None, np.inf
    for i in range(len(poly) - 1):
        q = _closest_point_on_segment(poly[i], poly[i + 1], p)
        d2 = np.sum((q - p) ** 2)
        if d2 < best_d2:
            best_d2, best = d2, q
    return best


def get_lengths_and_alpha(deformed_list, L0):
    pts = np.array(deformed_list)
    seg = pts[1:] - pts[:-1]
    L = np.linalg.norm(seg, axis=1)  # absolute lengths
    alpha = L / L0                   # normalized stretch factors
    return L, alpha


def get_LSM_quality(deformed_list, target_curve):
    """Compute mean squared error distance to target curve."""
    pts = np.array(deformed_list)
    tgt = np.array(target_curve)
    errs = []
    for p in pts:
        q = closest_point_on_polyline(tgt, p)
        errs.append(np.sum((q - p) ** 2))
    return float(np.mean(errs))


# ---------------- Relaxation ---------------- #

def relax_to_curve(target_curve, flat_line, max_iters=2000, step_size=0.01, tol=1e-6,
                   k=1.0, length_tol=0.1, lock_x=False, length_method='spring'):
    """
    Relax 3D polyline to match target curve, with link length constraints.
    """
    target = np.array(target_curve)
    flat_line = np.array(flat_line)
    current = flat_line.copy()

    n = len(current)
    L0 = np.linalg.norm(flat_line[-1] - flat_line[0]) / (n - 1)  # rest length
    Lmin = (1 - length_tol) * L0
    Lmax = (1 + length_tol) * L0

    for it in range(max_iters):
        prev = current.copy()

        if lock_x:
            current[1:-1, 0] = flat_line[1:-1, 0]

        # Step 1: move toward target curve
        for i in range(1, n - 1):
            t = closest_point_on_polyline(target, current[i])
            direction = t - current[i]
            current[i] += step_size * direction

        # Step 2: enforce length constraints
        if length_method == 'spring':
            current = _enforce_lengths_spring(current, L0, k)
        elif length_method == 'relu':
            current = _enforce_lengths_relu(current, Lmin, Lmax)
        elif length_method == 'hybrid':
            current = _enforce_lengths_spring(current, L0, k)
            current = _enforce_lengths_relu(current, Lmin, Lmax)

        # Convergence check
        avg_move = np.mean(np.linalg.norm(current - prev, axis=1))
        if avg_move < tol:
            print(f"Converged in {it+1} iterations")
            break

    return [tuple(p) for p in current], L0


def _enforce_lengths_spring(current, L0, k):
    n = len(current)
    forces = np.zeros_like(current)
    for i in range(n - 1):
        seg = current[i + 1] - current[i]
        length = np.linalg.norm(seg)
        if length > 1e-12:
            u = seg / length
            f = k * (length - L0) * u
            if i > 0:
                forces[i] += f
            if i < n - 2:
                forces[i + 1] -= f
    current[1:-1] += forces[1:-1]
    return current


def _enforce_lengths_relu(current, Lmin, Lmax):
    n = len(current)
    for i in range(n - 1):
        seg = current[i + 1] - current[i]
        length = np.linalg.norm(seg)
        if length > 1e-12:
            target_len = np.clip(length, Lmin, Lmax)
            delta = (target_len - length) / 2.0
            u = seg / length
            if i > 0:
                current[i] -= delta * u
            if i < n - 2:
                current[i + 1] += delta * u
            if i == 0:
                current[i + 1] = current[i] + target_len * u
            elif i == n - 2:
                current[i] = current[i + 1] - target_len * u
    return current
def ribbon_from_curve(x, y, z, sheet_width=0.4, n_strips=35, translate=(0.0, 0.0, 0.0), twist_per_width=0.0):
    """
    Build a ribbon 'sheet' by offsetting the curve along an in-plane normal.
    Returns X, Y, Z of shape (n_strips, N) to use with plot_wireframe/plot_surface.
    """
    x = np.asarray(x); y = np.asarray(y); z = np.asarray(z)

    # tangent in XY and corresponding in-plane normal
    tx = np.gradient(x);  ty = np.gradient(y)
    tlen = np.hypot(tx, ty) 
    nx = -ty / tlen;       ny = tx / tlen

    v = np.linspace(-sheet_width/2, sheet_width/2, n_strips)  # lateral offsets
    X = x[None, :] + np.outer(v, nx)
    Y = y[None, :] + np.outer(v, ny)
    Z = z[None, :] + twist_per_width * np.outer(v, np.ones_like(z))

    # translate whole sheet for nicer composition
    X += translate[0]; Y += translate[1]; Z += translate[2]
    return X, Y, Z

# ---------------- Example Usage ---------------- #

if __name__ == "__main__":
    # 3D target curve: simple arc in x-y, with small z twist
    num_points_target = 50
    theta = np.linspace(0, np.pi, num_points_target)
    x_t = np.sin(theta*0.7) 
    y_t = np.sin(theta) * 0.5
    z_t = 0.1 * np.sin(2*theta)  # some 3D variation
    target_curve = list(zip(x_t, y_t, z_t))

    # Flat line in 3D (straight along x-axis)
    num_points = 20
    x_f = np.linspace(0, 1, num_points)
    y_f = np.zeros(num_points)
    z_f = np.zeros(num_points)
    flat_line = list(zip(x_f, y_f, z_f))

    # Relaxation
    deformed, L0 = relax_to_curve(target_curve, flat_line, k=0.3,
                                  length_method='relu', length_tol=0.15)
    L, alpha = get_lengths_and_alpha(deformed, L0)
    mse = get_LSM_quality(deformed, target_curve)

    # --- Visualization --- #
    deformed_np = np.array(deformed)
    target_np = np.array(target_curve)
    flat_np = np.array(flat_line)

    fig = plt.figure(figsize=(10, 5))

    # 3D geometry
    ax3d = fig.add_subplot(121, projection='3d')

    # --- Target sheet (add this) ---
    Xs, Ys, Zs = ribbon_from_curve(
        target_np[:, 0], target_np[:, 1], target_np[:, 2],
        sheet_width=0.4,        # adjust width
        n_strips=35,            # adjust line density
        translate=(0.0, 0.0, 0.0),  # shift if needed
        twist_per_width=0.1     # 0 for flat sheet
    )
    ax3d.plot_wireframe(Xs, Ys, Zs, linewidth=0.6, rstride=1, cstride=2, alpha=0.9)

    # Existing curves
    ax3d.plot(target_np[:, 0], target_np[:, 1], target_np[:, 2], 'b--', label='Target Curve')
    ax3d.plot(flat_np[:, 0], flat_np[:, 1], flat_np[:, 2], 'k:', label='Flat Line')
    ax3d.plot(deformed_np[:, 0], deformed_np[:, 1], deformed_np[:, 2], 'go-', label='Deformed Line')
    ax3d.set_title("3D Geometry")
    ax3d.legend()

    # Alpha & Lengths
    ax2d = fig.add_subplot(122)
    idx = np.arange(len(L))
    ax2d.plot(idx, alpha, 'r.-', label='Alpha (L/L0)')
    ax2d.plot(idx, L, 'g.-', label='Segment Length (abs)')
    ax2d.axhline(L0 * (1 - 0.15), color='k', linestyle='--', alpha=0.6, label='Tolerance Band')
    ax2d.axhline(L0 * (1 + 0.15), color='k', linestyle='--', alpha=0.6)
    ax2d.set_title(f"Stretch/Length Tracking | MSE={mse:.3e}")
    ax2d.set_xlabel("Segment Index")
    ax2d.set_ylabel("Value")
    ax2d.legend()

    plt.tight_layout()
    plt.show()
