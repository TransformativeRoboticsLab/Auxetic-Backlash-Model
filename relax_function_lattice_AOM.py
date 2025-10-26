import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D projection)
# (uniform edge color requested) remove Line3DCollection coloring
import matplotlib as mpl



# ---------------------- Target surface helpers ---------------------- #

def make_target_surface(Nu=40, Nv=40):
    u = np.linspace(0, 1, Nu)
    v = np.linspace(0, 1, Nv)
    U, V = np.meshgrid(u, v, indexing='ij')
    X = U
    Y = V
    # Removed the additive bilinear term 0.1*(U-0.5)*(V-0.5)
    Z = 0.4 * np.sin(2*np.pi*U) * np.sin(2*np.pi*V)
    return X, Y, Z

def interp_surface_at_uv(X, Y, Z, u, v):

    Nu, Nv = X.shape
    fu = np.clip(u * (Nu - 1), 0.0, Nu - 1)
    fv = np.clip(v * (Nv - 1), 0.0, Nv - 1)
    iu = int(np.floor(fu)); iv = int(np.floor(fv))
    iu1 = min(iu + 1, Nu - 1); iv1 = min(iv + 1, Nv - 1)
    du = fu - iu; dv = fv - iv
    def bilerp(A):
        return ((1-du)*(1-dv)*A[iu, iv] +
                (1-du)*dv    *A[iu, iv1] +
                du    *(1-dv)*A[iu1, iv] +
                du    *dv    *A[iu1, iv1])
    x = bilerp(X); y = bilerp(Y); z = bilerp(Z)
    return np.array([x, y, z], dtype=float)

def surface_normals(X, Y, Z):

    Nu, Nv = X.shape
    # Finite differences in u and v
    def dfdc(A, axis):
        # central differences, forward/backward at borders
        if axis == 0:
            d = np.zeros_like(A)
            d[1:-1, :] = 0.5*(A[2:, :] - A[:-2, :])
            d[0,    :] = A[1, :] - A[0, :]
            d[-1,   :] = A[-1, :] - A[-2, :]
            return d
        else:
            d = np.zeros_like(A)
            d[:, 1:-1] = 0.5*(A[:, 2:] - A[:, :-2])
            d[:, 0   ] = A[:, 1] - A[:, 0]
            d[:, -1  ] = A[:, -1] - A[:, -2]
            return d
    Xu, Xv = dfdc(X,0), dfdc(X,1)
    Yu, Yv = dfdc(Y,0), dfdc(Y,1)
    Zu, Zv = dfdc(Z,0), dfdc(Z,1)
    Tu = np.stack([Xu, Yu, Zu], axis=-1)   # (Nu,Nv,3)
    Tv = np.stack([Xv, Yv, Zv], axis=-1)
    N = np.cross(Tu, Tv)
    nrm = np.linalg.norm(N, axis=-1, keepdims=True)
    nrm = np.where(nrm < 1e-12, 1e-12, nrm)
    N = N / nrm
    return N

# ---------------------- Lattice construction ---------------------- #

def make_initial_lattice(Nu=20, Nv=20, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0)):

    W, H = size
    x0, y0, z0 = origin
    xs = np.linspace(x0, x0 + W, Nu)
    ys = np.linspace(y0, y0 + H, Nv)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = np.full_like(X, z0)
    P = np.stack([X, Y, Z], axis=-1)  # (Nu, Nv, 3)
    return P


def rest_lengths(P0):

    # Horizontal edges: between (i,j) and (i, j+1)
    seg_x = P0[:, 1:, :] - P0[:, :-1, :]
    L0x = np.linalg.norm(seg_x, axis=-1)
    # Vertical edges: between (i,j) and (i+1, j)
    seg_y = P0[1:, :, :] - P0[:-1, :, :]
    L0y = np.linalg.norm(seg_y, axis=-1)
    return L0x, L0y


# ---------------------- Constraint enforcement ---------------------- #

def enforce_lengths_relu_grid(P, L0x, L0y, length_tol=0.1, fix_boundary=True, sweeps=1):

    Nu, Nv, _ = P.shape
    Lmin_x, Lmax_x = (1 - length_tol) * L0x, (1 + length_tol) * L0x
    Lmin_y, Lmax_y = (1 - length_tol) * L0y, (1 + length_tol) * L0y

    for _ in range(sweeps):
        # Accumulate node updates here (avoids double-adding and makes boundary policy easy)
        dP = np.zeros_like(P)

        # ---- Horizontal edges (between (i,j) and (i, j+1)) ----
        seg = P[:, 1:, :] - P[:, :-1, :]
        L = np.linalg.norm(seg, axis=-1)
        mask = L > 1e-12
        u = np.zeros_like(seg)
        u[mask] = seg[mask] / L[mask][..., None]
        target = np.clip(L, Lmin_x, Lmax_x)
        delta = 0.5 * (target - L)
        dvec = (delta[..., None]) * u  # (Nu, Nv-1, 3)

        # Left node (i, j): subtract dvec; skip j=0 if fixed
        j_left_start = 1 if fix_boundary else 0
        if j_left_start < Nv-1:
            dP[:, j_left_start: Nv-1, :] -= dvec[:, j_left_start:, :]

        # Right node (i, j+1): add dvec; skip j+1 = Nv-1 if fixed -> j <= Nv-3
        j_right_end = Nv-2 if fix_boundary else Nv-1
        if 0 < j_right_end:
            dP[:, 1: 1 + j_right_end, :] += dvec[:, 0:j_right_end, :]

        # ---- Vertical edges (between (i,j) and (i+1, j)) ----
        seg = P[1:, :, :] - P[:-1, :, :]
        L = np.linalg.norm(seg, axis=-1)
        mask = L > 1e-12
        u = np.zeros_like(seg)
        u[mask] = seg[mask] / L[mask][..., None]
        target = np.clip(L, Lmin_y, Lmax_y)
        delta = 0.5 * (target - L)
        dvec = (delta[..., None]) * u  # (Nu-1, Nv, 3)

        # Top node (i, j): subtract dvec; skip i=0 if fixed
        i_top_start = 1 if fix_boundary else 0
        if i_top_start < Nu-1:
            dP[i_top_start: Nu-1, :, :] -= dvec[i_top_start:, :, :]

        # Bottom node (i+1, j): add dvec; skip i+1 = Nu-1 if fixed -> i <= Nu-3
        i_bot_end = Nu-2 if fix_boundary else Nu-1
        if 0 < i_bot_end:
            dP[1: 1 + i_bot_end, :, :] += dvec[0:i_bot_end, :, :]

        # Apply accumulated updates
        P += dP

        # No need to "undo" boundary nodes; we've avoided touching them above.

    return P


# ---------------------- Relaxation loop ---------------------- #

def relax_lattice_to_surface(Xt, Yt, Zt,
                             Nu=20, Nv=20,
                             max_iters=1500, step_size=0.01, tol=1e-5,
                             length_tol=0.15, fix_boundary=True,
                             length_method='relu', length_sweeps=2,
                             move_mode='normal',  # 'normal' (along surface normal) or 'full' (full vector)
                             init_size=(1.0, 1.0), init_origin=(0.0, 0.0, 0.0)):

    assert Xt.shape == Yt.shape == Zt.shape
    Nu_t, Nv_t = Xt.shape
    assert Nu_t == Nu and Nv_t == Nv, "Target and lattice must share (Nu, Nv)."

    # Initial flat sheet
    P0 = make_initial_lattice(Nu, Nv, size=init_size, origin=init_origin)
    P  = P0.copy()

    # Rest lengths
    L0x, L0y = rest_lengths(P0)

    # Precompute surface normals (for 'normal' move mode)
    Nsurf = surface_normals(Xt, Yt, Zt)

    for it in range(max_iters):
        Prev = P.copy()

        # Pin boundary before attraction (if desired)
        if fix_boundary:
            P[0, :, :]  = P0[0, :, :]
            P[-1, :, :] = P0[-1, :, :]
            P[:, 0, :]  = P0[:, 0, :]
            P[:, -1, :] = P0[:, -1, :]

        # --- Attraction toward param-matched target surface ---
        i_start = 1 if fix_boundary else 0
        i_stop  = Nu - (1 if fix_boundary else 0)
        j_start = 1 if fix_boundary else 0
        j_stop  = Nv - (1 if fix_boundary else 0)

        for i in range(i_start, i_stop):
            u = i / (Nu - 1 if Nu > 1 else 1)
            for j in range(j_start, j_stop):
                v = j / (Nv - 1 if Nv > 1 else 1)
                t = interp_surface_at_uv(Xt, Yt, Zt, u, v)
                if move_mode == 'normal':
                    n = Nsurf[i, j, :]
                    d = t - P[i, j, :]
                    d_along_n = np.dot(d, n) * n
                    P[i, j, :] += step_size * d_along_n
                else:
                    P[i, j, :] += step_size * (t - P[i, j, :])

        # --- Enforce edge lengths (backlash window) ---
        P = enforce_lengths_relu_grid(P, L0x, L0y, length_tol=length_tol,
                                      fix_boundary=fix_boundary, sweeps=length_sweeps)

        # Convergence check
        avg_move = np.mean(np.linalg.norm(P - Prev, axis=-1))
        if avg_move < tol:
            print(f"Converged in {it+1} iterations, avg_move={avg_move:.3e}")
            break

    return P, L0x, L0y


# ---------------------- Diagnostics ---------------------- #

def edge_lengths_and_alpha(P, L0x, L0y):

    seg_x = P[:, 1:, :] - P[:, :-1, :]
    Lx = np.linalg.norm(seg_x, axis=-1)
    alpha_x = Lx / L0x

    seg_y = P[1:, :, :] - P[:-1, :, :]
    Ly = np.linalg.norm(seg_y, axis=-1)
    alpha_y = Ly / L0y
    return Lx, alpha_x, Ly, alpha_y


# ---------------------- Demo / Example ---------------------- #
if __name__ == "__main__":
    # ==== 15 x 15 visible nodes ====
    Nu, Nv = 15, 15

    # Target surface (wireframe plane)
    Xt, Yt, Zt = make_target_surface(Nu=Nu, Nv=Nv)

    # Relax the lattice onto the target surface, moving along the surface NORMAL only
    P, L0x, L0y = relax_lattice_to_surface(
        Xt, Yt, Zt,
        Nu=Nu, Nv=Nv,
        max_iters=2000,
        step_size=0.02,
        tol=1e-6,
        length_tol=0.15,
        fix_boundary=True,
        length_sweeps=4,
        move_mode='normal',  # <-- only along desired plane's local normal
        init_size=(1.0, 1.0),
        init_origin=(0.0, 0.0, 0.0),
    )

    # Diagnostics
    def edge_lengths_and_alpha(P, L0x, L0y):
        seg_x = P[:, 1:, :] - P[:, :-1, :]
        Lx = np.linalg.norm(seg_x, axis=-1)
        ax = Lx / L0x
        seg_y = P[1:, :, :] - P[:-1, :, :]
        Ly = np.linalg.norm(seg_y, axis=-1)
        ay = Ly / L0y
        return Lx, ax, Ly, ay

    Lx, ax, Ly, ay = edge_lengths_and_alpha(P, L0x, L0y)

    # ------------- Visualization ------------- #
    fig = plt.figure(figsize=(12, 5))

    # (1) 3D: target wireframe + deformed lattice with uniform edge color
ax3d = fig.add_subplot(1, 2, 1, projection='3d')
ax3d.plot_wireframe(Xt, Yt, Zt, rstride=1, cstride=1, linewidth=0.6, alpha=0.5)
# draw grid lines uniformly
edge_color = 'tab:blue'
for i in range(Nu):
    ax3d.plot(P[i, :, 0], P[i, :, 1], P[i, :, 2], '-', linewidth=1.4, alpha=0.95, color=edge_color)
for j in range(Nv):
    ax3d.plot(P[:, j, 0], P[:, j, 1], P[:, j, 2], '-', linewidth=1.4, alpha=0.95, color=edge_color)

# Visible nodes
ax3d.scatter(P[..., 0], P[..., 1], P[..., 2], s=18, c='k', alpha=0.8)

# Fit axes limits
all_x = np.concatenate([Xt.ravel(), P[...,0].ravel()])
all_y = np.concatenate([Yt.ravel(), P[...,1].ravel()])
all_z = np.concatenate([Zt.ravel(), P[...,2].ravel()])
ax3d.set_xlim(all_x.min(), all_x.max())
ax3d.set_ylim(all_y.min(), all_y.max())
ax3d.set_zlim(all_z.min(), all_z.max())
ax3d.set_title("15x15 lattice of cells")

    # (2) Alpha heatmap for horizontal edges (vertical similar if wanted)
ax2 = fig.add_subplot(1, 2, 2)
im1 = ax2.imshow(ax.T, origin='lower', aspect='equal', cmap='viridis')
ax2.set_title("Alpha_x (horizontal edges)")
plt.colorbar(im1, ax=ax2, fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()
