import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation


# ════════════════════════════════════════════════════════════════════════════
# HELPERS
# ════════════════════════════════════════════════════════════════════════════

def make_initial_lattice(Nu=15, Nv=15, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    W, H = size
    x0, y0, z0 = origin
    X, Y = np.meshgrid(np.linspace(x0, x0+W, Nu), np.linspace(y0, y0+H, Nv), indexing='ij')
    return np.stack([X, Y, np.full_like(X, z0)], axis=-1)


def rest_lengths(P0):
    L0x = np.linalg.norm(P0[:, 1:] - P0[:, :-1], axis=-1)
    L0y = np.linalg.norm(P0[1:]   - P0[:-1],   axis=-1)
    return L0x, L0y


def compute_spring_forces(P, L0x, L0y, k=50.0):
    F = np.zeros_like(P)
    seg = P[:, 1:] - P[:, :-1]
    L   = np.linalg.norm(seg, axis=-1, keepdims=True).clip(1e-12)
    Fm  = k * (L[..., 0] - L0x)[..., None] * (seg / L)
    F[:, :-1] += Fm;  F[:, 1:] -= Fm
    seg = P[1:]  - P[:-1]
    L   = np.linalg.norm(seg, axis=-1, keepdims=True).clip(1e-12)
    Fm  = k * (L[..., 0] - L0y)[..., None] * (seg / L)
    F[:-1] += Fm;  F[1:] -= Fm
    return F


def _to_zero_based(pairs_1based):
    if not len(pairs_1based):
        return np.empty((0, 2), dtype=int)
    return np.array([(r-1, c-1) for r, c in pairs_1based], dtype=int)


def _in_bounds(idx, Nu, Nv):
    if not idx.size:
        return idx
    return idx[np.all((idx >= 0) & (idx < np.array([Nu, Nv])), axis=1)]


# ════════════════════════════════════════════════════════════════════════════
# SINGLE SIMULATION
# ════════════════════════════════════════════════════════════════════════════

def simulate_dynamic_deformation(
    Nu=15, Nv=15, steps=600, dt=0.01,
    mass=0.05, k_spring=120.0, damping=0.85,
    gravity=np.array([0, 0, -9.81]),
    fix_boundary=True,
    bolt_cells_1based=(),
    servo_cells_1based=(
        (3,3),(3,8),(3,13),(8,3),(8,8),(8,13),(13,3),(13,8),(13,13),
    ),
    servo_profile=lambda t, ni, nj: 0.08 * np.sin(2*np.pi*ni - np.pi*t),
):
    P0 = make_initial_lattice(Nu, Nv)
    P, V = P0.copy(), np.zeros_like(P0)
    L0x, L0y = rest_lengths(P0)
    bolt_idx  = _in_bounds(_to_zero_based(bolt_cells_1based),  Nu, Nv)
    servo_idx = _in_bounds(_to_zero_based(servo_cells_1based), Nu, Nv)
    sni = servo_idx[:,0]/(Nu-1) if servo_idx.size else np.array([])
    snj = servo_idx[:,1]/(Nv-1) if servo_idx.size else np.array([])
    trajectory = []

    for step in range(steps):
        t = step * dt
        F = compute_spring_forces(P, L0x, L0y, k_spring) + mass * gravity
        V += (F/mass)*dt;  V *= damping;  P += V*dt
        if fix_boundary:
            P[0,:,:]=P0[0,:,:]; V[0,:,:]=0; P[-1,:,:]=P0[-1,:,:]; V[-1,:,:]=0
            P[:,0,:]=P0[:,0,:]; V[:,0,:]=0; P[:,-1,:]=P0[:,-1,:]; V[:,-1,:]=0
        if bolt_idx.size:
            P[bolt_idx[:,0],bolt_idx[:,1]] = P0[bolt_idx[:,0],bolt_idx[:,1]]
            V[bolt_idx[:,0],bolt_idx[:,1]] = 0
        if servo_idx.size:
            ii, jj = servo_idx[:,0], servo_idx[:,1]
            P[ii,jj,:2] = P0[ii,jj,:2]
            P[ii,jj, 2] = P0[ii,jj,2] + servo_profile(t, sni, snj)
            V[ii,jj]    = 0
        trajectory.append(P.copy())

    return np.array(trajectory)


# ════════════════════════════════════════════════════════════════════════════
# SHARED VISUALIZATION HELPERS
# ════════════════════════════════════════════════════════════════════════════

def _quads(P, Nu, Nv):
    return [[P[i,j], P[i+1,j], P[i+1,j+1], P[i,j+1]]
            for i in range(Nu-1) for j in range(Nv-1)]

def _fcolors(P, Nu, Nv, cmap, zlo=-0.3, zhi=0.3):
    return [cmap(np.clip((P[i:i+2,j:j+2,2].mean()-zlo)/(zhi-zlo), 0, 1))
            for i in range(Nu-1) for j in range(Nv-1)]

def _draw_surface(ax, P, Nu, Nv, title, cmap, servo_idx, bolt_idx):
    ax.cla()
    poly = Poly3DCollection(_quads(P, Nu, Nv), zsort='average', antialiased=False)
    poly.set_facecolor(_fcolors(P, Nu, Nv, cmap))
    poly.set_edgecolor('0.3');  poly.set_linewidth(0.3);  poly.set_alpha(0.88)
    ax.add_collection3d(poly)
    if servo_idx.size > 0:
        ax.scatter(*P[servo_idx[:,0], servo_idx[:,1]].T, s=45, c='red',    zorder=5, label='Servo')
    if bolt_idx.size > 0:
        ax.scatter(*P[bolt_idx[:,0],  bolt_idx[:,1] ].T, s=45, c='orange', zorder=5, label='Bolt')
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_zlim(-0.3, 0.3)
    ax.view_init(elev=25, azim=35)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')


# ════════════════════════════════════════════════════════════════════════════
# SINGLE ANIMATION (filled polytopes)
# ════════════════════════════════════════════════════════════════════════════

def animate_lattice(trajectory, bolt_cells_1based, servo_cells_1based, interval=30):
    steps, Nu, Nv, _ = trajectory.shape
    servo_idx = _to_zero_based(servo_cells_1based)
    bolt_idx  = _to_zero_based(bolt_cells_1based)
    cmap = plt.get_cmap('coolwarm')

    fig = plt.figure(figsize=(8, 6))
    ax  = fig.add_subplot(111, projection='3d')

    def update(frame):
        _draw_surface(ax, trajectory[frame], Nu, Nv,
                      f"Sine-wave driven lattice  |  t = {frame*0.01:.2f} s",
                      cmap, servo_idx, bolt_idx)
        ax.legend(loc='upper left', fontsize=8)
        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    plt.tight_layout();  plt.show()
    return ani


# ════════════════════════════════════════════════════════════════════════════
# MPC ─ EMPIRICAL COMPLIANCE MATRIX
# ════════════════════════════════════════════════════════════════════════════

def compute_empirical_compliance(
    Nu, Nv, k_spring, mass, gravity, servo_idx,
    fix_boundary=True, n_relax=400, dt_r=0.005, damp_r=0.65, eps=0.02,
):
    """
    Estimate C[n_free, n_servo] = dz_free / du_servo via finite differences.

    Each column s is computed by relaxing the lattice to static equilibrium
    with servo s perturbed by +eps and comparing against the zero-servo
    baseline.  This captures the full geometric nonlinearity of the mesh —
    no analytical stiffness matrix required.
    """
    P0 = make_initial_lattice(Nu, Nv)
    L0x, L0y = rest_lengths(P0)
    ii_s, jj_s = servo_idx[:,0], servo_idx[:,1]
    n_s = len(servo_idx)

    # Node classification
    fixed = np.zeros((Nu, Nv), dtype=bool)
    if fix_boundary:
        fixed[0,:] = fixed[-1,:] = fixed[:,0] = fixed[:,-1] = True
    srv = np.zeros((Nu, Nv), dtype=bool);  srv[ii_s, jj_s] = True
    free_nodes = np.argwhere(~fixed & ~srv)       # (n_free, 2)
    ii_f, jj_f = free_nodes[:,0], free_nodes[:,1]

    def relax(u_vec):
        P = P0.copy();  V = np.zeros_like(P)
        for _ in range(n_relax):
            F = compute_spring_forces(P, L0x, L0y, k_spring) + mass * gravity
            V += (F/mass)*dt_r;  V *= damp_r;  P += V*dt_r
            if fix_boundary:
                P[0,:,:]=P0[0,:,:]; V[0,:,:]=0; P[-1,:,:]=P0[-1,:,:]; V[-1,:,:]=0
                P[:,0,:]=P0[:,0,:]; V[:,0,:]=0; P[:,-1,:]=P0[:,-1,:]; V[:,-1,:]=0
            P[ii_s,jj_s,:2] = P0[ii_s,jj_s,:2]
            P[ii_s,jj_s, 2] = P0[ii_s,jj_s,2] + u_vec
            V[ii_s,jj_s]    = 0
        return P[ii_f, jj_f, 2].copy()

    print("  Building empirical compliance matrix...")
    z_base = relax(np.zeros(n_s))
    C = np.zeros((len(free_nodes), n_s))
    for s in range(n_s):
        u_p = np.zeros(n_s);  u_p[s] = eps
        C[:, s] = (relax(u_p) - z_base) / eps
        print(f"    column {s+1}/{n_s}")
    return C, free_nodes


def build_mpc_gain(C, free_nodes, Nu, Nv, lambda_reg=0.05):
    """
    Precompute MPC feedback gain  L = (HC)^{+λ}  (Tikhonov pseudoinverse).

    Observation operator H selects free nodes that lie on the two
    cross-section planes the controller 'looks at':
        X–Z plane  →  all free nodes at  j = Nv // 2
        Y–Z plane  →  all free nodes at  i = Nu // 2

    L maps an (n_obs,) cross-section z-error vector to a (n_servo,)
    servo correction.
    """
    cs_mask = (free_nodes[:,1] == Nv//2) | (free_nodes[:,0] == Nu//2)
    H  = np.eye(len(free_nodes))[cs_mask]     # (n_obs, n_free)
    HC = H @ C                                 # (n_obs, n_servo)
    n_s = C.shape[1]
    # Solve  (HC^T HC + λI) L = HC^T
    L = np.linalg.solve(HC.T @ HC + lambda_reg * np.eye(n_s), HC.T)
    return L, H, cs_mask


# ════════════════════════════════════════════════════════════════════════════
# COMPARISON SIMULATION — open-loop  vs  MPC closed-loop
# ════════════════════════════════════════════════════════════════════════════

def simulate_comparison(
    Nu=15, Nv=15, steps=600, dt=0.01,
    mass=0.05, k_spring=120.0, damping=0.85,
    gravity=np.array([0, 0, -9.81]),
    fix_boundary=True,
    bolt_cells_1based=(),
    servo_cells_1based=(
        (3,3),(3,8),(3,13),(8,3),(8,8),(8,13),(13,3),(13,8),(13,13),
    ),
    servo_profile=lambda t, ni, nj: 0.05*np.sin(2*np.pi*ni - np.pi*t),
    mpc_lambda=0.05,
    u_clip=0.15,
):
    """
    Run both controllers in lockstep.

    MPC law (feedforward + cross-section feedback):

        u_mpc(t) = u_open(t) + L · [y_target(t) − y_actual(t)]

    where  y = z at cross-section free nodes,
    and    L = Tikhonov pseudoinverse of  H·C  (compliance projected
                                                onto the observation planes).

    The static compliance matrix C is estimated empirically so that the
    linearisation captures geometry rather than assuming a flat membrane.
    """
    P0 = make_initial_lattice(Nu, Nv)
    L0x, L0y = rest_lengths(P0)
    bolt_idx  = _in_bounds(_to_zero_based(bolt_cells_1based),  Nu, Nv)
    servo_idx = _in_bounds(_to_zero_based(servo_cells_1based), Nu, Nv)
    ii_s, jj_s = servo_idx[:,0], servo_idx[:,1]
    sni = servo_idx[:,0]/(Nu-1)
    snj = servo_idx[:,1]/(Nv-1)

    # ── MPC precomputation ──────────────────────────────────────────────── #
    C, free_nodes = compute_empirical_compliance(
        Nu, Nv, k_spring, mass, gravity, servo_idx, fix_boundary
    )
    L_mpc, H_mpc, cs_mask = build_mpc_gain(C, free_nodes, Nu, Nv, lambda_reg=mpc_lambda)
    ii_f, jj_f = free_nodes[:,0], free_nodes[:,1]
    ni_f = free_nodes[:,0]/(Nu-1)
    nj_f = free_nodes[:,1]/(Nv-1)
    print(f"  MPC ready — observing {cs_mask.sum()} cross-section free nodes.")

    P_ol, V_ol = P0.copy(), np.zeros_like(P0)
    P_cl, V_cl = P0.copy(), np.zeros_like(P0)

    def step(P, V, z_off):
        """Single timestep: integrate → apply all constraints."""
        F = compute_spring_forces(P, L0x, L0y, k_spring) + mass * gravity
        V += (F/mass)*dt;  V *= damping;  P += V*dt
        if fix_boundary:
            P[0,:,:]=P0[0,:,:]; V[0,:,:]=0; P[-1,:,:]=P0[-1,:,:]; V[-1,:,:]=0
            P[:,0,:]=P0[:,0,:]; V[:,0,:]=0; P[:,-1,:]=P0[:,-1,:]; V[:,-1,:]=0
        if bolt_idx.size:
            P[bolt_idx[:,0],bolt_idx[:,1]] = P0[bolt_idx[:,0],bolt_idx[:,1]]
            V[bolt_idx[:,0],bolt_idx[:,1]] = 0
        P[ii_s,jj_s,:2] = P0[ii_s,jj_s,:2]
        P[ii_s,jj_s, 2] = P0[ii_s,jj_s,2] + z_off
        V[ii_s,jj_s]    = 0

    traj_ol, traj_cl = [], []
    for s in range(steps):
        t = s * dt
        u_open = servo_profile(t, sni, snj)

        # Open-loop: servos exactly follow the nominal profile
        step(P_ol, V_ol, u_open)

        # MPC: observe current cross-section z, compute correction, then step
        y_act = H_mpc @ P_cl[ii_f, jj_f, 2]
        y_tgt = H_mpc @ servo_profile(t, ni_f, nj_f)
        u_cl  = np.clip(u_open + L_mpc @ (y_tgt - y_act), -u_clip, u_clip)
        step(P_cl, V_cl, u_cl)

        traj_ol.append(P_ol.copy())
        traj_cl.append(P_cl.copy())

    return np.array(traj_ol), np.array(traj_cl), free_nodes, cs_mask


# ════════════════════════════════════════════════════════════════════════════
# COMPARISON ANIMATION  (2 × 2 layout)
# ════════════════════════════════════════════════════════════════════════════

def animate_comparison(
    traj_ol, traj_cl, free_nodes, cs_mask,
    servo_cells_1based, bolt_cells_1based,
    Nu, Nv, servo_profile, dt=0.01, interval=40,
):
    """
    2 × 2 figure:
      [Open-loop 3D surface]  |  [MPC 3D surface]
      [X–Z cross-section]     |  [Y–Z cross-section]

    Cross-section plots show the target sine wave (black dashed),
    open-loop (blue), and MPC (red) at cell-centre slices through
    j = Nv//2  and  i = Nu//2.
    """
    steps     = len(traj_ol)
    servo_idx = _to_zero_based(servo_cells_1based)
    bolt_idx  = _to_zero_based(bolt_cells_1based)
    cmap      = plt.get_cmap('coolwarm')

    xz_j  = Nv // 2;  yz_i  = Nu // 2
    x_arr = np.linspace(0, 1, Nu)
    y_arr = np.linspace(0, 1, Nv)

    fig = plt.figure(figsize=(14, 10))
    ax_ol = fig.add_subplot(2, 2, 1, projection='3d')
    ax_cl = fig.add_subplot(2, 2, 2, projection='3d')
    ax_xz = fig.add_subplot(2, 2, 3)
    ax_yz = fig.add_subplot(2, 2, 4)

    def update(frame):
        t     = frame * dt
        P_ol  = traj_ol[frame]
        P_cl  = traj_cl[frame]

        # ── 3D surfaces ──────────────────────────────────────────────────── #
        _draw_surface(ax_ol, P_ol, Nu, Nv, f"Open-Loop   t = {t:.2f} s",        cmap, servo_idx, bolt_idx)
        _draw_surface(ax_cl, P_cl, Nu, Nv, f"MPC Closed-Loop   t = {t:.2f} s",  cmap, servo_idx, bolt_idx)

        # Target z along each cross-section
        # X–Z: sine varies with x (ni = x_arr), nj fixed at mid-column
        z_tgt_xz = servo_profile(t, x_arr, xz_j/(Nv-1))          # (Nu,)
        # Y–Z: sine is constant in y for a wave travelling along x
        z_tgt_yz = float(servo_profile(t, yz_i/(Nu-1), 0.0))      # scalar

        # ── X–Z cross-section  (all i, j = xz_j) ────────────────────────── #
        ax_xz.cla()
        ax_xz.plot(x_arr, z_tgt_xz,       'k--', lw=1.8, label='Target',    zorder=3)
        ax_xz.plot(x_arr, P_ol[:,xz_j,2], color='steelblue', lw=1.5, label='Open-Loop')
        ax_xz.plot(x_arr, P_cl[:,xz_j,2], color='crimson',   lw=1.5, label='MPC')
        ax_xz.axhline(0, color='0.75', lw=0.8)
        ax_xz.set_xlim(0, 1);  ax_xz.set_ylim(-0.25, 0.25)
        ax_xz.set_xlabel('X  (i-axis)');  ax_xz.set_ylabel('Z displacement (m)')
        ax_xz.set_title(f'X–Z cross-section  (column j = {xz_j})', fontsize=10)
        ax_xz.legend(fontsize=8, loc='upper right');  ax_xz.grid(alpha=0.3)

        # ── Y–Z cross-section  (i = yz_i, all j) ────────────────────────── #
        ax_yz.cla()
        ax_yz.axhline(z_tgt_yz,            color='k',         ls='--', lw=1.8, label='Target', zorder=3)
        ax_yz.plot(y_arr, P_ol[yz_i,:,2],  color='steelblue', lw=1.5, label='Open-Loop')
        ax_yz.plot(y_arr, P_cl[yz_i,:,2],  color='crimson',   lw=1.5, label='MPC')
        ax_yz.axhline(0, color='0.75', lw=0.8)
        ax_yz.set_xlim(0, 1);  ax_yz.set_ylim(-0.25, 0.25)
        ax_yz.set_xlabel('Y  (j-axis)');  ax_yz.set_ylabel('Z displacement (m)')
        ax_yz.set_title(f'Y–Z cross-section  (row i = {yz_i})', fontsize=10)
        ax_yz.legend(fontsize=8, loc='upper right');  ax_yz.grid(alpha=0.3)

        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    fig.suptitle('Surface Shape Control: Open-Loop vs MPC', fontsize=12)
    plt.tight_layout()
    plt.show()
    return ani


# ════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    BOLT_CELLS  = ()
    SERVO_CELLS = (
        (3, 3),  (3, 8),  (3, 13),
        (8, 3),  (8, 8),  (8, 13),
        (13, 3), (13, 8), (13, 13),
    )

    A, omega = 0.05, np.pi
    # DC offset of 0.04 m is intentional: it gives the MPC a persistent
    # gravity-compensation error to correct, making the comparison legible.
    profile = lambda t, ni, nj: 0.04 + A * np.sin(2 * np.pi * ni - omega * t)

    # ── View 1: single animated surface ─────────────────────────────────── #
    traj = simulate_dynamic_deformation(
        Nu=15, Nv=15, steps=600, dt=0.01,
        mass=0.03, k_spring=120.0, damping=0.5,
        gravity=np.array([0, 0, -9.81]),
        fix_boundary=True,
        bolt_cells_1based=BOLT_CELLS,
        servo_cells_1based=SERVO_CELLS,
        servo_profile=profile,
    )
    animate_lattice(traj, BOLT_CELLS, SERVO_CELLS)

    # ── View 2: open-loop vs MPC  (compliance precompute ~5 s) ──────────── #
    traj_ol, traj_cl, free_nodes, cs_mask = simulate_comparison(
        Nu=15, Nv=15, steps=600, dt=0.01,
        mass=0.03, k_spring=120.0, damping=0.5,
        gravity=np.array([0, 0, -9.81]),
        fix_boundary=True,
        bolt_cells_1based=BOLT_CELLS,
        servo_cells_1based=SERVO_CELLS,
        servo_profile=profile,
        mpc_lambda=0.05,   # Tikhonov regularisation on the servo correction
        u_clip=0.15,       # actuator travel limit ±0.15 m
    )
    animate_comparison(
        traj_ol, traj_cl, free_nodes, cs_mask,
        SERVO_CELLS, BOLT_CELLS, 15, 15, profile, dt=0.01,
    )