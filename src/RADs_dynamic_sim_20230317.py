import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


# ---------- Helpers ---------- #

def make_initial_lattice(Nu=15, Nv=15, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    W, H = size
    x0, y0, z0 = origin
    xs = np.linspace(x0, x0 + W, Nu)
    ys = np.linspace(y0, y0 + H, Nv)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = np.full_like(X, z0)
    return np.stack([X, Y, Z], axis=-1)  # (Nu, Nv, 3)


def rest_lengths(P0):
    L0x = np.linalg.norm(P0[:, 1:, :] - P0[:, :-1, :], axis=-1)
    L0y = np.linalg.norm(P0[1:, :, :] - P0[:-1, :, :], axis=-1)
    return L0x, L0y


def compute_spring_forces(P, L0x, L0y, k=50.0):
    Nu, Nv, _ = P.shape
    F = np.zeros_like(P)

    # Horizontal springs
    seg = P[:, 1:, :] - P[:, :-1, :]          # (Nu, Nv-1, 3)
    L   = np.linalg.norm(seg, axis=-1, keepdims=True).clip(1e-12)
    d   = seg / L
    Fm  = k * (L[..., 0] - L0x)[..., None] * d
    F[:, :-1, :] += Fm
    F[:,  1:, :] -= Fm

    # Vertical springs
    seg = P[1:, :, :] - P[:-1, :, :]          # (Nu-1, Nv, 3)
    L   = np.linalg.norm(seg, axis=-1, keepdims=True).clip(1e-12)
    d   = seg / L
    Fm  = k * (L[..., 0] - L0y)[..., None] * d
    F[:-1, :, :] += Fm
    F[ 1:, :, :] -= Fm

    return F


# ---------- Dynamic simulation ---------- #

def _to_zero_based(pairs_1based):
    return np.array([(r - 1, c - 1) for (r, c) in pairs_1based], dtype=int)


def simulate_dynamic_deformation(
    Nu=15, Nv=15, steps=600, dt=0.01,
    mass=0.05, k_spring=120.0, damping=0.85,
    gravity=np.array([0, 0, -9.81]),
    fix_boundary=True,
    bolt_cells_1based=(),
    servo_cells_1based=(
        (3, 3),  (3, 8),  (3, 13),
        (8, 3),  (8, 8),  (8, 13),
        (13, 3), (13, 8), (13, 13),
    ),
    # --- Servo profile ---
    # Signature: servo_profile(t, ni, nj)
    #   t  – simulation time
    #   ni – normalised i-position of this node in [0, 1]
    #   nj – normalised j-position of this node in [0, 1]
    # Returns a z-offset (metres) for that node at that time.
    #
    # Default: a sine wave travelling along the i-axis.
    #   amplitude  A  = 0.08 m
    #   wavelength λ  = 1 lattice-unit (full width)
    #   wave speed  c  = 0.5 lattice-units per second  →  ω = 2π·c/λ
    servo_profile=lambda t, ni, nj: 0.08 * np.sin(2 * np.pi * ni - 2 * np.pi * 0.5 * t),
    ):
    """
    Run dynamic lattice simulation.

    Servo nodes are driven with a *position-dependent* z(t, ni, nj) profile,
    so each node receives a different offset — enough to approximate a spatial
    sine wave sampled at the servo locations.
    """
    P0 = make_initial_lattice(Nu, Nv)
    P  = P0.copy()
    V  = np.zeros_like(P)
    F  = np.zeros_like(P)

    L0x, L0y = rest_lengths(P0)

    in_bounds = lambda i, j: (0 <= i < Nu) and (0 <= j < Nv)
    bolt_idx  = _to_zero_based(bolt_cells_1based)
    servo_idx = _to_zero_based(servo_cells_1based)
    bolt_idx  = np.array([ij for ij in bolt_idx  if in_bounds(*ij)], dtype=int)
    servo_idx = np.array([ij for ij in servo_idx if in_bounds(*ij)], dtype=int)

    # Pre-compute normalised positions for every servo node once
    if servo_idx.size > 0:
        servo_ni = servo_idx[:, 0] / (Nu - 1)   # i / (Nu-1)  →  [0, 1]
        servo_nj = servo_idx[:, 1] / (Nv - 1)   # j / (Nv-1)  →  [0, 1]

    def clamp_bolts():
        if bolt_idx.size == 0:
            return
        ii, jj = bolt_idx[:, 0], bolt_idx[:, 1]
        P[ii, jj, :] = P0[ii, jj, :]
        V[ii, jj, :] = 0.0

    def drive_servos(t):
        if servo_idx.size == 0:
            return
        ii, jj = servo_idx[:, 0], servo_idx[:, 1]
        # Vectorised: evaluate profile for every servo node at once
        z_offsets = servo_profile(t, servo_ni, servo_nj)   # shape (n_servos,)
        P[ii, jj, 0:2] = P0[ii, jj, 0:2]                  # lock x, y
        P[ii, jj, 2]   = P0[ii, jj, 2] + z_offsets
        V[ii, jj, :]   = 0.0

    trajectory = []

    for step in range(steps):
        t = step * dt

        F[:] = 0.0
        F += compute_spring_forces(P, L0x, L0y, k=k_spring)
        F += mass * gravity

        A  = F / mass
        V += A * dt
        V *= damping
        P += V * dt

        # --- Constraints (order matters) ---
        if fix_boundary:
            P[0,  :, :] = P0[0,  :, :];  
            V[0,  :, :] = 0.0
            P[-1, :, :] = P0[-1, :, :];  
            V[-1, :, :] = 0.0
            P[:,  0, :] = P0[:,  0, :];  
            V[:,  0, :] = 0.0
            P[:, -1, :] = P0[:, -1, :];  
            V[:, -1, :] = 0.0

        clamp_bolts()
        drive_servos(t)

        trajectory.append(P.copy())

    return np.array(trajectory)   # (steps, Nu, Nv, 3)


# ---------- Visualization ---------- #

def animate_lattice(trajectory, bolt_cells_1based, servo_cells_1based, interval=30):
    steps, Nu, Nv, _ = trajectory.shape

    def to_zero_based(pairs):
        return np.array([(r - 1, c - 1) for (r, c) in pairs], dtype=int)

    bolt_idx  = to_zero_based(bolt_cells_1based)
    servo_idx = to_zero_based(servo_cells_1based)

    fig = plt.figure(figsize=(8, 6))
    ax  = fig.add_subplot(111, projection='3d')

    def update(frame):
        ax.cla()
        P = trajectory[frame]
        ax.set_title(f"Sine-wave driven lattice  |  t = {frame * 0.01:.2f} s")

        for i in range(Nu):
            ax.plot(P[i, :, 0], P[i, :, 1], P[i, :, 2], color='tab:blue', lw=1.2)
        for j in range(Nv):
            ax.plot(P[:, j, 0], P[:, j, 1], P[:, j, 2], color='tab:blue', lw=1.2)

        ax.scatter(P[..., 0], P[..., 1], P[..., 2], s=10, c='k', alpha=0.6)

        if servo_idx.size > 0:
            ax.scatter(P[servo_idx[:, 0], servo_idx[:, 1], 0],
                       P[servo_idx[:, 0], servo_idx[:, 1], 1],
                       P[servo_idx[:, 0], servo_idx[:, 1], 2],
                       s=60, c='red', label='Servo (Actuated)')
        if bolt_idx.size > 0:
            ax.scatter(P[bolt_idx[:, 0], bolt_idx[:, 1], 0],
                       P[bolt_idx[:, 0], bolt_idx[:, 1], 1],
                       P[bolt_idx[:, 0], bolt_idx[:, 1], 2],
                       s=60, c='orange', label='Bolt (Fixed)')

        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(-0.3, 0.3)
        ax.view_init(elev=25, azim=35)
        ax.legend(loc='upper left')
        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    plt.tight_layout()
    plt.show()


# ---------- Run ---------- #

if __name__ == "__main__":
    BOLT_CELLS  = ()
    SERVO_CELLS = (
        (3, 3),  (3, 8),  (3, 13),
        (8, 3),  (8, 8),  (8, 13),
        (13, 3), (13, 8), (13, 13),
    )

    # ------------------------------------------------------------------ #
    # Servo profile options — pick one (or write your own):
    #
    #  Travelling wave along i-axis:
    #    lambda t, ni, nj: A * sin(2π·ni - ω·t)
    #
    #  Standing wave:
    #    lambda t, ni, nj: A * sin(2π·ni) * cos(ω·t)
    #
    #  Diagonal travelling wave:
    #    lambda t, ni, nj: A * sin(2π*(ni + nj)/2 - ω·t)
    #
    # Parameters:
    #   A = amplitude (metres)   ω = angular frequency (rad/s)
    # ------------------------------------------------------------------ #

    A, omega = 0.05, np.pi   # 0.05 m amplitude, 0.5 Hz wave speed

    traj = simulate_dynamic_deformation(
        Nu=15, Nv=15,
        steps=600, dt=0.01,
        mass=0.03, k_spring=120.0, damping=0.5,
        gravity=np.array([0, 0, -9.81]),
        fix_boundary=True,
        bolt_cells_1based=BOLT_CELLS,
        servo_cells_1based=SERVO_CELLS,
        servo_profile=lambda t, ni, nj: A * np.sin(2 * np.pi * ni - omega * t),
    )

    animate_lattice(traj, BOLT_CELLS, SERVO_CELLS)