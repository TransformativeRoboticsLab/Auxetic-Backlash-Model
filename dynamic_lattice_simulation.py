import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


# This script extends the static lattice relaxation model into a dynamic simulation,
# introducing *time evolution* (t) and *physical deformation* under gravity.

# Each lattice node is modeled as a point mass connected to its neighbors via springs.
# The system evolves according to Newtonian dynamics (F = m * a), allowing the lattice
# to deform, oscillate, and settle over time.

# Core physics:
# --------------
# - Each node has:
#     * Mass (m)
#     * Velocity (v)
#     * Forces from springs, damping, and gravity

# - Springs connect nodes horizontally and vertically, maintaining rest lengths L0x and L0y.
# - Damping dissipates energy to simulate realistic settling.
# - Boundary nodes are fixed to anchor the structure.

# Mathematically:
#     F_total = F_spring + F_gravity + F_damping
#     a = F_total / m
#     v += a * dt
#     p += v * dt

#Increase gravity → faster, deeper sag

# Decrease gravity → slower, gentler deformation

# Outputs:
# ---------
# - A 3D animated lattice deformation under gravity.
# - (Optional) Z(x, y, t) displacement data for later comparison with experimental results.




# ---------- Helpers from lattice solver ---------- #

def make_initial_lattice(Nu=15, Nv=15, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    """Create flat NxN lattice of nodes in 3D."""
    W, H = size
    x0, y0, z0 = origin
    xs = np.linspace(x0, x0 + W, Nu)
    ys = np.linspace(y0, y0 + H, Nv)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = np.full_like(X, z0)
    P = np.stack([X, Y, Z], axis=-1)  # (Nu, Nv, 3)
    return P


def rest_lengths(P0):
    """Compute horizontal/vertical rest lengths for lattice."""
    seg_x = P0[:, 1:, :] - P0[:, :-1, :]
    seg_y = P0[1:, :, :] - P0[:-1, :, :]
    L0x = np.linalg.norm(seg_x, axis=-1)
    L0y = np.linalg.norm(seg_y, axis=-1)
    return L0x, L0y


def compute_spring_forces(P, L0x, L0y, k=50.0):
    """Compute spring forces on each node from neighboring edges."""
    Nu, Nv, _ = P.shape
    F = np.zeros_like(P)

    # Horizontal springs
    for i in range(Nu):
        for j in range(Nv - 1):
            p1, p2 = P[i, j], P[i, j + 1]
            seg = p2 - p1
            L = np.linalg.norm(seg)
            if L < 1e-12:
                continue
            dir = seg / L
            Fmag = k * (L - L0x[i, j])
            F[i, j] += Fmag * dir
            F[i, j + 1] -= Fmag * dir

    # Vertical springs
    for i in range(Nu - 1):
        for j in range(Nv):
            p1, p2 = P[i, j], P[i + 1, j]
            seg = p2 - p1
            L = np.linalg.norm(seg)
            if L < 1e-12:
                continue
            dir = seg / L
            Fmag = k * (L - L0y[i, j])
            F[i, j] += Fmag * dir
            F[i + 1, j] -= Fmag * dir

    return F


# ---------- Dynamic simulation ---------- #

def simulate_dynamic_deformation(Nu=15, Nv=15, steps=400, dt=0.01,
                                 mass=0.1, k_spring=50.0, damping=0.98,
                                 gravity=np.array([0, 0, -9.81]),
                                 fix_boundary=True):
    """Run dynamic lattice simulation under gravity."""

    # Initial flat lattice
    P0 = make_initial_lattice(Nu, Nv)
    P = P0.copy()
    V = np.zeros_like(P)  # velocity
    F = np.zeros_like(P)  # total force

    L0x, L0y = rest_lengths(P0)

    trajectory = []

    for step in range(steps):
        # Reset forces
        F[:] = 0.0

        # Spring forces
        F += compute_spring_forces(P, L0x, L0y, k=k_spring)

        # Gravity
        F += mass * gravity

        # Acceleration & integration
        A = F / mass
        V += A * dt
        V *= damping  # damping for stability
        P += V * dt

        # Fix boundary nodes
        if fix_boundary:
            P[0, :, :] = P0[0, :, :]
            P[-1, :, :] = P0[-1, :, :]
            P[:, 0, :] = P0[:, 0, :]
            P[:, -1, :] = P0[:, -1, :]

        # Record snapshot
        trajectory.append(P.copy())

    return np.array(trajectory)  # shape: (steps, Nu, Nv, 3)


# ---------- Visualization ---------- #

def animate_lattice(trajectory, interval=30):
    """Animate lattice deformation in 3D."""
    steps, Nu, Nv, _ = trajectory.shape
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    def update(frame):
        ax.cla()
        P = trajectory[frame]
        ax.set_title(f"Dynamic Deformation | t={frame}")
        for i in range(Nu):
            ax.plot(P[i, :, 0], P[i, :, 1], P[i, :, 2], color='tab:blue', lw=1.5)
        for j in range(Nv):
            ax.plot(P[:, j, 0], P[:, j, 1], P[:, j, 2], color='tab:blue', lw=1.5)
        ax.scatter(P[..., 0], P[..., 1], P[..., 2], s=8, c='k')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(-0.5, 0.5)
        ax.view_init(elev=25, azim=30)
        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    plt.tight_layout()
    plt.show()


# ---------- Run demo ---------- #

if __name__ == "__main__":
    traj = simulate_dynamic_deformation(
        Nu=15, Nv=15,
        steps=500, dt=0.01,
        mass=0.05, k_spring=80.0, damping=0.8,
        gravity=np.array([0, 0, -9.81]),
        fix_boundary=True
    )

    animate_lattice(traj)


# V += A * dt
# V *= damping
# P += V * dt
# A is acceleration from total forces (spring + gravity)

# V is velocity

# damping is a multiplier between 0 and 1, applied each timestep

# so each frame, the velocity is scaled by that damping factor.

# If damping = 0.99, it means:

# keep 99% of your velocity from the previous step (lose 1%).
# That’s a 1% velocity loss per frame, which is small but accumulates, so it slowly settles.

# If damping = 5.0, that means:

# multiply velocity by 5 every frame — so it actually amplifies energy violently,
# and the simulation explodes (NaNs, huge oscillations, instability).

# the correct damping range should always be between 0 and 1.