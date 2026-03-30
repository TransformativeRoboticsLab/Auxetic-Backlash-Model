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

def _to_zero_based(pairs_1based):
    """Convert 1-based (row,col) pairs to 0-based (i,j) numpy array of indices."""
    return np.array([(r-1, c-1) for (r, c) in pairs_1based], dtype=int)

def simulate_dynamic_deformation(
    Nu=15, Nv=15, steps=400, dt=0.01,
    mass=0.1, k_spring=50.0, damping=0.98,
    gravity=np.array([0, 0, -9.81]),
    fix_boundary=True,

 
    bolt_cells_1based=((1,2),(1,14),(8,2),(8,14),(15,2),(15,14)),
    servo_cells_1based=((2,3),(2,8),(2,13),
                        (6,3),(6,8),(7,13),
                        (10,3),(10,8),(10,13),
                        (14,3),(14,8),(14,13)),

    # --- Actuation profile for servo cells (return z-offset at time t) ---
    # Example: gentle downwards pulse with a smooth settle
    servo_profile=lambda t: -0.05*np.sin(np.pi*min(t,1.0)) if t <= 1.0 else -0.05
):
    """
    Run dynamic lattice simulation under gravity + Jacob's intermediate conditions:
      - bolt cells are clamped to their initial positions
      - servo cells receive a prescribed z(t) displacement (actuation)
    """

    # Initial flat lattice
    P0 = make_initial_lattice(Nu, Nv)
    P   = P0.copy()
    V   = np.zeros_like(P)  # velocity
    F   = np.zeros_like(P)  # total force

    L0x, L0y = rest_lengths(P0)

    # Build index sets (0-based) and sanity-filter for range
    bolt_idx  = _to_zero_based(bolt_cells_1based)
    servo_idx = _to_zero_based(servo_cells_1based)
    in_bounds = lambda i,j: (0 <= i < Nu) and (0 <= j < Nv)
    bolt_idx  = np.array([ij for ij in bolt_idx  if in_bounds(*ij)], dtype=int)
    servo_idx = np.array([ij for ij in servo_idx if in_bounds(*ij)], dtype=int)

    # Quick helpers to clamp sets of nodes
    def clamp_nodes_to_initial(idxs):
        """Set selected nodes back to P0 and zero their velocities."""
        if idxs.size == 0:
            return
        ii, jj = idxs[:,0], idxs[:,1]
        P[ii, jj, :] = P0[ii, jj, :]
        V[ii, jj, :] = 0.0

    def set_servo_nodes(t):
        """Drive selected nodes with z(t) offset relative to initial position."""
        if servo_idx.size == 0:
            return
        z = servo_profile(t)
        ii, jj = servo_idx[:,0], servo_idx[:,1]
        # keep x,y at initial; set z = initial z + z(t)
        P[ii, jj, 0:2] = P0[ii, jj, 0:2]
        P[ii, jj, 2]   = P0[ii, jj, 2] + z
        V[ii, jj, :]   = 0.0

    trajectory = []

    for step in range(steps):
        t = step * dt

        # Reset forces
        F[:] = 0.0

        # Springs + gravity
        F += compute_spring_forces(P, L0x, L0y, k=k_spring)
        F += mass * gravity

        # Integrate
        A = F / mass
        V += A * dt
        V *= damping  # 0 < damping <= 1 recommended
        P += V * dt

        # --- Apply constraints (order matters) ---
        # 1) Global boundary pins (if requested)
        if fix_boundary:
            P[0, :, :]  = P0[0, :, :];   V[0, :, :]  = 0.0
            P[-1, :, :] = P0[-1, :, :];  V[-1, :, :] = 0.0
            P[:, 0, :]  = P0[:, 0, :];   V[:, 0, :]  = 0.0
            P[:, -1, :] = P0[:, -1, :];  V[:, -1, :] = 0.0

        # 2) Frame bolts (orange) – fixed to initial
        clamp_nodes_to_initial(bolt_idx)

        # 3) Actuated servo cells (red) – prescribed z(t)
        set_servo_nodes(t)

        # Record
        trajectory.append(P.copy())

    return np.array(trajectory)  # (steps, Nu, Nv, 3)



# ---------- Visualization ---------- #

def animate_lattice(trajectory, bolt_cells_1based, servo_cells_1based, interval=30):
    """Animate lattice deformation in 3D, with bolt and servo cells highlighted."""
    steps, Nu, Nv, _ = trajectory.shape

    # Convert 1-based coordinates to 0-based
    def to_zero_based(pairs):
        return np.array([(r - 1, c - 1) for (r, c) in pairs], dtype=int)

    bolt_idx = to_zero_based(bolt_cells_1based)
    servo_idx = to_zero_based(servo_cells_1based)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    def update(frame):
        ax.cla()
        P = trajectory[frame]
        ax.set_title(f"Dynamic Lattice | t={frame}")
        
        # plot grid lines
        for i in range(Nu):
            ax.plot(P[i, :, 0], P[i, :, 1], P[i, :, 2], color='tab:blue', lw=1.2)
        for j in range(Nv):
            ax.plot(P[:, j, 0], P[:, j, 1], P[:, j, 2], color='tab:blue', lw=1.2)
        
        # scatter all nodes
        ax.scatter(P[..., 0], P[..., 1], P[..., 2], s=10, c='k', alpha=0.6)

        # overlay servo and bolt nodes
        if len(servo_idx) > 0:
            ax.scatter(P[servo_idx[:,0], servo_idx[:,1], 0],
                       P[servo_idx[:,0], servo_idx[:,1], 1],
                       P[servo_idx[:,0], servo_idx[:,1], 2],
                       s=50, c='red', label='Servo (Actuated)')
        if len(bolt_idx) > 0:
            ax.scatter(P[bolt_idx[:,0], bolt_idx[:,1], 0],
                       P[bolt_idx[:,0], bolt_idx[:,1], 1],
                       P[bolt_idx[:,0], bolt_idx[:,1], 2],
                       s=50, c='orange', label='Bolt (Fixed)')

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(-0.5, 0.5)
        ax.view_init(elev=25, azim=35)
        ax.legend(loc='upper left')
        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    plt.tight_layout()
    plt.show()



# ---------- Run demo ---------- #

if __name__ == "__main__":
    bolt_cells = ()
    servo_cells = ((3,3),(3,8),(3,13),
                   (8,3),(8,8),(8,13),
                   (13,3),(13,8),(13,13))
    traj = simulate_dynamic_deformation(
        Nu=15, 
        Nv=15,
        steps=600, 
        dt=0.01,
        mass=0.05, 
        k_spring=120.0, 
        damping=0.85,
        gravity=np.array([0, 0, -9.81]),
        fix_boundary=True,
        servo_profile=lambda t: -0.04*np.sin(np.pi*min(t,1.0)) if t <= 1.0 else -0.04
    )

    animate_lattice(traj,  bolt_cells, servo_cells)


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