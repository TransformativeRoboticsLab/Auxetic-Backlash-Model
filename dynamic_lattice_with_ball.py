import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# ---------- Helpers from lattice solver ---------- #

def make_initial_lattice(Nu=15, Nv=15, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    W, H = size
    x0, y0, z0 = origin
    xs = np.linspace(x0, x0 + W, Nu)
    ys = np.linspace(y0, y0 + H, Nv)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = np.full_like(X, z0)
    return np.stack([X, Y, Z], axis=-1)


def rest_lengths(P0):
    seg_x = P0[:, 1:, :] - P0[:, :-1, :]
    seg_y = P0[1:, :, :] - P0[:-1, :, :]
    return np.linalg.norm(seg_x, axis=-1), np.linalg.norm(seg_y, axis=-1)


def compute_spring_forces(P, L0x, L0y, k=50.0):
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


# ---------- Lattice & Servo Setup ---------- #

def _to_zero_based(pairs_1based):
    return np.array([(r-1, c-1) for (r, c) in pairs_1based], dtype=int)


# ---------- Dynamic Deformation Simulation ---------- #

def simulate_dynamic_deformation_with_sequential_servos(
    Nu=15, Nv=15, steps=1200, dt=0.01,
    mass=0.1, k_spring=50.0, damping=0.98,
    gravity=np.array([0, 0, -9.81]),
    fix_boundary=True,
    bolt_cells_1based=((1,2),(1,14),(8,2),(8,14),(15,2),(15,14)),
    servo_cells_1based=((2,8),(8,2),(8,14),(14,8)),  
    servo_T=1.0,
    servo_amp=-0.02
):
    """Sequential actuation of four servos (plus-shaped layout)."""
    P0 = make_initial_lattice(Nu, Nv)
    P, V = P0.copy(), np.zeros_like(P0)
    F = np.zeros_like(P)
    L0x, L0y = rest_lengths(P0)

    bolt_idx = _to_zero_based(bolt_cells_1based)
    servo_idx = _to_zero_based(servo_cells_1based)
    servo_sequence = list(map(tuple, servo_idx))
    servo_displacements = {i: [] for i in range(len(servo_sequence))}

    trajectory, active_log = [], []

    for step in range(steps):
        t = step * dt
        F[:] = 0
        F += compute_spring_forces(P, L0x, L0y, k=k_spring)
        F += mass * gravity

        # Integrate motion
        A = F / mass
        V += A * dt
        V *= damping
        P += V * dt

        # Clamp boundaries and bolts
        if fix_boundary:
            P[0,:,:]=P0[0,:,:]; P[-1,:,:]=P0[-1,:,:]
            P[:,0,:]=P0[:,0,:]; P[:,-1,:]=P0[:,-1,:]
            V[0,:,:]=V[-1,:,:]=V[:,0,:]=V[:,-1,:]=0
        if len(bolt_idx)>0:
            ii,jj=bolt_idx[:,0],bolt_idx[:,1]
            P[ii,jj,:]=P0[ii,jj,:]; V[ii,jj,:]=0

        # Sequential servo control
        total_servos = len(servo_sequence)
        idx = int((t // servo_T) % total_servos)
        ci, cj = servo_sequence[idx]
        local_t = (t % servo_T)
        z_offset = servo_amp * np.sin(np.pi * local_t / servo_T)
        P[ci, cj, 2] = P0[ci, cj, 2] + z_offset
        V[ci, cj, :] = 0
        for i in range(total_servos):
            if i == idx:
                servo_displacements[i].append(z_offset)
            else:
                servo_displacements[i].append(0)
        active_log.append((t, ci, cj))
        trajectory.append(P.copy())

    return np.array(trajectory), servo_displacements


# ---------- Ball Simulation ---------- #

def interp_bilinear(Z, x, y, Xmin=0.0, Xmax=1.0, Ymin=0.0, Ymax=1.0):
    Nu, Nv = Z.shape
    u = (x - Xmin) * (Nu - 1) / (Xmax - Xmin)
    v = (y - Ymin) * (Nv - 1) / (Ymax - Ymin)
    i = int(np.clip(np.floor(u), 0, Nu - 2))
    j = int(np.clip(np.floor(v), 0, Nv - 2))
    du, dv = u - i, v - j
    z00, z10 = Z[i, j], Z[i + 1, j]
    z01, z11 = Z[i, j + 1], Z[i + 1, j + 1]
    return (1 - du)*(1 - dv)*z00 + du*(1 - dv)*z10 + (1 - du)*dv*z01 + du*dv*z11


def surface_gradients(Z, x, y, size=(1.0, 1.0)):
    Nu, Nv = Z.shape
    W, H = size
    u = np.clip(x * (Nu - 1) / W, 1, Nu - 2)
    v = np.clip(y * (Nv - 1) / H, 1, Nv - 2)
    i, j = int(round(u)), int(round(v))
    dzdx = (Z[i + 1, j] - Z[i - 1, j]) * (Nu - 1) / (2 * W)
    dzdy = (Z[i, j + 1] - Z[i, j - 1]) * (Nv - 1) / (2 * H)
    return dzdx, dzdy


def simulate_ball_over_lattice(traj, dt, size=(1.0,1.0),
                               ball_radius=0.03, g=9.81, mu=0.03,
                               damping=0.995, slope_gain=0.6,
                               x0=0.5, y0=0.5):
    steps, Nu, Nv, _ = traj.shape
    W, H = size
    bx, by = x0, y0
    Z0 = traj[0,:,:,2]
    bz = interp_bilinear(Z0, bx*W, by*H, 0,W,0,H) + ball_radius + 1e-3
    vx = vy = vz = 0.0
    out = np.zeros((steps,3))
    for t in range(steps):
        Z = traj[t,:,:,2]
        z_surf = interp_bilinear(Z,bx*W,by*H,0,W,0,H)
        dzdx,dzdy = surface_gradients(Z,bx*W,by*H,size=size)
        ax, ay, az = 0.0, 0.0, -g
        contact = (bz - ball_radius) <= z_surf
        if contact:
            bz = z_surf + ball_radius
            if vz < 0: vz = 0
            ax += -slope_gain*g*dzdx
            ay += -slope_gain*g*dzdy
            vx *= (1 - mu); vy *= (1 - mu)
        vx += ax*dt; vy += ay*dt; vz += az*dt
        bx += vx*dt; by += vy*dt; bz += vz*dt
        vx *= damping; vy *= damping; vz *= damping
        bx = np.clip(bx,0,1); by = np.clip(by,0,1)
        out[t] = (bx, by, bz)
    return out


# ---------- Visualization ---------- #

def animate_with_ball(traj, ball_traj, bolt_cells_1based, servo_cells_1based, interval=30):
    steps, Nu, Nv, _ = traj.shape
    bolt_idx = _to_zero_based(bolt_cells_1based)
    servo_idx = _to_zero_based(servo_cells_1based)
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    def update(frame):
        ax.cla()
        P = traj[frame]
        bx, by, bz = ball_traj[frame]
        ax.set_title(f"Lattice + Ball | frame {frame}")
        for i in range(Nu): ax.plot(P[i,:,0], P[i,:,1], P[i,:,2], color='tab:blue', lw=1.2)
        for j in range(Nv): ax.plot(P[:,j,0], P[:,j,1], P[:,j,2], color='tab:blue', lw=1.2)
        ax.scatter(P[servo_idx[:,0],servo_idx[:,1],0],
                   P[servo_idx[:,0],servo_idx[:,1],1],
                   P[servo_idx[:,0],servo_idx[:,1],2],
                   s=40,c='red',label='Servos')
        ax.scatter(P[bolt_idx[:,0],bolt_idx[:,1],0],
                   P[bolt_idx[:,0],bolt_idx[:,1],1],
                   P[bolt_idx[:,0],bolt_idx[:,1],2],
                   s=40,c='orange',label='Bolts')
        ax.scatter([bx],[by],[bz], color='gold', s=100, edgecolor='black')
        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_zlim(-0.2,0.3)
        ax.view_init(elev=25, azim=35)
        return []

    ani = FuncAnimation(fig, update, frames=range(0, steps, 2), interval=interval)
    plt.tight_layout(); plt.show()


# ---------- Run Simulation ---------- #

if __name__ == "__main__":
    bolt_cells = ((1,2),(1,14),(8,2),(8,14),(15,2),(15,14))
    servo_cells = ((2,8),(8,4),(8,12),(14,8))  # plus-shape servos

    traj, servo_displacements = simulate_dynamic_deformation_with_sequential_servos(
        Nu=15, Nv=15,
        steps=1600, dt=0.01,
        mass=0.05, k_spring=80.0, damping=0.85,
        gravity=np.array([0,0,-9.81]),
        servo_T=1.0, servo_amp=-0.02,
        servo_cells_1based=servo_cells

    )

    ball_traj = simulate_ball_over_lattice(traj, dt=0.01, size=(1.0,1.0),
                                           ball_radius=0.03, g=9.81, mu=0.03,
                                           damping=0.995, x0=0.5, y0=0.5)

    animate_with_ball(traj, ball_traj, bolt_cells, servo_cells)

    # ---- Four plots (one per servo) ----
    t = np.arange(len(traj)) * 0.01
    labels = ["Top (2,8)", "Left (8,2)", "Right (8,14)", "Bottom (14,8)"]
    fig, axs = plt.subplots(4, 1, figsize=(6,6), sharex=True)
    for i, key in enumerate(servo_displacements.keys()):
        axs[i].plot(t, np.array(servo_displacements[key])*1000, lw=2)
        axs[i].set_ylabel(f"{labels[i]}\nZ [mm]")
        axs[i].grid(True)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("Servo Actuation Profiles (Plus Shape, 2 cm Up Each)")
    plt.tight_layout()
    plt.show()
