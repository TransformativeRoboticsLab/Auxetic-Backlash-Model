"""
simulate_closed_loop.py

Identify a reduced-order plant model from experimental OptiTrack data,
design a per-mode PID controller, simulate the closed-loop response,
and compare against the open-loop experimental measurements.

Modeling approach
─────────────────
1.  SVD of the (n_markers × n_frames) Z-displacement matrix produces
    spatial basis functions (columns of U) and temporal modal coordinates.
2.  Each modal coordinate is treated as the output of a second-order
    plant driven by a sine-wave command:
        G(s) = K ωₙ² / (s² + 2ζωₙs + ωₙ²)
    Parameters (K, ωₙ, ζ) are estimated from the experimental PSD.
3.  A discrete-time PID controller is designed per mode to track the
    desired sine reference with improved accuracy.
4.  The full surface is reconstructed from the controlled modal
    coordinates via the spatial basis.

Usage:
    python simulate_closed_loop.py data/RADs_*.csv -o results/
    python simulate_closed_loop.py test1.csv test2.csv test3.csv
"""

import sys
import glob
import re
import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy.optimize import minimize_scalar, minimize
from pathlib import Path


# ═════════════════════════════════════════════
#  OPTITRACK CSV PARSER (self-contained)
# ═════════════════════════════════════════════
def parse_optitrack_csv(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    meta_cells = lines[0].strip().replace("\r", "").split(",")
    meta = {}
    for i in range(0, len(meta_cells) - 1, 2):
        key = meta_cells[i].strip()
        val = meta_cells[i + 1].strip()
        if key:
            meta[key] = val

    frame_rate = float(
        meta.get("Export Frame Rate", meta.get("Capture Frame Rate", 120))
    )
    meta["frame_rate"] = frame_rate

    name_cells = lines[3].strip().replace("\r", "").split(",")
    raw_names = [c.strip().strip('"') for c in name_cells[2:]]

    marker_order = []
    seen = set()
    for i in range(0, len(raw_names), 3):
        name = raw_names[i]
        if name and name not in seen:
            marker_order.append(name)
            seen.add(name)

    data_rows = []
    for line in lines[7:]:
        parts = line.strip().replace("\r", "").split(",")
        row = []
        for p in parts:
            try:
                row.append(float(p))
            except ValueError:
                row.append(np.nan)
        data_rows.append(row)

    data = np.array(data_rows)
    time = data[:, 1]

    markers = {}
    for idx, name in enumerate(marker_order):
        col_base = 2 + 3 * idx
        markers[name] = {
            "x": data[:, col_base].copy(),
            "y": data[:, col_base + 1].copy(),
            "z": data[:, col_base + 2].copy(),
        }

    return meta, time, markers


def extract_test_label(filepath):
    stem = Path(filepath).stem
    period_match = re.search(r"(\d+\.?\d*)s[_\s]*period", stem, re.IGNORECASE)
    amp_match = re.search(r"(\d+\.?\d*)deg[_\s]*amplitude", stem, re.IGNORECASE)
    parts = []
    if period_match:
        parts.append(f"T={period_match.group(1)}s")
    if amp_match:
        parts.append(f"A={amp_match.group(1)}deg")
    if parts:
        return " / ".join(parts)
    label = stem.replace("_", " ").strip()
    return label[:40] if len(label) > 40 else label


def extract_sort_key(filepath):
    stem = Path(filepath).stem
    m = re.search(r"(\d+\.?\d*)s[_\s]*period", stem, re.IGNORECASE)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+\.?\d*)deg[_\s]*amplitude", stem, re.IGNORECASE)
    if m:
        return float(m.group(1))
    return 0.0


# ═════════════════════════════════════════════
#  SVD-BASED SURFACE DECOMPOSITION
# ═════════════════════════════════════════════
class SurfaceBasisModel:
    """
    Decomposes a marker-Z matrix into spatial modes (SVD)
    and fits a second-order transfer function to each retained mode.
    """

    def __init__(self, Z, time, fs, n_modes=None, energy_threshold=0.99):
        """
        Parameters
        ----------
        Z : ndarray (n_markers, n_frames)
            Mean-subtracted Z displacements.
        time : ndarray (n_frames,)
        fs : float
        n_modes : int or None
            If None, retain enough modes to capture energy_threshold.
        energy_threshold : float
            Fraction of total variance to retain (default 99%).
        """
        self.fs = fs
        self.dt = 1.0 / fs
        self.time = time
        self.n_frames = len(time)
        self.n_markers = Z.shape[0]

        # Per-marker means (for reconstruction)
        self.z_means = np.nanmean(Z, axis=1)
        Z_centered = Z - self.z_means[:, None]

        # Replace any remaining NaNs with 0 for SVD
        Z_clean = np.nan_to_num(Z_centered, nan=0.0)

        # SVD
        self.U, self.S, self.Vt = np.linalg.svd(Z_clean, full_matrices=False)

        # Determine number of modes
        energy_cumulative = np.cumsum(self.S ** 2) / np.sum(self.S ** 2)
        if n_modes is None:
            n_modes = int(np.searchsorted(energy_cumulative, energy_threshold) + 1)
            n_modes = max(n_modes, 1)
        self.n_modes = min(n_modes, len(self.S))

        # Truncate
        self.U_r = self.U[:, : self.n_modes]          # (n_markers, n_modes)
        self.S_r = self.S[: self.n_modes]              # (n_modes,)
        self.Vt_r = self.Vt[: self.n_modes, :]         # (n_modes, n_frames)

        # Modal time series: q_i(t) = S_i * V_i(t)
        self.modal_coords = self.S_r[:, None] * self.Vt_r  # (n_modes, n_frames)

        # Energy fractions
        total_energy = np.sum(self.S ** 2)
        self.mode_energy = self.S_r ** 2 / total_energy

        # Fit second-order plant per mode
        self.plants = []
        for i in range(self.n_modes):
            plant = self._fit_second_order(self.modal_coords[i])
            self.plants.append(plant)

    def _fit_second_order(self, q):
        """
        Estimate second-order parameters (wn, zeta, K) from the PSD
        of a modal time series q(t).

        G(s) = K * wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
        """
        # Welch PSD
        nperseg = min(len(q), int(8 * self.fs))
        freqs, psd = sig.welch(q, fs=self.fs, nperseg=nperseg,
                               noverlap=nperseg // 2)

        # Find dominant frequency (skip DC)
        mask = freqs > 0.02
        if not np.any(mask):
            return {"wn": 1.0, "zeta": 0.5, "K": 1.0, "peak_freq": 0.0}

        peak_idx = np.argmax(psd[mask])
        peak_freq = freqs[mask][peak_idx]
        wn = 2 * np.pi * peak_freq  # rad/s

        if wn < 1e-6:
            wn = 2 * np.pi * 0.1

        # Estimate damping from half-power bandwidth
        peak_power = psd[mask][peak_idx]
        half_power = peak_power / 2.0
        above_half = psd[mask] >= half_power
        if np.sum(above_half) >= 2:
            indices = np.where(above_half)[0]
            bw = freqs[mask][indices[-1]] - freqs[mask][indices[0]]
            zeta = bw / (2 * peak_freq) if peak_freq > 0 else 0.5
        else:
            zeta = 0.3

        zeta = np.clip(zeta, 0.05, 0.95)

        # Gain: match the RMS of the modal coordinate
        rms_q = np.std(q)
        K = rms_q if rms_q > 1e-10 else 1.0

        return {
            "wn": wn,
            "zeta": zeta,
            "K": K,
            "peak_freq": peak_freq,
        }

    def reconstruct(self, modal_coords):
        """
        Reconstruct full surface Z from modal coordinates.

        Parameters
        ----------
        modal_coords : ndarray (n_modes, n_frames)

        Returns
        -------
        Z_recon : ndarray (n_markers, n_frames)
        """
        Z_recon = self.U_r @ modal_coords + self.z_means[:, None]
        return Z_recon

    def summary(self):
        lines = []
        lines.append(f"  SVD Basis Model: {self.n_modes} modes, "
                      f"{self.n_markers} markers, {self.n_frames} frames")
        for i in range(self.n_modes):
            p = self.plants[i]
            lines.append(
                f"    Mode {i}: energy={self.mode_energy[i]*100:.1f}%  "
                f"f={p['peak_freq']:.4f} Hz  "
                f"wn={p['wn']:.3f} rad/s  "
                f"zeta={p['zeta']:.3f}  "
                f"K={p['K']:.3f} mm"
            )
        return "\n".join(lines)


# ═════════════════════════════════════════════
#  DISCRETE PID CONTROLLER
# ═════════════════════════════════════════════
class PIDController:
    """
    Discrete-time PID with anti-windup and derivative filtering.
    """

    def __init__(self, Kp, Ki, Kd, dt, output_limit=None, d_filter_coeff=0.1):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.dt = dt
        self.output_limit = output_limit
        self.d_filter_coeff = d_filter_coeff

        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_derivative = 0.0

    def reset(self):
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_derivative = 0.0

    def step(self, error):
        # Proportional
        P = self.Kp * error

        # Integral with anti-windup
        self.integral += error * self.dt
        I = self.Ki * self.integral

        # Filtered derivative
        raw_deriv = (error - self.prev_error) / self.dt
        alpha = self.d_filter_coeff
        filtered_deriv = alpha * raw_deriv + (1 - alpha) * self.prev_derivative
        D = self.Kd * filtered_deriv

        output = P + I + D

        # Clamp and anti-windup
        if self.output_limit is not None:
            if abs(output) > self.output_limit:
                output = np.sign(output) * self.output_limit
                # Back-calculate integral to prevent windup
                self.integral = (output - P - D) / self.Ki if abs(self.Ki) > 1e-12 else 0

        self.prev_error = error
        self.prev_derivative = filtered_deriv

        return output


# ═════════════════════════════════════════════
#  SECOND-ORDER PLANT SIMULATOR
# ═════════════════════════════════════════════
class SecondOrderPlant:
    """
    Simulates G(s) = K * wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
    in discrete time via Tustin (bilinear) transform.
    """

    def __init__(self, wn, zeta, K, dt):
        self.wn = wn
        self.zeta = zeta
        self.K = K
        self.dt = dt

        # Convert to discrete-time transfer function via bilinear transform
        num_c = [K * wn ** 2]
        den_c = [1.0, 2 * zeta * wn, wn ** 2]
        self.sys_c = sig.TransferFunction(num_c, den_c)
        self.sys_d = self.sys_c.to_discrete(dt, method="bilinear")

        # State: store last two inputs and outputs for difference equation
        self.b = self.sys_d.num
        self.a = self.sys_d.den

        # Normalize
        self.b = self.b / self.a[0]
        self.a = self.a / self.a[0]

        self.reset()

    def reset(self):
        self.x_hist = [0.0, 0.0]  # past inputs
        self.y_hist = [0.0, 0.0]  # past outputs

    def step(self, u):
        """Advance one time step with input u, return output y."""
        b, a = self.b, self.a
        y = (b[0] * u + b[1] * self.x_hist[0] + b[2] * self.x_hist[1]
             - a[1] * self.y_hist[0] - a[2] * self.y_hist[1])

        # Shift history
        self.x_hist[1] = self.x_hist[0]
        self.x_hist[0] = u
        self.y_hist[1] = self.y_hist[0]
        self.y_hist[0] = y

        return y

    def simulate_open_loop(self, u_array):
        """Simulate the plant response to input array u."""
        self.reset()
        y = np.zeros(len(u_array))
        for i, u in enumerate(u_array):
            y[i] = self.step(u)
        return y


# ═════════════════════════════════════════════
#  PID AUTO-TUNING (per mode)
# ═════════════════════════════════════════════
def design_pid_for_mode(plant_params, dt, bandwidth_mult=3.0):
    """
    Design a PID controller for one modal plant using
    a frequency-domain approach.

    Strategy: set the crossover frequency at bandwidth_mult × wn
    to get faster tracking, with enough phase margin for stability.
    """
    wn = plant_params["wn"]
    zeta = plant_params["zeta"]
    K = plant_params["K"]

    if wn < 1e-6 or K < 1e-10:
        return PIDController(1.0, 0.1, 0.0, dt)

    # Target crossover frequency
    wc = bandwidth_mult * wn

    # At crossover, the plant magnitude is approximately:
    # |G(jwc)| ≈ K * wn^2 / sqrt((wn^2 - wc^2)^2 + (2*zeta*wn*wc)^2)
    plant_mag_at_wc = K * wn ** 2 / np.sqrt(
        (wn ** 2 - wc ** 2) ** 2 + (2 * zeta * wn * wc) ** 2
    )

    # Kp to achieve unity gain at crossover
    Kp = 1.0 / (plant_mag_at_wc + 1e-12)

    # Ki for zero steady-state error (integral action below wn)
    Ki = Kp * wn * 0.3

    # Kd for phase lead (derivative action around crossover)
    Kd = Kp * zeta / (wc + 1e-12) * 0.5

    # Limit output to prevent unreasonable commands
    output_limit = K * 20.0

    return PIDController(Kp, Ki, Kd, dt, output_limit=output_limit)


# ═════════════════════════════════════════════
#  CLOSED-LOOP SIMULATION
# ═════════════════════════════════════════════
def simulate_closed_loop_modes(basis_model, reference_signal_func, time):
    """
    Simulate closed-loop control of each SVD mode independently.

    Parameters
    ----------
    basis_model : SurfaceBasisModel
    reference_signal_func : callable(t) -> ndarray (n_modes,)
        Returns the desired modal coordinates at time t.
    time : ndarray (n_frames,)

    Returns
    -------
    q_cl : ndarray (n_modes, n_frames)   — closed-loop modal coordinates
    q_ol : ndarray (n_modes, n_frames)   — open-loop (plant-only) modal coords
    u_cl : ndarray (n_modes, n_frames)   — control signals
    q_ref : ndarray (n_modes, n_frames)  — reference signals
    """
    dt = basis_model.dt
    n_modes = basis_model.n_modes
    n_frames = len(time)

    q_cl = np.zeros((n_modes, n_frames))
    q_ol = np.zeros((n_modes, n_frames))
    u_cl = np.zeros((n_modes, n_frames))
    q_ref = np.zeros((n_modes, n_frames))

    for m in range(n_modes):
        p = basis_model.plants[m]
        plant_cl = SecondOrderPlant(p["wn"], p["zeta"], p["K"], dt)
        plant_ol = SecondOrderPlant(p["wn"], p["zeta"], p["K"], dt)
        controller = design_pid_for_mode(p, dt)

        for k in range(n_frames):
            ref = reference_signal_func(time[k])[m]
            q_ref[m, k] = ref

            # Closed loop: controller drives plant
            error = ref - q_cl[m, k - 1] if k > 0 else ref
            u = controller.step(error)
            u_cl[m, k] = u
            q_cl[m, k] = plant_cl.step(u)

            # Open loop: same reference fed directly as plant input
            q_ol[m, k] = plant_ol.step(ref)

    return q_cl, q_ol, u_cl, q_ref


def build_reference_from_experiment(basis_model):
    """
    Construct a reference signal generator that produces a clean sine
    wave matching the dominant frequency and amplitude of each
    experimental mode. This is the 'ideal' target the controller tracks.
    """
    modal_refs = []
    for i in range(basis_model.n_modes):
        q_exp = basis_model.modal_coords[i]
        p = basis_model.plants[i]
        freq = p["peak_freq"]
        amplitude = np.std(q_exp) * np.sqrt(2)  # sine amplitude from RMS

        # Estimate phase offset by correlating with sine/cosine at peak_freq
        t = basis_model.time
        if freq > 0.01:
            cos_ref = np.cos(2 * np.pi * freq * t)
            sin_ref = np.sin(2 * np.pi * freq * t)
            a_cos = np.mean(q_exp * cos_ref) * 2
            a_sin = np.mean(q_exp * sin_ref) * 2
            phase = np.arctan2(a_sin, a_cos)
        else:
            phase = 0.0

        modal_refs.append({
            "freq": freq,
            "amplitude": amplitude,
            "phase": phase,
        })

    def reference_func(t):
        ref = np.zeros(basis_model.n_modes)
        for i, mr in enumerate(modal_refs):
            if mr["freq"] > 0.01:
                ref[i] = mr["amplitude"] * np.cos(
                    2 * np.pi * mr["freq"] * t + mr["phase"]
                )
        return ref

    return reference_func, modal_refs


# ═════════════════════════════════════════════
#  PROCESS ONE EXPERIMENTAL FILE
# ═════════════════════════════════════════════
def process_file(csv_path):
    """
    Load experimental data, build basis model, simulate closed-loop,
    return all results for plotting.
    """
    label = extract_test_label(csv_path)
    meta, time, markers = parse_optitrack_csv(csv_path)
    fs = meta["frame_rate"]
    dt = 1.0 / fs
    marker_names = list(markers.keys())
    n_markers = len(marker_names)

    # Build Z matrix
    Z = np.zeros((n_markers, len(time)))
    for idx, name in enumerate(marker_names):
        Z[idx, :] = markers[name]["z"]

    # Handle NaNs: forward-fill then backward-fill per marker
    for i in range(n_markers):
        z = Z[i, :]
        nans = np.isnan(z)
        if nans.any():
            valid = ~nans
            if valid.sum() > 1:
                z[nans] = np.interp(
                    np.flatnonzero(nans), np.flatnonzero(valid), z[valid]
                )
            Z[i, :] = z

    # Build basis model
    basis = SurfaceBasisModel(Z, time, fs, energy_threshold=0.995)

    # Build reference signal from experimental data
    ref_func, modal_refs = build_reference_from_experiment(basis)

    # Simulate
    q_cl, q_ol, u_cl, q_ref = simulate_closed_loop_modes(basis, ref_func, time)

    # Reconstruct surfaces
    Z_exp = Z  # experimental (open loop, real)
    Z_ol_sim = basis.reconstruct(q_ol)  # simulated open loop
    Z_cl_sim = basis.reconstruct(q_cl)  # simulated closed loop
    Z_ref = basis.reconstruct(q_ref)     # ideal reference surface

    # Error metrics
    exp_error = Z_exp - Z_ref
    ol_error = Z_ol_sim - Z_ref
    cl_error = Z_cl_sim - Z_ref

    rms_exp = np.sqrt(np.mean(exp_error ** 2))
    rms_ol = np.sqrt(np.mean(ol_error ** 2))
    rms_cl = np.sqrt(np.mean(cl_error ** 2))

    return {
        "filepath": csv_path,
        "label": label,
        "sort_key": extract_sort_key(csv_path),
        "time": time,
        "fs": fs,
        "marker_names": marker_names,
        "markers": markers,
        "basis": basis,
        "modal_refs": modal_refs,
        "q_exp": basis.modal_coords,
        "q_cl": q_cl,
        "q_ol": q_ol,
        "q_ref": q_ref,
        "u_cl": u_cl,
        "Z_exp": Z_exp,
        "Z_ol_sim": Z_ol_sim,
        "Z_cl_sim": Z_cl_sim,
        "Z_ref": Z_ref,
        "rms_exp_error": rms_exp,
        "rms_ol_error": rms_ol,
        "rms_cl_error": rms_cl,
    }


# ═════════════════════════════════════════════
#  PER-FILE PLOTS
# ═════════════════════════════════════════════
def plot_single_results(result, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    t = result["time"]
    fs = result["fs"]
    label = result["label"]
    basis = result["basis"]
    n_modes = basis.n_modes

    # Determine a good window: 3 periods of the dominant frequency
    dom_freq = basis.plants[0]["peak_freq"]
    if dom_freq > 0.01:
        t_window = 3.0 / dom_freq
    else:
        t_window = min(30.0, t[-1])
    n_window = min(len(t), int(t_window * fs))
    tw = t[:n_window]

    colors = {"exp": "#2196F3", "ol_sim": "#FF9800", "cl_sim": "#4CAF50", "ref": "#9E9E9E"}

    # ── Fig 1: Modal coordinates comparison (first 3 modes) ──
    n_plot_modes = min(n_modes, 4)
    fig, axes = plt.subplots(n_plot_modes, 1, figsize=(14, 3.5 * n_plot_modes), sharex=True)
    if n_plot_modes == 1:
        axes = [axes]

    for i in range(n_plot_modes):
        ax = axes[i]
        ax.plot(tw, result["q_ref"][i, :n_window], color=colors["ref"],
                lw=1.5, ls="--", label="Reference", alpha=0.7)
        ax.plot(tw, result["q_exp"][i, :n_window], color=colors["exp"],
                lw=1.0, label="Experiment (OL)", alpha=0.8)
        ax.plot(tw, result["q_ol"][i, :n_window], color=colors["ol_sim"],
                lw=1.0, label="Sim Open-Loop", alpha=0.8)
        ax.plot(tw, result["q_cl"][i, :n_window], color=colors["cl_sim"],
                lw=1.2, label="Sim Closed-Loop", alpha=0.9)

        energy_pct = basis.mode_energy[i] * 100
        ax.set_ylabel(f"Mode {i}\n({energy_pct:.1f}% energy)")
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(loc="upper right", fontsize=8, ncol=4)
            ax.set_title(f"[{label}] Modal Coordinates: Experiment vs Simulation")

    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    fig.savefig(out_dir / "01_modal_comparison.png", dpi=150)
    plt.close(fig)

    # ── Fig 2: Surface Z at a representative marker ──
    mid_marker = len(result["marker_names"]) // 2
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    ax = axes[0]
    ax.plot(tw, result["Z_ref"][mid_marker, :n_window], color=colors["ref"],
            lw=1.5, ls="--", label="Reference", alpha=0.7)
    ax.plot(tw, result["Z_exp"][mid_marker, :n_window], color=colors["exp"],
            lw=1.0, label="Experiment (OL)", alpha=0.8)
    ax.plot(tw, result["Z_ol_sim"][mid_marker, :n_window], color=colors["ol_sim"],
            lw=1.0, label="Sim Open-Loop", alpha=0.8)
    ax.plot(tw, result["Z_cl_sim"][mid_marker, :n_window], color=colors["cl_sim"],
            lw=1.2, label="Sim Closed-Loop", alpha=0.9)
    ax.set_ylabel("Z Position (mm)")
    ax.set_title(f"[{label}] Marker {result['marker_names'][mid_marker]} — Z Tracking")
    ax.legend(fontsize=8, ncol=4)
    ax.grid(True, alpha=0.3)

    # Tracking error
    ax = axes[1]
    err_exp = result["Z_exp"][mid_marker, :n_window] - result["Z_ref"][mid_marker, :n_window]
    err_ol = result["Z_ol_sim"][mid_marker, :n_window] - result["Z_ref"][mid_marker, :n_window]
    err_cl = result["Z_cl_sim"][mid_marker, :n_window] - result["Z_ref"][mid_marker, :n_window]
    ax.plot(tw, err_exp, color=colors["exp"], lw=0.8, label="Experiment Error", alpha=0.7)
    ax.plot(tw, err_ol, color=colors["ol_sim"], lw=0.8, label="Sim OL Error", alpha=0.7)
    ax.plot(tw, err_cl, color=colors["cl_sim"], lw=1.0, label="Sim CL Error", alpha=0.8)
    ax.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax.set_ylabel("Tracking Error (mm)")
    ax.set_xlabel("Time (s)")
    ax.legend(fontsize=8, ncol=3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_dir / "02_marker_tracking.png", dpi=150)
    plt.close(fig)

    # ── Fig 3: Control effort per mode ──
    fig, axes = plt.subplots(n_plot_modes, 1, figsize=(14, 3.0 * n_plot_modes), sharex=True)
    if n_plot_modes == 1:
        axes = [axes]

    for i in range(n_plot_modes):
        ax = axes[i]
        ax.plot(tw, result["u_cl"][i, :n_window], color="#E91E63", lw=0.7, alpha=0.8)
        ax.set_ylabel(f"Mode {i}\nControl u")
        ax.grid(True, alpha=0.3)
    axes[0].set_title(f"[{label}] Control Effort per Mode")
    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    fig.savefig(out_dir / "03_control_effort.png", dpi=150)
    plt.close(fig)

    # ── Fig 4: RMS error across all markers ──
    n_markers = result["Z_exp"].shape[0]
    rms_per_marker_exp = np.sqrt(np.mean(
        (result["Z_exp"] - result["Z_ref"]) ** 2, axis=1))
    rms_per_marker_ol = np.sqrt(np.mean(
        (result["Z_ol_sim"] - result["Z_ref"]) ** 2, axis=1))
    rms_per_marker_cl = np.sqrt(np.mean(
        (result["Z_cl_sim"] - result["Z_ref"]) ** 2, axis=1))

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(n_markers)
    w = 0.25
    ax.bar(x - w, rms_per_marker_exp, w, color=colors["exp"],
           alpha=0.7, label="Experiment (OL)")
    ax.bar(x, rms_per_marker_ol, w, color=colors["ol_sim"],
           alpha=0.7, label="Sim Open-Loop")
    ax.bar(x + w, rms_per_marker_cl, w, color=colors["cl_sim"],
           alpha=0.7, label="Sim Closed-Loop")

    short_names = [n.replace("Unlabeled ", "") for n in result["marker_names"]]
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=90, fontsize=7)
    ax.set_ylabel("RMS Tracking Error (mm)")
    ax.set_title(f"[{label}] Per-Marker RMS Error vs Reference")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "04_rms_per_marker.png", dpi=150)
    plt.close(fig)

    # ── Fig 5: PSD comparison for dominant mode ──
    fig, ax = plt.subplots(figsize=(12, 5))
    nperseg = min(len(t), int(8 * fs))
    for name, q, color in [
        ("Experiment (OL)", result["q_exp"][0], colors["exp"]),
        ("Sim Open-Loop", result["q_ol"][0], colors["ol_sim"]),
        ("Sim Closed-Loop", result["q_cl"][0], colors["cl_sim"]),
        ("Reference", result["q_ref"][0], colors["ref"]),
    ]:
        freqs, psd = sig.welch(q, fs=fs, nperseg=nperseg)
        ax.semilogy(freqs, psd, color=color, lw=1.5, label=name, alpha=0.8)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("PSD (mm^2/Hz)")
    ax.set_title(f"[{label}] Mode 0 PSD — Experiment vs Simulation")
    ax.set_xlim(0, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "05_psd_comparison.png", dpi=150)
    plt.close(fig)


# ═════════════════════════════════════════════
#  CROSS-FILE COMPARISON PLOTS
# ═════════════════════════════════════════════
def plot_cross_comparison(all_results, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = len(all_results)
    labels = [r["label"] for r in all_results]
    x = np.arange(n)
    colors_bar = plt.cm.tab10(np.linspace(0, 1, max(n, 3)))

    # ── C1: RMS error reduction bar chart ──
    fig, ax = plt.subplots(figsize=(max(8, 2.5 * n), 6))
    w = 0.25
    rms_exp = [r["rms_exp_error"] for r in all_results]
    rms_ol = [r["rms_ol_error"] for r in all_results]
    rms_cl = [r["rms_cl_error"] for r in all_results]

    ax.bar(x - w, rms_exp, w, color="#2196F3", alpha=0.7, label="Experiment (OL)")
    ax.bar(x, rms_ol, w, color="#FF9800", alpha=0.7, label="Sim Open-Loop")
    ax.bar(x + w, rms_cl, w, color="#4CAF50", alpha=0.7, label="Sim Closed-Loop")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Surface RMS Error (mm)")
    ax.set_title("Tracking Error: Open-Loop Experiment vs Closed-Loop Simulation")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "C1_rms_error_comparison.png", dpi=150)
    plt.close(fig)

    # ── C2: Error reduction percentage ──
    fig, ax = plt.subplots(figsize=(max(8, 2.5 * n), 5))
    reduction_vs_exp = [
        100 * (1 - r["rms_cl_error"] / r["rms_exp_error"])
        if r["rms_exp_error"] > 1e-10 else 0
        for r in all_results
    ]
    reduction_vs_ol = [
        100 * (1 - r["rms_cl_error"] / r["rms_ol_error"])
        if r["rms_ol_error"] > 1e-10 else 0
        for r in all_results
    ]

    ax.bar(x - 0.18, reduction_vs_exp, 0.35, color="#2196F3", alpha=0.7,
           label="CL vs Experiment OL")
    ax.bar(x + 0.18, reduction_vs_ol, 0.35, color="#FF9800", alpha=0.7,
           label="CL vs Sim OL")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("RMS Error Reduction (%)")
    ax.set_title("Closed-Loop Error Reduction Across Tests")
    ax.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "C2_error_reduction.png", dpi=150)
    plt.close(fig)

    # ── C3: Model parameters across tests ──
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Natural frequency of mode 0
    wn_vals = [r["basis"].plants[0]["wn"] / (2 * np.pi) for r in all_results]
    axes[0].bar(x, wn_vals, color=colors_bar[:n], alpha=0.7, edgecolor="k")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[0].set_ylabel("Natural Freq (Hz)")
    axes[0].set_title("Mode 0 — Identified fn")
    axes[0].grid(True, alpha=0.3, axis="y")

    # Damping ratio
    zeta_vals = [r["basis"].plants[0]["zeta"] for r in all_results]
    axes[1].bar(x, zeta_vals, color=colors_bar[:n], alpha=0.7, edgecolor="k")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[1].set_ylabel("Damping Ratio")
    axes[1].set_title("Mode 0 — Identified zeta")
    axes[1].grid(True, alpha=0.3, axis="y")

    # Gain
    K_vals = [r["basis"].plants[0]["K"] for r in all_results]
    axes[2].bar(x, K_vals, color=colors_bar[:n], alpha=0.7, edgecolor="k")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[2].set_ylabel("Modal Gain K (mm)")
    axes[2].set_title("Mode 0 — Identified K")
    axes[2].grid(True, alpha=0.3, axis="y")

    plt.suptitle("Identified Plant Parameters Across Tests", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "C3_plant_parameters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── C4: Mode energy distribution per test ──
    fig, ax = plt.subplots(figsize=(max(8, 2.5 * n), 6))
    max_modes = max(r["basis"].n_modes for r in all_results)
    bottom = np.zeros(n)
    mode_colors = plt.cm.Set2(np.linspace(0, 1, max_modes))

    for m in range(max_modes):
        energies = []
        for r in all_results:
            if m < r["basis"].n_modes:
                energies.append(r["basis"].mode_energy[m] * 100)
            else:
                energies.append(0)
        ax.bar(x, energies, bottom=bottom, color=mode_colors[m],
               label=f"Mode {m}", edgecolor="white", lw=0.5)
        bottom += np.array(energies)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Energy (%)")
    ax.set_title("SVD Mode Energy Distribution Across Tests")
    ax.legend(fontsize=8, loc="upper right", ncol=min(max_modes, 5))
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "C4_mode_energy.png", dpi=150)
    plt.close(fig)

    # ── C5: PSD overlay — mode 0 closed-loop vs experiment ──
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    for i, r in enumerate(all_results):
        nperseg = min(len(r["time"]), int(8 * r["fs"]))
        f_exp, p_exp = sig.welch(r["q_exp"][0], fs=r["fs"], nperseg=nperseg)
        f_cl, p_cl = sig.welch(r["q_cl"][0], fs=r["fs"], nperseg=nperseg)
        axes[0].semilogy(f_exp, p_exp, color=colors_bar[i], lw=1.5,
                          label=r["label"], alpha=0.8)
        axes[1].semilogy(f_cl, p_cl, color=colors_bar[i], lw=1.5,
                          label=r["label"], alpha=0.8)

    axes[0].set_title("Experiment (Open-Loop) — Mode 0 PSD")
    axes[1].set_title("Sim Closed-Loop — Mode 0 PSD")
    for ax in axes:
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("PSD (mm^2/Hz)")
        ax.set_xlim(0, 2.0)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "C5_psd_ol_vs_cl.png", dpi=150)
    plt.close(fig)

    print(f"\nComparison plots saved to {out_dir.resolve()}/")


# ═════════════════════════════════════════════
#  CONSOLE REPORT
# ═════════════════════════════════════════════
def print_report(all_results):
    sep = "=" * 80

    for r in all_results:
        print(f"\n{sep}")
        print(f"  {r['label']}  ({Path(r['filepath']).name})")
        print(sep)
        print(r["basis"].summary())
        print()
        print(f"  RMS Surface Tracking Error:")
        print(f"    Experiment (open-loop) :  {r['rms_exp_error']:.4f} mm")
        print(f"    Simulated open-loop    :  {r['rms_ol_error']:.4f} mm")
        print(f"    Simulated closed-loop  :  {r['rms_cl_error']:.4f} mm")

        if r["rms_exp_error"] > 1e-10:
            pct = 100 * (1 - r["rms_cl_error"] / r["rms_exp_error"])
            print(f"    CL reduction vs exp    :  {pct:.1f}%")
        if r["rms_ol_error"] > 1e-10:
            pct = 100 * (1 - r["rms_cl_error"] / r["rms_ol_error"])
            print(f"    CL reduction vs sim OL :  {pct:.1f}%")

    if len(all_results) > 1:
        print(f"\n{'─'*80}")
        print("  CROSS-TEST COMPARISON")
        print(f"{'─'*80}")
        print(f"  {'Test':<25s} {'RMS Exp':>10s} {'RMS OL':>10s} {'RMS CL':>10s} "
              f"{'CL Reduction':>14s} {'Modes':>6s}")
        print("  " + "-" * 77)
        for r in all_results:
            pct = 100 * (1 - r["rms_cl_error"] / r["rms_exp_error"]) \
                if r["rms_exp_error"] > 1e-10 else 0
            print(
                f"  {r['label']:<25s} {r['rms_exp_error']:10.4f} "
                f"{r['rms_ol_error']:10.4f} {r['rms_cl_error']:10.4f} "
                f"{pct:13.1f}% {r['basis'].n_modes:6d}"
            )


# ═════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="Simulate closed-loop surface control and compare with experiments."
    )
    parser.add_argument(
        "files", nargs="+",
        help="One or more OptiTrack CSV files"
    )
    parser.add_argument(
        "-o", "--output", default="closed_loop_analysis",
        help="Output directory (default: closed_loop_analysis)"
    )
    parser.add_argument(
        "-m", "--modes", type=int, default=None,
        help="Number of SVD modes to retain (default: auto via 99.5%% energy)"
    )
    args = parser.parse_args()

    # Expand globs
    csv_paths = []
    for pattern in args.files:
        expanded = sorted(glob.glob(pattern))
        if expanded:
            csv_paths.extend(expanded)
        else:
            csv_paths.append(pattern)

    csv_paths = sorted(set(csv_paths), key=extract_sort_key)
    print(f"Processing {len(csv_paths)} file(s)...\n")

    out_root = Path(args.output)
    all_results = []

    for csv_path in csv_paths:
        print(f"{'─'*60}")
        print(f"  {csv_path}")
        result = process_file(csv_path)
        all_results.append(result)

        # Per-file plots
        safe_name = Path(csv_path).stem
        plot_single_results(result, out_root / safe_name)
        print(f"  Plots -> {out_root / safe_name}/")

    # Console report
    print_report(all_results)

    # Cross-file comparison
    if len(all_results) >= 1:
        plot_cross_comparison(all_results, out_root / "comparison")

    print(f"\nDone. All output in: {out_root.resolve()}/")


if __name__ == "__main__":
    main()