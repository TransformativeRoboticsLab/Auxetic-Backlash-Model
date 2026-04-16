"""
analyze_surface_dynamics.py

Analyze and compare the dynamics of OptiTrack marker grids oscillating in Z
across one or more test files at different excitation frequencies.

Handles momentary marker dropouts via NaN-aware cubic interpolation.

Usage:
    # Single file
    python analyze_surface_dynamics.py data/test_12s.csv

    # Multiple files
    python analyze_surface_dynamics.py data/test_4s.csv data/test_8s.csv data/test_12s.csv

    # Glob pattern
    python analyze_surface_dynamics.py data/RADs_*.csv

    # Custom output directory
    python analyze_surface_dynamics.py data/*.csv -o results/
"""

import sys
import glob
import re
import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import signal, interpolate
from pathlib import Path


# ─────────────────────────────────────────────
#  1. PARSE OPTITRACK CSV
# ─────────────────────────────────────────────
def parse_optitrack_csv(filepath):
    """
    Read a Motive CSV export and return structured marker data.

    Returns
    -------
    meta : dict          — capture metadata (frame rate, units, etc.)
    time : ndarray (N,)  — timestamps in seconds
    markers : dict        — {name: {"x": (N,), "y": (N,), "z": (N,)}}
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Row 0: metadata key-value pairs
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
    meta["filepath"] = str(filepath)

    # Row 3: marker names (each repeated 3x for X, Y, Z)
    name_cells = lines[3].strip().replace("\r", "").split(",")
    raw_names = [c.strip().strip('"') for c in name_cells[2:]]

    # Build ordered unique marker list
    marker_order = []
    seen = set()
    for i in range(0, len(raw_names), 3):
        name = raw_names[i]
        if name and name not in seen:
            marker_order.append(name)
            seen.add(name)

    # Parse numeric data (row 7 onward)
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


# ─────────────────────────────────────────────
#  2. EXTRACT A HUMAN-READABLE TEST LABEL
# ─────────────────────────────────────────────
def extract_test_label(filepath):
    """
    Pull a short label from the filename for plot legends.
    Tries to parse period/amplitude from naming conventions like
    'RADs_12s_period_1440deg_amplitude.csv' -> 'T=12s / A=1440deg'
    Falls back to the stem if no pattern matches.
    """
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

    # Fallback: clean up the stem
    label = stem.replace("_", " ").strip()
    if len(label) > 40:
        label = label[:37] + "..."
    return label


def extract_sort_key(filepath):
    """Numeric sort key so files appear in order of period or amplitude."""
    stem = Path(filepath).stem
    period_match = re.search(r"(\d+\.?\d*)s[_\s]*period", stem, re.IGNORECASE)
    if period_match:
        return float(period_match.group(1))
    amp_match = re.search(r"(\d+\.?\d*)deg[_\s]*amplitude", stem, re.IGNORECASE)
    if amp_match:
        return float(amp_match.group(1))
    return 0.0


# ─────────────────────────────────────────────
#  3. GAP DETECTION AND INTERPOLATION
# ─────────────────────────────────────────────
def diagnose_gaps(signal_array):
    """Return list of (start_index, length) for each NaN gap."""
    is_nan = np.isnan(signal_array)
    gaps = []
    if not np.any(is_nan):
        return gaps

    in_gap = False
    start = 0
    for i, v in enumerate(is_nan):
        if v and not in_gap:
            in_gap = True
            start = i
        elif not v and in_gap:
            gaps.append((start, i - start))
            in_gap = False
    if in_gap:
        gaps.append((start, len(signal_array) - start))
    return gaps


def interpolate_gaps(time, signal_array, max_gap_samples=60):
    """
    Fill NaN gaps with cubic interpolation when the gap is short enough.
    Gaps longer than max_gap_samples remain NaN — they represent genuine
    tracking loss, not momentary flickers.
    """
    out = signal_array.copy()
    valid = ~np.isnan(out)
    if valid.sum() < 4:
        return out

    gaps = diagnose_gaps(out)
    for start, length in gaps:
        if length <= max_gap_samples:
            interp_fn = interpolate.interp1d(
                time[valid], out[valid], kind="cubic", fill_value="extrapolate"
            )
            out[start : start + length] = interp_fn(time[start : start + length])

    return out


# ─────────────────────────────────────────────
#  4. PER-MARKER SPECTRAL + KINEMATIC ANALYSIS
# ─────────────────────────────────────────────
def analyze_marker_z(time, z_raw, fs, marker_name):
    """
    For a single marker's Z trace:
      - interpolate short gaps
      - remove mean (DC offset)
      - estimate dominant frequency via Welch PSD
      - compute amplitude, velocity, acceleration

    Returns a result dict or None if data is too sparse.
    """
    z_interp = interpolate_gaps(time, z_raw, max_gap_samples=int(0.5 * fs))
    valid_frac = np.sum(~np.isnan(z_interp)) / len(z_interp)
    if valid_frac < 0.5:
        return None

    z_clean = z_interp.copy()
    mean_z = np.nanmean(z_clean)
    z_clean[np.isnan(z_clean)] = mean_z
    z_detrend = z_clean - mean_z

    # Welch PSD
    nperseg = min(len(z_detrend), int(8 * fs))
    freqs, psd = signal.welch(
        z_detrend, fs=fs, nperseg=nperseg, noverlap=nperseg // 2
    )

    # Dominant frequency (ignore DC bin)
    freq_mask = freqs > 0.01
    peak_idx = np.argmax(psd[freq_mask])
    dominant_freq = freqs[freq_mask][peak_idx]
    dominant_period = 1.0 / dominant_freq if dominant_freq > 0 else np.inf

    # Amplitude from RMS (sine wave: amplitude = sqrt(2) * RMS)
    rms = np.sqrt(np.nanmean(z_detrend**2))
    amplitude = rms * np.sqrt(2)
    peak_to_peak = np.nanmax(z_detrend) - np.nanmin(z_detrend)

    # Velocity, acceleration via central differences
    dt = 1.0 / fs
    velocity = np.gradient(z_clean, dt)
    acceleration = np.gradient(velocity, dt)

    return {
        "name": marker_name,
        "mean_z": mean_z,
        "amplitude": amplitude,
        "peak_to_peak": peak_to_peak,
        "dominant_freq_hz": dominant_freq,
        "dominant_period_s": dominant_period,
        "rms": rms,
        "max_velocity": np.max(np.abs(velocity)),
        "max_acceleration": np.max(np.abs(acceleration)),
        "valid_fraction": valid_frac,
        "freqs": freqs,
        "psd": psd,
        "z_detrend": z_detrend,
        "z_interp": z_interp,
        "velocity": velocity,
        "acceleration": acceleration,
    }


# ─────────────────────────────────────────────
#  5. SURFACE-LEVEL PHASE ANALYSIS
# ─────────────────────────────────────────────
def analyze_phase_coherence(results, fs):
    """
    Cross-correlate all marker pairs at the dominant frequency to
    detect phase differences across the surface.
    """
    names = [r["name"] for r in results]
    n = len(results)
    phase_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            zi = results[i]["z_detrend"]
            zj = results[j]["z_detrend"]
            f, Pxy = signal.csd(
                zi, zj, fs=fs, nperseg=min(len(zi), int(8 * fs))
            )
            dom_f = results[i]["dominant_freq_hz"]
            fidx = np.argmin(np.abs(f - dom_f))
            phase_diff = np.angle(Pxy[fidx], deg=True)
            phase_matrix[i, j] = phase_diff
            phase_matrix[j, i] = -phase_diff

    return names, phase_matrix


# ─────────────────────────────────────────────
#  6. ANALYZE ONE FILE -> RETURN SUMMARY DICT
# ─────────────────────────────────────────────
def analyze_single_file(csv_path):
    """
    Run the full per-marker analysis on one CSV.
    Returns a summary dict with everything needed for
    per-file plots and cross-file comparison.
    """
    label = extract_test_label(csv_path)
    print(f"\n{'─'*60}")
    print(f"  Loading: {csv_path}")
    print(f"  Label:   {label}")

    meta, time, markers = parse_optitrack_csv(csv_path)
    fs = meta["frame_rate"]

    print(f"  Rate={fs} Hz  Duration={time[-1]-time[0]:.1f}s  "
          f"Frames={len(time)}  Markers={len(markers)}")

    results = []
    skipped = []
    for name in markers:
        r = analyze_marker_z(time, markers[name]["z"], fs, name)
        if r is not None:
            results.append(r)
        else:
            skipped.append(name)

    if skipped:
        print(f"  Skipped (too sparse): {', '.join(skipped)}")

    if not results:
        print("  No usable markers — skipping file entirely.")
        return None

    # Aggregate stats
    freqs = [r["dominant_freq_hz"] for r in results]
    amps = [r["amplitude"] for r in results]
    p2ps = [r["peak_to_peak"] for r in results]
    max_vels = [r["max_velocity"] for r in results]
    max_accs = [r["max_acceleration"] for r in results]

    # Phase coherence
    names_list, phase_mat = analyze_phase_coherence(results, fs)
    upper = phase_mat[np.triu_indices(len(names_list), k=1)]

    # Gap stats
    total_gaps = 0
    total_gap_samples = 0
    for name in markers:
        gaps = diagnose_gaps(markers[name]["z"])
        total_gaps += len(gaps)
        total_gap_samples += sum(g[1] for g in gaps)

    summary = {
        "filepath": csv_path,
        "label": label,
        "sort_key": extract_sort_key(csv_path),
        "meta": meta,
        "time": time,
        "markers": markers,
        "results": results,
        "fs": fs,
        "dominant_freq_mean": np.mean(freqs),
        "dominant_freq_std": np.std(freqs),
        "dominant_period_mean": 1.0 / np.mean(freqs) if np.mean(freqs) > 0 else np.inf,
        "amplitude_mean": np.mean(amps),
        "amplitude_std": np.std(amps),
        "amplitude_all": np.array(amps),
        "p2p_mean": np.mean(p2ps),
        "p2p_std": np.std(p2ps),
        "max_vel_mean": np.mean(max_vels),
        "max_vel_std": np.std(max_vels),
        "max_acc_mean": np.mean(max_accs),
        "max_acc_std": np.std(max_accs),
        "phase_mean_abs": np.mean(np.abs(upper)) if len(upper) > 0 else 0,
        "phase_max_abs": np.max(np.abs(upper)) if len(upper) > 0 else 0,
        "phase_matrix": phase_mat,
        "phase_names": names_list,
        "total_gaps": total_gaps,
        "total_gap_samples": total_gap_samples,
        "total_possible_samples": len(time) * len(markers),
        "n_markers": len(markers),
        "n_markers_valid": len(results),
    }

    return summary


# ─────────────────────────────────────────────
#  7. PER-FILE DIAGNOSTIC PLOTS
# ─────────────────────────────────────────────
def plot_single_file(summary, out_dir):
    """Generate the 6 diagnostic figures for one test file."""
    s = summary
    time = s["time"]
    results = s["results"]
    markers = s["markers"]
    fs = s["fs"]
    label = s["label"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmap = plt.cm.turbo(np.linspace(0, 1, len(results)))

    # --- Z traces ---
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for i, r in enumerate(results):
        axes[0].plot(time, markers[r["name"]]["z"], color=cmap[i], alpha=0.5, lw=0.4)
        axes[1].plot(time, r["z_detrend"], color=cmap[i], alpha=0.6, lw=0.4)
    axes[0].set_ylabel("Raw Z (mm)")
    axes[0].set_title(f"[{label}] Raw Marker Z Positions")
    axes[0].grid(True, alpha=0.3)
    axes[1].set_ylabel("Z - Mean (mm)")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_title(f"[{label}] Mean-Subtracted Z")
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "01_z_traces.png", dpi=150)
    plt.close(fig)

    # --- PSD ---
    fig, ax = plt.subplots(figsize=(12, 5))
    for i, r in enumerate(results):
        ax.semilogy(r["freqs"], r["psd"], color=cmap[i], alpha=0.6, lw=0.8)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("PSD (mm^2/Hz)")
    ax.set_title(f"[{label}] Z-Axis Power Spectral Density")
    ax.set_xlim(0, 2.0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "02_psd.png", dpi=150)
    plt.close(fig)

    # --- Spatial amplitude ---
    fig, ax = plt.subplots(figsize=(8, 7))
    xs = [np.nanmean(markers[r["name"]]["x"]) for r in results]
    ys = [np.nanmean(markers[r["name"]]["y"]) for r in results]
    amps = [r["amplitude"] for r in results]
    sc = ax.scatter(xs, ys, c=amps, s=200, cmap="plasma", edgecolors="k", lw=0.8)
    for r, x, y in zip(results, xs, ys):
        ax.annotate(f'{r["amplitude"]:.2f}', (x, y),
                    textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=7)
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(f"[{label}] Z-Oscillation Amplitude (mm)")
    plt.colorbar(sc, label="Amplitude (mm)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "03_spatial_amplitude.png", dpi=150)
    plt.close(fig)

    # --- Kinematics (first 3 periods of representative marker) ---
    rep = results[0]
    n_show = min(len(time), int(3 * s["dominant_period_mean"] * fs))
    t_short = time[:n_show]
    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
    axes[0].plot(t_short, rep["z_detrend"][:n_show], "k", lw=0.8)
    axes[0].set_ylabel("Z - Mean (mm)")
    axes[0].set_title(f"[{label}] Kinematics — {rep['name']}")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(t_short, rep["velocity"][:n_show], "tab:blue", lw=0.8)
    axes[1].set_ylabel("Velocity (mm/s)")
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(t_short, rep["acceleration"][:n_show], "tab:red", lw=0.8)
    axes[2].set_ylabel("Acceleration (mm/s^2)")
    axes[2].set_xlabel("Time (s)")
    axes[2].grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "04_kinematics.png", dpi=150)
    plt.close(fig)

    # --- Phase coherence ---
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(s["phase_matrix"], cmap="RdBu", vmin=-180, vmax=180, aspect="equal")
    short = [n.replace("Unlabeled ", "") for n in s["phase_names"]]
    ax.set_xticks(range(len(short)))
    ax.set_yticks(range(len(short)))
    ax.set_xticklabels(short, rotation=90, fontsize=7)
    ax.set_yticklabels(short, fontsize=7)
    ax.set_title(f"[{label}] Phase Difference at Dominant Freq (deg)")
    plt.colorbar(im, label="Phase (deg)")
    plt.tight_layout()
    fig.savefig(out_dir / "05_phase_coherence.png", dpi=150)
    plt.close(fig)

    # --- Gap report ---
    fig, ax = plt.subplots(figsize=(14, 5))
    for i, name in enumerate(markers):
        z = markers[name]["z"]
        gaps = diagnose_gaps(z)
        for start, length in gaps:
            ax.barh(i, length / fs, left=time[start], height=0.8, color="red", alpha=0.7)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels([n.replace("Unlabeled ", "") for n in markers], fontsize=7)
    ax.set_xlabel("Time (s)")
    ax.set_title(f"[{label}] Marker Dropout Gaps")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(out_dir / "06_gap_report.png", dpi=150)
    plt.close(fig)


# ─────────────────────────────────────────────
#  8. CROSS-TEST COMPARISON PLOTS
# ─────────────────────────────────────────────
def plot_comparison(summaries, out_dir):
    """
    Generate comparison figures that overlay or juxtapose metrics
    from all test files side by side.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_tests = len(summaries)
    labels = [s["label"] for s in summaries]
    colors = plt.cm.tab10(np.linspace(0, 1, max(n_tests, 3)))

    # ── Fig C1: Amplitude box plot across tests ──
    fig, ax = plt.subplots(figsize=(max(6, 2 * n_tests), 6))
    amp_data = [s["amplitude_all"] for s in summaries]
    bp = ax.boxplot(amp_data, patch_artist=True, widths=0.6)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel("Z Amplitude (mm)")
    ax.set_title("Amplitude Distribution — All Markers per Test")
    ax.grid(True, alpha=0.3, axis="y")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    fig.savefig(out_dir / "C1_amplitude_comparison.png", dpi=150)
    plt.close(fig)

    # ── Fig C2: PSD overlay (mean across markers per test) ──
    fig, ax = plt.subplots(figsize=(12, 6))
    for i, s in enumerate(summaries):
        all_psds = np.array([r["psd"] for r in s["results"]])
        mean_psd = np.mean(all_psds, axis=0)
        freqs = s["results"][0]["freqs"]
        ax.semilogy(freqs, mean_psd, color=colors[i], lw=1.8,
                     label=s["label"], alpha=0.85)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Mean PSD (mm^2/Hz)")
    ax.set_title("Mean Z-Axis PSD — Comparison Across Tests")
    ax.set_xlim(0, 3.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "C2_psd_overlay.png", dpi=150)
    plt.close(fig)

    # ── Fig C3: Summary bar charts (freq, amplitude, velocity, acceleration) ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    x_pos = np.arange(n_tests)

    ax = axes[0, 0]
    vals = [s["dominant_freq_mean"] for s in summaries]
    errs = [s["dominant_freq_std"] for s in summaries]
    ax.bar(x_pos, vals, yerr=errs, color=colors[:n_tests], alpha=0.7,
           edgecolor="k", capsize=4)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Dominant Frequency (Hz)")
    ax.set_title("Dominant Frequency")
    ax.grid(True, alpha=0.3, axis="y")

    ax = axes[0, 1]
    vals = [s["amplitude_mean"] for s in summaries]
    errs = [s["amplitude_std"] for s in summaries]
    ax.bar(x_pos, vals, yerr=errs, color=colors[:n_tests], alpha=0.7,
           edgecolor="k", capsize=4)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Amplitude (mm)")
    ax.set_title("Mean Z Amplitude")
    ax.grid(True, alpha=0.3, axis="y")

    ax = axes[1, 0]
    vals = [s["max_vel_mean"] for s in summaries]
    errs = [s["max_vel_std"] for s in summaries]
    ax.bar(x_pos, vals, yerr=errs, color=colors[:n_tests], alpha=0.7,
           edgecolor="k", capsize=4)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Max Velocity (mm/s)")
    ax.set_title("Mean Peak Velocity")
    ax.grid(True, alpha=0.3, axis="y")

    ax = axes[1, 1]
    vals = [s["max_acc_mean"] for s in summaries]
    errs = [s["max_acc_std"] for s in summaries]
    ax.bar(x_pos, vals, yerr=errs, color=colors[:n_tests], alpha=0.7,
           edgecolor="k", capsize=4)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Max Acceleration (mm/s^2)")
    ax.set_title("Mean Peak Acceleration")
    ax.grid(True, alpha=0.3, axis="y")

    plt.suptitle("Cross-Test Dynamics Comparison", fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(out_dir / "C3_summary_bars.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── Fig C4: Phase coherence comparison ──
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(6, 2 * n_tests), 5))
    means = [s["phase_mean_abs"] for s in summaries]
    maxes = [s["phase_max_abs"] for s in summaries]
    ax.bar(x_pos - width / 2, means, width, label="Mean |dphi|",
           color=colors[:n_tests], alpha=0.6, edgecolor="k")
    ax.bar(x_pos + width / 2, maxes, width, label="Max |dphi|",
           color=colors[:n_tests], alpha=0.3, edgecolor="k", hatch="//")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Phase Difference (deg)")
    ax.set_title("Surface Phase Coherence — Comparison Across Tests")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "C4_phase_comparison.png", dpi=150)
    plt.close(fig)

    # ── Fig C5: Time-domain overlay (mean Z of all markers per test) ──
    fig, ax = plt.subplots(figsize=(14, 5))
    for i, s in enumerate(summaries):
        all_z = np.array([r["z_detrend"] for r in s["results"]])
        mean_z = np.mean(all_z, axis=0)
        max_period = max(ss["dominant_period_mean"] for ss in summaries)
        n_show = min(len(s["time"]), int(3 * max_period * s["fs"]))
        ax.plot(s["time"][:n_show], mean_z[:n_show],
                color=colors[i], lw=1.2, label=s["label"], alpha=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean Z - offset (mm)")
    ax.set_title("Mean Surface Z Motion — First 3 Periods of Slowest Test")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_dir / "C5_time_overlay.png", dpi=150)
    plt.close(fig)

    # ── Fig C6: Gap comparison ──
    fig, ax = plt.subplots(figsize=(max(6, 2 * n_tests), 5))
    gap_pcts = [
        100.0 * s["total_gap_samples"] / s["total_possible_samples"]
        if s["total_possible_samples"] > 0 else 0
        for s in summaries
    ]
    gap_counts = [s["total_gaps"] for s in summaries]
    ax.bar(x_pos - width / 2, gap_counts, width, label="Dropout Events",
           color="tab:orange", alpha=0.7, edgecolor="k")
    ax2 = ax.twinx()
    ax2.bar(x_pos + width / 2, gap_pcts, width, label="% Frames Lost",
            color="tab:red", alpha=0.5, edgecolor="k", hatch="//")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Dropout Events")
    ax2.set_ylabel("% Frames Lost")
    ax.set_title("Data Quality — Dropout Comparison Across Tests")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(out_dir / "C6_gap_comparison.png", dpi=150)
    plt.close(fig)

    print(f"\nComparison plots saved to {out_dir.resolve()}/")


# ─────────────────────────────────────────────
#  9. CONSOLE REPORT
# ─────────────────────────────────────────────
def print_summary(summaries):
    """Print a combined console report for all tests."""
    sep = "=" * 78

    for s in summaries:
        print(f"\n{sep}")
        print(f"  {s['label']}  ({Path(s['filepath']).name})")
        print(sep)
        print(f"  Markers: {s['n_markers']} total, {s['n_markers_valid']} valid")
        print(f"  Duration: {s['time'][-1] - s['time'][0]:.1f}s  @ {s['fs']} Hz")
        print(f"  Dominant Freq:  {s['dominant_freq_mean']:.4f} +/- {s['dominant_freq_std']:.4f} Hz  "
              f"(T = {s['dominant_period_mean']:.2f} s)")
        print(f"  Amplitude:      {s['amplitude_mean']:.3f} +/- {s['amplitude_std']:.3f} mm")
        print(f"  Peak-to-Peak:   {s['p2p_mean']:.3f} +/- {s['p2p_std']:.3f} mm")
        print(f"  Max Velocity:   {s['max_vel_mean']:.2f} +/- {s['max_vel_std']:.2f} mm/s")
        print(f"  Max Accel:      {s['max_acc_mean']:.1f} +/- {s['max_acc_std']:.1f} mm/s^2")
        print(f"  Phase Spread:   mean |dphi| = {s['phase_mean_abs']:.1f} deg, "
              f"max |dphi| = {s['phase_max_abs']:.1f} deg")
        gap_pct = 100.0 * s["total_gap_samples"] / s["total_possible_samples"] \
            if s["total_possible_samples"] > 0 else 0
        print(f"  Dropouts:       {s['total_gaps']} events, "
              f"{s['total_gap_samples']} frames ({gap_pct:.2f}%)")

    # Cross-test comparison table
    if len(summaries) > 1:
        print(f"\n{'─'*78}")
        print("  CROSS-TEST COMPARISON TABLE")
        print(f"{'─'*78}")
        header = (f"  {'Test':<22s} {'Freq(Hz)':>9s} {'Period(s)':>10s} "
                  f"{'Amp(mm)':>9s} {'P2P(mm)':>9s} {'Vmax':>9s} "
                  f"{'Amax':>10s} {'|dphi|':>7s} {'Gaps':>6s}")
        print(header)
        print("  " + "-" * 76)
        for s in summaries:
            gap_pct = 100.0 * s["total_gap_samples"] / s["total_possible_samples"] \
                if s["total_possible_samples"] > 0 else 0
            print(
                f"  {s['label']:<22s} {s['dominant_freq_mean']:9.4f} "
                f"{s['dominant_period_mean']:10.2f} "
                f"{s['amplitude_mean']:9.3f} {s['p2p_mean']:9.3f} "
                f"{s['max_vel_mean']:9.2f} {s['max_acc_mean']:10.1f} "
                f"{s['phase_mean_abs']:6.1f}  {gap_pct:5.1f}%"
            )

    # Per-marker detail table for each test
    for s in summaries:
        print(f"\n  Per-marker detail: {s['label']}")
        print(f"  {'Marker':<20s} {'Freq(Hz)':>9s} {'Amp(mm)':>9s} {'P2P(mm)':>9s} "
              f"{'Vmax(mm/s)':>11s} {'Amax(mm/s2)':>12s} {'Valid%':>7s}")
        print("  " + "-" * 70)
        for r in sorted(s["results"], key=lambda x: x["name"]):
            print(
                f"  {r['name']:<20s} {r['dominant_freq_hz']:9.4f} "
                f"{r['amplitude']:9.3f} {r['peak_to_peak']:9.3f} "
                f"{r['max_velocity']:11.2f} {r['max_acceleration']:12.1f} "
                f"{r['valid_fraction']*100:6.1f}%"
            )


# ─────────────────────────────────────────────
#  10. MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Analyze and compare OptiTrack surface dynamics across tests."
    )
    parser.add_argument(
        "files", nargs="+",
        help="One or more OptiTrack CSV files (globs supported by shell)"
    )
    parser.add_argument(
        "-o", "--output", default="surface_analysis",
        help="Output directory (default: surface_analysis)"
    )
    args = parser.parse_args()

    # Expand any remaining globs (Windows compatibility)
    csv_paths = []
    for pattern in args.files:
        expanded = sorted(glob.glob(pattern))
        if expanded:
            csv_paths.extend(expanded)
        else:
            csv_paths.append(pattern)

    # Sort by extracted numeric key (period or amplitude)
    csv_paths = sorted(set(csv_paths), key=extract_sort_key)

    print(f"Found {len(csv_paths)} file(s) to analyze.")
    out_root = Path(args.output)

    # ── Analyze each file ──
    summaries = []
    for csv_path in csv_paths:
        s = analyze_single_file(csv_path)
        if s is not None:
            summaries.append(s)
            # Per-file diagnostic plots in a subdirectory
            safe_name = Path(csv_path).stem
            plot_single_file(s, out_root / safe_name)

    if not summaries:
        print("\nNo files produced usable results.")
        sys.exit(1)

    # ── Console report ──
    print_summary(summaries)

    # ── Comparison plots (only meaningful with 2+, but safe with 1) ──
    plot_comparison(summaries, out_root / "comparison")

    print(f"\nDone. All output in: {out_root.resolve()}/")


if __name__ == "__main__":
    main()