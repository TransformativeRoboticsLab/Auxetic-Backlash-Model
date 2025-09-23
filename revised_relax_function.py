import numpy as np
import matplotlib.pyplot as plt

def _closest_point_on_segment(a, b, p):
    ab = b - a
    t = np.dot(p - a, ab) / (np.dot(ab, ab) + 1e-12)
    t = np.clip(t, 0.0, 1.0)
    return a + t * ab


def closest_point_on_polyline(poly, p):
    best, best_d2 = None, np.inf
    for i in range(len(poly)-1):
        q = _closest_point_on_segment(poly[i], poly[i+1], p)
        d2 = np.sum((q - p)**2)
        if d2 < best_d2:
            best_d2, best = d2, q
    return best

def get_lengths_and_alpha(deformed_list, L0):
    pts = np.array(deformed_list)
    seg = pts[1:] - pts[:-1]
    L = np.linalg.norm(seg, axis=1)       # segment lengths
    alpha = L / L0
    return L, alpha

def get_LSM_quality(deformed_list, target_curve):
    """
    Simple MSE distance-to-target: for each deformed node,
    measure squared distance to nearest point on the target curve.
    Compute how close the deformed curve is to the target curve.
    """
    pts = np.array(deformed_list)
    tgt = np.array(target_curve)
    errs = []
    for p in pts:
        q = closest_point_on_polyline(tgt, p)
        errs.append(np.sum((q - p) ** 2))
    return float(np.mean(errs))

def relax_to_curve(target_curve, flat_line, max_iters=2000, step_size=0.01, tol=1e-6, 
                   use_smoothing=False, smooth_factor=0.1, smooth_tol=0.001, k=1.0, 
                   length_tol=0.001, lock_x=False, length_method='spring'):
    """
    Iteratively adjust internal points of a flat line to approximate a target curve.

    Parameters:
        target_curve (list of (x, y)): Target curve as linked (x, y) tuples.
        flat_line (list of (x, y)): Initial flat line as linked (x, y) tuples.
        max_iters (int): Maximum number of iterations.
        step_size (float): Gradient descent step size.
        tol (float): Convergence tolerance (average movement).
        k (float): spring constant for length correction.
        length_tol (float): tolerance for segment length (as fraction of L0).
        lock_x (bool): if True, nodes keep original x-coordinates.
        length_method (str): 'spring', 'relu', or 'hybrid' for length enforcement.

    Returns:
        tuple: (deformed_line, L0) where deformed_line approximates the target curve.
    """
    target = np.array(target_curve)
    flat_line = np.array(flat_line)  
    current = flat_line.copy()     

    assert np.allclose(current[0], [0, 0]) and np.allclose(current[-1], [1, 0]), \
           "First and last points must be fixed at (0,0) and (1,0)."

    n = len(current)
    L0 = 1.0 / (n - 1)  # rest length of each identical link
    Lmin = -length_tol * L0 + L0# (1.0 - length_tol) * L0
    Lmax = length_tol * L0 + L0# (1.0 + length_tol) * L0

    for iteration in range(max_iters):
        prev = current.copy()

        # Lock x-coordinates if requested
        if lock_x:
            current[1:-1, 0] = flat_line[1:-1, 0]

        # Step 1: Move towards target curve
        for i in range(1, n - 1):  # Skip fixed endpoints
            target_point = closest_point_on_polyline(target, current[i])
            direction = target_point - current[i]
            current[i] += step_size * direction

        # Step 2: Enforce length constraints
        if length_method == 'spring':
            current = _enforce_lengths_spring(current, L0, k)
        elif length_method == 'relu':
            current = _enforce_lengths_relu(current, Lmin, Lmax)
        elif length_method == 'hybrid':
            current = _enforce_lengths_spring(current, L0, k)
            current = _enforce_lengths_relu(current, Lmin, Lmax)

        # Step 3: Optional smoothing
        if use_smoothing:
            current = _apply_smoothing(current, smooth_factor, smooth_tol)

        # Check convergence
        avg_move = np.mean(np.linalg.norm(current - prev, axis=1))
        if avg_move < tol:
            print(f"Converged after {iteration + 1} iterations")
            break

    return [tuple(p) for p in current], L0


def _enforce_lengths_spring(current, L0, k):
    """Apply spring forces to maintain segment lengths near L0"""
    n = len(current)
    forces = np.zeros_like(current)
    
    # Calculate spring forces for each segment
    for i in range(n - 1):
        segment = current[i + 1] - current[i]
        length = np.linalg.norm(segment)
        
        if length > 1e-12:  # Avoid division by zero
            unit_vector = segment / length
            force_magnitude = k * (length - L0)
            force = force_magnitude * unit_vector
            
            # Apply equal and opposite forces to segment endpoints
            if i > 0:  # Don't move first endpoint
                forces[i] += force
            if i < n - 2:  # Don't move last endpoint
                forces[i + 1] -= force
    
    # Apply forces (but keep endpoints fixed)
    current[1:-1] += forces[1:-1]
    return current


def _enforce_lengths_relu(current, Lmin, Lmax):
    """Hard ReLU segment lengths to be within [Lmin, Lmax]"""
    n = len(current)
    
    # Process each segment and adjust both endpoints
    for i in range(n - 1):
        segment = current[i + 1] - current[i]
        length = np.linalg.norm(segment)
        
        if length > 1e-12:  # Some small number, avoid division by zero
            target_length = np.clip(length, Lmin, Lmax)
            
            # if length and target length are not basically exact, keep modifying
            if abs(length - target_length) > 1e-12:
                unit_vector = segment / length
                correction = (target_length - length) / 2.0
                
                # Move both endpoints toward/away from each other
                # But respect the fixed endpoints constraint
                if i > 0:  # Don't move first point
                    current[i] -= correction * unit_vector
                if i < n - 2:  # Don't move last point  
                    current[i + 1] += correction * unit_vector
                
                # Special case: if one endpoint is fixed, move the other one fully
                if i == 0:  # First segment - only move second point
                    current[i + 1] = current[i] + target_length * unit_vector
                elif i == n - 2:  # Last segment - only move second-to-last point
                    current[i] = current[i + 1] - target_length * unit_vector
    return current


def _apply_smoothing(current, smooth_factor, smooth_tol):
    """Apply smoothing to reduce sharp variations in segment lengths"""
    n = len(current)
    
    # Calculate all segment lengths
    segment_lengths = np.zeros(n - 1)
    for i in range(n - 1):
        segment_lengths[i] = np.linalg.norm(current[i + 1] - current[i])
    
    # Apply smoothing corrections
    for i in range(1, n - 1):  # Only modify internal points
        left_len = segment_lengths[i - 1] if i > 0 else 0
        right_len = segment_lengths[i] if i < n - 1 else 0
        
        if i > 1:  # Check smoothness with left neighbor
            two_left_len = segment_lengths[i - 2]
            avg_len = (left_len + right_len) / 2
            
            if abs(two_left_len - avg_len) > smooth_tol:
                # Apply small correction
                left_segment = current[i] - current[i - 1]
                left_length = np.linalg.norm(left_segment)
                if left_length > 1e-12:
                    correction = smooth_factor * (two_left_len - avg_len) / avg_len
                    current[i] += correction * (left_segment / left_length)
        
        if i < n - 2:  # Check smoothness with right neighbor
            two_right_len = segment_lengths[i + 1]
            avg_len = (left_len + right_len) / 2
            
            if abs(two_right_len - avg_len) > smooth_tol:
                # Apply small correction
                right_segment = current[i + 1] - current[i]
                right_length = np.linalg.norm(right_segment)
                if right_length > 1e-12:
                    correction = smooth_factor * (two_right_len - avg_len) / avg_len
                    current[i] += correction * (right_segment / right_length)
    return current


def get_alpha_from_deformed(deformed_list, L0):
    """
    Compute the stretch factor (alpha) for each segment.
    """
    points = np.array(deformed_list)
    diffs = points[1:] - points[:-1]
    dists = np.linalg.norm(diffs, axis=1)
    alpha_list = dists / L0
    return alpha_list




num_points = 25
num_points_target = 30
x_vals = np.linspace(0, 1, num_points)
x_vals_target = np.linspace(0, 1, num_points_target)
# y_vals = np.sqrt(0.25 - (x_vals - 0.5)**2)  # Semicircle radius 0.5
y_vals_target = 1/0.2 * (0.2969*np.sqrt(x_vals_target) - 0.1260*x_vals_target - 0.3516*x_vals_target**2 + 0.2843*x_vals_target**3 - 0.1015*x_vals_target**4) 
# y_vals_target = np.sqrt(0.25 - (x_vals_target - 0.5)**2 ) /10  # Semicircle radius 0.5
target_curve = list(zip(x_vals_target, y_vals_target)) # removed the 0.1 scaling
print(target_curve)

# Flat line from (0,0) to (1,0)
flat_line = list(zip(x_vals, np.zeros_like(x_vals)))

# Run the relaxation
deformed, L0 = relax_to_curve(
    target_curve, flat_line,
    k=0.3, lock_x=False, length_method='relu'
)
# --- Compute alpha + segment lengths ---
L, alpha = get_lengths_and_alpha(deformed, L0)

# x positions for segments (midpoints)
deformed_np = np.array(deformed)
segment_x = 0.5 * (deformed_np[:-1, 0] + deformed_np[1:, 0])

# tolerance bands (normalized around 1.0)
# in real structures, links can’t stretch/shrink arbitrarily; they only tolerate a certain percentage
tol = 0.15
alpha_min = 1 - tol
alpha_max = 1 + tol

mse = get_LSM_quality(deformed, target_curve)

# --- Plot ---
fig, ax1 = plt.subplots(figsize=(8,5))

# Left axis: geometry
ax1.plot(*zip(*target_curve), label='NACA 2412 Target Curve', linestyle='--', linewidth=2)
ax1.plot(*zip(*flat_line),    label='Flat Line', linestyle=':')
ax1.plot(*zip(*deformed),     label='Deformed Line', marker='o')
ax1.set_xlabel("Arbitrary Length")
ax1.set_ylabel("Arbitrary Height")

# Right axis (first): alpha (stretch factor)
ax2 = ax1.twinx()
x_vals = np.linspace(0, 1, len(alpha))
ax2.plot(x_vals, alpha, 'r.-', label='Alpha (L/L0)')
ax2.axhline(alpha_min, color='k', linestyle='--', linewidth=1, alpha=0.6, label='Tolerance Band')
ax2.axhline(alpha_max, color='k', linestyle='--', linewidth=1, alpha=0.6)
ax2.set_ylabel("Alpha (stretch factor)", color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Far-right axis: absolute segment lengths
ax3 = ax1.twinx()
ax3.spines['right'].set_position(("outward", 60))  # offset the axis
L, _ = get_lengths_and_alpha(deformed, L0)
ax3.plot(x_vals, L, 'g.-', label='Segment Length')
ax3.set_ylabel("Segment Length (absolute)", color='g')
ax3.tick_params(axis='y', labelcolor='g')

# Legends (merged)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines3, labels3 = ax3.get_legend_handles_labels()
ax1.legend(lines1 + lines2 + lines3, labels1 + labels2 + labels3,
           loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

# Title
plt.title(f"RADs Iterative Alpha Solver | MSE={mse:.3e}")
plt.tight_layout()
plt.show()