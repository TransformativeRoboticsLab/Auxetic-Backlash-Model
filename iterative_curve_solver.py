import numpy as np

def relax_to_curve(target_curve, flat_line, max_iters=1000, step_size=0.01, tol=1e-6):
    """
    Iteratively adjust internal points of a flat line to approximate a target curve.

    Parameters:
        target_curve (list of (x, y)): Target curve as linked (x, y) tuples.
        flat_line (list of (x, y)): Initial flat line as linked (x, y) tuples.
        max_iters (int): Maximum number of iterations.
        step_size (float): Gradient descent step size.
        tol (float): Convergence tolerance (average movement).

    Returns:
        list of (x, y): Deformed line approximating the target curve.
    """
    target = np.array(target_curve)
    current = np.array(flat_line)

    # assert len(target) == len(current), "Target and flat line must have the same number of points."
    assert np.allclose(current[0], [0, 0]) and np.allclose(current[-1], [1, 0]), "First and last points must be fixed."

    n = len(current)

    for _ in range(max_iters):
        prev = current.copy()

        for i in range(1, n - 1):  # Skip endpoints
            # Move towards the target curve point
            direction = target[i] - current[i]
            current[i] += step_size * direction

        # Optional: enforce constant segment lengths
        # This can be used to mimic physical link lengths more realistically
        for _ in range(1):  # One pass per iteration to preserve link lengths
            for i in range(1, n - 2):
    

                two_left = current[i-1] - current[i - 2]
                left = current[i] - current[i - 1]
                two_right = current[i + 2] - current[i+1]
                right = current[i + 1] - current[i]
                "Lengths of neighboring line segments"
                left_len = np.linalg.norm(left)
                two_left_len = np.linalg.norm(two_left)
                right_len = np.linalg.norm(right)
                two_right_len = np.linalg.norm(two_right)

                # bounds on segment length differences
                upper_bound = ((left_len+right_len)/2)+0.001
                lower_bound = ((left_len+right_len)/2)-0.001

                desired_left_len = np.linalg.norm(target[i] - target[i - 1])
                desired_right_len = np.linalg.norm(target[i + 1] - target[i])

                if left_len > 1e-6:
                    current[i] -= 0.5 * (left_len - desired_left_len) * (left / left_len)
                if right_len > 1e-6:
                    current[i] += 0.5 * (right_len - desired_right_len) * (right / right_len)
                # print(two_left_len)
                # if two_left_len > upper_bound:
                #     current[i] -= 0.01
                # if two_left_len < lower_bound:
                #     current[i] += 0.01
                # if two_right_len > upper_bound:
                #     current[i] -= 0.01
                # if two_right_len < lower_bound:
                #     current[i] += 0.01
                
        # Convergence check
        avg_move = np.mean(np.linalg.norm(current - prev, axis=1))
        if avg_move < tol:
            break

    return [tuple(p) for p in current]


def get_alpha_from_deformed(deformed_list):
    """
    Compute the distance between each pair of adjacent points in a list.
    Parameters:
        points (list of (x, y)): A list of 2D points.
    Returns:
        list of float: Distances between each pair of adjacent points and resulting alpha function values
    """
    points = np.array(deformed_list)
    diffs = points[1:] - points[:-1]
    dists = np.linalg.norm(diffs, axis=1)
    N_cells = len(deformed_list)
    alpha_list = [i*(N_cells) for i in dists.tolist()]
    return alpha_list


# Use a semicircle as target curve
import matplotlib.pyplot as plt

num_points = 15
num_points_target = 15
x_vals = np.linspace(0, 1, num_points)
x_vals_target = np.linspace(0, 1, num_points_target)
# y_vals = np.sqrt(0.25 - (x_vals - 0.5)**2)  # Semicircle radius 0.5
y_vals_target = 0.12/0.2 * (0.2969*np.sqrt(x_vals_target) - 0.1260*x_vals_target - 0.3516*x_vals_target**2 + 0.2843*x_vals_target**3 - 0.1015*x_vals_target**4) 
# y_vals_target = np.sqrt(0.25 - (x_vals_target - 0.5)**2 ) /10  # Semicircle radius 0.5
target_curve = list(zip(x_vals_target, y_vals_target*0.1))
print(target_curve)

# Flat line from (0,0) to (1,0)
flat_line = list(zip(x_vals, np.zeros_like(x_vals)))

# Run the relaxation
deformed = relax_to_curve(target_curve, flat_line)
alpha = get_alpha_from_deformed(deformed)
print(alpha)

# Plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis
ax1.plot(*zip(*target_curve), label='NACA 2412 Target Curve', linestyle='--', linewidth=3)
ax1.plot(*zip(*flat_line), label='Flat Line', linestyle=':')
ax1.plot(*zip(*deformed), label='Deformed Line', marker='o')
ax2.plot(x_vals[1:], alpha, label='Alpha Line', marker='.', c='red')
plt.plot()
plt.legend()
ax1.set_xlabel("Arbitrary Length")
ax1.set_ylabel("Arbitrary Height")
ax2.set_ylabel("Alpha")
plt.title("RADs Iterative Alpha Solver")
plt.axis('equal')
plt.show()