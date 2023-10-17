# TRL 2023
# Jacob Miske

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb


class auxetic_cell:
    # Class to describe an individual auxetic bilayer cell
    def __init__(self, x, y):
        self.arm_length = 2                                       # [cm]
        self.arms = 4
        self.layers = 2
        self.center_x = x
        self.center_y = y
        self.top_layer_angle = 0                                  # [radians]
        self.bottom_layer_angle = 45                              # [radians]
        self.arm_angle = 2*(np.pi)/(float(self.arms)) # [radians]
        self.arm_positions = []
    def get_position_list(self):
        cell_positions = []
        cell_center = [self.center_x, self.center_y]
        cell_positions.append(cell_center)
        # assume equally spaced arms in angle
        # First four cell positions are for top layer
        for i in range(4):
            cell_i_position = [self.arm_length * np.sin(i * self.arm_angle + self.top_layer_angle),
                               self.arm_length * np.cos(i * self.arm_angle + self.top_layer_angle)]
            cell_positions.append(cell_i_position)
        # Second four cell positions are for bottom layer
        for j in range(4):
            cell_j_position = [self.arm_length * np.sin(j * self.arm_angle + self.bottom_layer_angle),
                               self.arm_length * np.cos(j * self.arm_angle + self.bottom_layer_angle)]
            cell_positions.append(cell_j_position)
        # return list of size = 1 (center) + arm ends (layers*arms)
        return cell_positions
    def update_arm_positions(self):
        # rewrites arm positions based on connections
        print('not ready')
    def connect(self, other_cell_):
        # connects a cell arm to another cell arm
        print('not ready')


def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """
    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i


def get_bezier_parameters(X, Y, degree=3):
    """ Least square qbezier fit using penrose pseudoinverse.

    Parameters:

    X: array of x data.
    Y: array of y data. Y[0] is the y point for X[0].
    degree: degree of the Bézier curve. 2 for quadratic, 3 for cubic.

    Based on https://stackoverflow.com/questions/12643079/b%C3%A9zier-curve-fitting-with-scipy
    and probably on the 1998 thesis by Tim Andrew Pastva, "Bézier Curve Fitting".
    """
    if degree < 1:
        raise ValueError('degree must be 1 or greater.')

    if len(X) != len(Y):
        raise ValueError('X and Y must be of the same length.')

    if len(X) < degree + 1:
        raise ValueError(f'There must be at least {degree + 1} points to '
                         f'determine the parameters of a degree {degree} curve. '
                         f'Got only {len(X)} points.')

    def bpoly(n, t, k):
        """ Bernstein polynomial when a = 0 and b = 1. """
        return t ** k * (1 - t) ** (n - k) * comb(n, k)
        # return comb(n, i) * ( t**(n-i) ) * (1 - t)**i

    def bmatrix(T):
        """ Bernstein matrix for Bézier curves. """
        return np.matrix([[bpoly(degree, t, k) for k in range(degree + 1)] for t in T])

    def least_square_fit(points, M):
        M_ = np.linalg.pinv(M)
        return M_ * points

    T = np.linspace(0, 1, len(X))
    M = bmatrix(T)
    points = np.array(list(zip(X, Y)))

    final = least_square_fit(points, M).tolist()
    final[0] = [X[0], Y[0]]
    final[len(final) - 1] = [X[len(X) - 1], Y[len(Y) - 1]]
    return final


def bezier_curve(points, nTimes=50):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.

       points should be a list of lists, or list of tuples
       such as [ [1,1],
                 [2,3],
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000

        See http://processingjs.nihongoresources.com/bezierinfo/
    """

    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)

    return xvals, yvals


def naca0018(x):
    # Thickness distribution for NACA 0018
    t = 0.18
    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    return yt


def compute_endpoint(x, y, angle, length):
    """
    Compute the endpoint of a line segment.
    Starting point (x, y), given angle and length.
    """
    end_x = x + length * np.cos(angle)
    end_y = y + length * np.sin(angle)
    return end_x, end_y


def plot_linkage_system(angle1, angle2, lengths):
    """
    Plot a linkage system with given angles and lengths.
    """
    x = [0, 0]  # Initial x coordinates
    y = [0, 0]  # Initial y coordinates

    # Compute coordinates for the second linkage
    x[1], y[1] = compute_endpoint(x[0], y[0], angle1, lengths[0])

    # Compute coordinates for the third linkage
    x2, y2 = compute_endpoint(x[1], y[1], angle1 + angle2, lengths[1])

    # Plot the linkages
    plt.figure(figsize=(10, 10))
    plt.plot([x[0], x[1]], [y[0], y[1]], '-o', label="Link 1")
    plt.plot([x[1], x2], [y[1], y2], '-o', label="Link 2")
    plt.legend()
    plt.grid(True)
    plt.axis("equal")
    plt.title("Linkage System")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()





if __name__ == '__main__':
    print("This model assumes links can flex")

    # Show NACA 0018 shape
    # Create x values from 0 to 1, 200 points
    x = np.linspace(0, 1, 200)
    # Calculate upper and lower surface points
    upper_y = naca0018(x)
    lower_y = -naca0018(x)
    # Plot the airfoil shape
    plt.plot(x, upper_y, 'b-')
    plt.plot(x, lower_y, 'b-')
    plt.title('NACA 0018 Airfoil Shape')
    plt.xlabel('Chord')
    plt.ylabel('Thickness')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

    # Linkages plot
    # Lengths of the linkages
    link_lengths = [5, 3]

    # Angles (in radians) between the linkages
    angle_1 = np.radians(45)  # 45 degrees
    angle_2 = np.radians(-30)  # -30 degrees

    plot_linkage_system(angle_1, angle_2, link_lengths)

    # bezier example
    points = []
    xpoints = [2.70, 2.69, 2.68, 2.66, 2.64, 2.63, 2.61, 2.61, 2.64, 2.68,
               2.74, 2.82, 2.90, 2.99, 3.07, 3.16, 3.24, 3.33, 3.42]
    ypoints = [8.95, 8.85, 8.75, 8.65, 8.55, 8.47, 8.40, 8.32,
               8.27, 8.23, 8.18, 8.18, 8.18, 8.18, 8.19, 8.19,
               8.19, 8.20, 8.20]

    for i in range(len(xpoints)):
        points.append([xpoints[i], ypoints[i]])

    # Plot the original points
    plt.plot(xpoints, ypoints, "ro", label='Original Points')
    # Get the Bezier parameters based on a degree.
    data = get_bezier_parameters(xpoints, ypoints, degree=4)
    x_val = [x[0] for x in data]
    y_val = [x[1] for x in data]
    print(data)
    # Plot the control points
    plt.plot(x_val, y_val, 'k--o', label='Control Points')
    # Plot the resulting Bezier curve
    xvals, yvals = bezier_curve(data, nTimes=1000)
    plt.plot(xvals, yvals, 'b-', label='B Curve')
    plt.xlabel("X Dimension (m)")

    plt.legend()
    plt.show()