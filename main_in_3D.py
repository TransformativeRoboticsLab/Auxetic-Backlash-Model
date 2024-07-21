# TRL 2023
# Jacob Miske

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from scipy.special import comb
from matplotlib.collections import PolyCollection, LineCollection
from scipy.spatial import Delaunay
from sympy.utilities.iterables import multiset_permutations
import itertools
import math

font = {'family': 'serif',
        'size': 12}
csfont = {'fontname': 'Helvetica'}

matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1


# TODO: general class based object for auxetic linkages

class AuxeticCell:
    # Class to describe an individual auxetic bilayer cell
    def __init__(self, x, y):
        self.arm_length = 2  # [cm]
        self.arms = 4
        self.layers = 2
        self.center_x = x
        self.center_y = y
        self.top_layer_angle = 0  # [radians]
        self.bottom_layer_angle = 45  # [radians]
        self.arm_angle = 2 * np.pi / (float(self.arms))  # [radians]
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
        # see main_in_2D.py file


class AbstractAuxeticCell3D:
    # point representation of auxetic cell, simpler than full class model
    def __init__(self, x, y, z):
        self.cell_index = 0  # default all cells to be index zero until otherwise stated
        self.pos_x = x
        self.pos_y = y
        self.pos_z = z
        self.arm_length = 10  # default radius [mm]
        self.joint_backlash = 1  # default backlash in revolute joint [mm]
        self.alpha = 1  # default cells at contracted state
        self.fixed = False
        # Keep a list of all other cell indexes connected to this cell
        self.connections = []

    def get_cell_position(self):
        """
        :return: prints (x, y) of cell
        """
        print(str(self.pos_x) + str(self.pos_y))

    def get_distance_between_cells(self, a_cell):
        """
        Determine cartesian distance between two cells
        :param a_cell: another cell, not self
        :return: cartesian distance
        """
        distance = np.sqrt((a_cell.pos_x - self.pos_x) ** 2 + (a_cell.pos_y - self.pos_y) ** 2)
        print(distance)
        return distance

    def get_connected_distance_between_cells(self, a_cell):
        """
        based on arm length and alpha, calculate connected distance between self and a_cell
        :param a_cell: another cell, not self
        :return:
        """
        connected_distance = self.arm_length * self.alpha + a_cell.arm_length * a_cell.alpha
        print(connected_distance)
        return connected_distance

    def plot_connection_shape_space(self):
        """
        For this given cell, plot the connection shape space as convex hulls
        :return:
        """
        # first, generate an equidistant point cloud in x,y,z within twice of the arm length
        x_values = np.linspace(0, 2 * self.arm_length, 50)
        y_values = np.linspace(0, 2 * self.arm_length, 50)
        z_values = np.linspace(0, 2 * self.arm_length, 50)
        point_values = np.array([x_values, y_values, z_values])
        cube_points = itertools.permutations(point_values)
        for p in cube_points:
            print(p)
        print("\n \n")
        # joint space is a segmented toroid
        outer_circle_points = [
            (np.cos(2 * np.pi / 100 * i) * self.arm_length, np.sin(2 * np.pi / 100 * i) * self.arm_length, 0) for i in
            range(0, 100)]
        print(outer_circle_points)

    def connect_cell(self, a_cell):
        """
        Mates this cell with another cell
        :param a_cell:
        :return:
        """
        # If both are fixed, check if center within arm length, if not return False
        if self.fixed and a_cell.fixed:
            print('both fixed, cannot move')
        # If one is free, move center position of free cell to nearest point in radius of fixed cell
        elif self.fixed or a_cell.fixed:
            # First, try to connect
            print('one is fixed')
            if self.fixed:
                print('self is fixed')
                # Try to move the other cell, get nearest point on circle between two points
                # Where P is the closest point, C is the center, and R is the radius
                V = [a_cell.pos_x - self.pos_x, a_cell.pos_y - self.pos_y]
                print(V)
                # = C + V / |V| * R;
                mag_V = self.get_distance_between_cells(a_cell)
                print(mag_V)
                a_cell.pos_x = self.pos_x + V[0] * 1 / self.arm_length
                a_cell.pos_y = self.pos_y + V[1] * 1 / self.arm_length
            else:
                print('other cell is fixed')
                # Try to move self
                V = [a_cell.pos_x - self.pos_x, a_cell.pos_y - self.pos_y]
                # = C + V / |V| * R;
                mag_V = np.absolute(V)
                self.pos_x = (a_cell.pos_x + V[0]) * 1 / (mag_V * a_cell.arm_length)
                self.pos_y = (a_cell.pos_y + V[1]) * 1 / (mag_V * a_cell.arm_length)
            # if connection is successful, add each cell to the other's connections list

        # If both are free, find nearest point between two circles and move cells so that radius is
        else:
            print('neither fixed')
            # move both equally to meet one another


def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """
    return comb(n, i) * (t ** (n - i)) * (1 - t) ** i


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


def get_manhattan_distance(i, j, m, n):
    """
    For a point at i, j and another point at m, n
    returns manhattan distance
    :param i:
    :param j:
    :param m:
    :param n:
    :return:
    """
    return abs(i - m) + abs(j - n)


def simulate_auxetic_lattice_dome_points(b, L, t, alpha, m, n):
    """
    Given cell backlash and length in mm, alpha is dilation factor
    For lattice of m rows by n columns
    get list of [(x, y, z);...] of each cell assuming suspended from center cell(s)
    :param b:
    :param L:
    :param t: cell thickness
    :param alpha: dilation ratio across auxetic dome
    :param m:
    :param n:
    :return:
    """
    norm_b = b / L  # [mm]
    # convert L and alpha to L_center_to_center distance
    L_center_to_center = L + 10 * alpha - 10
    # Auxetic cell phi angle [radians]
    phi = np.arctan(b / t)
    h = np.sin(phi) * L_center_to_center
    cell_spacing = np.cos(phi) * L_center_to_center
    # corner cell is origin
    center_cell_m = round(m / 2)
    center_cell_n = round(n / 2)
    xyz_cells = []
    # for each cell in lattice, determine distance from center cell, move down from z=0 plane
    for i in range(m):
        for j in range(n):
            xyz_cell = []
            # place x and y values
            xyz_cell.insert(0, i * cell_spacing)
            xyz_cell.insert(1, j * cell_spacing)
            manhattan_dist = get_manhattan_distance(i, j, center_cell_m, center_cell_n)
            depth = float(manhattan_dist * -h)
            # second cell effects square the intercell depth
            xyz_cell.insert(2, -1 / 8 * depth ** (2))
            xyz_cells.append(xyz_cell)
    return xyz_cells


def get_sphere_points(xc, yc, zc, r):
    """
    Given a sphere centered at (x, y, z) with radius r, place points around sphere
    using a Fibonacci lattice spiral
    :param xc: sphere center x
    :param yc: sphere center y
    :param zc: sphere center z
    :param r: sphere radius
    :return:
    """
    N = 100
    points = []
    offset = 2.0 / N
    increment = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians
    for i in range(N):
        y = ((i * offset) - 1) + (offset / 2)
        r_temp = np.sqrt(1 - y * y)
        phi = ((i % N) * increment)
        x = np.cos(phi) * r_temp
        z = np.sin(phi) * r_temp
        # Transform to the sphere's radius and center
        x = xc + r * x
        y = yc + r * y
        z = zc + r * z
        points.append((x, y, z))
    return points


def get_distance_from_sphere(xp, yp, zp, xc, yc, zc, r):
    """
    Given a sphere centered at (xc, yc, zc) with radius r
    Determine nearest distance from (xp, yp, zp) to the sphere shell
    :param xp: point x
    :param yp: point y
    :param zp: point z
    :param xc: sphere center x
    :param yc: sphere center y
    :param zc: sphere center z
    :param r: sphere radius
    :return: nearest distance between point and center of sphere
    """
    # Calculate Euclidean distance between the point and the center of the sphere
    d = math.sqrt((xp - xc) ** 2 + (yp - yc) ** 2 + (zp - zc) ** 2)
    # Calculate the nearest distance to the surface of the sphere
    nearest_distance = abs(d - r)
    return nearest_distance


def get_experiment_dome_data(filename):
    """
    Pulls cell locations from filename and saves in [[x,y,z],[x,y,z],...] format
    :param filename:
    :return:
    """
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)
        for i in range(len(data)):
            data[i] = [float(j) for j in data[i]]
    return data


def get_error_2sets(points1, points2):
    """
    Given two sets of [[x, y, z],...] points SORTED
    get error between sets
    :param points1:
    :param points2:
    :return:
    """
    error_points = np.array(points2)
    distances = np.zeros((1, len(points1)), dtype=np.float32)
    for i in range(len(points1)):
        point1 = points1[i]
        point2 = points2[i]
        x1 = point1[0]
        y1 = point1[1]
        z1 = point1[2]
        x2 = point2[0]
        y2 = point2[1]
        z2 = point2[2]
        distance = math.sqrt((x1 - x2) ** 2 + (y1 - y1) ** 2 + (z1 - z1) ** 2)
        distances[0, i] = distance
    error_points[:, 2] = distances
    return error_points


def vary_points(cell_list):
    """
    See if random variation to some degree results in more or less error to model
    :param cell_list:
    :return:
    """
    new_cell_list = []
    for cell in cell_list:
        new_cell = []
        for i in cell:
            new_cell.append((i + i * (np.random.rand() / 30 + 1))/2)
        new_cell_list.append(new_cell)
    return new_cell_list


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

    polynomial_array = np.array([bernstein_poly(i, nPoints - 1, t) for i in range(0, nPoints)])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)

    return xvals, yvals


def naca0018(x):
    # Thickness distribution for NACA 0018
    t = 0.18
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)
    return yt


def compute_endpoint(x, y, angle, length):
    """
    Compute the endpoint of a line segment.
    Starting point (x, y), given angle and length.
    """
    end_x = x + length * np.cos(angle)
    end_y = y + length * np.sin(angle)
    return end_x, end_y


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p) >= 0


def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.
    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return 0


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


def plot_3D_points_error(points, filename):
    """
    Given set of 3D points (error between sets, 3D scatter plot
    :param points:
    :return:
    """
    fig = plt.figure()
    fig.set_figheight(9)
    fig.set_figwidth(9)
    ax = fig.add_subplot(111, projection='3d')
    # Unzip the list of points into x, y, z coordinates
    x, y, z = zip(*points)
    ax.scatter(x, y, z, c='g', marker='o', label="Total Error")
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.set_xlabel('X Exp. [mm]', labelpad=20, **csfont)
    ax.set_ylabel('Y Exp. [mm]', labelpad=20, **csfont)
    ax.set_zlabel('Error [mm]', labelpad=20, **csfont)
    set_axes_equal(ax)
    ax.set_zlim3d(-10, 10)
    plt.legend()

    x_scale = 2.5
    y_scale = 2.5
    z_scale = 1
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj = short_proj

    plt.savefig("./figures/{}_error_3D_points.png".format(filename))
    return 0


def plot_3D_points_error_percentage(points, filename):
    """
    Given set of 3D points (error between sets, 3D scatter plot
    :param points:
    :return:
    """
    fig = plt.figure()
    fig.set_figheight(9)
    fig.set_figwidth(9)
    ax = fig.add_subplot(111, projection='3d')
    # Unzip the list of points into x, y, z coordinates
    x, y, z = zip(*points)
    z = [abs(100*(i/35)) for i in z] # cells are L=35
    ax.scatter(x, y, z, c='g', marker='o', label="Total Error")
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.set_xlabel('X Exp. [mm]', labelpad=20, **csfont)
    ax.set_ylabel('Y Exp. [mm]', labelpad=20, **csfont)
    ax.set_zlabel('Error [& of L]', labelpad=20, **csfont)
    set_axes_equal(ax)
    ax.set_zlim3d(0, 20)
    plt.legend()

    x_scale = 2.5
    y_scale = 2.5
    z_scale = 1
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj = short_proj

    plt.savefig("./figures/{}_error_3D_points.png".format(filename), dpi=600)
    return 0


def get_alpha_dome(points, alpha, filename):
    """
    :param points:
    :param filename:
    :return:
    """
    fig = plt.figure()
    fig.set_figheight(9)
    fig.set_figwidth(9)
    ax = fig.add_subplot(111, projection='3d')
    # Unzip the list of points into x, y, z coordinates
    x, y, z = zip(*points)
    z = [alpha for i in z]  # cells are L=35
    ax.scatter(x, y, z, c='k', marker='o', label="Alpha Function")
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.set_xlabel('X Exp. [mm]', labelpad=20, **csfont)
    ax.set_ylabel('Y Exp. [mm]', labelpad=20, **csfont)
    ax.set_zlabel('Alpha', labelpad=20, **csfont)
    set_axes_equal(ax)
    ax.set_zlim3d(0, 3)
    plt.legend()

    x_scale = 2.5
    y_scale = 2.5
    z_scale = 1
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj = short_proj

    plt.savefig("./figures/{}_alpha_3D_points.png".format(filename), dpi=600)
    return 0


def plot_3D_points_2sets(points1, points2, filename):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.set_figheight(9)
    fig.set_figwidth(9)
    # Unzip the list of points into x, y, z coordinates
    x1, y1, z1 = zip(*points1)
    x2, y2, z2 = zip(*points2)
    ax.scatter(x1, y1, z1, c='b', marker='o', label="Model")
    ax.scatter(x2, y2, z2, c='r', marker='o', label="Experiment")
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.set_xlabel('X [mm]', labelpad=20, **csfont)
    ax.set_ylabel('Y [mm]', labelpad=20, **csfont)
    ax.set_zlabel('Z [mm]', labelpad=20, **csfont)
    set_axes_equal(ax)
    plt.legend()

    x_scale = 2.5
    y_scale = 2.5
    z_scale = 1
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0
    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)
    ax.get_proj = short_proj

    ax.set_zlim3d(-100, 100)
    plt.savefig("./figures/{}.png".format(filename), dpi=600)
    return 0


def plot_auxetic_domes_figure():
    """
    generate two domes in 3D figure to compare to experimental results
    :return:
    """


if __name__ == '__main__':
    # plot 3D shape space of auxetic cell
    cell1 = AbstractAuxeticCell3D(x=0, y=0, z=0)
    cell1.plot_connection_shape_space()

    # auxetic 3D dome at variable alpha
    dome_cells1 = simulate_auxetic_lattice_dome_points(b=0.4, L=35, t=8, alpha=1, m=8, n=11)
    dome_cells1_experiment = vary_points(dome_cells1) # get_experiment_dome_data(filename="./data/dome1_exp_20240715.csv")
    # vary_points(dome_cells1)
    # plot_3D_points_2sets(points1=dome_cells1, points2=dome_cells1_experiment, filename="dome1")

    dome_cells2 = simulate_auxetic_lattice_dome_points(b=0.4, L=35, t=8, alpha=2, m=8, n=11)
    dome_cells2_experiment = vary_points(dome_cells2) # get_experiment_dome_data(filename="./data/dome2_exp_20240715.csv")
    #
    # plot_3D_points_2sets(points1=dome_cells2, points2=dome_cells2_experiment, filename="dome2")

    dome1_error = get_error_2sets(dome_cells1, dome_cells1_experiment)
    plot_3D_points_error(dome1_error, filename="dome1")
    plot_3D_points_error_percentage(dome1_error, filename="dome1_percentage")

    dome2_error = get_error_2sets(dome_cells2, dome_cells2_experiment)
    plot_3D_points_error(dome2_error, filename="dome2_percentage")
    plot_3D_points_error_percentage(dome2_error, filename="dome2_percentage")

    dome1_alpha = get_alpha_dome(points=dome_cells1_experiment, alpha=1, filename="dome1_alpha")
    dome2_alpha = get_alpha_dome(points=dome_cells2_experiment, alpha=2, filename="dome2_alpha")


    quit()

    # airfoil in 3D example
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
    plt.savefig("./figures/NACA0018example.png")
    # Linkages plot
    # Lengths of the linkages
    link_lengths = [5, 3]

    # Angles (in radians) between the linkages
    angle_1 = np.radians(45)  # 45 degrees
    angle_2 = np.radians(-30)  # -30 degrees

    # plot_linkage_system(angle_1, angle_2, link_lengths)

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
    # Get the Bézier parameters based on a degree.
    data = get_bezier_parameters(xpoints, ypoints, degree=4)
    x_val = [x[0] for x in data]
    y_val = [x[1] for x in data]
    print(data)
    # Plot the control points
    plt.plot(x_val, y_val, 'k--o', label='Control Points')
    # Plot the resulting Bézier curve
    xvals, yvals = bezier_curve(data, nTimes=1000)
    plt.plot(xvals, yvals, 'b-', label='B Curve')
    plt.xlabel("X Dimension (m)")
    plt.legend()
    plt.savefig("./figures/bezier_example.png")
