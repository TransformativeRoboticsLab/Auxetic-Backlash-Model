# TRL
# Jacob Miske

from math import *
from itertools import *
import numpy as np
import matplotlib.pyplot as plt

# References
# even spaced points on a convex hull
# https://stackoverflow.com/questions/75416188/get-evenly-spaced-points-from-a-curved-shape

class AuxeticCell:
    # Class to describe an individual auxetic bilayer cell
    def __init__(self, x, y):
        self.cell_index = 0  # default all cells to be index zero until otherwise stated
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


class AbstractAuxeticCell:
    # point representation of auxetic cell, simpler than full class model
    def __init__(self, x, y):
        self.cell_index = 0  # default all cells to be index zero until otherwise stated
        self.pos_x = x
        self.pos_y = y
        self.arm_length = 10  # default radius [mm]
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


def plot_abstract_cells(list_cells):
    """
    Given a list of AbstractAuxeticCells, plot centers and connected links
    :param list_cells:
    :return:
    """
    plt.figure(0)
    for cell in list_cells:
        print(vars(cell))
        plt.scatter(cell.pos_x, cell.pos_y, label=cell.cell_index)
        plt.legend()
        plt.savefig("TRL_abstract_cells_plot.png")
    plt.show()
    plt.close()
    return 0

def get_relu(a, b):
    # returns ReLU function response
    return lambda x: max(0, a * (x - b))


def plot_auxetic_cell(list_of_cell_positions):
    # for a cell, plot it's position in space
    plt.figure(0)
    center_x = list_of_cell_positions[0][0]
    center_y = list_of_cell_positions[0][1]
    xs = [i[0] for i in list_of_cell_positions]
    ys = [j[1] for j in list_of_cell_positions]
    for point in list_of_cell_positions[1:]:
        plt.plot([center_x, point[0]], [center_y, point[1]], 'b')
    plt.scatter(xs, ys)
    plt.show()
    plt.close()


def scissor_position(length, angle, number_of_links):
    """
    Calculate the positions of the scissor mechanism with backlash
    """
    # ReLU strength
    a = 0.1

    x_positions = [0]
    y_positions = [0]

    for i in range(number_of_links):
        # Calculate position of the next hinge based on current hinge
        x_next = x_positions[-1] + length * np.sin(angle)
        y_next = y_positions[-1] + length * np.cos(angle)

        x_positions.append(x_next)
        y_positions.append(y_next)

        # Reverse the direction of the angle for the next linkage
        angle = -angle

        x_next = x_positions[-1] + length * np.sin(angle)
        y_next = y_positions[-1] + length * np.cos(angle)

        x_positions.append(x_next)
        y_positions.append(y_next)

        # Apply decrement to neutral pi/4 angle with ReLU relationship
        angle = (angle - a * (angle - 0.785))
        print(angle)
    print(x_positions)
    print(y_positions)
    return x_positions, y_positions


def plot_scissor_mechanism(xs, ys, number_of_series):
    plt.figure(figsize=(10, 6))
    for i in range(number_of_series):
        plt.plot(xs, ys, '-o')

    plt.title("Scissor Mechanisms in Series")
    plt.xlabel("Distance in X (cm)")
    plt.ylabel("Distance in Y (cm)")
    plt.xlim([-2, 2])
    plt.ylim([0, 80])

    plt.legend()
    plt.grid(True)
    plt.show()


def draw_lattice_with_variable_dot_size(rows, cols, distance, dot_sizes):
    # Create a meshgrid for the lattice points
    x = 3.5 * np.arange(0, cols * distance, distance)
    y = 3.5 * np.arange(0, rows * distance, distance)
    X, Y = np.meshgrid(x, y)
    # Experimental results from lab notebook
    x_experimental = [0, 0.02, 0.04, 0.05, 0.05,
                      1.01, 1.02, 1.02, 1.03, 1.06,
                      2.01, 2.03, 2.04, 2.05, 2.06,
                      3.03, 3.04, 3.05, 3.08, 3.08,
                      4.05, 4.06, 4.07, 4.09, 4.11]
    y_experimental = [0, 1.02, 2.03, 3.04, 4.06,
                      0.02, 1.05, 2.04, 3.04, 4.06,
                      0.03, 1.03, 2.04, 3.05, 4.09,
                      0.03, 1.05, 2.05, 3.06, 4.09,
                      0.04, 1.05, 2.06, 3.09, 4.10]
    # scaling and linear set
    x_experimental = [3.5 * i for i in x_experimental]
    y_experimental = [3.5 * j for j in y_experimental]
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    # Plot the lattice points
    plt.scatter(X, Y, s=dot_sizes, c='blue', label="Configuration Space", alpha=0.3)
    plt.scatter(x_experimental, y_experimental, s=15, c='red', facecolors='r', label="Experiment", alpha=0.5)
    plt.scatter([0], [0], facecolors='none', edgecolors='g', s=50, label="Locked Cell")
    plt.xlim([-1, 15])
    plt.ylim([-1, 15])
    plt.xlabel("Distance in X [cm]")
    plt.ylabel("Distance in Y [cm]")
    plt.title("Flat Lattice with Backlash Pulled from Corner")
    # Setting equal scaling and showing the plot
    leg = plt.legend(loc='lower right')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)
    leg.get_frame().set_linewidth(2)
    leg.get_frame().set_edgecolor('k')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("./figures/2D_lattice_experiment.png", dpi=400)
    plt.close()


def plot_mass_on_spring(m, k, w):
    """
    For a single mass m on spring k with force=1 at frequency w, plot distance from equilibrium over time
    :param m:
    :param k:
    :return:
    """
    t = np.linspace(0, 100, 100)
    y = []
    for i in t:
        y = np.cos(w * t)
    plt.figure(1)
    plt.scatter(t, y)
    plt.show()
    plt.close()


def plot_N_masses_with_stiffness_in_between(m, k, w):
    """
    Assuming all masses are m=1, spring constant k, freq of first mass at w
    Dead zone between masses of length 0.1, length = 1
    :param m:
    :param k:
    :param w:
    :return:
    """
    t = np.linspace(0, 200, 2000)
    # first mass starts at position 0
    y1 = [0]
    y2 = [1]
    y3 = [2]
    y4 = [3]
    y5 = [4]
    # initially, all masses are unstressed
    y1_accel = [0]
    y2_accel = [0]
    y3_accel = [0]
    y4_accel = [0]
    y5_accel = [0]
    # for each timestep, determine the accel on each mass, y1 thru y5
    for count, i in enumerate(t, 0):
        y1_accel.append(0.01 * cos(w * t[count]) + k * ((y2[count] - y1[count]) - 1))
        y2_accel.append(- k * (y2[count] - y1[count]) + k * (y3[count] - y2[count]))
        y3_accel.append(- k * (y3[count] - y2[count]) + k * (y4[count] - y3[count]))
        y4_accel.append(- k * (y4[count] - y3[count]) + k * (y5[count] - y4[count]))
        y5_accel.append(- k * (y5[count] - y4[count]))
        # Now adjust the position of individual masses based on the change in acceleration over the time step
        # last position plus movement 2(a)(dt)
        y1.append(y1[-1] + 2 * y1_accel[count] * 1)
        y2.append(y2[-1] + 2 * y2_accel[count] * 1)
        y3.append(y3[-1] + 2 * y3_accel[count] * 1)
        y4.append(y4[-1] + 2 * y4_accel[count] * 1)
        y5.append(y5[-1] + 2 * y3_accel[count] * 1)
    # TODO: this is a hard code line solution for t
    t = np.linspace(0, 200, 2001)
    plt.figure(1)
    plt.scatter(t, y1, label="Joint 1 (Driving)")
    plt.scatter(t, y2, label="Joint 2")
    plt.scatter(t, y3, label="Joint 3")
    plt.scatter(t, y4, label="Joint 4")
    plt.scatter(t, y5, label="Joint 5")
    plt.title("N Joints in Series with Linear Coupling")
    plt.ylabel("Distance along x-Direction (cm)")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.show()
    plt.close()


def plot_N_masses_with_stiffness_in_between_with_RELU_gap(m, k, w, l_gap):
    """
    Assuming all masses are m=1, spring constant k, freq of first mass at w
    Dead zone between masses of length l_gap
    :param m:
    :param k:
    :param w:
    :return:
    """
    t = np.linspace(0, 200, 2000)
    # first mass starts at position 0, equilibrium length of 1 length unit
    y1 = [0]
    y2 = [1]
    y3 = [2]
    y4 = [3]
    y5 = [4]
    # initially, all masses are unstressed
    y1_accel = [0]
    y2_accel = [0]
    y3_accel = [0]
    y4_accel = [0]
    y5_accel = [0]
    # for each timestep, determine the accel on each mass, y1 thru y5
    for count, i in enumerate(t, 0):
        # Find distance between each set of masses
        y1_to_y2_distance = abs(y2[count] - y1[count])
        y2_to_y3_distance = abs(y3[count] - y2[count])
        y3_to_y4_distance = abs(y4[count] - y3[count])
        y4_to_y5_distance = abs(y5[count] - y4[count])
        # Apply ReLU relationship
        if y1_to_y2_distance > 1 + l_gap or y1_to_y2_distance < 1 - l_gap:
            y1_to_y2_distance = y2[count] - y1[count]
        else:
            y1_to_y2_distance = 1
        if y2_to_y3_distance > 1 + l_gap or y2_to_y3_distance < 1 - l_gap:
            y2_to_y3_distance = y3[count] - y2[count]
        else:
            y2_to_y3_distance = 1
        if y3_to_y4_distance > 1 + l_gap or y3_to_y4_distance < 1 - l_gap:
            y3_to_y4_distance = y4[count] - y3[count]
        else:
            y3_to_y4_distance = 1
        if y4_to_y5_distance > 1 + l_gap or y4_to_y5_distance < 1 - l_gap:
            y4_to_y5_distance = y5[count] - y4[count]
        else:
            y4_to_y5_distance = 1
        print("y1toy2 dist")
        print(y1_to_y2_distance)
        # append accelerations
        y1_accel.append(0.01 * cos(w * t[count]) + k * (y1_to_y2_distance - 1))
        y2_accel.append(- k * (y1_to_y2_distance - 1) + k * (y2_to_y3_distance - 1))
        y3_accel.append(- k * (y2_to_y3_distance - 1) + k * (y3_to_y4_distance - 1))
        y4_accel.append(- k * (y3_to_y4_distance - 1) + k * (y4_to_y5_distance - 1))
        y5_accel.append(- k * (y4_to_y5_distance - 1))
        print("y1 accel")
        print(y1_accel[-1])
        # Now adjust the position of individual masses based on the change in acceleration over the time step
        # last position plus movement 2(a)(dt)
        y1.append(y1[-1] + 2 * y1_accel[count] * 1)
        y2.append(y2[-1] + 2 * y2_accel[count] * 1)
        y3.append(y3[-1] + 2 * y3_accel[count] * 1)
        y4.append(y4[-1] + 2 * y4_accel[count] * 1)
        y5.append(y5[-1] + 2 * y3_accel[count] * 1)
    # TODO: this is a hard code line solution for t
    t = np.linspace(0, 200, 2001)

    plt.figure(1)
    plt.scatter(t, y1, label="Joint 1 (Driving)")
    plt.scatter(t, y2, label="Joint 2")
    plt.scatter(t, y3, label="Joint 3")
    plt.scatter(t, y4, label="Joint 4")
    plt.scatter(t, y5, label="Joint 5")
    plt.title("N Joints in Series with ReLU Coupling m={} k={} w={} l_gap={}".format(m, k, w, l_gap))
    plt.ylabel("Distance along x-Direction (cm)")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.xlim([0, 200])
    plt.savefig("TRL_ReLU_coupling m={} k={} w={} l_gap={}.png".format(m, k, w, l_gap))
    plt.close()


def plot_N_masses_with_stiffness_in_between_with_RELU_gap_fixed_BC(m, k, w, l_gap):
    """
    Assuming all masses are m=1, spring constant k, freq of first mass at w
    Dead zone between masses of length l_gap
    :param m:
    :param k:
    :param w:
    :return:
    """
    t = np.linspace(0, 200, 2000)
    # first mass starts at position 0, equilibrium length of 1 length unit
    y1 = [0]
    y2 = [1]
    y3 = [2]
    y4 = [3]
    y5 = [4]
    # initially, all masses are unstressed
    y1_accel = [0]
    y2_accel = [0]
    y3_accel = [0]
    y4_accel = [0]
    y5_accel = [0]
    # for each timestep, determine the accel on each mass, y1 thru y5
    for count, i in enumerate(t, 0):
        # Find distance between each set of masses
        y1_to_y2_distance = abs(y2[count] - y1[count])
        y2_to_y3_distance = abs(y3[count] - y2[count])
        y3_to_y4_distance = abs(y4[count] - y3[count])
        y4_to_y5_distance = abs(y5[count] - y4[count])
        # Apply ReLU relationship
        if y1_to_y2_distance > 1 + l_gap or y1_to_y2_distance < 1 - l_gap:
            y1_to_y2_distance = y2[count] - y1[count]
        else:
            y1_to_y2_distance = 1
        if y2_to_y3_distance > 1 + l_gap or y2_to_y3_distance < 1 - l_gap:
            y2_to_y3_distance = y3[count] - y2[count]
        else:
            y2_to_y3_distance = 1
        if y3_to_y4_distance > 1 + l_gap or y3_to_y4_distance < 1 - l_gap:
            y3_to_y4_distance = y4[count] - y3[count]
        else:
            y3_to_y4_distance = 1
        if y4_to_y5_distance > 1 + l_gap or y4_to_y5_distance < 1 - l_gap:
            y4_to_y5_distance = y5[count] - y4[count]
        else:
            y4_to_y5_distance = 1
        print("y1toy2 dist")
        print(y1_to_y2_distance)
        # append accelerations
        y1_accel.append(0.01 * cos(w * t[count]) + k * (y1_to_y2_distance - 1))
        y2_accel.append(- k * (y1_to_y2_distance - 1) + k * (y2_to_y3_distance - 1))
        y3_accel.append(- k * (y2_to_y3_distance - 1) + k * (y3_to_y4_distance - 1))
        y4_accel.append(- k * (y3_to_y4_distance - 1) + k * (y4_to_y5_distance - 1))
        y5_accel.append(- k * (y4_to_y5_distance - 1))
        print("y1 accel")
        print(y1_accel[-1])
        # Now adjust the position of individual masses based on the change in acceleration over the time step
        # last position plus movement 2(a)(dt)
        y1.append(y1[-1] + 2 * y1_accel[count] * 1)
        y2.append(y2[-1] + 2 * y2_accel[count] * 1)
        y3.append(y3[-1] + 2 * y3_accel[count] * 1)
        y4.append(y4[-1] + 2 * y4_accel[count] * 1)
        # Boundary condition, fifth mass stuck to x=5
        y5.append(5)
    # TODO: this is a hard code line solution for t
    t = np.linspace(0, 200, 2001)

    plt.figure(1)
    plt.scatter(t, y1, label="Joint 1 (Driving)")
    plt.scatter(t, y2, label="Joint 2")
    plt.scatter(t, y3, label="Joint 3")
    plt.scatter(t, y4, label="Joint 4")
    plt.scatter(t, y5, label="Joint 5")
    plt.title("N Joints in Series with ReLU Modeled Coupling, Fixed BC")
    plt.ylabel("Distance along x-Direction (cm)")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.xlim([0, 200])
    plt.show()
    plt.close()


def plot_N_masses_with_stiffness_in_between_with_RELU_gap_2D(m, k, w, l_gap):
    """
    Assuming all masses are m=1, spring constant k, freq of first mass at w
    Dead zone between masses of length l_gap
    :param m:
    :param k:
    :param w:
    :return:
    """
    t = np.linspace(0, 200, 2000)
    # first mass starts at position 0, equilibrium length of 1 length unit
    y1 = [0]
    y2 = [1]
    y3 = [2]
    y4 = [3]
    y5 = [4]
    # initially, all masses are unstressed
    y1_accel = [0]
    y2_accel = [0]
    y3_accel = [0]
    y4_accel = [0]
    y5_accel = [0]
    # for each timestep, determine the accel on each mass, y1 thru y5
    for count, i in enumerate(t, 0):
        # Find distance between each set of masses
        y1_to_y2_distance = abs(y2[count] - y1[count])
        y2_to_y3_distance = abs(y3[count] - y2[count])
        y3_to_y4_distance = abs(y4[count] - y3[count])
        y4_to_y5_distance = abs(y5[count] - y4[count])
        # Apply ReLU relationship
        if y1_to_y2_distance > 1 + l_gap or y1_to_y2_distance < 1 - l_gap:
            y1_to_y2_distance = y2[count] - y1[count]
        else:
            y1_to_y2_distance = 1
        if y2_to_y3_distance > 1 + l_gap or y2_to_y3_distance < 1 - l_gap:
            y2_to_y3_distance = y3[count] - y2[count]
        else:
            y2_to_y3_distance = 1
        if y3_to_y4_distance > 1 + l_gap or y3_to_y4_distance < 1 - l_gap:
            y3_to_y4_distance = y4[count] - y3[count]
        else:
            y3_to_y4_distance = 1
        if y4_to_y5_distance > 1 + l_gap or y4_to_y5_distance < 1 - l_gap:
            y4_to_y5_distance = y5[count] - y4[count]
        else:
            y4_to_y5_distance = 1
        print("y1toy2 dist")
        print(y1_to_y2_distance)
        # append accelerations
        y1_accel.append(0.01 * cos(w * t[count]) + k * (y1_to_y2_distance - 1))
        y2_accel.append(- k * (y1_to_y2_distance - 1) + k * (y2_to_y3_distance - 1))
        y3_accel.append(- k * (y2_to_y3_distance - 1) + k * (y3_to_y4_distance - 1))
        y4_accel.append(- k * (y3_to_y4_distance - 1) + k * (y4_to_y5_distance - 1))
        y5_accel.append(- k * (y4_to_y5_distance - 1))
        print("y1 accel")
        print(y1_accel[-1])
        # Now adjust the position of individual masses based on the change in acceleration over the time step
        # last position plus movement 2(a)(dt)
        y1.append(y1[-1] + 2 * y1_accel[count] * 1)
        y2.append(y2[-1] + 2 * y2_accel[count] * 1)
        y3.append(y3[-1] + 2 * y3_accel[count] * 1)
        y4.append(y4[-1] + 2 * y4_accel[count] * 1)
        y5.append(y5[-1] + 2 * y3_accel[count] * 1)
    # TODO: this is a hard code line solution for t
    t = np.linspace(0, 200, 2001)

    plt.figure(1)
    plt.scatter(t, y1, label="Joint 1 (Driving)")
    plt.scatter(t, y2, label="Joint 2")
    plt.scatter(t, y3, label="Joint 3")
    plt.scatter(t, y4, label="Joint 4")
    plt.scatter(t, y5, label="Joint 5")
    plt.title("N Joints in Series with ReLU Coupling m={} k={} w={} l_gap={}".format(m, k, w, l_gap))
    plt.ylabel("Distance along x-Direction (cm)")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.xlim([0, 200])
    plt.savefig("TRL_ReLU_coupling m={} k={} w={} l_gap={}.png".format(m, k, w, l_gap))
    plt.close()


if __name__ == '__main__':
    # AbstractAuxeticCell tests
    cell1 = AbstractAuxeticCell(x=0, y=0)
    cell2 = AbstractAuxeticCell(x=50, y=80)
    cell1.fixed = True
    cell1.get_distance_between_cells(cell2)
    cell1.get_connected_distance_between_cells(cell2)
    cell1.connect_cell(cell2)

    plot_abstract_cells([cell1, cell2])

    # Parameters for the scissor mechanism
    # LENGTH = 2.0  # Length of each linkage
    # ANGLE = np.pi / 4  # Angle between the links
    # NUMBER_OF_SERIES = 5  # Number of scissor mechanisms in series
    # xs, ys = scissor_position(LENGTH, ANGLE, number_of_links=20)

    # plot_scissor_mechanism(xs, ys, NUMBER_OF_SERIES)
    # 2 mm gap

    # Figure 2.5 of RADs
    # dot_sizes = np.array(
    #     [[1, 3, 5, 7, 9], [3, 5, 7, 9, 11], [5, 7, 9, 11, 13], [7, 9, 11, 13, 15], [9, 11, 13, 15, 17]])
    # # scale up for visualization
    # dot_sizes = np.multiply(dot_sizes, 35)
    # draw_lattice_with_variable_dot_size(5, 5, 1,
    #                                     dot_sizes=dot_sizes)  # 5x5 lattice with points separated by 1 unit distance

    # dot_sizes = np.array([[30, 17, 12, 17, 30], [17, 12, 5, 12, 17], [12, 5, 1, 5, 12], [5, 12, 5, 12, 5], [1, 5, 12, 5, 1]])
    # draw_lattice_with_variable_dot_size(5, 5, 1, dot_sizes=dot_sizes)  # 5x5 lattice with points separated by 1 unit distance
    # plot_mass_on_spring(m=1, k=1, w=0.1)
    # plot_N_masses_with_stiffness_in_between(m=13, k=0.1, w=0.2)
    # plot_N_masses_with_stiffness_in_between_with_RELU_gap(m=13, k=0.1, w=0.2, l_gap=0.35)
    # plot_N_masses_with_stiffness_in_between_with_RELU_gap_fixed_BC(m=13, k=0.1, w=0.2, l_gap=0.35)
