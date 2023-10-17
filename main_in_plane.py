# TRL
# Jacob Miske

from math import *
from itertools import *
import numpy as np
import matplotlib.pyplot as plt


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


class abstract_auxetic_cell:
    # point representation of auxetic cell
    def __init__(self, x, y):
        self.position_x = x
        self.position_y = y



def get_relu(a, b):
    # returns ReLU function response
    return lambda x: max(0, a*(x-b))


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
        angle = (angle - a*(angle - 0.785))
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
    x = np.arange(0, cols * distance, distance)
    y = np.arange(0, rows * distance, distance)
    X, Y = np.meshgrid(x, y)
    # Plot the lattice points
    plt.scatter(X, Y, s=dot_sizes, c='blue')
    plt.xlim([-1, 5])
    plt.ylim([-1, 5])
    plt.xlabel("Distance in X [cm]")
    plt.ylabel("Distance in Y [cm]")
    plt.title("2D Lattice Cell Variance from Backlash")
    # Setting equal scaling and showing the plot
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()





if __name__ == '__main__':
    # Parameters for the scissor mechanism
    LENGTH = 2.0  # Length of each linkage
    ANGLE = np.pi / 4  # Angle between the links
    NUMBER_OF_SERIES = 5  # Number of scissor mechanisms in series
    xs, ys = scissor_position(LENGTH, ANGLE, number_of_links=20)

    plot_scissor_mechanism(xs, ys, NUMBER_OF_SERIES)

    dot_sizes = np.array([[1, 5, 12, 5, 1],[5, 12, 17, 12, 5],[12, 17, 30, 17, 12],[5, 12, 17, 12, 5], [1, 5, 12, 5, 1]])
    draw_lattice_with_variable_dot_size(5, 5, 1, dot_sizes=dot_sizes)  # 5x5 lattice with points separated by 1 unit distance

    dot_sizes = np.array([[30, 17, 12, 17, 30], [17, 12, 5, 12, 17], [12, 5, 1, 5, 12], [5, 12, 5, 12, 5], [1, 5, 12, 5, 1]])
    draw_lattice_with_variable_dot_size(5, 5, 1, dot_sizes=dot_sizes)  # 5x5 lattice with points separated by 1 unit distance

    # cell1 = auxetic_cell(x=0, y=0)
    # cell1_positions = cell1.get_position_list()
    #
    # plot_auxetic_cell(list_of_cell_positions=cell1_positions)
    #
    # # Set fixed variable of the structures shape
    # link_length = 2 # [cm]
    #
    # ReLU_example = get_relu(1, 1)
    # ReLU_example(5)
    #
    # a11 = lambda x, y: get_relu(1, 0.564)(x)
    # print(list(takewhile(lambda x: x > 0, accumulate(range(10, 2000), a11))))
    # target = [0.6]