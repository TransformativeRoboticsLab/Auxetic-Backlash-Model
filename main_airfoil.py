# Jacob Miske
# MIT License
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


def get_airfoil_positions():
    """
    Get x and y of airfoil shape
    :return:
    """
    t_steps = 30
    t = []
    x = []
    y = []
    # Get NACA0018 profile and set to x and y
    with open("./naca0018.dat", 'rb') as file:
        for line in file:
            x_value = str(line[0:6].decode("utf-8"))
            y_value = str(line[12:19].decode("utf-8"))
            x.append(float(x_value))
            y.append(float(y_value))
    return x, y, t


def get_airfoil_slope_function(x, y):
    """
    Given a set of airfoil points, generates a slope function across airfoil
    :return:
    """
    a = []
    first = True
    for count, point in enumerate(x, 0):
        if first:
            first = False
        else:
            slope = (y[count] - y[count-1])/(x[count] - x[count-1])
            a.append(slope)
    print(a)
    return a, t


def get_airfoil_alpha_function(x, y):
    """
    Given a set of airfoil points, generates a slope function across airfoil
    :return:
    """
    points = len(x)
    alpha_1_length = float(1/points)
    alpha = []
    first = True
    for count, point in enumerate(x, 0):
        if first:
            first = False
        else:
            alpha_value = np.cos(np.arctan(abs((y[count] - y[count-1])/(x[count] - x[count-1]))))*alpha_1_length*20+1
            alpha.append(alpha_value)
    print(alpha)
    print(type(alpha))
    return alpha, t


def plot_airfoil_slope(x, a, y, alpha):
    """
    Generates plot of slope across NACA0018 airfoil
    :param x:
    :param a:
    :return:
    """
    plt.figure(1)
    plt.scatter(x=x, y=y, label="NACA0018 Airfoil")
    plt.scatter(x=x[1:], y=a, label="Slope")
    plt.scatter(x=x[1:], y=alpha, label="Alpha")
    plt.title("NACA0018 Airfoil a(x) (Tip on left, tail on right)")
    plt.xlabel("Length of wing (a.u.)")
    plt.ylabel("Height of wing (a.u.)")
    plt.grid()
    plt.legend()
    plt.savefig("NACA0018 Airfoil a(x)")
    plt.show()
    plt.close()

if __name__ == '__main__':
    x_values, y_values, t = get_airfoil_positions()
    a_values, t = get_airfoil_slope_function(x=x_values, y=y_values)
    alpha_values, t = get_airfoil_alpha_function(x=x_values, y=y_values)
    plot_airfoil_slope(x_values, a_values, y_values, alpha_values)


