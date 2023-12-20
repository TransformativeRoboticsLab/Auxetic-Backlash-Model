# Jacob Miske
# MIT License
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import re

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
    with open("./data/naca{}.dat".format(NACA_number), 'rb') as file:
        for line in file:
            line_as_list_of_strings = re.findall(r'\S+', line.decode("utf-8"))
            # x_value = str(line[0:6].decode("utf-8"))
            # y_value = str(line[12:19].decode("utf-8"))
            x_value = line_as_list_of_strings[0]
            y_value = line_as_list_of_strings[1]
            x.append(float(x_value))
            y.append(float(y_value))
    return x, y, t


def get_airfoil_from_NACA_values(m, p, th, c=1):
    """
    :param m: Max chamber in percentage
    :param p: Position of max chamber in tenths of chord, less than c
    :param t: Max thickness in percentage of chord
    :param c: chord length (typically x=0 to x=1, i.e. 1)
    Reference:
    https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf
    :return:
    """
    # convert to decimal
    m = float(m/100)
    p = float(p/10)
    th = float(th/100)
    # Generic x values from 0 to chord length
    x = np.linspace(0, c, 1000)
    # y value of camber
    y_c = [0]
    # y_value of upper profile
    y_u = [0]
    # y_value of lower profile
    y_l = [0]
    # thickness dist
    y_t = [0]
    for count, i in enumerate(x, 0):
        theta = np.arctan((y_c[count] - y_c[count-1])/0.001)
        y_t.append((th/0.2)*(0.2969*np.sqrt(i) - 0.126*i - 0.3516*i**2 + 0.2843*i**3 - 0.1015*i**4) )
        if i < p:
            y_c.append((m/p**2)*(2*p*i - i**2))
        if i >= p:
            y_c.append((m/(1-p)**2)*((1-2*p)+2*p*i - i**2))
        y_u.append(y_c[count] + y_t[count]*np.cos(theta))
        y_l.append(y_c[count] + y_t[count]*np.cos(theta))
    # TODO: return other values than only x and upper profile
    return x, y_u


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
    return a, t


def get_airfoil_alpha_function(x, y):
    """
    Given a set of airfoil points, generates an alpha function across airfoil
    :return:
    """
    points = len(x)
    alpha_1_length = float(2/points)
    alpha = []
    first = True
    print(x)
    print(y)
    for count, point in enumerate(x, 0):
        if first:
            first = False
        else:
            # New alpha calculation based on length change prop to expansion ratio
            alpha_value = (np.sqrt((y[count] - y[count-1])**2+(x[count] - x[count-1])**2))/alpha_1_length + 1
            alpha.append(alpha_value)
    return alpha, t


def plot_airfoil_slope(x, a, y, alpha):
    """
    Generates plot of slope and alpha across an airfoil
    :param x:
    :param a:
    :return:
    """
    plt.figure(1)
    plt.scatter(x=x, y=y, label="NACA{} Airfoil".format(NACA_number))
    plt.scatter(x=x[1:], y=a, label="Slope")
    plt.scatter(x=x[1:], y=alpha, label="Alpha")
    plt.title("NACA{} Airfoil a(x) (Tip on left, tail on right)".format(NACA_number))
    plt.xlabel("Length of wing (a.u.)")
    plt.ylabel("Height of wing (a.u.)")
    plt.ylim((-3, 3))
    plt.grid()
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil a(x)".format(NACA_number))
    # plt.show()
    plt.close()


if __name__ == '__main__':
    NACA_numbers = ["0006", "0009", "0018", "0021", "0024", "1408", "1410", "2408", "4412",
                    "4415", "6409", "6412"]

    data = False
    if data:
        for number in NACA_numbers:
            NACA_number = number
            x_values, y_values, t = get_airfoil_positions()
            a_values, t = get_airfoil_slope_function(x=x_values, y=y_values)
            alpha_values, t = get_airfoil_alpha_function(x=x_values, y=y_values)
            plot_airfoil_slope(x_values, a_values, y_values, alpha_values)
    generated = True
    if generated:
        for number in NACA_numbers:
            NACA_number = number
            m1 = float(NACA_number[:1])
            p1 = float(NACA_number[1:2])
            th1 = float(NACA_number[2:4])
            t=1
            x_values, y_values = get_airfoil_from_NACA_values(m=m1, p=p1, th=th1, c=1)
            y_values = y_values[:1000]
            a_values, t = get_airfoil_slope_function(x=x_values, y=y_values)
            alpha_values, t = get_airfoil_alpha_function(x=x_values, y=y_values)
            plot_airfoil_slope(x_values, a_values, y_values, alpha_values)


