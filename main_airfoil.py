# Jacob Miske
# MIT License
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
import re

font = {'family' : 'serif',
        'size'   : 16}
csfont = {'fontname':'Helvetica'}

matplotlib.rc('font', **font)


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
    for count, point in enumerate(x, 0):
        if first:
            first = False
        else:
            # Constrained alpha function, gives a(x) = 1 at tip and tail of airfoil
            alpha_value = np.cos(np.arctan(abs((y[count] - y[count - 1]) / (x[count] - x[count - 1])))) * alpha_1_length + 1
            # Non-constrained, alpha function based on length change prop to expansion ratio (projective)
            # alpha_value = (np.sqrt((y[count] - y[count-1])**2+(x[count] - x[count-1])**2))/alpha_1_length + 1
            alpha.append(alpha_value)
    return alpha, t


def get_camber_length(x, y):
    """
    Returns float of the airfoil profiles camber length given x and y points along the airfoil
    :param x: x values of airfoil points
    :param y: y values of airfoil points
    :return:
    """
    # Sum up floating point number
    L_camber = 0.0
    for count, i in enumerate(x[1:], 0):
        L_i = np.sqrt((x[count] - x[count-1])**2 + (y[count] - y[count-1])**2)
        L_camber += L_i
    return L_camber


def get_normalized_backlash(b, L):
    """
    For a revolute robotic joint, return normalized backlash value
    :param b: outer radius minus inner radius of joint, backlash
    :param L: Float for cell size
    :return: Float between 0 and 1, typically far less than 0.1
    """
    # Init b_norm
    b_norm = b/L
    return b_norm


def get_max_curvature(d, b, L):
    """
    Given a set of cells of thickness d, backlash b, and cell link length L, get max curvature.
    All metrics need to be in the same units of length
    :param d: Float for cell thickness
    :param b: Float for cell backlash
    :param L: Float for cell size
    :return:
    """
    theta = np.arctan(b/d)
    max_curvature = np.sin(theta)/L
    return max_curvature


def get_die_off(d, b, L, angle_range):
    """
    Given a set of cells of thickness d, backlash b, and cell link length L, get die off distance.
    All metrics need to be in the same units of length
    :param d: Float for cell thickness
    :param b: Float for cell backlash
    :param L: Float for cell size
    :param angle_range: degrees of rotation on cell
    :return:
    """
    delta_theta = np.arcsin(2*b/L)
    res = np.floor(angle_range/delta_theta)
    return res


def get_airfoil_error(m, p, th, c, cell_size):
    """
    From the airfoil points and a cell size,
    :param m: Max chamber in percentage
    :param p: Position of max chamber in tenths of chord, less than c
    :param t: Max thickness in percentage of chord
    :param c: chord length (typically x=0 to x=1, i.e. 1)
    :param cell_size: length of cell mechanism making the airfoil
    Reference:
    https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf
    :return:
    """
    y_error = []
    # convert to decimal
    m = float(m / 100)
    p = float(p / 10)
    th = float(th / 100)
    # Set of 1000 x values from 0 to chord length
    x = np.linspace(0, c, 1000)
    # the x location of each cell
    x_cells = np.linspace(0, c, int(round(cell_size)))
    # y value of camber, start at origin
    y_c = [0]
    # y_value of upper profile
    y_u = [0]
    # y_value of lower profile
    y_l = [0]
    # thickness dist
    y_t = [0]
    for count, i in enumerate(x, 0):
        theta = np.arctan((y_c[count] - y_c[count - 1]) / 0.001)
        y_t.append((th / 0.2) * (0.2969 * np.sqrt(i) - 0.126 * i - 0.3516 * i ** 2 + 0.2843 * i ** 3 - 0.1015 * i ** 4))
        if i < p:
            y_c.append((m / p ** 2) * (2 * p * i - i ** 2))
        if i >= p:
            y_c.append((m / (1 - p) ** 2) * ((1 - 2 * p) + 2 * p * i - i ** 2))
        y_upper_value = y_c[count] + y_t[count] * np.cos(theta)
        y_u.append(y_upper_value)
        y_lower_value = y_c[count] + y_t[count] * np.cos(theta)
        y_l.append(y_lower_value)
    # interpolation for error calculation
    y_interp = np.interp(x_cells, x, y_u[:1000])
    # plt.figure(1)
    # plt.plot(x_cells, y_interp, 'o')
    # plt.show()
    # plt.close()
    y_interp_approx = np.interp(x, x_cells, y_interp)
    for count, i in enumerate(y_interp_approx, 0):
        y_error.append(y_interp_approx[count] - y_u[count])
    return x, y_error


def plot_airfoil_slope(x, a, y, alpha):
    """
    Generates plot of slope and alpha across an airfoil
    :param x:
    :param a:
    :return:
    """
    plt.figure(1, figsize=(8, 6), dpi=80)
    plt.scatter(x=x, y=y, label="NACA{} Airfoil".format(NACA_number))
    plt.scatter(x=x[1:], y=a, label="Slope")
    plt.scatter(x=x[1:], y=alpha, label="Alpha")
    plt.title("NACA{} Airfoil a(x) (Tip on left, tail on right)".format(NACA_number), **csfont)
    plt.xlabel("Length of wing (a.u.)", **csfont)
    plt.ylabel("Height of wing (a.u.)", **csfont)
    plt.ylim((-3, 3))
    plt.grid()
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil a(x)".format(NACA_number))
    # plt.show()
    plt.close()
    plt.figure(2, figsize=(8, 6), dpi=80)
    plt.scatter(x=x, y=y, label="NACA{} Airfoil".format(NACA_number))
    plt.title("NACA{} Airfoil a(x) (Tip on left, tail on right)".format(NACA_number), **csfont)
    plt.xlabel("Length of wing (a.u.)", **csfont)
    plt.ylabel("Height of wing (a.u.)", **csfont)
    plt.ylim((0, 0.2))
    plt.grid()
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil".format(NACA_number))
    # plt.show()
    plt.close()


def plot_b_L_maxK(b_list, L_list, max_K_list):
    """
    Takes in three lists, generated 3D plot
    :param b_list:
    :param L_list:
    :param max_K_list:
    :return:
    """
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(b_list, L_list, max_K_list)
    ax.set_title("Max Curvature as Function of Cell Size and Normalized Backlash")
    ax.set_xlabel('Normalized Backlash [unitless]')
    ax.set_ylabel('Cell Size [mm]')
    ax.set_zlabel('Max Curvature [1/mm]')
    ax.xaxis.labelpad = 30
    ax.yaxis.labelpad = 30
    ax.zaxis.labelpad = 30
    ax.view_init(elev=30, azim=60, roll=0)
    plt.savefig("./figures/max_curvature.png")
    plt.show()
    plt.close()
    return 0


def plot_b_L_dieoff(b_list, L_list, do_list):
    """
    Takes in three lists, generated 3D plot of die-off distance
    :param b_list:
    :param L_list:
    :param do_list:
    :return:
    """
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(b_list, L_list, do_list)
    ax.set_title("Die-off Distance as Function of Cell Size and Normalized Backlash")
    ax.set_xlabel('Normalized Backlash')
    ax.set_ylabel('Cell Size')
    ax.set_zlabel('Die-Off Distance')
    ax.xaxis.labelpad = 30
    ax.yaxis.labelpad = 30
    ax.zaxis.labelpad = 30
    plt.savefig("./figures/DO_distance.png")
    plt.show()
    plt.close()
    return


def plot_airfoil_error(x, y):
    """
    Plots the error between rigid cells and true airfoil shape
    :param x:
    :param y:
    :return:
    """
    plt.figure(1, figsize=(8, 6), dpi=80)
    plt.scatter(x=x, y=y, label="NACA{} Airfoil Error".format(NACA_number))
    plt.title("NACA{} Airfoil Error with Rigid Cells".format(NACA_number), **csfont)
    plt.xlabel("Length of wing (a.u.)", **csfont)
    plt.ylabel("Error from True Airfoil (a.u.)", **csfont)
    plt.grid()
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil error (5 cells)".format(NACA_number))
    plt.show()
    plt.close()
    return 0


if __name__ == '__main__':
    print("TRL Airfoil Functions")
    NACA_numbers = ["0018", "0024", "1408", "1410", "2408", "4412"]

    # Determine characteristic relationships
    # Take ranges for thickness, backlash, cell size [in millimeters]
    d_values = np.linspace(0.1, 1.1, 10)
    b_values = np.linspace(0.1, 2.1, 10)
    L_values = np.linspace(10, 50, 100)
    # 5 mm default cell size
    default_cell_size = 5
    data = []
    # Run through range of reasonable values for b, d, and L
    for thickness in d_values:
        for backlash in b_values:
            for cell_size in L_values:
                b_norm = get_normalized_backlash(b= backlash, L=cell_size)
                max_k = get_max_curvature(d=thickness, b=backlash, L=cell_size)
                do = get_die_off(d=thickness, b=backlash, L=cell_size, angle_range=60)
                data.append([thickness, backlash, cell_size, b_norm, max_k, do])
    # Change data type for input array
    data_array = np.array(data)
    # Generate plots from the data
    # plot_b_L_maxK(b_list=data_array[:, 3], L_list=data_array[:, 2], max_K_list=data_array[:, 4])
    # plot_b_L_dieoff(b_list=data_array[:, 3], L_list=data_array[:, 2], do_list=data_array[:, 5], )

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
        # Gather plots of NACA profile based on mathematical model and gather alpha funcs
        for number in NACA_numbers:
            NACA_number = number
            m1 = float(NACA_number[:1])
            p1 = float(NACA_number[1:2])
            th1 = float(NACA_number[2:4])
            t=1
            x_values, y_values = get_airfoil_from_NACA_values(m=m1, p=p1, th=th1, c=1)
            y_values = y_values[:1000]
            # a_values, t = get_airfoil_slope_function(x=x_values, y=y_values)
            # alpha_values, t = get_airfoil_alpha_function(x=x_values, y=y_values)
            # plot_airfoil_slope(x_values, a_values, y_values, alpha_values)

            # for each NACA number of interest with chord line length = 1, set cell size and determine error
            L_c = get_camber_length(x=x_values, y=y_values)
            print("Camber line length")
            print(L_c)
            x, y_err = get_airfoil_error(m=m1, p=p1, th=th1, c=1, cell_size=default_cell_size)
            plot_airfoil_error(x, y_err)


