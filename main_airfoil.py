# Jacob Miske
# MIT License
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
import re
import pandas as pd

font = {'family': 'serif',
        'size': 16}
csfont = {'fontname':'Helvetica'}

matplotlib.rc('font', **font)


def get_airfoil_positions(NACA_number):
    """
    Get x and y of airfoil shape
    :return:
    """
    t_steps = 30
    t = []
    x = []
    y = []
    # Get NACAXXXX profile and set to x and y
    with open("./data/naca_profiles/naca{}.dat".format(NACA_number), 'rb') as file:
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

def get_top_airfoil_contour_from_NACA_values(m, p, th, c=1):
    """
    Similar to get_airfoil_from_NACA_values, but only returns top contour values
    :param m: Max chamber in percentage
    :param p: Position of max chamber in tenths of chord, less than c
    :param t: Max thickness in percentage of chord
    :param c: chord length (typically x=0 to x=1, i.e. 1)
    Reference:
    https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf
    :return:
    """
    # convert percentages to floating, decimal values
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
            alpha_value = np.cos(np.arctan(abs((y[count] - y[count - 1]) / (x[count] - x[count - 1])))) * (1-0.2*x[count]) * alpha_1_length + 1
            # Non-constrained, alpha function based on length change prop to expansion ratio (projective)
            # alpha_value = (np.sqrt((y[count] - y[count-1])**2+(x[count] - x[count-1])**2))/alpha_1_length + 1
            alpha.append(alpha_value)
    return alpha, t


def get_airfoil_alpha_function_slope(x, y):
    """
    Given a set of airfoil points, generates an alpha function across airfoil based on slope across the foil
    :return:
    """
    points = len(x)
    chord_width = float(max(x) - min(x))
    flat_cell_width = chord_width/(2*float(points))
    alpha = []
    # for every x value except the last
    for count, point in enumerate(x[:-1], 0):
        slope = abs((y[count+1] - y[count])/(x[count+1] - x[count]))
        # to find distance between cells when extended
        arctan_angle = np.arctan(slope)
        distance_to_next_cell = flat_cell_width/(np.cos(arctan_angle))
        alpha.append(distance_to_next_cell/flat_cell_width)
    return alpha, t


def get_airfoil_mocap():
    """
    From optitrack file, get the right x,y, and z values for each cell
    :return:
    """
    df = pd.read_csv('./data/motion_capture/MLS 45s on column 4 20240319.csv', sep=',', header=None)
    first_points = list(df.values[0])
    # convert to normalized airfoil dims
    xs = [(i/300)-0.385 for i in first_points[4::3]]
    zs = [j/700 for j in first_points[3::3]]
    plt.figure(0)
    plt.scatter(xs, zs)
    plt.savefig("./figures/motion_capture_test_20240319.png")
    plt.close()
    # data = zip(xs, zs)
    # data = sorted(data, key=lambda k: [k[1], k[0]])
    # xs = data[0:]
    # zs = data[1:]
    # print(xs)
    return xs, zs


def get_airfoil_mocap_multitest():
    """
    From optitrack files, get the right x, y, and z values for each cell
    :return:
    """
    path = './data/motion_capture/multiform/'
    extension = 'csv'
    files = glob.glob(path+'*.{}'.format(extension))
    print(files)
    count = 1
    xs_lists = []
    zs_lists = []
    filenames = []
    for file in files:
        df = pd.read_csv(str(file), sep=',', header=None)
        # remove side points
        # df = df.drop(df.columns[[2, 3, 4, 14, 15, 16, 35, 36, 37]], axis=1)
        # take first row
        first_points = list(df.values[0])
        # convert to normalized airfoil dims
        xs = [(i/300)-0.46 for i in first_points[2::3]]
        zs = [(j/700)+0.011 for j in first_points[3::3]]
        xs_measured = [i for i in first_points[2::3]]
        zs_measured = [j for j in first_points[3::3]]
        fig, ax = plt.figure(0)
        plt.scatter(xs_measured, zs_measured)
        ax.grid(False)
        plt.savefig("./figures/airfoils/motion_capture_multitest_{}.png".format(count))
        plt.close()
        # data = zip(xs, zs)
        # data = sorted(data, key=lambda k: [k[1], k[0]])
        # xs = data[0:]
        # zs = data[1:]
        # print(xs)
        count += 1
        xs_lists.append(xs)
        zs_lists.append(zs)
        filenames.append(str(file))
    return xs_lists, zs_lists, filenames


def get_airfoil_photo_capture():
    """
    Grabs x and z data from photo capture
    :return:
    """
    path = './data/photo_capture/'
    extension = 'csv'
    files = glob.glob(path + '*.{}'.format(extension))
    print(files)
    count = 1
    xs_lists = []
    zs_lists = []
    filenames = []
    for file in files:
        df = pd.read_csv(str(file), sep=',', header=None)
        # convert to normalized airfoil dims if necessary
        xs = [(-i/400) + 1 for i in df.iloc[:,0]]
        zs = [j/2000 for j in df.iloc[:,1]]
        # measured values
        xs_measured = [i for i in xs]
        zs_measured = [j for j in zs]
        plt.figure(count)
        plt.scatter(xs_measured, zs_measured)
        plt.savefig("./figures/airfoils/photo_capture_multitest_{}.png".format(count))
        plt.close()
        # data = zip(xs, zs)
        # data = sorted(data, key=lambda k: [k[1], k[0]])
        # xs = data[0:]
        # zs = data[1:]
        # print(xs)
        count += 1
        xs_lists.append(xs)
        zs_lists.append(zs)
        filenames.append(str(file))
    return xs_lists, zs_lists, filenames


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


def validate_alpha_function(x, y, alpha_list, number):
    """
    Given a list of alpha across a curve,
    compare the resulting distance to true airfoil circumference
    :param x: x positions of each cell
    :param y: respective y positions of each cell, one to one lined up with x
    :param alpha_list: alpha at each cell corresponding to x and y
    :return:
    """
    points = len(x)
    chord_width = float(max(x) - min(x))
    flat_cell_width = 2*chord_width / (float(points))
    print(flat_cell_width)
    distance_contour = []
    distance_alpha = []
    length_of_airfoil_to_alpha = 0
    length_of_airfoil_contour = 0
    for count, i in enumerate(x[:-1], 0):
        # for each point in ideal airfoil, calculate distance between two points and sum
        distance_i = np.sqrt((y[count] - y[count + 1])**2 + (x[count] - x[count + 1])**2)
        distance_contour.append(distance_i)
        length_of_airfoil_contour += distance_i
        # for each alpha, determine distance and sum
        distance_i_alpha = flat_cell_width*alpha_list[count]
        distance_alpha.append(distance_i_alpha)
        length_of_airfoil_to_alpha += distance_i_alpha
    # calculate error
    flat_total_distace = len(x) * flat_cell_width
    print(flat_total_distace)
    error = length_of_airfoil_to_alpha - length_of_airfoil_contour
    error_percentage = (error/length_of_airfoil_contour)*100
    # Check difference visually
    fig1, ax1 = plt.subplots()
    ax1.scatter(x[:-1], distance_contour, c='orange', label='Contour Distances')
    ax1.scatter(x[:-1], distance_alpha, c='r', label='Alpha Distances')
    ax2 = ax1.twinx()
    ax2.set_ylabel('alpha')
    ax2.scatter(x[:-1], alpha_list)
    fig1.set_figheight(9)
    fig1.set_figwidth(9)
    plt.title("Alpha Validation Function")
    ax2.set_ylabel("Alpha(x) Function")
    plt.xlabel("Chord Length of Foil (a.u.)")
    plt.savefig("./figures/validate_alpha_{}.png".format(number))
    plt.legend()
    plt.close()

    # linear shift to alpha
    distance_contour = [i*18 for i in distance_contour]

    # Curve Validation figure visually
    x_cells = [0, 0.090, 0.180, 0.272, 0.360, 0.45, 0.54, 0.63, 0.72, 0.81, 0.9, 1]
    y_cells = [0, 0.0455, 0.071, 0.089, 0.085, 0.071, 0.061, 0.050, 0.0412, 0.032, 0.024, 0]
    fig1, ax1 = plt.subplots()
    ax1.scatter(x, y, label='Airfoil Points')
    ax1.scatter(x_cells, y_cells, label='Measured Curve')
    ax1.set_ylim([0, 0.2])
    ax1.set_xlabel("Width of Airfoil (Arbitrary)")
    ax1.set_ylabel("Height of Airfoil (Arbitrary)")
    ax2 = ax1.twinx()
    ax2.scatter(x[:-1], distance_contour, label="Conformal Alpha Function", c='r')
    ax2.set_ylim([1, 2])
    fig1.set_figheight(9)
    fig1.set_figwidth(9)
    plt.title("Alpha Validation - NACA 0018 Experiment")
    ax2.set_ylabel("Alpha Function")
    plt.xlabel("Chord Length of Foil (a.u.)")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc=1)
    plt.savefig("./figures/Experimental_alpha_{}.png".format(number))
    plt.close()

    # linear shift to angle
    distance_contour = [round(i*35.5, 1) for i in distance_contour]

    # Curve Validation figure visually
    x_cells = [0, 0.090, 0.180, 0.272, 0.360, 0.45, 0.54, 0.63, 0.72, 0.81, 0.9, 1]
    y_cells = [0, 0.0455, 0.071, 0.089, 0.085, 0.071, 0.061, 0.050, 0.0412, 0.032, 0.024, 0]
    fig1, ax1 = plt.subplots()
    ax1.scatter(x, y, label='Airfoil Points')
    ax1.scatter(x_cells, y_cells, label='Measured Curve')
    ax1.set_ylim([0, 0.2])
    ax1.set_xlabel("Width of Airfoil (Arbitrary)")
    ax1.set_ylabel("Height of Airfoil (Arbitrary)")
    ax2 = ax1.twinx()
    ax2.scatter(x[1:20], distance_contour[1:20], label="Conformal Alpha Angle", c='r')
    ax2.set_ylim([15, 75])
    fig1.set_figheight(9)
    fig1.set_figwidth(9)
    plt.title("Alpha Validation - NACA 0018 Experiment")
    ax2.set_ylabel("Alpha Angle")
    plt.xlabel("Chord Length of Foil (a.u.)")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc=1)
    plt.savefig("./figures/angle_alpha_{}.png".format(number))
    plt.close()

    return error, error_percentage


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


def plot_airfoil_slope(x, a, y, alpha, x_exp, z_exp):
    """
    Generates plot of slope and alpha across an airfoil
    :param x:
    :param a:
    :return:
    """
    plt.figure(1, figsize=(9, 8), dpi=80)
    plt.scatter(x=x, y=y, label="NACA{} Airfoil".format(NACA_number))
    plt.scatter(x=x[1:], y=a, label="Slope")
    plt.scatter(x=x[1:], y=alpha, label=r"$\alpha$ Function")
    plt.title("NACA{} Airfoil a(x) (Tip on left, tail on right)".format(NACA_number), **csfont)
    plt.xlabel("Length of wing (a.u.)", **csfont)
    plt.ylabel("Height of wing (a.u.)", **csfont)
    plt.ylim((-3, 3))
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil a(x)".format(NACA_number))
    # plt.show()
    plt.close()
    return 0


def plot_airfoil_slope_experimental(x_exp, z_exp, x_a, alpha, x_exp_2, z_exp_2, x_a_2, alpha_2):
    """
    Generates plot of slope and alpha across an airfoil, two profiles for blended wing
    :param x_exp:
    :param z_exp:
    :param x_a:
    :param alpha:
    :return:
    """
    blended_wing = True
    # get perfect airfoil points from equation directly
    m1 = float(NACA_number[:1])
    p1 = float(NACA_number[1:2])
    th1 = float(NACA_number[2:4])
    x_values, z_values = get_airfoil_from_NACA_values(m=m1, p=p1, th=th1, c=1)
    z_values = z_values[:1000] # cut z values down to only 1000
    # Get 0018 airfoil
    x_values_0018, z_values_0018 = get_airfoil_from_NACA_values(m=0, p=0, th=18, c=1)
    z_values_0018 = z_values_0018[:1000]  # cut z values down to only 1000
    # Plot to compare to experimental results
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.set_figheight(9)
    fig.set_figwidth(9)
    # FIRST ROW PLOT use data set on angles in /data/
    ax1.plot(x_values, z_values, '-', c='m', label="NACA{} Airfoil".format(NACA_number))
    ax1.plot(x_values_0018, z_values_0018, '-', c='r', label="NACA{} Airfoil".format("0018"))
    ax1.scatter(x_exp, z_exp, marker="^", s=50, c='b', label="Experimental NACA2408")
    if blended_wing:
        ax1.scatter(x_exp_2, z_exp_2, marker="^", s=50, c='orange', label="Experimental NACA0018")
    # sort the experimental values from smallest x to largest x similar to x, z
    x_exp, z_exp = (list(x) for x in zip(*sorted(zip(x_exp, z_exp))))
    # for each x_exp value, get closest true z_value
    z_true = []
    for i in x_exp:
        interp_value = np.interp(i, x_values, z_values)
        z_true.append(interp_value)
    print(z_exp)
    print(z_true)
    z_error = [abs(a_i - b_i) for a_i, b_i in zip(z_exp, z_true)]
    z_error_2 = [abs(a_i - b_i)*15 for a_i, b_i in zip(z_exp_2, z_true)]

    # Set first subplot y_label
    ax1.set_ylabel("Height of wing \n (Arbitrary Length)", **csfont)
    fig.suptitle(r"Blended Body NACA Foil - Experiment and $\alpha$(x) Comparison".format(NACA_number), **csfont)
    # Used to generate dual Y axis plot # ax2 = ax1.twinx()
    # SECOND ROW PLOT error between experiment and real NACA profile
    # z_error is in decimal, so we convert to percentage here:
    z_error = [i*100 for i in z_error]
    # linear leveling
    ax2.plot(x_exp, z_error, '-x', c='m')
    if blended_wing:
        ax2.plot(x_exp_2, z_error_2, '-o', c='m')
    # 0 to 2 percent of full length error
    ax2.set_ylim([0,2])
    ax2.set_ylabel("Experimental Error \n (% of Full Length)", **csfont)

    # THIRD ROW PLOT use data set on angles in /data/
    ax3.scatter(x_a, alpha, marker='o', c='#000000', label=r"$\alpha$ - 0018")
    ax3.set_ylabel('Alpha Function \n (degrees)', **csfont)
    if blended_wing:
        ax3.scatter(x_a_2, alpha_2, marker='x', c='#222222', label=r"$\alpha$ - 2408")
    # ax4 = ax3.twinx()
    # alpha_expansion_ratio = alpha
    # ax4.scatter(x_alpha, alpha_expansion_ratio, c='#000000')
    # ax4.set_ylabel("")

    plt.xlabel("Length of wing \n (Arbitrary Length)", **csfont)
    # generate limits and legend
    ax1.set_xlim((0, 1))
    ax2.set_xlim((0, 1))
    ax3.set_xlim((0, 1))
    ax1.set_ylim((0, 0.2))
    ax3.set_ylim((20, 60))
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    h3, l3 = ax3.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc=1)
    ax3.legend(h3, l3, loc=1)


    plt.savefig("./figures/NACA{} Airfoil Experimental Blended Wing {}.png".format(NACA_number, count))
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
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    sc = ax.scatter(b_list, L_list, c=max_K_list)
    ax.set_title("Max Curvature given Cell Size and Backlash")
    ax.set_xlabel('Normalized Backlash [n.d.]')
    ax.set_ylabel('Cell Size [mm]')
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.tick_params(axis='both', which='major', pad=10)
    plt.colorbar(sc, label="Max Curvature [1/mm]")
    plt.savefig("./figures/max_curvature_2d.png")
    # plt.show()
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
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    sc = ax.scatter(b_list, L_list, c=do_list)
    ax.set_title("Die-off Distance given Cell Size and Backlash")
    ax.set_xlabel('Normalized Backlash [n.d.]')
    ax.set_ylabel('Cell Size [mm]')
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.tick_params(axis='both', which='major', pad=10)
    plt.colorbar(sc, label="Die-off [cells]")
    plt.savefig("./figures/DO_distance_2d.png")
    # plt.show()
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
    plt.ylim([-0.03, 0])
    plt.grid()
    plt.legend()
    plt.savefig("./figures/NACA{} Airfoil error.png".format(NACA_number))
    plt.show()
    plt.close()
    return 0


if __name__ == '__main__':
    print("TRL Airfoil Study Model")
    NACA_numbers = ["0018", "0024", "1408", "1410", "2408", "4412"]
    print("Examining airfoils: {}".format(str(NACA_numbers)))

    study_cell_geometry = False
    if study_cell_geometry:
        # Determine characteristic relationships
        # Take ranges for thickness, backlash, cell size [in millimeters]
        d_values = np.linspace(0.1, 1.1, 10)
        b_values = np.linspace(0.2, 2.1, 10)
        L_values = np.linspace(10, 50, 100)
        # 5 mm default cell size
        default_cell_size = 5
        data = []
        # Run through range of reasonable values for b, d, and L
        print("For each thickness, backlash, and cell size - determining b_norm, DO, and dDO/dx")
        for thickness in d_values:
            for backlash in b_values:
                for cell_size in L_values:
                    b_norm = get_normalized_backlash(b= backlash, L=cell_size)
                    max_k = get_max_curvature(d=thickness, b=backlash, L=cell_size)
                    do = get_die_off(d=thickness, b=backlash, L=cell_size, angle_range=60)/cell_size
                    data.append([thickness, backlash, cell_size, b_norm, max_k, do])

        # Change data type for input array
        data_array = np.array(data)
        # Generate plots from the data
        plot_b_L_maxK(b_list=data_array[:, 3], L_list=data_array[:, 2], max_K_list=data_array[:, 4])
        plot_b_L_dieoff(b_list=data_array[:, 3], L_list=data_array[:, 2], do_list=data_array[:, 5])

    airfoil_data = True
    if airfoil_data:
        # for each NACA airfoil profile under study
        for number in NACA_numbers:
            NACA_number = number
            print("Calculating Alpha(x) for airfoil: {}".format(number))
            x_values, y_values, t = get_airfoil_positions(number)
            a_values, t = get_airfoil_slope_function(x=x_values, y=y_values)

            # Original alpha function
            # alpha_values, t = get_airfoil_alpha_function(x=x_values, y=y_values)

            # Slope based alpha method
            alpha_values, t = get_airfoil_alpha_function_slope(x=x_values, y=y_values)

            # validate alpha function
            e, ep = validate_alpha_function(x=x_values, y=y_values, alpha_list=alpha_values, number=number)

            # For two particular profiles, get slope and alpha plots against experimental data
            # pulling angles from data directly, TODO: code in data pull
            x_alpha = list(np.linspace(0, 0.95, 22))
            # angles set manually on RADs lattice and recorded here
            alpha_1 = [42, 42, 41, 38, 35, 32, 31, 30, 29, 27, 26, 25, 25, 27, 29, 31, 34, 37, 39, 41, 42, 44]  # Degrees
            alpha_2 = [45, 42, 39, 38, 37, 35, 34, 33, 32, 31, 30, 30, 31, 32, 34, 35, 37, 38, 40, 42, 44, 45]  # Degrees

            if NACA_number == "0018":
                # x_lists, z_lists, fnames = get_airfoil_mocap_multitest()
                x_lists_0018, z_lists_0018, fnames = get_airfoil_photo_capture()
                for count, sample in enumerate(x_lists_0018, 0):
                    plot_airfoil_slope_experimental(x_exp=x_lists_0018[count], z_exp=z_lists_0018[count], x_a=x_alpha, alpha=alpha_1,
                                                    x_exp_2=x_lists_0018[count], z_exp_2=x_lists_0018[count], x_a_2=x_alpha, alpha_2=alpha_2)
            if NACA_number == "2408":
                # x_lists, z_lists, fnames = get_airfoil_mocap_multitest()
                x_lists_2408, z_lists_2408, fnames = get_airfoil_photo_capture()
                for count, sample in enumerate(x_lists_2408, 0):
                    plot_airfoil_slope_experimental(x_exp=x_lists_0018[0], z_exp=z_lists_0018[0], x_a=x_alpha, alpha=alpha_1,
                                                    x_exp_2=x_lists_2408[2], z_exp_2=z_lists_2408[2], x_a_2=x_alpha, alpha_2=alpha_2)

    # Look at generate NACA profile from model
    generated = False
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




