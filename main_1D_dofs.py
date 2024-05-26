# Jacob Miske
# MIT License
from numpy import genfromtxt
import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import itertools
import ast
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker

# grabs closest number in myList to myNumber
# min(myList, key=lambda x:abs(x-myNumber))


font = {'family': 'serif',
        'size': 20}
csfont = {'fontname': 'Helvetica'}

matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1


def main():
    """
    Function to set up definitions
    :return:
    """
    # degrees of rotation range left to still be considered free dof
    dof_range = 20  # degrees
    return 0


def set_1D_chain(N, L_list, angle, fixed_dict, gap):
    """
    Set a 1D chain of revolute joints with lengths L_list and initial angles angle_list
    :param N: Number of cells
    :param L_list: lengths between cell i and cell i+1 (N-1 long)
    :param angle: initial angle between each cell in chain (scalar value), in degrees
    :param fixed_dict: indexs of fixed cells in chain (N long), keys are angles
    :return: x and y of each cell in chain, angle_list
    """
    # Assume fixed gap leading to 5 deg range between cells, full range is 60 deg
    # b = 2  # degrees
    full_range = 60
    # angle list is modified once the fixed_dict is parsed, starts fully collapsed at full_range
    angle_list = [full_range for _ in range(N)]
    top_angle_list = []
    bottom_angle_list = []
    # Begin by setting full free chain, first cell at x=0, y=0
    x = [0]
    y = [0]
    for i in range(N - 1):
        L_i = L_list[i]
        x.append((float(np.cos(angle)) * L_i) + x[i - 1])
        y.append(0)
    # then begin applying fixed_dict constraints and reducing angle_list elements
    for i in list(fixed_dict.keys()):
        angle_list[i] = fixed_dict[i]
    # for each element, determine the closest fixed angle and set to that angle
    for count, i in enumerate(angle_list, 0):
        # get index of closest fixed angle
        closest_fixed_cell = min(list(fixed_dict.keys()), key=lambda x: abs(x - count))
        # get index of second-closest fixed angle along chain, set angle between these two based on respective distance
        # TODO: set intermediate cells based on nearest neighbors (not just closest)
        # fixed_dict_without_closest = fixed_dict
        # fixed_dict_without_closest = fixed_dict_without_closest.pop(closest_fixed_cell, None)
        # second_closest_fixed_cell = min(list(fixed_dict_without_closest.keys()), key=lambda x: abs(x - count))
        # print(second_closest_fixed_cell)
        # set cell to the angle determined by nearest fixed cell
        angle_list[count] = fixed_dict[closest_fixed_cell]
    angles = angle_list
    for count2, i in enumerate(angles, 0):
        closest_fixed_cell = min(list(fixed_dict.keys()), key=lambda x: abs(x - count2))
        # Determine distance from each cell to fixed cell and add backlash according to ReLU relationship
        distance_to_fixed = L_list[0] * abs(count2 - closest_fixed_cell)

        # set upper and lower angle bounds for each element
        if count2 not in fixed_dict.keys():
            top_angle_list.append(np.min((60, fixed_dict[closest_fixed_cell] + gap * distance_to_fixed)))
            bottom_angle_list.append(np.max((0, fixed_dict[closest_fixed_cell] - gap * distance_to_fixed)))
        else:
            top_angle_list.append(fixed_dict[closest_fixed_cell])
            bottom_angle_list.append(fixed_dict[closest_fixed_cell])
    return x, y, angle_list, top_angle_list, bottom_angle_list


def get_DoF_from_chain(top_angle_list, bottom_angle_list):
    """
    Simply determines how many "DoF" are present based on each cell's range compared to the threshold
    :param top_angle_list:
    :param bottom_angle_list:
    :return:
    """
    dof_range = 20  # degrees threshold for DOF
    count_dof = 0
    for count, cell in enumerate(top_angle_list, 0):
        if (top_angle_list[count] - bottom_angle_list[count]) > dof_range:
            count_dof += 1
    return count_dof


def plot_chain(x, y, top_angle_list, bottom_angle_list):
    """
    Generates plot of each
    :return:
    """
    # Angular range view
    # plt.figure(0)
    # plt.xlabel("X Coordinate")
    # plt.ylabel("Y Coordinate")
    # Top view
    fig1 = plt.figure(1, figsize=(8, 8))
    ax1 = fig1.add_subplot(121, projection='3d')
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    ax1.bar3d(x, y, bottom_angle_list, dx=0.1, dy=0.1, dz=top_angle_list, shade=True)
    plt.ylim((-1, 1))
    plt.xlim((-10, 0))
    ax1.set_xlabel("X Dimension")
    ax1.set_ylabel("Y Dimension")
    ax1.set_zlabel("Angle Range")
    ax1.set_zlim(0, 90)
    plt.savefig("representation_of_angle_range_per_cell.png")
    plt.close()
    # Side view
    # plt.figure(2)
    # plt.xlabel("X Coordinate")
    # plt.ylabel("Z Coordinate")


def dual_relu(x, b):
    """ReLU returns 1 if x>0, else 0."""
    dual_relu = np.maximum(0, x - b) + np.minimum(0, x + b)
    dual_relu[0:9] = -100
    dual_relu[71:80] = 100
    return dual_relu


if __name__ == '__main__':
    # set global variables in main function
    main()
    # Set specific case variables
    N = 20
    Ls = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    Ls_2 = [2 for i in range(20)]
    L2_3 = [3 for i in range(20)]
    angle = [10]
    b = [0.75, 1, 1.5, 2, 2.5, 3]  # degrees
    fixs_dict = {0: 30, 19: 30}

    # change to rerun 1D analysis, don't use when plotting
    DOF_analysis = False

    if DOF_analysis:
        for backlash in b:
            x, y, a, ta, ba = set_1D_chain(N, Ls, angle, fixed_dict=fixs_dict, gap=backlash)
            plot_chain(x=x, y=y, top_angle_list=ta, bottom_angle_list=ba)
        # Generate angle amplitude to DoF plot
        dof_count = get_DoF_from_chain(top_angle_list=ta, bottom_angle_list=ba)
        print(dof_count)
        angle_settings = list(np.linspace(5, 55, 51))

        dof_result = [[], [], [], [], [], []]
        DO_number_at_b = []
        # run only one cell fixed
        for count, backlash in enumerate(b, 0):
            for angle in angle_settings:
                fixs_dict = {10: angle}
                x, y, a, ta, ba = set_1D_chain(N, Ls_2, angle, fixed_dict=fixs_dict, gap=backlash)
                plot_chain(x=x, y=y, top_angle_list=ta, bottom_angle_list=ba)
                # Generate angle amplitude to DoF plot
                print("Angle of single, center cell in 20 cell chain: {}".format(angle))
                dof_count = get_DoF_from_chain(top_angle_list=ta, bottom_angle_list=ba)
                print("DoF Count: {}".format(dof_count))
                dof_result[count].append(dof_count)
                DO_number_at_b.append([backlash, angle, dof_count])
        # DO distance
        a1 = 20 - max(dof_result[0])
        a2 = 20 - max(dof_result[1])
        a3 = 20 - max(dof_result[2])
        print(a1, a2, a3)

        # DO at backlash and angle figure
        backlash_list = [i[0] for i in DO_number_at_b]
        angle_list = [i[1] for i in DO_number_at_b]
        DO_list = [i[2] for i in DO_number_at_b]

        DO_number_array = np.array(DO_number_at_b)
        DO_number_array_ddx = np.gradient(DO_number_array)
        print(DO_number_array_ddx)

        # angle_list = [i[:, 0] for i in DO_number_array_ddx]
        # backlash_list = [i[:, 1] for i in DO_number_array_ddx]
        dDO_dx = [i[:, 2] for i in DO_number_array_ddx]
        print(dDO_dx)
        # zip data in form of many lists to rows

        rows = []
        # rows.append(angle_settings)
        # rows.append(dof_result)
        # rows.append(angle_list)
        # rows.append(backlash_list)
        # rows.append(DO_list)
        # collect all modeling data into one
        count = 0
        for DOF_list in dof_result:
            rows.append([angle_settings,
                         DOF_list,
                         angle_list[51 * (count):51 * (count + 1)],
                         backlash_list[51 * (count):51 * (count + 1)],
                         DO_list[51 * (count):51 * (count + 1)]])
            count += 1

        data_path = "./data_from_1D_DOF.csv"
        with open(data_path, "w") as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow(row)

    # Use data from data_from_1D_DOF.csv, not directly calculated variables
    data = []
    with open('data_from_1D_DOF.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            data.append(row)
            print('\n')
    print(type(data[0][0]))

    dof_result_0p75deg = data[0][1]
    dof_result_0p75deg = [float(i) for i in dof_result_0p75deg.strip('][').split(', ')]

    dof_result_1deg = data[1][1]
    dof_result_1deg = [float(i) for i in dof_result_1deg.strip('][').split(', ')]

    dof_result_2deg = data[3][1]
    dof_result_2deg = [float(i) for i in dof_result_2deg.strip('][').split(', ')]

    dof_result_3deg = data[5][1]
    dof_result_3deg = [float(i) for i in dof_result_3deg.strip('][').split(', ')]

    angle_list = data[0][2] + ',' + data[1][2] + ',' + data[2][2] + ',' + data[3][2] + ',' + data[4][2] + ',' + data[5][
        2]

    backlash_list = data[0][3] + ',' + data[1][3] + ',' + data[2][3] + ',' + data[3][3] + ',' + data[4][3] + ',' + \
                    data[5][3]
    DO_list = data[0][4] + ',' + data[1][4] + ',' + data[2][4] + ',' + data[3][4] + ',' + data[4][4] + ',' + data[5][4]
    angle_list = [float(i) for i in angle_list.replace('[', '').replace(']', '').split(',')]
    backlash_list = [float(i) for i in backlash_list.replace('[', '').replace(']', '').split(',')]
    DO_list = [float(i) for i in DO_list.replace('[', '').replace(']', '').split(',')]

    angle_settings = list(np.linspace(5, 55, 51))

    colors = itertools.cycle(["r", "b", "g", "k"])
    # DOF result figures
    fig = plt.figure(4, figsize=(10, 10))
    ax = fig.add_subplot()
    plt.scatter(angle_settings, dof_result_0p75deg, s=30, color=next(colors), label="dθ = 0.75 deg")
    plt.scatter(angle_settings, dof_result_1deg, s=30, color=next(colors), label="dθ = 1 deg")
    plt.scatter(angle_settings, dof_result_2deg, s=30, color=next(colors), label="dθ = 2 deg")
    plt.scatter(angle_settings, dof_result_3deg, s=30, color=next(colors), label="dθ = 3 deg")
    plt.title("Varying DoF of Revolute Joints in Series")
    plt.xlabel("Fixed Angle of Cell i=0 in 1D Chain ")
    plt.ylabel("Number of Free Revolute Cells")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # Force matplotlib to only use integers on axis markings
    plt.xlim([5, 30])
    plt.ylim((0, 20))
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 6
    plt.legend()
    ax.tick_params(axis='both', which='major', pad=10)
    plt.grid()
    plt.savefig("./figures/dof_to_fixed_cell_angle_relation.png")

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    # for pcolormesh --> , shading='flat', vmin=DO_list.min(), vmax=DO_list.max()
    sc = ax.scatter(angle_list, backlash_list, c=DO_list)
    # z_for_plot = np.array([[i*i + j*j for j in backlash_list for i in angle_list]])
    ax.set_title("Max DO Distance - Varying Backlash and Angle")
    ax.set_xlabel('Fixed Angle of Cell i=0')
    ax.set_ylabel(r"$\Delta \Theta_i$")
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 6
    ax.tick_params(axis='both', which='major', pad=10)
    plt.colorbar(sc, label="DO Distance")
    plt.savefig("./figures/DO_mapping_to_b_and_angle.png")
    plt.close()

    # alternative contour plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot()
    sc = plt.tricontourf(angle_list, backlash_list, DO_list)
    ax.set_title("Max DO Distance - Varying Backlash and Angle", pad=30)
    ax.set_xlabel('Fixed Angle of Cell i=0')
    ax.set_ylabel(r"$\Delta \Theta_i$")
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 6
    ax.tick_params(axis='both', which='major', pad=10)
    cbar = plt.colorbar(sc,
                        ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16],
                        extend='both',
                        label="Die-Off Distance"
                        )
    # plt.colorbar(sc, label="DO Distance")
    plt.savefig("./figures/DO_mapping_to_b_and_angle_contour.png")
    plt.close()

    # ReLU plot
    b1 = 10
    b2 = 20
    X = np.arange(-40, 40, 1)
    Y1 = dual_relu(X, 0)
    Y2 = dual_relu(X, b1)
    Y3 = dual_relu(X, b2)
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot()
    sc = ax.plot(X, Y1, label="Linear")
    sc = ax.plot(X, Y2, label="ReLU Model 1")
    sc = ax.plot(X, Y3, label="ReLU Model 2")
    plt.xlim([-30, 30])
    plt.ylim([-40, 40])
    plt.grid()
    ax.set_title("Unit Cell - Non-Linear Model Using ReLU")
    ax.set_xlabel('Unit Cell Angle [degrees]')
    ax.set_ylabel('Rotational Stiffness [(rev)(N/mm)]')
    plt.legend()
    plt.savefig("./figures/ReLU_demo_figure.png")
    plt.close()


    quit()
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot()
    sc = ax.scatter(angle_list, backlash_list, c=dDO_dx[1])
    ax.set_title("dDO/dx for b and drive angle", pad=30)
    ax.set_xlabel('Fixed Angle of Cell i=0')
    ax.set_ylabel(r"$\Delta \Theta_i$")
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 6
    ax.tick_params(axis='both', which='major', pad=10)
    plt.colorbar(sc, label="dDO/dx")
    plt.savefig("./figures/DO_mapping_to_b_and_angle_ddx.png")
    plt.close()



    quit()
    # example
    fig = plt.figure(5, figsize=(10, 10))
    ax = fig.add_subplot()
    y1 = np.polyfit(angle_settings, [(9 - i) / 9 for i in dof_result[0]], 3)
    y2 = np.polyfit(angle_settings, [(9 - i) / 9 for i in dof_result[1]], 3)
    y3 = np.polyfit(angle_settings, [(9 - i) / 9 for i in dof_result[2]], 3)
    plt.plot(angle_settings, y1, 'r', label="dθ = 1 deg")
    plt.plot(angle_settings, y2, 'b', label="dθ = 2 deg")
    plt.plot(angle_settings, y3, 'g', label="dθ = 3 deg")
    plt.title("Effect on Range for Joints in Series")
    plt.xlabel(r"Cell Index $i$")
    plt.ylabel("Range of Motion Remaining")
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    plt.legend()
    plt.xlim([0, 30])
    ax.tick_params(axis='both', which='major', pad=10)
    plt.grid()
    plt.savefig("./figures/dof_to_angle die off.png")
