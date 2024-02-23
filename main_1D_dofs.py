# Jacob Miske
# MIT License
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# grabs closest number in myList to myNumber
# min(myList, key=lambda x:abs(x-myNumber))

font = {'family': 'serif',
        'size': 18}
csfont = {'fontname': 'Helvetica'}

matplotlib.rc('font', **font)


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
    dof_range = 20  # degrees
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
    dual_relu = np.maximum(0, x-b) + np.minimum(0, x+b)
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
    b = [0.5, 1, 1.5, 2, 2.5, 3]  # degrees
    fixs_dict = {0: 30, 19: 30}

    DOF_analysis = True

    if DOF_analysis:
        for backlash in b:
            x, y, a, ta, ba = set_1D_chain(N, Ls, angle, fixed_dict=fixs_dict, gap=backlash)
            plot_chain(x=x, y=y, top_angle_list=ta, bottom_angle_list=ba)
        # Generate angle amplitude to DoF plot
        dof_count = get_DoF_from_chain(top_angle_list=ta, bottom_angle_list=ba)
        print(dof_count)
        angle_settings = list(np.linspace(5, 55, 50))

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

        # DOF result figures
        fig = plt.figure(4, figsize=(10, 10))
        ax = fig.add_subplot()
        plt.stem(angle_settings, dof_result[1], 'r', label="dθ = 1 deg")
        plt.stem(angle_settings, dof_result[3], 'b', label="dθ = 2 deg")
        plt.stem(angle_settings, dof_result[5], 'g', label="dθ = 3 deg")
        plt.title("# of DOF of Revolute Joints in Series")
        plt.xlabel("Fixed Angle of 1st Cell in Chain ")
        plt.ylabel("# of DoF")
        plt.xlim((0, 60))
        plt.ylim((0, 20))
        ax.xaxis.labelpad = 6
        ax.yaxis.labelpad = 6
        plt.legend()
        plt.xlim([0, 30])
        ax.tick_params(axis='both', which='major', pad=10)
        plt.grid()
        plt.savefig("dof_to_angle_relation 1 Ls.png")

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot()
        sc = ax.scatter(angle_list, backlash_list, c=DO_list)
        # z_for_plot = np.array([[i*i + j*j for j in backlash_list for i in angle_list]])
        ax.set_title("Max DO Distance - Varying Backlash and Driven Angle")
        ax.set_xlabel('Driven cell angle')
        ax.set_ylabel('dθ [degrees]')
        ax.xaxis.labelpad = 6
        ax.yaxis.labelpad = 6
        ax.tick_params(axis='both', which='major', pad=10)
        #plt.colorbar(sc, label="DO Distance")
        plt.savefig("./figures/DO_mapping_to_b_and_angle.png")
        plt.show()
        plt.close()

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot()
        sc = ax.scatter(angle_list, backlash_list, c=dDO_dx[1])
        ax.set_title("dDO/dx for b and drive angle")
        ax.set_xlabel('Driven cell angle')
        ax.set_ylabel('Normalized Backlash [n.d.]')
        ax.xaxis.labelpad = 6
        ax.yaxis.labelpad = 6
        ax.tick_params(axis='both', which='major', pad=10)
        plt.colorbar(sc, label="dDO/dx")
        plt.savefig("./figures/DO_mapping_to_b_and_angle_ddx.png")
        plt.show()
        plt.close()

    # ReLU plot
    b = 10
    X = np.arange(-60, 60, 1)
    Y1 = dual_relu(X, 0)
    Y2 = dual_relu(X, b)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot()
    sc = ax.plot(X, Y1, label="Linear")
    sc = ax.plot(X, Y2, label="ReLU Stiffness")
    ax.set_title("ReLU Function as Non-Linear Stiffness Model")
    ax.set_xlabel('Driven Cell Angle')
    ax.set_ylabel('Stiffness [N/mm]')
    plt.legend()
    plt.savefig("./figures/ReLU_demo_figure.png")
    plt.show()
    plt.close()


    # hardcoded example
    # fig = plt.figure(5, figsize=(10, 10))
    # ax = fig.add_subplot()
    # y1 = np.polyfit(angle_settings, [(9-i)/9 for i in dof_result[0]], 3)
    # y2 = np.polyfit(angle_settings, [(9-i)/9 for i in dof_result[1]], 3)
    # y3 = np.polyfit(angle_settings, [(9-i)/9 for i in dof_result[2]], 3)
    # plt.plot(angle_settings, y1, 'r', label="dθ = 1 deg")
    # plt.plot(angle_settings, y2, 'b', label="dθ = 2 deg")
    # plt.plot(angle_settings, y3, 'g', label="dθ = 3 deg")
    # plt.title("Effect on Range for Joints in Series")
    # plt.xlabel("Distance from i=0")
    # plt.ylabel("Range of Motion Remaining")
    # ax.xaxis.labelpad = 10
    # ax.yaxis.labelpad = 10
    # plt.legend()
    # plt.xlim([0,30])
    # ax.tick_params(axis='both', which='major', pad=10)
    # plt.grid()
    # plt.savefig("dof_to_angle die off.png")
