# Jacob Miske
# MIT License
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# grabs closest number in myList to myNumber
# min(myList, key=lambda x:abs(x-myNumber))

font = {'family': 'serif',
        'size': 16}
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


def set_1D_chain(N, L_list, angle, fixed_dict):
    """
    Set a 1D chain of revolute joints with lengths L_list and initial angles angle_list
    :param N: Number of cells
    :param L_list: lengths between cell i and cell i+1 (N-1 long)
    :param angle: initial angle between each cell in chain (scalar value), in degrees
    :param fixed_dict: indexs of fixed cells in chain (N long), keys are angles
    :return: x and y of each cell in chain, angle_list
    """
    # Assume fixed gap leading to 5 deg range between cells, full range is 60 deg
    gap = 2  # degrees
    full_range = 60
    # angle list is modified once the fixed_dict is parsed, starts fully collapsed at full_range
    angle_list = [full_range for i in range(N)]
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
        distance_to_fixed = abs(count2 - closest_fixed_cell)

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


if __name__ == '__main__':
    # set global variables in main function
    main()
    # Set specific case variables
    N = 20
    Ls = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    angle = [10]
    fixs_dict = {0: 30, 19: 30}
    x, y, a, ta, ba = set_1D_chain(N, Ls, angle, fixed_dict=fixs_dict)
    plot_chain(x=x, y=y, top_angle_list=ta, bottom_angle_list=ba)
    # Generate angle amplitude to DoF plot
    dof_count = get_DoF_from_chain(top_angle_list=ta, bottom_angle_list=ba)
    print(dof_count)
    angle_settings = list(np.linspace(5, 55, 50))
    dof_result = []
    # run only one cell fixed
    for angle in angle_settings:
        fixs_dict = {10: angle}
        x, y, a, ta, ba = set_1D_chain(N, Ls, angle, fixed_dict=fixs_dict)
        plot_chain(x=x, y=y, top_angle_list=ta, bottom_angle_list=ba)
        # Generate angle amplitude to DoF plot
        print("Angle of single, center cell in 20 cell chain: {}".format(angle))
        dof_count = get_DoF_from_chain(top_angle_list=ta, bottom_angle_list=ba)
        print("DoF Count: {}".format(dof_count))
        dof_result.append(dof_count)

    fig = plt.figure(4, figsize=(8, 8))
    ax = fig.add_subplot()
    plt.stem(angle_settings, dof_result)
    plt.title("# of DOF of Revolute Joints in Series")
    plt.xlabel("Fixed Angle of 1st Cell in Chain ")
    plt.ylabel("# of DoF")
    plt.xlim((0, 60))
    plt.ylim((0, 12))
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.tick_params(axis='both', which='major', pad=10)
    plt.grid()
    plt.savefig("dof_to_angle_relation 20 deg 20 cells.png")
