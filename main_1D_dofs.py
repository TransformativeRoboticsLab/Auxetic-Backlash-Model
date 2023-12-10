# Jacob Miske
# MIT License
import numpy as np
import matplotlib.pyplot as plt

# grabs closest number in myList to myNumber
#min(myList, key=lambda x:abs(x-myNumber))


def main():
    """
    Function to set up definitions
    :return:
    """
    # Percentage of rotation range left to still be considered free dof
    dof_range = 1
    return 0


def set_1D_chain(N, L_list, angle, fixed_dict):
    """
    Set a 1D chain of revolute joints with lengths L_list and initial angles angle_list
    :param N: Number of cells
    :param L_list: lengths between cell i and cell i+1 (N-1 long)
    :param angle: initial angle between each cell in chain (scalar value), in degrees
    :param fixed_dict: indexs of fixed cells in chain (N long), keys are angles
    :return: x and y of each cell in chain, range_list
    """
    # Assume fixed gap leading to 5 deg range between cells, full range is 45 deg
    gap = 5 #degrees
    full_range = 45
    # range list is modified once the fixed_dict is parsed
    range_list = [full_range for i in range(N)]
    # Begin by setting full free chain, first cell at x=0, y=0
    x = [0]
    y = [0]
    for i in range(N-1):
        L_i = L_list[i]
        x.append(np.cos(angle)*L_i)
        y.append(0)
    # then begin applying fixed_dict constraints and reducing range_list elements
    for i in list(fixed_dict.keys()):
        range_list[i] = fixed_dict[i]
    print(range_list)
    return x, y, range_list


def plot_chain(x, y, ):
    """
    Generates plot of each
    :return:
    """
    # Angular range view
    plt.figure(0)
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")

    # Top view
    plt.figure(1)
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")

    # Side view
    plt.figure(2)

    plt.xlabel("X Coordinate")
    plt.ylabel("Z Coordinate")


if __name__ == '__main__':
    main()
    N = 10
    Ls = [1,1,1,1,1,1,1,1,1]
    angle = [10]
    fixs_dict = {2:20, 4:25}
    set_1D_chain(N, Ls, angle, fixed_dict=fixs_dict)
