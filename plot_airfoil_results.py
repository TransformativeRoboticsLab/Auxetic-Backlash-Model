# Jacob Miske 2024
# Airfoil from auxetic lattice
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

# default
NACA_number = '0018'

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
