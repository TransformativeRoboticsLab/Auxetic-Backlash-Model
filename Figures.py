import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math



def figure_1():
    mu = 0
    variance = 10
    sigma = math.sqrt(variance)
    x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
    x_floor = [np.floor(i) for i in x]

    y = stats.norm.pdf(x, mu, sigma)
    y = [100*i for i in y]
    y = [np.floor(i)+1 for i in y]

    variance_2 = 4
    sigma_2 = math.sqrt(variance_2)
    x_2 = np.linspace(mu - 3 * sigma_2, mu + 3 * sigma, 100)
    x_floor_2 = [np.floor(i) for i in x_2]
    y_2 = stats.norm.pdf(x_2, mu, sigma_2)
    y_2 = [65 * i for i in y_2]
    y_2 = [np.floor(i) + 1 for i in y_2]

    plt.plot(x, y, 'r', label="L=1, l_gap=0.2L")
    plt.plot(x_2, y_2, 'b', label="L=1, l_gap=0.3L")
    plt.title("Degree of Freedom in 5x5 grid of revolute joints")
    plt.legend()
    plt.ylabel("Degrees of Freedom")
    plt.xlabel("Angle Change in Driven Revolute Joint")
    plt.show()
    return 0


def figure_3():
    x_3 = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    y_3 = [0, 0.1, 0.3, 0.4, 0.44, 0.43, 0.4, 0.2, 0.0]
    plt.plot(x_3, y_3, 'b', label="Frustrated cells")
    plt.title("Frustated N Joints to Curve")
    plt.legend()
    plt.ylabel("Distance (y)")
    plt.xlabel("Distance (x)")
    plt.ylim([0,8])
    plt.show()
    return 0

if __name__ == '__main__':
    # figure_1()
    figure_3()