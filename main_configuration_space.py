# Resources
# https://stackoverflow.com/questions/27270477/3d-convex-hull-from-point-cloud

# Plan
# I want to create an array of cell objects (e.g. 5x5)
# Each cell has a 3D position and orientation (6 numbers)
# x,y,z,theta,phi,psi

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from matplotlib.collections import PolyCollection, LineCollection
from scipy.spatial import Delaunay


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p) >= 0


def plot_in_hull_2D(p, hull):
    """
    Plot `in_hull` for 2d data
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)
    # plot triangulation
    poly = PolyCollection(hull.points[hull.simplices], facecolors='w', edgecolors='b')
    plt.clf()
    plt.title('in hull')
    plt.gca().add_collection(poly)
    plt.plot(hull.points[:, 0], hull.points[:, 1], 'o')
    # plot the convex hull
    edges = set()
    edge_points = []

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(hull.points[[i, j]])

    for ia, ib in hull.convex_hull:
        add_edge(ia, ib)

    lines = LineCollection(edge_points, color='g')
    plt.gca().add_collection(lines)
    plt.show()

    # plot tested points `p` - black are inside hull, red outside
    inside = in_hull(p, hull)
    plt.plot(p[inside, 0], p[inside, 1], '.k')
    plt.plot(p[-inside, 0], p[-inside, 1], '.r')


if __name__ == '__main__':
    # Create a single cell, fixed at the (0, 0, 0) origin point
    cell_center = np.array([0,0,0])
    tested = np.random.rand(20, 3)
    cloud = np.random.rand(50, 3)
    print(cloud[0])
    print(type(cloud[0]))

    print(in_hull(cell_center, cloud))

    # tested_2d = np.random.rand(20, 2)
    # cloud_2d = np.random.rand(50, 2)
    # print(plot_in_hull_2D(tested_2d, cloud_2d))
