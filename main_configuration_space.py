# Resources
# https://stackoverflow.com/questions/27270477/3d-convex-hull-from-point-cloud

# Plan
# I want to create an array of cell objects (e.g. 5x5)
# Each cell has a 3D position and orientation (6 numbers)
# x,y,z,theta,phi,psi

import numpy as np
import matplotlib.pyplot as plt

def plot_in_hull(p, hull):
    """
    plot relative to `in_hull` for 2d data
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection

    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    # plot triangulation
    poly = PolyCollection(hull.points[hull.vertices], facecolors='w', edgecolors='b')
    plt.clf()
    plt.title('in hull')
    plt.gca().add_collection(poly)
    plt.plot(hull.points[:,0], hull.points[:,1], 'o', hold=1)


    # plot the convex hull
    edges = set()
    edge_points = []

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(hull.points[ [i, j] ])

    for ia, ib in hull.convex_hull:
        add_edge(ia, ib)

if __name__ == '__main__':
    tested = np.random.rand(20, 3)
    cloud = np.random.rand(50, 3)

    print(in_hull(tested, cloud))