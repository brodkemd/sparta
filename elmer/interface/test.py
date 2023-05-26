import numpy as np

from mayavi import mlab


file = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh/mesh.nodes"

with open(file, 'r') as f:
    lines = [[eval(val) for val in line.strip().split()] for line in f.readlines()]

lines = np.array(lines)

mlab.points3d(lines[:, 2], lines[:, 3], lines[:, 4], scale_factor=1.)

mlab.show()