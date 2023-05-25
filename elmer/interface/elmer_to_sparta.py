import os
import numpy as np

#
# use 303 for type in elements file
#

file_exts = [
    "boundary",
    "elements",
    "header",
    "nodes"
]

file_stem = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh/mesh"

cur_max = 0
for ext in ["elements"]:
    file = f"{file_stem}.{ext}"
    
    lines = []
    with open(file, 'r') as f:
        lines = np.array([
            [eval(val) for val in line.strip().split()] for line in f.readlines()
        ])
    
    for i in range(4, lines.shape[1]):
        new_max = max(lines[:, i])
        if new_max > cur_max: cur_max = new_max


    