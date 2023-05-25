import os
import numpy as np
from mayavi import mlab

sparta_surf_file  = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh/out.surf"
out_node_file     = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh_out/mesh.nodes"
out_element_file = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh_out/mesh.elements" 

with open(sparta_surf_file, 'r') as f:
    lines = [line.strip().split() for line in f.readlines()]
    
start = 0
points = 0
triangles = 0

i = -1
while i < len(lines):
    i+=1
    line = lines[i]
    if not len(line): continue

    if line[-1] == "points":
        points = eval(line[0])
    
    if line[-1] == "triangles":
        triangles = eval(line[0])
    
    if line[0] == "Points":
        start = i+2
        break
    
    i+=1


nodes = lines[start:start+points]
triangles = lines[start+points+3:]

with open(out_node_file, 'w') as f:
    for i in range(len(nodes)):
        nodes[i].insert(1, "-1")
        f.write(" ".join(nodes[i]) + "\n")

with open(out_element_file, 'w') as f:
    for i in range(len(triangles)):
        triangles[i].insert(1, "1")
        triangles[i].insert(2, "303")
        f.write(" ".join(triangles[i]) + "\n")
