import numpy as np
from collections import OrderedDict

def read_file(_file_path:str) -> list[list[str]]:
    with open(_file_path, 'r') as f:
        return [
            [eval(val) for val in line.strip().split()] for line in f.readlines()
        ]

file = "/home/marekbrodke/Downloads/tutorials-GUI-files/TemperatureGeneric/case.dat"
with open(file, 'r') as f:
    lines = [[val for val in line.strip().split()] for line in f.readlines()]
    
    for i in range(len(lines)):
        if lines[i][0] == "Perm:":
            lines = [[eval(val) for val in line] for line in lines[i+1:]]
            break

    start_vals = 0
    for i, item in enumerate(lines):
        if len(item) == 1:
            start_vals = i
            break

perms = np.array(lines[:start_vals])
vals = np.array([val[0] for val in lines[start_vals:]])

_lines = read_file("/home/marekbrodke/Downloads/tutorials-GUI-files/TemperatureGeneric/mesh.nodes")
node_data = OrderedDict()
for line in _lines: node_data[line[0]] = line[1:]

lines = np.array([[*node_data[item][1:]] for item in node_data])