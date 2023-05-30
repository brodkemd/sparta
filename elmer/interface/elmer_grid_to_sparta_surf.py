import os
from collections import OrderedDict
import numpy as np
from my import error

from mayavi import mlab

#
# use 303 for type in elements file
#

def read_file(_file_path:str) -> list[list[str]]:
    with open(_file_path, 'r') as f:
        return [
            [eval(val) for val in line.strip().split()] for line in f.readlines()
        ]

class elmer_to_sparta:
    file_stem = ""
    boundary_data = []
    node_data = OrderedDict()

    def __init__(self, _file_stem:str) -> None:
        self.file_stem = _file_stem
        
        self.load_boundary()
        self.load_nodes()


    def load_boundary(self) -> None:
        _lines = read_file(f"{self.file_stem}.boundary")
        
        for i, line in enumerate(_lines):
            if line[4] != 303:
                error(f"line {i} in boundary file, boundary elements must be triangles")

        self.boundary_data = _lines


    def load_nodes(self) -> None:
        _lines = read_file(f"{self.file_stem}.nodes")
        
        for line in _lines: self.node_data[line[0]] = line[1:]


    def update_boundary_file(self) -> None:
        with open(f"{self.file_stem}.boundary", 'w') as f:
            for i, line in enumerate(self.boundary_data):
                line = 2*[line[0]] + line[2:]
                f.write(" ".join([str(val) for val in line]) + "\n")
            

    def show_nodes(self) -> None:
        lines = np.array([[*self.node_data[item][1:]] for item in self.node_data])
        mlab.points3d(lines[:, 0], lines[:, 1], lines[:, 2], scale_factor=1.)
        mlab.show()


    def make_sparta_surf(self):
        tris = OrderedDict()
        boundary_element_ids = OrderedDict()
        boundary_element_start = 5
        for i in range(len(self.boundary_data)):
            tris[self.boundary_data[i][0]] = self.boundary_data[i][boundary_element_start:]

            for j, item in enumerate(self.boundary_data[i][boundary_element_start:]):
                if item not in boundary_element_ids:
                    boundary_element_ids[item] = []
                
                boundary_element_ids[item].append((i, j+boundary_element_start))
        
        boundary_element_ids = OrderedDict(sorted(boundary_element_ids.items(), key=lambda x: x))
        boundary_element_id_list = list(boundary_element_ids.keys())
        for i in range(len(boundary_element_id_list)-1):
            if boundary_element_id_list[i+1] - boundary_element_id_list[i] != 1:
                error(f"Got: {i+1} {i}")
        
        # with open("out.txt", 'w') as f:
        #     for item in list(boundary_element_id_list):
        #         f.write(f"{item}\n")
                
        with open(f"{self.file_stem}.surf", 'w') as f:
            print(f"{self.file_stem}.surf")            
            f.write(
                f"# Surface element file written by SPARTA\n\n{len(boundary_element_ids)} points\n{len(tris)} triangles\n\nPoints\n\n"
            )

            for item in self.node_data:
                f.write(f"{item} " + " ".join([str(val) for val in self.node_data[item][1:]]) + "\n")
                # print(item, self.node_data[item])
                # exit()

            f.write("\nTriangles\n\n")
            
            for item in tris: f.write(f"{item} " + " ".join([str(val) for val in tris[item]]) + "\n")

if __name__ == "__main__":
    file_stem = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh/mesh"

    inst = elmer_to_sparta(file_stem)
    # inst.update_boundary_file()
    # inst.show_nodes()
    inst.make_sparta_surf()

    