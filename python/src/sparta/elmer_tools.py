from collections import OrderedDict
from ._src import error
import numpy as np
import os, shutil
from sys import float_info

def read_file(_file_path:str) -> list[list[str]]:
    with open(_file_path, 'r') as f:
        return [
            [eval(val) for val in line.strip().split()] for line in f.readlines()
        ]

class Elmer_Grid_to_Sparta_Surf:
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
            

    # def show_nodes(self) -> None:
    #     lines = np.array([[*self.node_data[item][1:]] for item in self.node_data])
    #     mlab.points3d(lines[:, 0], lines[:, 1], lines[:, 2], scale_factor=1.)
    #     mlab.show()


    def make_sparta_surf(self, out_file:str=None) -> str:
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
        
        if out_file is None:
            out_file = f"{self.file_stem}.surf"

        # if os.path.exists(out_file):
        #     error("Output file already exists, " + out_file)
                
        with open(out_file, 'w') as f:
            # print(out_file)
            f.write(
                f"# Surface element file written by SPARTA\n\n{len(self.node_data)} points\n{len(tris)} triangles\n\nPoints\n\n"
            )

            for item in self.node_data:
                f.write(f"{item} " + " ".join([f"{{:.{float_info.dig}f}}".format(float(val)) for val in self.node_data[item][1:]]) + "\n")
                # print(item, self.node_data[item])
                # exit()

            f.write("\nTriangles\n\n")
            
            for item in tris: f.write(f"{item} " + " ".join([str(val) for val in tris[item]]) + "\n")

        return out_file


class Sparta_Surf_Temperature_file(Elmer_Grid_to_Sparta_Surf):
    data_file = ""
    data = np.array([])
    def __init__(self, _data_file:str, _file_stem:str) -> None:
        self.file_stem = _file_stem
        self.data_file = _data_file

        super().__init__(self.file_stem)
        self.load_boundary()
        self.load_nodes()
        self.load_data()


    def load_data(self) -> None:
        with open(self.data_file, 'r') as f:
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

        # perms = np.array(lines[:start_vals])
        self.data = np.array([val[0] for val in lines[start_vals:]])
    
    def make_sparta_surf_data(self, data_file:str=None) -> str:
        tris = OrderedDict()
        temps = OrderedDict()
        if data_file is None:
            data_file = f"{self.file_stem}.surf.data"

        if os.path.exists(data_file):
            error("Data file already exists, " + data_file)
        
        boundary_element_start = 5
        for i in range(len(self.boundary_data)):
            temps[self.boundary_data[i][0]] = np.average([self.data[item - 1] for item in self.boundary_data[i][boundary_element_start:]])
            tris[self.boundary_data[i][0]] = self.boundary_data[i][boundary_element_start:]

        header = "\n".join([
            "ITEM: TIMESTEP", 
            "NAN",
            "ITEM: NUMBER OF SURFS",
            str(len(tris)),
            "ITEM: BOX BOUNDS NAN NAN NAN",
            "-NAN NAN",
            "-NAN NAN",
            "-NAN NAN",
            "ITEM: SURFS id c"
        ])

        with open(data_file, 'w') as f:
            f.write(header + "\n")
            for item in temps:
                f.write(f"{item} {temps[item]}\n")
        
        return data_file


def itemize_boundary_file(_source_file:str, _output_file:str):
    # loading the lines
    with open(_source_file, 'r') as f:
        lines = [[eval(val) for val in line.strip().split()] for line in f.readlines()]

    # setting the boundary the element id to the index
    for i in range(len(lines)): lines[i][1] = lines[i][0]

    # writing back to source file
    with open(_output_file, 'w') as f:
        for line in lines:
            f.write(" ".join([str(item) for item in line]) + "\n")


def itemize_elements_file(_source_file:str, _output_file:str):
    # loading the lines
    with open(_source_file, 'r') as f:
        lines = [[eval(val) for val in line.strip().split()] for line in f.readlines()]

    # setting the boundary the element id to the index
    for i in range(len(lines)): lines[i][1] = lines[i][0]

    # writing back to source file
    with open(_output_file, 'w') as f:
        for line in lines:
            f.write(" ".join([str(item) for item in line]) + "\n")