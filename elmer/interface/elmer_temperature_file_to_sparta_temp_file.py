from mayavi import mlab
from collections import OrderedDict
import numpy as np
from my import error

from elmer_grid_to_sparta_surf import elmer_to_sparta

def read_file(_file_path:str) -> list[list[str]]:
    with open(_file_path, 'r') as f:
        return [
            [eval(val) for val in line.strip().split()] for line in f.readlines()
        ]


class Sparta_Surf_Temperature_file(elmer_to_sparta):
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
    
    def make_sparta_surf_data(self) -> None:
        tris = OrderedDict()
        temps = OrderedDict()
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

        with open(f"{self.file_stem}.surf.data", 'w') as f:
            f.write(header + "\n")
            for item in temps:
                f.write(f"{item} {temps[item]}\n")

if __name__ == "__main__":
    base = "/home/marekbrodke/Downloads/tutorials-GUI-files/TemperatureGeneric"
    inst = Sparta_Surf_Temperature_file(f"{base}/case.dat", f"{base}/mesh")
    inst.make_sparta_surf_data()
    inst.make_sparta_surf()
# print(inst.boundary_data[:4])


# lines = np.array([[*node_data[item][1:]] for item in node_data])
# plot = mlab.points3d(lines[:, 0], lines[:, 1], lines[:, 2], vals, scale_factor=0.0075)
# mlab.scalarbar(plot)
# mlab.show()




    