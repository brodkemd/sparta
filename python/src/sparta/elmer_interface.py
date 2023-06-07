from ._src import error
from types import FunctionType
import os, sys
from shutil import copy2

cwd:str = os.path.sep.join(__file__.split(os.path.sep)[:-1])

_materials = {
    "Air (room temperature)" : {
        "Density" : 1.205,
        "Heat conductivity" : 0.0257,
        "Heat capacity" :1005.0,
        "Heat expansion coeff." :3.43e-3,
        "Viscosity" :1.983e-5,
        "Turbulent Prandtl number" :0.713,
        "Sound speed" :343.0,
        "Relative Permittivity" :1.00059,
        "Relative Permeability" :1.00000037
    },
    "Water (room temperature)" : {
        "Density" :998.3,
        "Heat conductivity" :0.58,
        "Heat capacity" :4183.0,
        "Heat expansion coeff." :0.207e-3,
        "Viscosity" :1.002e-3,
        "Turbulent Prandtl number" :7.01,
        "Sound speed" :1497.0,
        "Relative Permittivity" :80.1,
        "Relative Permeability" :0.999992,
    },
    "Glycerol (room temperature)" : {
        "Density" :1261.0,
        "Heat conductivity" :0.28,
        "Heat capacity" :2400.0,
        "Heat expansion coeff." :0.610e-3,
        "Viscosity" :1.49,
        "Sound speed" :1904.0,
        "Relative Permittivity" :42.5
    },
    "Ethanol (room temperature)" : {
        "Density" :789.0,
        "Heat conductivity" :0.171,
        "Heat capacity" :2440.0,
        "Heat expansion coeff." :0.25e-3,
        "Viscosity" :1.074e-3,
        "Sound speed" :1144.0,
        "Relative Permittivity" :24.3
    },
    "Oil, olive (25 C)" : {
        "Density" :915.0,
        "Heat conductivity" :0.17,
        "Heat capacity" :1970.0,
        "Heat expansion coeff." :0.72e-3,
        "Viscosity" :0.081,
        "Sound speed" :1430.0,
        "Relative Permittivity" :3.1
    },
    "Aluminium" : {
        "Density" :2700.0,
        "Youngs modulus" :70.0e9,
        "Poisson ratio" :0.35,
        "Shear modulus" :26.0e9,
        "Bulk modulus" :76.0e9,
        "Heat expansion coeff." :23.1e-6,
        "Heat conductivity" :237.0,
        "Heat capacity" :897.0,
        "Electric resistivity" :26.50e-9,
        "Electric conductivity" :37.73e6,
        "Sound speed" :5000.0,
        "Relative Permeability" :1.000022
    },
    "Carbon Steel" : {
        "Density" :7850.0,
        "Youngs modulus" :200.0e9,
        "Poisson ratio" :0.285,
        "Heat expansion coeff." :13.8e-6,
        "Heat conductivity" :44.8,
        "Heat capacity" :1265.0,
        "Electric resistivity" :690.0e-9,
        "Electric conductivity" :1.449e6,
        "Tensile strength" :1079.0e6,
        "Yield strength" :472.0e6,
        "Sound speed" :5100.0
    },
    "Alloy Steel" : {
        "Density" :7850.0,
        "Youngs modulus" :200.0e9,
        "Poisson ratio" :0.285,
        "Heat expansion coeff." :12.0e-6,
        "Heat conductivity" :37.2,
        "Heat capacity" :976.0,
        "Electric resistivity" :731.0e-9,
        "Electric conductivity" :1.367e6,
        "Tensile strength" :1320.0e6,
        "Yield strength" :1080.0e6,
        "Sound speed" :5100.0
    },
    "Stainless Steel" : {
        "Density" :7925.0,
        "Youngs modulus" :200.0e9,
        "Poisson ratio" :0.285,
        "Heat expansion coeff." :14.9e-6,
        "Heat conductivity" :24.0,
        "Heat capacity" :460.0,
        "Electric resistivity" :548.0e-9,
        "Electric conductivity" :1.824e6,
        "Tensile strength" :671.0e6,
        "Yield strength" :380.0e6,
        "Sound speed" :5100.0
    },
    "Austenitic stainless steel" : {
        "Density" :7810.0,
        "Youngs modulus" :197.0e9,
        "Poisson ratio" :0.3,
        "Heat conductivity" :16.2,
        "Heat capacity" :500.0,
        "Heat expansion coeff." :15.7e-6,
        "Electric resistivity" :6.85e-7,
        "Electric conductivity" :14.60e3,
        "Magnetic permeability" :1.02,
        "Tensile strength" :379.0e6
    },
    "Copper" : {
        "Density" :8960.0,
        "Youngs modulus" :115.0e9,
        "Poisson ratio" :0.34,
        "Shear modulus" :48.0e9,
        "Bulk modulus" :140.0e9,
        "Heat expansion coeff." :16.5e-6,
        "Heat conductivity" :401.0,
        "Heat capacity" :385.0,
        "Electric resistivity" :16.78e-9,
        "Electric conductivity" :59.59e6,
        "Sound speed" :3810.0,
        "Relative Permeability" :0.999994
    },
    "Lead" : {
        "Density" :11340.0,
        "Youngs modulus" :16.0e9,
        "Poisson ratio" :0.44,
        "Shear modulus" :5.6e9,
        "Bulk modulus" :46.0e9,
        "Heat expansion coeff." :28.9e-6,
        "Heat conductivity" :35.3,
        "Heat capacity" :128.0,
        "Electric resistivity" :208.0e-9,
        "Electric conductivity" :4.808e6,
        "Sound speed" :1190.0
    },
    "Glass" : {
        "Density" :2235.0,
        "Heat expansion coeff." :3.5e-6,
        "Poisson ratio" :0.15,
        "Youngs modulus" :65.0e9,
        "Shear modulus" :28.2e9,
        "Heat conductivity" :1.14,
        "Heat capacity" :710.0
    },
    "Polycarbonate" : {
        "Density" :1220.0,
        "Youngs modulus" :2.2e9,
        "Poisson ratio" :0.37,
        "Tensile strength" :60.0e6,
        "Heat expansion coeff." :67.0e-6,
        "Heat capacity" :1250.0,
        "Heat conductivity" :0.205
    },
    "Polyvinyl Chloride" : {
        "Density" :1380.0,
        "Youngs modulus" :3100.0e6,
        "Poisson ratio" :0.41,
        "Tensile strength" :65.0e6,
        "Heat expansion coeff." :80.0e-6,
        "Heat capacity" :900.0,
        "Heat conductivity" :0.16
    },
    "Gold" : {
        "Density" :19300.0,
        "Youngs modulus" :78.0e9,
        "Poisson ratio" :0.44,
        "Shear modulus" :27.0e9,
        "Bulk modulus" :180.0e9,
        "Heat expansion coeff." :14.2e-6,
        "Heat conductivity" :318.0,
        "Heat capacity" :129.0,
        "Electric resistivity" :22.14e-9,
        "Electric conductivity" :45.17e6,
        "Sound speed" :2030.0
    },
    "Silver" : {
        "Density" :10490.0,
        "Youngs modulus" :83.0e9,
        "Poisson ratio" :0.37,
        "Shear modulus" :30.0e9,
        "Bulk modulus" :100.0e9,
        "Heat expansion coeff." :18.9e-6,
        "Heat conductivity" :429.0,
        "Heat capacity" :235.0,
        "Electric resistivity" :15.87e-9,
        "Electric conductivity" :63.01e6,
        "Sound speed" :2680.0
    },
    "Platinum" : {
        "Density" :21450.0,
        "Youngs modulus" :144.79e9,
        "Poisson ratio" :0.38    ,
        "Heat expansion coeff." :8.8e-6,
        "Heat conductivity" :71.6,
        "Heat capacity" :133.0,
        "Electric resistivity" :105.0e-9,
        "Electric conductivity" :9.523e6,
        "Sound speed" :2800.0,
    },
    "Iron" : {
        "Density" :7870.0,
        "Youngs modulus" :193.053e9    ,
        "Poisson ratio" :0.29    ,
        "Heat expansion coeff." :11.8e-6,
        "Heat conductivity" :80.2,
        "Heat capacity" :449.0,
        "Electric resistivity" :97.1e-9,
        "Electric conductivity" :10.30e6,
        "Sound speed" :5000.0
    },
    "Water" : {
        "Density" :910.0,
        "Heat conductivity" :"Variable Temperature; Real MATC \"9.828*exp(-5.7E-03*(tx+273.15))\"",
        "Heat capacity" :"Variable Temperature; Real MATC \"146.3+(7.253*(tx+273.15))\"",
        "Viscosity Model" :"Power law",
        "Viscosity Exponent" :"$ (1.0/3.0)",
        "Critical Shear Rate" :1.0E-6,
        "Viscosity" :"Variable Temperature; Real MATC \"(2.0*3.0*1.916E03 * exp( -139.0E03/(8.314 *(tx+273.15))))^(-1.0/3.0)\""
    },
    "Solid Silicon" : {
        "Density" :2330.0,
        "Youngs modulus" :185.0e9    ,
        "Poisson ratio" :0.28    ,
        "Heat expansion coeff." :4.68e-6,
        "Heat conductivity" :"Variable Temperature; Real; 0 156; 300 156; 550 72; 800 43; 1050 29; 1300 25; 1550 23; 1800 21; 3000 21; End",
        "Heat capacity" :555.8,
        "Emissivity" :0.7,
        "Electric resistivity" :1.0e3,
        "Electric conductivity" :1.0e-3,
        "Sound speed" :8433.0,
        "Melting Point" :1683.0,
        "Latent heat" :1.8e6
    },
    "Liquid Silicon" : {
        "Liquid" :"Logical True",
        "Viscosity" :8.0e-4,
        "Density" :2570,
        "Heat conductivity" :50.0,
        "Emissivity" :0.3,
        "Heat expansion coeff." :1.08e-4
    },
    "Fused Silica (25 C)" : {
        "Density" :2200.0,
        "Youngs modulus" :72.0e9    ,
        "Poisson ratio" :0.17    ,
        "Heat expansion coeff." :5.4e-7,
        "Heat conductivity" :1.46,
        "Heat capacity" :670.0,
        "Sound speed" :5900.0,
        "Relative Permittivity" :3.75,
        "Melting Point" :1956.0
    },
    "Liquid CO2": {
        "Liquid" :"Logical True",
        "Density" :"Variable Temperature; Real MATC \"44.01*2.768/(0.26212^(1+ (1-tx/304.21)^0.2908))\"",
        "Viscosity" :"Variable Temperature; Real MATC \"1.E-3*exp(-3.097 + 48.86/tx + 0.02381*tx -0.0000784*tx^2)\"",
        "Heat conductivity" :"Variable Temperature; Real MATC \"0.407 - 0.0008438*tx -9.626E-07*tx^2\"",
        "Heat capacity" :"Variable Temperature; Real MATC \"1E3*(-3553.844 + 46.88128*tx - 0.2017221*tx^2 + 2.897028E-04*tx^3)/44\""
    },
    "Gas CO2": {
        "Liquid" :"Logical False",
        "Density" :1.98,
        "Viscosity" :"Variable Temperature; Real MATC \"1.E-3*(0.002545 + 4.549E-05*tx - 8.649E-09*tx^2)\"",
        "Heat conductivity" :"Variable Temperature; Real MATC \"-0.007215 + 8.015001E-05*tx + 5.477E-09*tx^2 - 1.053E-11*tx^3\"",
        "Heat capacity" :"Variable Temperature; Real MATC \"1E3*(2.926801E+01 - 2.236208E-02*tx + 2.652535E-04*tx^2 - 4.153087E-07*tx^3 + 2.005667E-10*tx^4)/44\""
    }
}


class Elmer:
    sif_file:str             = ""
    exe:str                  = None
    tab:str                  = 2*" "
    commands:list[list[str]] = []


    def __init__(self, output_file:str, elmer_exe:str=None) -> None:
        self.sif_file = output_file
        self.exe      = elmer_exe


    def run(self):
        if self.exe is not None:
            self.write()
            os.system(f"{self.exe} {self.sif_file}")
        else: error("Can not run without setting elmer_exe in constructor (init method)")


    def write(self):
        with open(self.sif_file, 'w') as f:
            for section in self.commands:
                f.write(f"{section[0]}\n")
                for i in range(1, len(section)-1):
                    f.write(self.tab + section[i] + "\n")
                f.write(f"{section[-1]}\n\n")


    def _add_section(self, _name:str, _n:str|int, args:list[str]) -> None:
        self.commands.append([])
        _n = str(_n)
        if len(_n):
            self.commands[-1].append(f"{_name} {_n}")
        else:
            self.commands[-1].append(f"{_name}")

        if isinstance(args, list):
            for arg in args: self.commands[-1].append(arg)
        elif isinstance(args, str):
            self.commands[-1].append(args)
        else:
            raise TypeError("invalid input type for args")
        
        self.commands[-1].append("End")


    def Header(self, args:list[str]) -> None:
        self._add_section("Header", "", args)


    def default_Header(self) -> None:
        self.Header([
            "Mesh DB \".\" \".\"",
            "Include Path \"\"",
            "Results Directory \"\""
        ])


    def Simulation(self, args:list[str]) -> None:
        self._add_section("Simulation", "", args)


    def Constants(self, args:list[str]) -> None:
        self._add_section("Constants", "", args)


    def default_Constants(self) -> None:
        self.Constants([
            "Gravity(4) = 0 -1 0 9.82",
            "Stefan Boltzmann = 5.670374419e-08",
            "Permittivity of Vacuum = 8.85418781e-12",
            "Permeability of Vacuum = 1.25663706e-6",
            "Boltzmann Constant = 1.380649e-23",
            "Unit Charge = 1.6021766e-19"
        ])


    def Body(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Body", _n, args)


    def Material(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Material", _n, args)


    # if material name is none, return available materials
    def MaterialDB(self, _n:int|str=None, material_name:str=None) -> list[str]|None:
        if material_name is None:
            return list(_materials.keys())
        else:
            
            _material_name = material_name.lower()
            for item in list(_materials.keys()):
                new_item = item.lower()
                _materials[new_item] = _materials[item]
                del _materials[item]
                
            if _material_name in _materials:
                _materials_list = []
                _materials[_material_name]["Name"] = f"\"{material_name}\""
                for item in _materials[_material_name]:
                    if not isinstance(_materials[_material_name][item], str):
                        _materials[_material_name][item] = f"Real {{:.{sys.float_info[6]}f}}".format(_materials[_material_name][item])
                
                for item in _materials[_material_name]:
                    _materials_list.append(f"{item} = {_materials[_material_name][item]}")
                
                self.Material(_n, _materials_list)
            else:
                error("no material named: " + material_name)


    def Body_Force(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Body Force", _n, args)


    def Equation(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Equation", _n, args)


    def Solver(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Solver", _n, args)


    def Boundary_Condition(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Boundary Condition", _n, args)


    def Initial_Condition(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Initial Condition", _n, args)


    def Component(self, _n:str|int, args:list[str]) -> None:
        self._add_section("Component", _n, args)


"""
makes a python file from a sif file

"""
def py_from_sif(sif:str, exe:str=None, backup_sif:bool=False):
    # getting the python file path from the sif file
    py  = sif[:sif.rfind(".")+1] + "py"

    if backup_sif: copy2(sif, sif+".backup")

    print("Making:", py, "from", sif)

    # adding the required header to the code
    data = [ "from sparta.elmer_interface import Elmer, man\n" ]
    if exe is None: data.append(f"e = Elmer(\"{sif}\")\n")
    else:           data.append(f"e = Elmer(\"{sif}\", \"{exe}\")\n")

    # loading the lines from the sif file, ignoring empty lines
    lines = []
    with open(sif, 'r') as f: lines = [line for line in f.readlines() if len(line.strip())]

    # parsing lines in file
    for line in lines:
        # if the line starts with a space then is a section command
        if not line.startswith(" "):
            split_line = line.strip().split()
            if len(split_line) > 1:
                while "" in split_line: split_line.remove("")
                
                func = split_line[:-1]
                id = split_line[-1]

                for j in range(len(func)):
                    func[j] = func[j][0].upper() + func[j][1:].lower()
                
                data.append("e." + "_".join(func) + "(" + id + ", [")

            else:
                if split_line[0].lower() == "end":
                    data.append("])\n")
                else:
                    data.append("e."+split_line[0][0].upper() + split_line[0][1:].lower() + "([")

        # if here then it is contents of the command
        else:
            # adding the line to the new file
            line = line.replace("\"", "\\\"")
            data.append(f"  \"{line.strip()}\",")
    
    if exe is not None:
        data.append("e.run()\n")
    else:
        data.append("e.write()\n")
    
    # writing python file
    with open(py, "w") as f:
        for line in data: f.write(line + "\n")


"""
Pass the function to get info on it in the browser

"""
def man(func:FunctionType=None, exit_after_call:bool=True, browser:str = "xdg-open") -> None:
    if func is None:
        _doc_path = os.path.join(cwd, "elmer_doc", "ElmerSolverManual.pdf")
        print(f"opening {_doc_path} with {browser}")
        os.system(f"{browser} {_doc_path} &")
        if exit_after_call:
            print("Opening in browser and exiting")
            exit()
    else:
        error("Function method not implement in this function yet")
        # if func.__name__ in _docs:
        #     print(f"opening {_docs[func.__name__]} with {browser}")
        #     os.system(f"{browser} {_docs[func.__name__]} &")
        #     if exit_after_call:
        #         print("Opening in browser and exiting")
        #         exit()
        # else:
        #     error(f"{func.__name__} is not a valid function.")