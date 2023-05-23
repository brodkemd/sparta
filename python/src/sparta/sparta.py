from types import FunctionType
import os
from collections import OrderedDict
from .molecule import Molecule, Species, VSS
from ._src import error

cwd:str = os.path.sep.join(__file__.split(os.path.sep)[:-1])


_docs:dict[str, str] = {
    "Adapt_grid" : os.path.join(cwd, "doc", "adapt_grid.html"),
    "Balance_grid" : os.path.join(cwd, "doc", "balance_grid.html"),
    "Boundary" : os.path.join(cwd, "doc", "boundary.html"),
    "Bound_modify" : os.path.join(cwd, "doc", "bound_modify.html"),
    "Clear" : os.path.join(cwd, "doc", "clear.html"),
    "Collide" : os.path.join(cwd, "doc", "collide.html"),
    "Collide_modify" : os.path.join(cwd, "doc", "collide_modify.html"),
    "Compute" : os.path.join(cwd, "doc", "compute.html"),
    "Create_box" : os.path.join(cwd, "doc", "create_box.html"),
    "Create_grid" : os.path.join(cwd, "doc", "create_grid.html"),
    "Create_particles" : os.path.join(cwd, "doc", "create_particles.html"),
    "Dimension" : os.path.join(cwd, "doc", "dimension.html"),
    "Dump" : os.path.join(cwd, "doc", "dump.html"),
    "Dump_modify" : os.path.join(cwd, "doc", "dump_modify.html"),
    "Echo" : os.path.join(cwd, "doc", "echo.html"),
    "Fix" : os.path.join(cwd, "doc", "fix.html"),
    "Global" : os.path.join(cwd, "doc", "global.html"),
    "Group" : os.path.join(cwd, "doc", "group.html"),
    "If" : os.path.join(cwd, "doc", "if.html"),
    "Include" : os.path.join(cwd, "doc", "include.html"),
    "Jump" : os.path.join(cwd, "doc", "jump.html"),
    "Label" : os.path.join(cwd, "doc", "label.html"),
    "Log" : os.path.join(cwd, "doc", "log.html"),
    "Mixture" : os.path.join(cwd, "doc", "mixture.html"),
    "Move_surf" : os.path.join(cwd, "doc", "move_surf.html"),
    "Next" : os.path.join(cwd, "doc", "next.html"),
    "Package" : os.path.join(cwd, "doc", "package.html"),
    "Partition" : os.path.join(cwd, "doc", "partition.html"),
    "Print" : os.path.join(cwd, "doc", "print.html"),
    "Quit" : os.path.join(cwd, "doc", "quit.html"),
    "React" : os.path.join(cwd, "doc", "react.html"),
    "React_modify" : os.path.join(cwd, "doc", "react_modify.html"),
    "Read_grid" : os.path.join(cwd, "doc", "read_grid.html"),
    "Read_isurf" : os.path.join(cwd, "doc", "read_isurf.html"),
    "Read_particles" : os.path.join(cwd, "doc", "read_particles.html"),
    "Read_restart" : os.path.join(cwd, "doc", "read_restart.html"),
    "Read_surf" : os.path.join(cwd, "doc", "read_surf.html"),
    "Region" : os.path.join(cwd, "doc", "region.html"),
    "Remove_surf" : os.path.join(cwd, "doc", "remove_surf.html"),
    "Reset_timestep" : os.path.join(cwd, "doc", "reset_timestep.html"),
    "Restart" : os.path.join(cwd, "doc", "restart.html"),
    "Run" : os.path.join(cwd, "doc", "run.html"),
    "Scale_particles" : os.path.join(cwd, "doc", "scale_particles.html"),
    "Seed" : os.path.join(cwd, "doc", "seed.html"),
    "Shell" : os.path.join(cwd, "doc", "shell.html"),
    "Species" : os.path.join(cwd, "doc", "species.html"),
    "Species_modify" : os.path.join(cwd, "doc", "species_modify.html"),
    "Stats" : os.path.join(cwd, "doc", "stats.html"),
    "Stats_modify" : os.path.join(cwd, "doc", "stats_modify.html"),
    "Stats_style" : os.path.join(cwd, "doc", "stats_style.html"),
    "Suffix" : os.path.join(cwd, "doc", "suffix.html"),
    "Surf_collide" : os.path.join(cwd, "doc", "surf_collide.html"),
    "Surf_react" : os.path.join(cwd, "doc", "surf_react.html"),
    "Surf_modify" : os.path.join(cwd, "doc", "surf_modify.html"),
    "Timestep" : os.path.join(cwd, "doc", "timestep.html"),
    "Uncompute" : os.path.join(cwd, "doc", "uncompute.html"),
    "Undump" : os.path.join(cwd, "doc", "undump.html"),
    "Unfix" : os.path.join(cwd, "doc", "unfix.html"),
    "Units" : os.path.join(cwd, "doc", "units.html"),
    "Variable" : os.path.join(cwd, "doc", "variable.html"),
    "Write_grid" : os.path.join(cwd, "doc", "write_grid.html"),
    "Write_isurf" : os.path.join(cwd, "doc", "write_isurf.html"),
    "Write_restart" : os.path.join(cwd, "doc", "write_restart.html"),
    "Write_surf" : os.path.join(cwd, "doc", "write_surf.html")
}

"""
Pass the function to get info on it in the browser

"""
def man(func:FunctionType=None, exit_after_call:bool=True, browser:str = "google-chrome") -> None:
    if func is None:
        _doc_path = os.path.join(cwd, "doc", "Manual.html")
        print(f"opening {_doc_path} with {browser}")
        os.system(f"{browser} {_doc_path} &")
        if exit_after_call:
            print("Opening in browser and exiting")
            exit()
    else:
        if func.__name__ in _docs:
            print(f"opening {_docs[func.__name__]} with {browser}")
            os.system(f"{browser} {_docs[func.__name__]} &")
            if exit_after_call:
                print("Opening in browser and exiting")
                exit()
        else:
            error(f"{func.__name__} is not a valid function.")

class Sparta:
    commands:list = []
    _docs:dict[str, str] = {}

    def __init__(self,
        exe_path:str, # path to sparta exe
        command_file:str, # file to unload the commands to
        species_file:str = "", # file to unload the species to
        vss_file:str = "", # file to unload the vss data to
    ) -> None:
        super().__init__()

        if not len(exe_path): raise ValueError("exe_path must have nonzero length")
        if not len(command_file): raise ValueError("command_file must have nonzero length")

        self.exe_path:str = exe_path # path to sparta binary
        self.command_file:str =  command_file  # file to unload the commands to
        self.species_file:str = species_file # file to unload the species to
        self.vss_file:str =  vss_file # file to unload the vss data to
        
        # class variables
        self.commands = []
        self.species = []
        self.colors = OrderedDict()
        self.vss = []
        self.database:dict[str, Molecule] = {}

        self.db_file = os.path.join(cwd, "db.txt")
        with open(self.db_file, 'r') as f:
            lines = [line.strip().split(',') for line in f.readlines()]
        
        for line in lines:
            for i in range(len(line)):
                line[i] = eval(line[i].strip())
            self.database[line[0]] = Molecule(*line)
    
    
    def add_molecule_from_db(self, id:str, color:str="") -> None:
        if id not in self.database: raise IndexError(f"\"{id}\" is not available in the database of molecules")
        else:
            self.add_species(*self.database[id].species(), color=color)
            self.add_vss(*self.database[id].vss())


    """
    adds newline to file

    """
    def new_line(self) -> None: self.commands.append([])

    """
    Adds a comment to file

    """
    def comment(self, comment:str) -> None: self.commands.append([f"# {comment}"])

    '''
    
    '''
    def execute(self, 
        args:list[str]|str="", # command line args to sparta
        run_simulation:bool = True # whether or not to run the sparta exe
    ) -> None:
        if isinstance(args, list): args = " ".join(args)
        elif not isinstance(args, str): raise TypeError("args must list of strings or string")

        if len(self.species_file):
            print(f"saving species to  : \"{self.species_file}\"")
            with open(self.species_file, 'w') as f:
                for item in self.species:
                    f.write(str(item) + "\n")
        
        if len(self.vss_file):
            print(f"saving vss to      : \"{self.vss_file}\"")
            with open(self.vss_file, 'w') as f:
                for item in self.vss:
                    f.write(str(item) + "\n")
        
        with open(self.command_file, 'w') as f:
            print(f"saving commands to : \"{self.command_file}\"")
            for item in self.commands:
                f.write(" ".join(str(element) for element in item) + "\n")
        
        cmd = f"{self.exe_path} {args} < {self.command_file}"
        
        if run_simulation:
            print(f"running \"{cmd}\"")
            if os.system(cmd):
                print("\nError: Running sparta, check its output")
                exit()
        else:
            print(f"would have run \"{cmd}\"")


    def add_species(self,
        id:str, # exp. for oxygen this is O2 
        molweight:str|float, # molecular weight in amu (atomic mass units, e.g. 16 for oxygen)
        molmass:str|float, # molecular mass (mass units)
        rotational_dof:str|int, # rotational degrees of freedom (integer, unitless)
        inv_rotational_relaxation:str|float, # inverse rotational relaxtion number (unitless)
        vibrational_dof:str|int, # vibrational degrees of freedom (integer, unitless)
        inv_vibrational_relaxation:str|float, # inverse vibrational relaxation number
        vibrational_temp:str|float, # vibrational temperature (temperature units)
        species_wt:str|float, # species weight (unitless)
        charge:str|float, # multiple of electon charge (1 for a proton)
        color:str = "" # color to add to the molecule
    ):
        self.species.append(Species(id, molweight, molmass, rotational_dof, inv_rotational_relaxation, vibrational_dof, inv_vibrational_relaxation, vibrational_temp, species_wt, charge))
        
        if len(color): self.colors[id] = color


    def add_vss(self,
            id:str, # exp. for oxygen this is O2
            diameter:str|float, #  VHS or VSS diameter of particle (distance units)
            omega:str|float, # temperature-dependence of viscosity (unitless)
            tref:str|float, # reference temperature (temperature units)
            alpha:str|float # angular scattering parameter (unitless)
        ):           
            self.vss.append(VSS(id, diameter, omega, tref, alpha))
    
    
    def set_Species(self):
        self.Species([self.species_file, " ".join([item.id for item in self.species])])
    
    
    """
    run  
        'class_instance.man(class_instance.Adapt_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Adapt_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["adapt_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["adapt_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Balance_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Balance_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["balance_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["balance_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Boundary)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Boundary(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["boundary"]+args)
        elif isinstance(args, str):
            self.commands.append(["boundary "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Bound_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Bound_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["bound_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["bound_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Clear)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Clear(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["clear"]+args)
        elif isinstance(args, str):
            self.commands.append(["clear "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Collide)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Collide(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["collide"]+args)
        elif isinstance(args, str):
            self.commands.append(["collide "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Collide_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Collide_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["collide_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["collide_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Compute)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Compute(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["compute"]+args)
        elif isinstance(args, str):
            self.commands.append(["compute "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_box)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_box(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_box"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_box "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dimension)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dimension(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dimension"]+args)
        elif isinstance(args, str):
            self.commands.append(["dimension "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dump"]+args)
        elif isinstance(args, str):
            self.commands.append(["dump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dump_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dump_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dump_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["dump_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Echo)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Echo(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["echo"]+args)
        elif isinstance(args, str):
            self.commands.append(["echo "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Fix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Fix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["fix"]+args)
        elif isinstance(args, str):
            self.commands.append(["fix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Global)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Global(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["global"]+args)
        elif isinstance(args, str):
            self.commands.append(["global "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Group)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Group(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["group"]+args)
        elif isinstance(args, str):
            self.commands.append(["group "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.If)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def If(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["if"]+args)
        elif isinstance(args, str):
            self.commands.append(["if "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Include)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Include(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["include"]+args)
        elif isinstance(args, str):
            self.commands.append(["include "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Jump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Jump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["jump"]+args)
        elif isinstance(args, str):
            self.commands.append(["jump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Label)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Label(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["label"]+args)
        elif isinstance(args, str):
            self.commands.append(["label "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Log)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Log(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["log"]+args)
        elif isinstance(args, str):
            self.commands.append(["log "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Mixture)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Mixture(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["mixture"]+args)
        elif isinstance(args, str):
            self.commands.append(["mixture "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Move_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Move_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["move_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["move_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Next)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Next(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["next"]+args)
        elif isinstance(args, str):
            self.commands.append(["next "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Package)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Package(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["package"]+args)
        elif isinstance(args, str):
            self.commands.append(["package "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Partition)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Partition(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["partition"]+args)
        elif isinstance(args, str):
            self.commands.append(["partition "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Print)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Print(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["print"]+args)
        elif isinstance(args, str):
            self.commands.append(["print "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Quit)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Quit(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["quit"]+args)
        elif isinstance(args, str):
            self.commands.append(["quit "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.React)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def React(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["react"]+args)
        elif isinstance(args, str):
            self.commands.append(["react "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.React_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def React_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["react_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["react_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_isurf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_isurf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_isurf"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_isurf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Region)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Region(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["region"]+args)
        elif isinstance(args, str):
            self.commands.append(["region "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Remove_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Remove_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["remove_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["remove_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Reset_timestep)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Reset_timestep(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["reset_timestep"]+args)
        elif isinstance(args, str):
            self.commands.append(["reset_timestep "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Run)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Run(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["run"]+args)
        elif isinstance(args, str):
            self.commands.append(["run "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Scale_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Scale_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["scale_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["scale_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Seed)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Seed(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["seed"]+args)
        elif isinstance(args, str):
            self.commands.append(["seed "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Shell)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Shell(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["shell"]+args)
        elif isinstance(args, str):
            self.commands.append(["shell "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Species)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Species(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["species"]+args)
        elif isinstance(args, str):
            self.commands.append(["species "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Species_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Species_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["species_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["species_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats_style)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats_style(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats_style"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats_style "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Suffix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Suffix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["suffix"]+args)
        elif isinstance(args, str):
            self.commands.append(["suffix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_collide)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_collide(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_collide"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_collide "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_react)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_react(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_react"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_react "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Timestep)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Timestep(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["timestep"]+args)
        elif isinstance(args, str):
            self.commands.append(["timestep "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Uncompute)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Uncompute(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["uncompute"]+args)
        elif isinstance(args, str):
            self.commands.append(["uncompute "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Undump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Undump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["undump"]+args)
        elif isinstance(args, str):
            self.commands.append(["undump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Unfix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Unfix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["unfix"]+args)
        elif isinstance(args, str):
            self.commands.append(["unfix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Units)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Units(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["units"]+args)
        elif isinstance(args, str):
            self.commands.append(["units "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Variable)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Variable(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["variable"]+args)
        elif isinstance(args, str):
            self.commands.append(["variable "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_isurf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_isurf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_isurf"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_isurf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Adapt_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Adapt_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["adapt_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["adapt_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Balance_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Balance_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["balance_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["balance_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Boundary)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Boundary(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["boundary"]+args)
        elif isinstance(args, str):
            self.commands.append(["boundary "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Bound_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Bound_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["bound_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["bound_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Clear)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Clear(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["clear"]+args)
        elif isinstance(args, str):
            self.commands.append(["clear "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Collide)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Collide(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["collide"]+args)
        elif isinstance(args, str):
            self.commands.append(["collide "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Collide_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Collide_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["collide_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["collide_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Compute)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Compute(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["compute"]+args)
        elif isinstance(args, str):
            self.commands.append(["compute "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_box)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_box(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_box"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_box "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Create_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Create_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["create_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["create_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dimension)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dimension(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dimension"]+args)
        elif isinstance(args, str):
            self.commands.append(["dimension "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dump"]+args)
        elif isinstance(args, str):
            self.commands.append(["dump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Dump_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Dump_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["dump_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["dump_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Echo)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Echo(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["echo"]+args)
        elif isinstance(args, str):
            self.commands.append(["echo "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Fix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Fix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["fix"]+args)
        elif isinstance(args, str):
            self.commands.append(["fix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Global)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Global(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["global"]+args)
        elif isinstance(args, str):
            self.commands.append(["global "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Group)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Group(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["group"]+args)
        elif isinstance(args, str):
            self.commands.append(["group "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.If)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def If(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["if"]+args)
        elif isinstance(args, str):
            self.commands.append(["if "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Include)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Include(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["include"]+args)
        elif isinstance(args, str):
            self.commands.append(["include "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Jump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Jump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["jump"]+args)
        elif isinstance(args, str):
            self.commands.append(["jump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Label)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Label(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["label"]+args)
        elif isinstance(args, str):
            self.commands.append(["label "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Log)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Log(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["log"]+args)
        elif isinstance(args, str):
            self.commands.append(["log "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Mixture)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Mixture(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["mixture"]+args)
        elif isinstance(args, str):
            self.commands.append(["mixture "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Move_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Move_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["move_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["move_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Next)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Next(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["next"]+args)
        elif isinstance(args, str):
            self.commands.append(["next "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Package)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Package(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["package"]+args)
        elif isinstance(args, str):
            self.commands.append(["package "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Partition)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Partition(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["partition"]+args)
        elif isinstance(args, str):
            self.commands.append(["partition "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Print)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Print(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["print"]+args)
        elif isinstance(args, str):
            self.commands.append(["print "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Quit)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Quit(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["quit"]+args)
        elif isinstance(args, str):
            self.commands.append(["quit "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.React)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def React(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["react"]+args)
        elif isinstance(args, str):
            self.commands.append(["react "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.React_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def React_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["react_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["react_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_isurf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_isurf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_isurf"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_isurf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Read_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Read_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["read_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["read_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Region)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Region(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["region"]+args)
        elif isinstance(args, str):
            self.commands.append(["region "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Remove_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Remove_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["remove_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["remove_surf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Reset_timestep)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Reset_timestep(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["reset_timestep"]+args)
        elif isinstance(args, str):
            self.commands.append(["reset_timestep "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Run)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Run(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["run"]+args)
        elif isinstance(args, str):
            self.commands.append(["run "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Scale_particles)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Scale_particles(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["scale_particles"]+args)
        elif isinstance(args, str):
            self.commands.append(["scale_particles "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Seed)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Seed(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["seed"]+args)
        elif isinstance(args, str):
            self.commands.append(["seed "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Shell)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Shell(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["shell"]+args)
        elif isinstance(args, str):
            self.commands.append(["shell "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Species)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Species(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["species"]+args)
        elif isinstance(args, str):
            self.commands.append(["species "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Species_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Species_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["species_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["species_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Stats_style)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Stats_style(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["stats_style"]+args)
        elif isinstance(args, str):
            self.commands.append(["stats_style "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Suffix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Suffix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["suffix"]+args)
        elif isinstance(args, str):
            self.commands.append(["suffix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_collide)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_collide(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_collide"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_collide "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_react)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_react(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_react"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_react "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Surf_modify)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Surf_modify(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["surf_modify"]+args)
        elif isinstance(args, str):
            self.commands.append(["surf_modify "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Timestep)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Timestep(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["timestep"]+args)
        elif isinstance(args, str):
            self.commands.append(["timestep "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Uncompute)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Uncompute(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["uncompute"]+args)
        elif isinstance(args, str):
            self.commands.append(["uncompute "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Undump)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Undump(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["undump"]+args)
        elif isinstance(args, str):
            self.commands.append(["undump "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Unfix)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Unfix(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["unfix"]+args)
        elif isinstance(args, str):
            self.commands.append(["unfix "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Units)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Units(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["units"]+args)
        elif isinstance(args, str):
            self.commands.append(["units "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Variable)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Variable(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["variable"]+args)
        elif isinstance(args, str):
            self.commands.append(["variable "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_grid)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_grid(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_grid"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_grid "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_isurf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_isurf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_isurf"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_isurf "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_restart)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_restart(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_restart"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_restart "+args])
        else:
            raise TypeError("invalid input type for args")

    """
    run  
        'class_instance.man(class_instance.Write_surf)'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    """
    def Write_surf(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append(["write_surf"]+args)
        elif isinstance(args, str):
            self.commands.append(["write_surf "+args])
        else:
            raise TypeError("invalid input type for args")
