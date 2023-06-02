class Elmer:
    sif_file:str = ""
    tab = 2*" "
    commands:list[list[str]] = []


    def __init__(self, output_file:str) -> None:
        self.sif_file = output_file


    def write(self):
        with open(self.sif_file, 'w') as f:
            for section in self.commands:
                f.write(f"{section[0]}\n")
                for i in range(1, len(section)-1):
                    f.write(self.tab + section[i] + "\n")
                f.write(f"{section[-1]}\n\n")


    def _add_section(self, _name:str, _n:str, args) -> None:
        self.commands.append([])
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
        

    def Header(self, args) -> None:
        self._add_section("Header", "", args)


    def Simulation(self, args) -> None:
        self._add_section("Simulation", "", args)


    def Constants(self, args) -> None:
        self._add_section("Constants", "", args)


    def default_Constants(self) -> None:
        self.Constants([
            "Gravity(4) = 0 -1 0 9.82",
            "Stefan Boltzmann = 5.67e-08",
            "Permittivity of Vacuum = 8.8542e-12",
            "Permeability of Vacuum = 1.25663706e-6",
            "Boltzmann Constant = 1.3807e-23",
            "Unit Charge = 1.602e-19"
        ])


    def Body(self, _n:str, args) -> None:
        self._add_section("Body", _n, args)


    def Material(self, _n:str, args) -> None:
        self._add_section("Material", _n, args)


    def Body_Force(self, _n:str, args) -> None:
        self._add_section("Body Force", _n, args)


    def Equation(self, _n:str, args) -> None:
        self._add_section("Equation", _n, args)


    def Solver(self, _n:str, args) -> None:
        self._add_section("Solver", _n, args)


    def Boundary_Condition(self, _n:str, args) -> None:
        self._add_section("Boundary Condition", _n, args)


    def Initial_Condition(self, _n:str, args) -> None:
        self._add_section("Initial Condition", _n, args)


    def Component(self, _n:str, args) -> None:
        self._add_section("Component", _n, args)