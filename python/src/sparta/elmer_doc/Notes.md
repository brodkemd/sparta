# Input File
## Sections

Each section started with a row containing the name of the section, sometimes a number, followed by a number of keyword commands, and ended with a row containing the word End. The names for starting new sections are:

- Header
- Simulation
- Constants
- Body $n$
- Material $n$
- Body Force $n$
- Equation $n$
- Solver $n$
- Boundary Condition $n$
- Initial Condition $n$
- Component $n$

Here $n$ associated with the section name represents an integer identifier needed for distinguishing between
sections of the same type.


The description of the above is given in the subsections below

### Header
The location of mesh files is given in the header section. If the Elmer mesh files are located in the directory `./mymesh`, the
header section may simply be
```
Header
  Mesh DB "." "mymesh"
End
```
Note that separate equations can nevertheless be discretized using different meshes if the location of mesh
files is given in the solver section described below.

### Simulation
The simulation section is used for giving general information that is not specific to a particular PDE model involved in the simulation. This information describes the coordinate system used, indicates whether the problem is stationary or evolutionary, defines the file names for outputting, etc.
```
Simulation
  Coordinate System = "Cartesian 2D"
  Coordinate Mapping(3) = 1 2 3
  Coordinate Scaling = 0.001
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
  Post File = "case.vtu"
  Output File = "case.dat"
  Simulation Timing = True
End
```
Currently the preferred method for visualization is to use VTU output and compatible visualization tool,
such as Paraview. Hence the suffix .vtu is used in the value of the Post File.

### Constants
The constants section is used for defining certain physical constants. For example the gravity vector and the Stefan-Boltzmann constant may be defined using the commands
```
Constants
  Gravity(4) = 0 -1 0 9.82
  Stefan Boltzmann = 5.67e-08
End
```
If the constants are not actually needed in the simulation, this section can also be left empty.

