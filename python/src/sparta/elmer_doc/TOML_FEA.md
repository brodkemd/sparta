# Expressions
This config file uses toml, a very popular configuration file format. I added the ability to access other variables (possibly other things in the future as well) by using an expression with the format below:
```
"$(SOMETHING)"
```
Replace ```SOMETHING``` with a path to a variable or a command. For example, if you want to get the emi used in sparta, use the command ```"$(sparta.emi)"```. The are also commands available, these are currently:
- HOME : gets the home directory
- PWD : gets the present working directory

Here is a full example:
```toml
simulation_directory = "fea_data"

[sparta]
emi = 0.005
nevery = 1
connect = true
run_every = 50
energy_threshold = 0.0
pressure_threshold = 0.0
shear_threshold = 0.0
groupID = "all"
mixID = "all"
customID = "temperature"

compute = ["CCCC", "surf", "$(sparta.groupID)", "$(sparta.mixID)", "etot"]
surf_collide = ["SSSS", "diffuse", "s_$(sparta.customID)", "$(sparta.emi)"]
surf_modify = ["all", "collide", "$(sparta.surf_collide[0])"]

[elmer]
exe = "$(HOME)/Documents/Apps/elmer/install/bin/ElmerSolver"
sif = "$(simulation_directory)/case.sif"
base_temp = 250.0
meshDB = "sphere_mesh"
node_temperature_file_ext = "node_temperatures"
node_position_file_ext = "node_positions"
gravity_on = true
print_intensity = "limited"
Solver_Input_File = "$(elmer.sif)"
```


# Variable Documentation
- simulation_directory
    - required type = string
    - Description: 
- sparta
    - emi
        - required type = float
        - Description: 
    - nevery
        - required type = int
        - Description: 
    - connect
        - required type = bool
        - Description: 
    - run_every
        - required type = int
        - Description: 
    - energy_threshold
        - required type = float
        - Description: 
    - pressure_threshold
        - required type = float
        - Description: 
    - shear_threshold
        - required type = float
        - Description: 
    - groupID
        - required type = string
        - Description: 
    - mixID
        - required type = string
        - Description: 
    - customID
        - required type = string
        - Description: 
    - compute
        - required type = array[str]
        - Description: 
    - surf_collide
        - required type = array[str]
        - Description: 
    - surf_modify
        - required type = array[str]
        - Description: 
- elmer
    - exe
        - required type = string
        - Description: 
    - sif
        - required type = string
        - Description: 
    - base_temp
        - required type = float
        - Description: 
    - meshDB
        - required type = string
        - Description: 
    - node_temperature_file_ext
        - required type = string
        - Description: 
    - node_position_file_ext
        - required type = string
        - Description: 
    - gravity_on
        - required type = bool
        - Description: 
    - print_intensity
        - required type = string
        - Description: 
    - simulation
        - Max_Output_Level
            - required type = int
            - Description: 
        - Coordinate_System
            - required type = string
            - Description: 
        - Coordinate_Mapping
            - required type = array[int]
            - Description: 
        - Simulation_Type
            - required type = string
            - Description: 
        - Steady_State_Max_Iterations
            - required type = int
            - Description: 
        - Output_Intervals
            - required type = array[int]
            - Description: 
        - Timestep_intervals
            - required type = array[int]
            - Description: 
        - Timestep_Sizes
            - required type = array[float]
            - Description: 
        - Timestepping_Method
            - required type = string
            - Description: 
        - BDF_Order
            - required type = int
            - Description: 
        - Solver_Input_File
            - required type = string
            - Description: 
        - Output_File
            - required type = string
            - Description: 
    - thermal_solver
        - Equation
            - required type = string
            - Description: 
        - Procedure
            - required type = array[str]
            - Description: 
        - Variable
            - required type = string
            - Description: 
        - Exec_Solver
            - required type = string
            - Description: 
        - Stabilize
            - required type = bool
            - Description: 
        - Optimize_Bandwidth
            - required type = bool
            - Description: 
        - Steady_State_Convergence_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Convergence_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Max_Iterations
            - required type = int
            - Description: 
        - Nonlinear_System_Newton_After_Iterations
            - required type = int
            - Description: 
        - Nonlinear_System_Newton_After_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Relaxation_Factor
            - required type = int
            - Description: 
        - Linear_System_Solver
            - required type = string
            - Description: 
        - Linear_System_Iterative_Method
            - required type = string
            - Description: 
        - Linear_System_Max_Iterations
            - required type = int
            - Description: 
        - Linear_System_Convergence_Tolerance
            - required type = float
            - Description: 
        - BiCGstabl_polynomial_degree
            - required type = int
            - Description: 
        - Linear_System_Preconditioning
            - required type = string
            - Description: 
        - Linear_System_ILUT_Tolerance
            - required type = float
            - Description: 
        - Linear_System_Abort_Not_Converged
            - required type = bool
            - Description: 
        - Linear_System_Residual_Output
            - required type = int
            - Description: 
    - elastic_solver
        - Equation
            - required type = string
            - Description: 
        - Variable
            - required type = string
            - Description: 
        - Calculate_Principal
            - required type = bool
            - Description: 
        - Calculate_Stresses
            - required type = bool
            - Description: 
        - Procedure
            - required type = array[str]
            - Description: 
        - Exec_Solver
            - required type = string
            - Description: 
        - Stabilize
            - required type = bool
            - Description: 
        - Optimize_Bandwidth
            - required type = bool
            - Description: 
        - Steady_State_Convergence_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Convergence_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Max_Iterations
            - required type = int
            - Description: 
        - Nonlinear_System_Newton_After_Iterations
            - required type = int
            - Description: 
        - Nonlinear_System_Newton_After_Tolerance
            - required type = float
            - Description: 
        - Nonlinear_System_Relaxation_Factor
            - required type = int
            - Description: 
        - Linear_System_Solver
            - required type = string
            - Description: 
        - Linear_System_Iterative_Method
            - required type = string
            - Description: 
        - Linear_System_Max_Iterations
            - required type = int
            - Description: 
        - Linear_System_Convergence_Tolerance
            - required type = float
            - Description: 
        - BiCGstabl_polynomial_degree
            - required type = int
            - Description: 
        - Linear_System_Preconditioning
            - required type = string
            - Description: 
        - Linear_System_ILUT_Tolerance
            - required type = float
            - Description: 
        - Linear_System_Abort_Not_Converged
            - required type = bool
            - Description: 
        - Linear_System_Residual_Output
            - required type = int
            - Description: 
        - Linear_System_Precondition_Recompute
            - required type = int
            - Description: 
    - equation
        - Name
            - required type = string
            - Description: 
        - Active_Solvers
            - required type = array[int]
            - Description: 
    - material
        - Name
            - required type = string
            - Description: 
        - Poisson_ratio
            - required type = float
            - Description: 
        - Heat_Capacity
            - required type = float
            - Description: 
        - Density
            - required type = float
            - Description: 
        - Youngs_modulus
            - required type = float
            - Description: 
        - Heat_expansion_Coefficient
            - required type = float
            - Description: 
        - Sound_speed
            - required type = float
            - Description: 
        - Heat_Conductivity
            - required type = float
            - Description: 
