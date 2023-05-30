#include "materials.h"
#include "elmer.h"

void print() {
    std::cout << '\n';
}

template <typename T, typename... TAIL>
void print(const T &t, TAIL... tail) {
    std::cout << t << " ";
    print(tail...);
}

int main() {
    std::cout << "\n";

    Data data;
    std::vector<std::string> temp;
    data.Aluminum_Generic.to_string(temp);

    std::string file = "test_case.sif";
    std::string file_stem = file.substr(0, file.find("."));

    Elmer inst(file);
    
    inst.Header();
    
    inst.Constants();
    
    inst.Material(1, temp);

    inst.Simulation({
        "Max Output Level = 4",
        "Coordinate System = Cartesian",
        "Coordinate Mapping(3) = 1 2 3",
        "Simulation Type = Steady state",
        "Steady State Max Iterations = 100",
        "Output Intervals(1) = 1",
        "Solver Input File = " + file,
        "Post File = " + file_stem + ".vtu"
    });

    inst.Body(1, {
        " Target Bodies(1) = 1",
        "Name = \"Body Property 1\"",
        "Equation = 1",
        "Material = 1",
        "Body Force = 1"
    });

    inst.Solver(1, {
        "Equation = Heat Equation",
        "Procedure = \"HeatSolve\" \"HeatSolver\"",
        "Variable = Temperature",
        "Exec Solver = Always",
        "Stabilize = True",
        "Optimize Bandwidth = True",
        "Steady State Convergence Tolerance = 1.0e-5",
        "Nonlinear System Convergence Tolerance = 1.0e-7",
        "Nonlinear System Max Iterations = 20",
        "Nonlinear System Newton After Iterations = 3",
        "Nonlinear System Newton After Tolerance = 1.0e-3",
        "Nonlinear System Relaxation Factor = 1",
        "Linear System Solver = Iterative",
        "Linear System Iterative Method = BiCGStab",
        "Linear System Max Iterations = 500",
        "Linear System Convergence Tolerance = 1.0e-10",
        "BiCGstabl polynomial degree = 2",
        "Linear System Preconditioning = ILU0",
        "Linear System ILUT Tolerance = 1.0e-3",
        "Linear System Abort Not Converged = False",
        "Linear System Residual Output = 10",
        "Linear System Precondition Recompute = 1"
    });

    inst.Equation(1, {
        "Name = \"Heat Equation\"",
        "Active Solvers(1) = 1"
    });

    inst.Body_Force(1, {
        "Name = \"Heating\"",
        "Heat Source = 0.01"
    });

    inst.Boundary_Conditions(1, {
        "Target Boundaries(1) = 57",
        "Name = \"RoomTemp\"",
        "Temperature = 293.0"
    });

    inst.write();
    return 0;
}