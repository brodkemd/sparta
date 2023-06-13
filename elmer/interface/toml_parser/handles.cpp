#include "out.h"


void FixFea::handle_both(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("emi", &FixFea::handle_emi),
        std::make_pair("tsurf_file", &FixFea::handle_tsurf_file)
    };

    this->run_table(_caller, "both", tbl, options);
}


void FixFea::handle_emi(std::string _caller, toml::node_t val) {
    toml::set(_caller, "emi", toml::double_t, val, this->emi);
    if (emi <= 0.0 || emi > 1.0)
        error(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");
    this->print("emi = " + std::to_string(this->emi));
}


void FixFea::handle_tsurf_file(std::string _caller, toml::node_t val) {
    toml::set(_caller, "tsurf_file", toml::string_t, val, this->tsurf_file);
    this->print("tsurf_file = " + this->tsurf_file);
}


void FixFea::handle_sparta(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("nevery", &FixFea::handle_nevery),
        std::make_pair("groupID", &FixFea::handle_groupID),
        std::make_pair("mixID", &FixFea::handle_mixID),
        std::make_pair("customID", &FixFea::handle_customID)
    };

    this->run_table(_caller, "sparta", tbl, options);
}


void FixFea::handle_nevery(std::string _caller, toml::node_t val) {
    toml::set(_caller, "nevery", toml::integer_t, val, this->nevery);
    this->print("nevery = " + std::to_string(this->nevery));
}


void FixFea::handle_groupID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "groupID", toml::string_t, val, this->groupID);
    this->print("groupID = " + this->groupID);
}


void FixFea::handle_mixID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "mixID", toml::string_t, val, this->mixID);
    this->print("mixID = " + this->mixID);
}


void FixFea::handle_customID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "customId", toml::string_t, val, this->customID);
    this->print("customId = " + this->customID);
}


void FixFea::handle_elmer(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("exe", &FixFea::handle_exe),
        std::make_pair("sif", &FixFea::handle_sif),
        std::make_pair("meshDBstem", &FixFea::handle_meshDBstem),
        std::make_pair("header", &FixFea::handle_header),
        std::make_pair("simulation", &FixFea::handle_simulation),
        std::make_pair("constants", &FixFea::handle_constants),
        std::make_pair("solver", &FixFea::handle_solver),
        std::make_pair("equation", &FixFea::handle_equation),
        std::make_pair("material", &FixFea::handle_material),
        std::make_pair("body", &FixFea::handle_body),
        std::make_pair("initial_condition", &FixFea::handle_initial_condition),
        std::make_pair("boundary_condition", &FixFea::handle_boundary_condition)
    };

    this->run_table(_caller, "elmer", tbl, options);
}

void FixFea::handle_exe(std::string _caller, toml::node_t val) {
    toml::set(_caller, "exe", toml::string_t, val, this->elmer.exe);
    // Calls the function with path as argument
    // If the file/directory exists at the path returns 0
    // If block executes if path exists
    if (!(stat(this->elmer.exe.c_str(),  &sb) == 0))
        error(FLERR, "Illegal fix fea command, exe path does not exist");
    
    this->print("exe = " + this->elmer.exe);
}

void FixFea::handle_sif(std::string _caller, toml::node_t val) {
    toml::set(_caller, "sif", toml::string_t, val, this->elmer.sif);
    if (!(stat(this->elmer.sif.c_str(),  &sb) == 0))
        error(FLERR, "Illegal fix fea command, sif path does not exist");
    this->print("sif = " + this->elmer.sif);
}

void FixFea::handle_meshDBstem(std::string _caller, toml::node_t val) {
    toml::set(_caller, "meshDBstem", toml::string_t, val, this->elmer.meshDBstem);

    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer.meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error(FLERR, "Illegal fix fea command, mesh database incomplete, " + (this->elmer.meshDBstem + "." + exts[i]) + " does not exist");
    }
    this->print("meshDBstem = " + this->elmer.meshDBstem);
}

void FixFea::handle_header(std::string _caller, toml::node_t tbl) {
    // dict_t options = {
    //     std::make_pair("Mesh_DB", &FixFea::handle_Mesh_DB),
    //     std::make_pair("Include_Path", &FixFea::handle_Include_Path),
    //     std::make_pair("Results_Directory", &FixFea::handle_Results_Directory)
    // };



    // this->run_table(_caller, "header", tbl, options);
}



void FixFea::handle_simulation(std::string _caller, toml::node_t tbl) {


    this->run_table(_caller, "simulation", tbl, options);
}

void FixFea::handle_constants(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("Gravity", &FixFea::handle_Gravity),
        std::make_pair("Stefan_Boltzmann", &FixFea::handle_Stefan_Boltzmann),
        std::make_pair("Permittivity_of_Vacuum", &FixFea::handle_Permittivity_of_Vacuum),
        std::make_pair("Permeability_of_Vacuum", &FixFea::handle_Permeability_of_Vacuum),
        std::make_pair("Boltzmann_Constant", &FixFea::handle_Boltzmann_Constant),
        std::make_pair("Unit_Charge", &FixFea::handle_Unit_Charge)
    };

    this->run_table(_caller, "constants", tbl, options);
}

void FixFea::handle_solver(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "solver", tbl, options);
}

void FixFea::handle_equation(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "equation", tbl, options);
}

void FixFea::handle_material(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "material", tbl, options);
}


void FixFea::handle_body(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "body", tbl, options);
}

void FixFea::handle_initial_condition(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "initial_condition", tbl, options);
}

void FixFea::handle_boundary_condition(std::string _caller, toml::node_t tbl) {
    dict_t options = {
        std::make_pair("1", &FixFea::handle_1)
    };

    this->run_table(_caller, "boundary_condition", tbl, options);
}

// void FixFea::handle_Mesh_DB(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Mesh_DB", toml::list_t, val, this->Mesh_DB);
//     this->print("Mesh_DB = " + std::to_string(this->Mesh_DB));
// }

// void FixFea::handle_Include_Path(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Include_Path", toml::string_t, val, this->Include_Path);
//     this->print("Include_Path = " + std::to_string(this->Include_Path));
// }

// void FixFea::handle_Results_Directory(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Results_Directory", toml::string_t, val, this->Results_Directory);
//     this->print("Results_Directory = " + std::to_string(this->Results_Directory));
// }



// void FixFea::handle_Max_Output_Level(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Max_Output_Level", toml::integer_t, val, this->Max_Output_Level);
//     this->print("Max_Output_Level = " + std::to_string(this->Max_Output_Level));
// }

// void FixFea::handle_Coordinate_System(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Coordinate_System", toml::string_t, val, this->Coordinate_System);
//     this->print("Coordinate_System = " + std::to_string(this->Coordinate_System));
// }

// void FixFea::handle_Coordinate_Mapping(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Coordinate_Mapping", toml::list_t, val, this->Coordinate_Mapping);
//     this->print("Coordinate_Mapping = " + std::to_string(this->Coordinate_Mapping));
// }

// void FixFea::handle_Simulation_Type(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Simulation_Type", toml::string_t, val, this->Simulation_Type);
//     this->print("Simulation_Type = " + std::to_string(this->Simulation_Type));
// }

// void FixFea::handle_Steady_State_Max_Iterations(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Steady_State_Max_Iterations", toml::integer_t, val, this->Steady_State_Max_Iterations);
//     this->print("Steady_State_Max_Iterations = " + std::to_string(this->Steady_State_Max_Iterations));
// }

// void FixFea::handle_Output_Intervals(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Output_Intervals", toml::list_t, val, this->Output_Intervals);
//     this->print("Output_Intervals = " + std::to_string(this->Output_Intervals));
// }

// void FixFea::handle_Timestep_intervals(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Timestep_intervals", toml::list_t, val, this->Timestep_intervals);
//     this->print("Timestep_intervals = " + std::to_string(this->Timestep_intervals));
// }

// void FixFea::handle_Timestep_Sizes(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Timestep_Sizes", toml::list_t, val, this->Timestep_Sizes);
//     this->print("Timestep_Sizes = " + std::to_string(this->Timestep_Sizes));
// }

// void FixFea::handle_Timestepping_Method(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Timestepping_Method", toml::string_t, val, this->Timestepping_Method);
//     this->print("Timestepping_Method = " + std::to_string(this->Timestepping_Method));
// }

// void FixFea::handle_BDF_Order(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "BDF_Order", toml::integer_t, val, this->BDF_Order);
//     this->print("BDF_Order = " + std::to_string(this->BDF_Order));
// }

// void FixFea::handle_Solver_Input_File(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Solver_Input_File", toml::string_t, val, this->Solver_Input_File);
//     this->print("Solver_Input_File = " + std::to_string(this->Solver_Input_File));
// }

// void FixFea::handle_Output_File(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Output_File", toml::string_t, val, this->Output_File);
//     this->print("Output_File = " + std::to_string(this->Output_File));
// }

// void FixFea::handle_Post_File(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Post_File", toml::string_t, val, this->Post_File);
//     this->print("Post_File = " + std::to_string(this->Post_File));
// }



// void FixFea::handle_Gravity(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Gravity", toml::list_t, val, this->Gravity);
//     this->print("Gravity = " + std::to_string(this->Gravity));
// }

// void FixFea::handle_Stefan_Boltzmann(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Stefan_Boltzmann", toml::double_t, val, this->Stefan_Boltzmann);
//     this->print("Stefan_Boltzmann = " + std::to_string(this->Stefan_Boltzmann));
// }

// void FixFea::handle_Permittivity_of_Vacuum(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Permittivity_of_Vacuum", toml::double_t, val, this->Permittivity_of_Vacuum);
//     this->print("Permittivity_of_Vacuum = " + std::to_string(this->Permittivity_of_Vacuum));
// }

// void FixFea::handle_Permeability_of_Vacuum(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Permeability_of_Vacuum", toml::double_t, val, this->Permeability_of_Vacuum);
//     this->print("Permeability_of_Vacuum = " + std::to_string(this->Permeability_of_Vacuum));
// }

// void FixFea::handle_Boltzmann_Constant(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Boltzmann_Constant", toml::double_t, val, this->Boltzmann_Constant);
//     this->print("Boltzmann_Constant = " + std::to_string(this->Boltzmann_Constant));
// }

// void FixFea::handle_Unit_Charge(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Unit_Charge", toml::double_t, val, this->Unit_Charge);
//     this->print("Unit_Charge = " + std::to_string(this->Unit_Charge));
// }



// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Equation", &FixFea::handle_Equation),
//         std::make_pair("Procedure", &FixFea::handle_Procedure),
//         std::make_pair("Variable", &FixFea::handle_Variable),
//         std::make_pair("Exec_Solver", &FixFea::handle_Exec_Solver),
//         std::make_pair("Stabilize", &FixFea::handle_Stabilize),
//         std::make_pair("Optimize_Bandwidth", &FixFea::handle_Optimize_Bandwidth),
//         std::make_pair("Steady_State_Convergence_Tolerance", &FixFea::handle_Steady_State_Convergence_Tolerance),
//         std::make_pair("Nonlinear_System_Convergence_Tolerance", &FixFea::handle_Nonlinear_System_Convergence_Tolerance),
//         std::make_pair("Nonlinear_System_Max_Iterations", &FixFea::handle_Nonlinear_System_Max_Iterations),
//         std::make_pair("Nonlinear_System_Newton_After_Iterations", &FixFea::handle_Nonlinear_System_Newton_After_Iterations),
//         std::make_pair("Nonlinear_System_Newton_After_Tolerance", &FixFea::handle_Nonlinear_System_Newton_After_Tolerance),
//         std::make_pair("Nonlinear_System_Relaxation_Factor", &FixFea::handle_Nonlinear_System_Relaxation_Factor),
//         std::make_pair("Linear_System_Solver", &FixFea::handle_Linear_System_Solver),
//         std::make_pair("Linear_System_Iterative_Method", &FixFea::handle_Linear_System_Iterative_Method),
//         std::make_pair("Linear_System_Max_Iterations", &FixFea::handle_Linear_System_Max_Iterations),
//         std::make_pair("Linear_System_Convergence_Tolerance", &FixFea::handle_Linear_System_Convergence_Tolerance),
//         std::make_pair("BiCGstabl_polynomial_degree", &FixFea::handle_BiCGstabl_polynomial_degree),
//         std::make_pair("Linear_System_Preconditioning", &FixFea::handle_Linear_System_Preconditioning),
//         std::make_pair("Linear_System_ILUT_Tolerance", &FixFea::handle_Linear_System_ILUT_Tolerance),
//         std::make_pair("Linear_System_Abort_Not_Converged", &FixFea::handle_Linear_System_Abort_Not_Converged),
//         std::make_pair("Linear_System_Residual_Output", &FixFea::handle_Linear_System_Residual_Output)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Equation(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Equation", toml::string_t, val, this->Equation);
//     this->print("Equation = " + std::to_string(this->Equation));
// }

// void FixFea::handle_Procedure(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Procedure", toml::list_t, val, this->Procedure);
//     this->print("Procedure = " + std::to_string(this->Procedure));
// }

// void FixFea::handle_Variable(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Variable", toml::string_t, val, this->Variable);
//     this->print("Variable = " + std::to_string(this->Variable));
// }

// void FixFea::handle_Exec_Solver(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Exec_Solver", toml::string_t, val, this->Exec_Solver);
//     this->print("Exec_Solver = " + std::to_string(this->Exec_Solver));
// }

// void FixFea::handle_Stabilize(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Stabilize", toml::bool_t, val, this->Stabilize);
//     this->print("Stabilize = " + std::to_string(this->Stabilize));
// }

// void FixFea::handle_Optimize_Bandwidth(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Optimize_Bandwidth", toml::bool_t, val, this->Optimize_Bandwidth);
//     this->print("Optimize_Bandwidth = " + std::to_string(this->Optimize_Bandwidth));
// }

// void FixFea::handle_Steady_State_Convergence_Tolerance(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Steady_State_Convergence_Tolerance", toml::double_t, val, this->Steady_State_Convergence_Tolerance);
//     this->print("Steady_State_Convergence_Tolerance = " + std::to_string(this->Steady_State_Convergence_Tolerance));
// }

// void FixFea::handle_Nonlinear_System_Convergence_Tolerance(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Nonlinear_System_Convergence_Tolerance", toml::double_t, val, this->Nonlinear_System_Convergence_Tolerance);
//     this->print("Nonlinear_System_Convergence_Tolerance = " + std::to_string(this->Nonlinear_System_Convergence_Tolerance));
// }

// void FixFea::handle_Nonlinear_System_Max_Iterations(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Nonlinear_System_Max_Iterations", toml::integer_t, val, this->Nonlinear_System_Max_Iterations);
//     this->print("Nonlinear_System_Max_Iterations = " + std::to_string(this->Nonlinear_System_Max_Iterations));
// }

// void FixFea::handle_Nonlinear_System_Newton_After_Iterations(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Nonlinear_System_Newton_After_Iterations", toml::integer_t, val, this->Nonlinear_System_Newton_After_Iterations);
//     this->print("Nonlinear_System_Newton_After_Iterations = " + std::to_string(this->Nonlinear_System_Newton_After_Iterations));
// }

// void FixFea::handle_Nonlinear_System_Newton_After_Tolerance(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Nonlinear_System_Newton_After_Tolerance", toml::double_t, val, this->Nonlinear_System_Newton_After_Tolerance);
//     this->print("Nonlinear_System_Newton_After_Tolerance = " + std::to_string(this->Nonlinear_System_Newton_After_Tolerance));
// }

// void FixFea::handle_Nonlinear_System_Relaxation_Factor(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Nonlinear_System_Relaxation_Factor", toml::integer_t, val, this->Nonlinear_System_Relaxation_Factor);
//     this->print("Nonlinear_System_Relaxation_Factor = " + std::to_string(this->Nonlinear_System_Relaxation_Factor));
// }

// void FixFea::handle_Linear_System_Solver(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Solver", toml::string_t, val, this->Linear_System_Solver);
//     this->print("Linear_System_Solver = " + std::to_string(this->Linear_System_Solver));
// }

// void FixFea::handle_Linear_System_Iterative_Method(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Iterative_Method", toml::string_t, val, this->Linear_System_Iterative_Method);
//     this->print("Linear_System_Iterative_Method = " + std::to_string(this->Linear_System_Iterative_Method));
// }

// void FixFea::handle_Linear_System_Max_Iterations(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Max_Iterations", toml::integer_t, val, this->Linear_System_Max_Iterations);
//     this->print("Linear_System_Max_Iterations = " + std::to_string(this->Linear_System_Max_Iterations));
// }

// void FixFea::handle_Linear_System_Convergence_Tolerance(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Convergence_Tolerance", toml::double_t, val, this->Linear_System_Convergence_Tolerance);
//     this->print("Linear_System_Convergence_Tolerance = " + std::to_string(this->Linear_System_Convergence_Tolerance));
// }

// void FixFea::handle_BiCGstabl_polynomial_degree(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "BiCGstabl_polynomial_degree", toml::integer_t, val, this->BiCGstabl_polynomial_degree);
//     this->print("BiCGstabl_polynomial_degree = " + std::to_string(this->BiCGstabl_polynomial_degree));
// }

// void FixFea::handle_Linear_System_Preconditioning(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Preconditioning", toml::string_t, val, this->Linear_System_Preconditioning);
//     this->print("Linear_System_Preconditioning = " + std::to_string(this->Linear_System_Preconditioning));
// }

// void FixFea::handle_Linear_System_ILUT_Tolerance(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_ILUT_Tolerance", toml::double_t, val, this->Linear_System_ILUT_Tolerance);
//     this->print("Linear_System_ILUT_Tolerance = " + std::to_string(this->Linear_System_ILUT_Tolerance));
// }

// void FixFea::handle_Linear_System_Abort_Not_Converged(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Abort_Not_Converged", toml::bool_t, val, this->Linear_System_Abort_Not_Converged);
//     this->print("Linear_System_Abort_Not_Converged = " + std::to_string(this->Linear_System_Abort_Not_Converged));
// }

// void FixFea::handle_Linear_System_Residual_Output(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Linear_System_Residual_Output", toml::integer_t, val, this->Linear_System_Residual_Output);
//     this->print("Linear_System_Residual_Output = " + std::to_string(this->Linear_System_Residual_Output));
// }

// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Name", &FixFea::handle_Name),
//         std::make_pair("Active_Solvers", &FixFea::handle_Active_Solvers)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Name(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Name", toml::string_t, val, this->Name);
//     this->print("Name = " + std::to_string(this->Name));
// }

// void FixFea::handle_Active_Solvers(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Active_Solvers", toml::list_t, val, this->Active_Solvers);
//     this->print("Active_Solvers = " + std::to_string(this->Active_Solvers));
// }

// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Name", &FixFea::handle_Name),
//         std::make_pair("Poisson-ratio", &FixFea::handle_Poisson-ratio),
//         std::make_pair("Heat_Capacity", &FixFea::handle_Heat_Capacity),
//         std::make_pair("Density", &FixFea::handle_Density),
//         std::make_pair("Youngs_modulus", &FixFea::handle_Youngs_modulus),
//         std::make_pair("Heat_expansion_Coefficient", &FixFea::handle_Heat_expansion_Coefficient),
//         std::make_pair("Sound_speed", &FixFea::handle_Sound_speed),
//         std::make_pair("Heat_Conductivity", &FixFea::handle_Heat_Conductivity)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Name(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Name", toml::string_t, val, this->Name);
//     this->print("Name = " + std::to_string(this->Name));
// }

// void FixFea::handle_Poisson-ratio(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Poisson-ratio", toml::double_t, val, this->Poisson-ratio);
//     this->print("Poisson-ratio = " + std::to_string(this->Poisson-ratio));
// }

// void FixFea::handle_Heat_Capacity(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Heat_Capacity", toml::double_t, val, this->Heat_Capacity);
//     this->print("Heat_Capacity = " + std::to_string(this->Heat_Capacity));
// }

// void FixFea::handle_Density(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Density", toml::double_t, val, this->Density);
//     this->print("Density = " + std::to_string(this->Density));
// }

// void FixFea::handle_Youngs_modulus(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Youngs_modulus", toml::double_t, val, this->Youngs_modulus);
//     this->print("Youngs_modulus = " + std::to_string(this->Youngs_modulus));
// }

// void FixFea::handle_Heat_expansion_Coefficient(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Heat_expansion_Coefficient", toml::double_t, val, this->Heat_expansion_Coefficient);
//     this->print("Heat_expansion_Coefficient = " + std::to_string(this->Heat_expansion_Coefficient));
// }

// void FixFea::handle_Sound_speed(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Sound_speed", toml::double_t, val, this->Sound_speed);
//     this->print("Sound_speed = " + std::to_string(this->Sound_speed));
// }

// void FixFea::handle_Heat_Conductivity(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Heat_Conductivity", toml::double_t, val, this->Heat_Conductivity);
//     this->print("Heat_Conductivity = " + std::to_string(this->Heat_Conductivity));
// }



// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Target_Bodies", &FixFea::handle_Target_Bodies),
//         std::make_pair("Name", &FixFea::handle_Name),
//         std::make_pair("Equation", &FixFea::handle_Equation),
//         std::make_pair("Material", &FixFea::handle_Material),
//         std::make_pair("Initial_condition", &FixFea::handle_Initial_condition)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Target_Bodies(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Target_Bodies", toml::list_t, val, this->Target_Bodies);
//     this->print("Target_Bodies = " + std::to_string(this->Target_Bodies));
// }

// void FixFea::handle_Name(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Name", toml::string_t, val, this->Name);
//     this->print("Name = " + std::to_string(this->Name));
// }

// void FixFea::handle_Equation(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Equation", toml::integer_t, val, this->Equation);
//     this->print("Equation = " + std::to_string(this->Equation));
// }

// void FixFea::handle_Material(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Material", toml::integer_t, val, this->Material);
//     this->print("Material = " + std::to_string(this->Material));
// }

// void FixFea::handle_Initial_condition(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Initial_condition", toml::integer_t, val, this->Initial_condition);
//     this->print("Initial_condition = " + std::to_string(this->Initial_condition));
// }



// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Name", &FixFea::handle_Name),
//         std::make_pair("Temperature", &FixFea::handle_Temperature)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Name(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Name", toml::string_t, val, this->Name);
//     this->print("Name = " + std::to_string(this->Name));
// }

// void FixFea::handle_Temperature(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Temperature", toml::integer_t, val, this->Temperature);
//     this->print("Temperature = " + std::to_string(this->Temperature));
// }



// void FixFea::handle_1(std::string _caller, toml::node_t tbl) {
//     dict_t options = {
//         std::make_pair("Target_Boundaries", &FixFea::handle_Target_Boundaries),
//         std::make_pair("Name", &FixFea::handle_Name),
//         std::make_pair("Heat_Flux_BC", &FixFea::handle_Heat_Flux_BC),
//         std::make_pair("Heat_Flux", &FixFea::handle_Heat_Flux)
//     };

//     this->run_table(_caller, "1", tbl, options);
// }

// void FixFea::handle_Target_Boundaries(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Target_Boundaries", toml::list_t, val, this->Target_Boundaries);
//     this->print("Target_Boundaries = " + std::to_string(this->Target_Boundaries));
// }

// void FixFea::handle_Name(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Name", toml::string_t, val, this->Name);
//     this->print("Name = " + std::to_string(this->Name));
// }

// void FixFea::handle_Heat_Flux_BC(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Heat_Flux_BC", toml::bool_t, val, this->Heat_Flux_BC);
//     this->print("Heat_Flux_BC = " + std::to_string(this->Heat_Flux_BC));
// }

// void FixFea::handle_Heat_Flux(std::string _caller, toml::node_t val) {
//     toml::set(_caller, "Heat_Flux", toml::integer_t, val, this->Heat_Flux);
//     this->print("Heat_Flux = " + std::to_string(this->Heat_Flux));
// }

