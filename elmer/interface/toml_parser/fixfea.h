#ifndef FIXFEA_H
#define FIXFEA_H

#include "toml.h"

#include <sys/stat.h>
#include <iostream>

class Base {
    protected:
        std::string name;

    public:
        std::vector<std::string> contents;

        void write_contents(std::string& buffer) {
            for (auto it : this->contents) {
                std::cout << it << "\n";
            }
        }
        
        // void set_contents(std::vector<std::string> in) {
        //     this->contents = in
        // }

        // void set_contents(toml::table _tbl) {
        //     this->contents = _tbl;
        // }

        // std::size_t length() {
        //     return this->contents.size();
        // }
};


class Header : public Base {
    public:
        Header() {
            this->name = "Header";
        }
};


class Constants : public Base {
    public:
        Constants() {
            this->name = "Constants";
        }
};


class Simulation : public Base {
    public:
        Simulation() {
            this->name = "Simulation";
        }
};


class Elmer {
    public:
        Elmer() {
            this->name = "elmer";
            this->header = Header();
            this->simulation = Simulation();
            this->constants = Constants();
        }

        std::string name;

        Header     header;
        Simulation simulation;
        Constants  constants;


        std::string exe;
        std::string sif;
        std::string meshDBstem;

};


class FixFea {
    protected:
        // under "both"
        double emi;
        std::string tsurf_file;

        // under "sparta"
        int nevery;
        std::string groupID;
        std::string mixID;
        std::string customID;

        // under elmer
        Elmer elmer = Elmer();

        // Structure which would store the metadata
        struct stat sb;     

    public:
        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

        void print(std::string _msg) { std::cout << _msg << "\n"; }

        // /**
        //  * executes a table
        // */
        // template<typename dict>
        // void FixFea::eval_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        // /**
        //  * executes a table that wrapped as a node_t type
        // */
        // template<typename dict>
        // void FixFea::eval_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

    public:
        void handle_emi(std::string _caller, toml::node_t val);
        void handle_tsurf_file(std::string _caller, toml::node_t val);
        void handle_both(std::string _caller, toml::node_t tbl);
        void handle_nevery(std::string _caller, toml::node_t val);
        void handle_groupID(std::string _caller, toml::node_t val);
        void handle_mixID(std::string _caller, toml::node_t val);
        void handle_customID(std::string _caller, toml::node_t val);
        void handle_sparta(std::string _caller, toml::node_t tbl);
        void handle_exe(std::string _caller, toml::node_t val);
        void handle_sif(std::string _caller, toml::node_t val);
        void handle_meshDBstem(std::string _caller, toml::node_t val);
        // void handle_Mesh_DB(std::string _caller, toml::node_t val);
        // void handle_Include_Path(std::string _caller, toml::node_t val);
        // void handle_Results_Directory(std::string _caller, toml::node_t val);
        void handle_header(std::string _caller, toml::node_t tbl);
        // void handle_Max_Output_Level(std::string _caller, toml::node_t val);
        // void handle_Coordinate_System(std::string _caller, toml::node_t val);
        // void handle_Coordinate_Mapping(std::string _caller, toml::node_t val);
        // void handle_Simulation_Type(std::string _caller, toml::node_t val);
        // void handle_Steady_State_Max_Iterations(std::string _caller, toml::node_t val);
        // void handle_Output_Intervals(std::string _caller, toml::node_t val);
        // void handle_Timestep_intervals(std::string _caller, toml::node_t val);
        // void handle_Timestep_Sizes(std::string _caller, toml::node_t val);
        // void handle_Timestepping_Method(std::string _caller, toml::node_t val);
        // void handle_BDF_Order(std::string _caller, toml::node_t val);
        // void handle_Solver_Input_File(std::string _caller, toml::node_t val);
        // void handle_Output_File(std::string _caller, toml::node_t val);
        // void handle_Post_File(std::string _caller, toml::node_t val);
        void handle_simulation(std::string _caller, toml::node_t tbl);
        // void handle_Gravity(std::string _caller, toml::node_t val);
        // void handle_Stefan_Boltzmann(std::string _caller, toml::node_t val);
        // void handle_Permittivity_of_Vacuum(std::string _caller, toml::node_t val);
        // void handle_Permeability_of_Vacuum(std::string _caller, toml::node_t val);
        // void handle_Boltzmann_Constant(std::string _caller, toml::node_t val);
        // void handle_Unit_Charge(std::string _caller, toml::node_t val);
        void handle_constants(std::string _caller, toml::node_t tbl);
        // void handle_Equation(std::string _caller, toml::node_t val);
        // void handle_Procedure(std::string _caller, toml::node_t val);
        // void handle_Variable(std::string _caller, toml::node_t val);
        // void handle_Exec_Solver(std::string _caller, toml::node_t val);
        // void handle_Stabilize(std::string _caller, toml::node_t val);
        // void handle_Optimize_Bandwidth(std::string _caller, toml::node_t val);
        // void handle_Steady_State_Convergence_Tolerance(std::string _caller, toml::node_t val);
        // void handle_Nonlinear_System_Convergence_Tolerance(std::string _caller, toml::node_t val);
        // void handle_Nonlinear_System_Max_Iterations(std::string _caller, toml::node_t val);
        // void handle_Nonlinear_System_Newton_After_Iterations(std::string _caller, toml::node_t val);
        // void handle_Nonlinear_System_Newton_After_Tolerance(std::string _caller, toml::node_t val);
        // void handle_Nonlinear_System_Relaxation_Factor(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Solver(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Iterative_Method(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Max_Iterations(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Convergence_Tolerance(std::string _caller, toml::node_t val);
        // void handle_BiCGstabl_polynomial_degree(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Preconditioning(std::string _caller, toml::node_t val);
        // void handle_Linear_System_ILUT_Tolerance(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Abort_Not_Converged(std::string _caller, toml::node_t val);
        // void handle_Linear_System_Residual_Output(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_solver(std::string _caller, toml::node_t tbl);
        // void handle_Name(std::string _caller, toml::node_t val);
        // void handle_Active_Solvers(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_equation(std::string _caller, toml::node_t tbl);
        // void handle_Name(std::string _caller, toml::node_t val);
        // void handle_Poisson_ratio(std::string _caller, toml::node_t val);
        // void handle_Heat_Capacity(std::string _caller, toml::node_t val);
        // void handle_Density(std::string _caller, toml::node_t val);
        // void handle_Youngs_modulus(std::string _caller, toml::node_t val);
        // void handle_Heat_expansion_Coefficient(std::string _caller, toml::node_t val);
        // void handle_Sound_speed(std::string _caller, toml::node_t val);
        // void handle_Heat_Conductivity(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_material(std::string _caller, toml::node_t tbl);
        // void handle_Target_Bodies(std::string _caller, toml::node_t val);
        // void handle_Name(std::string _caller, toml::node_t val);
        // void handle_Equation(std::string _caller, toml::node_t val);
        // void handle_Material(std::string _caller, toml::node_t val);
        // void handle_Initial_condition(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_body(std::string _caller, toml::node_t tbl);
        // void handle_Name(std::string _caller, toml::node_t val);
        // void handle_Temperature(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_initial_condition(std::string _caller, toml::node_t tbl);
        // void handle_Target_Boundaries(std::string _caller, toml::node_t val);
        // void handle_Name(std::string _caller, toml::node_t val);
        // void handle_Heat_Flux_BC(std::string _caller, toml::node_t val);
        // void handle_Heat_Flux(std::string _caller, toml::node_t val);
        // void handle_1(std::string _caller, toml::node_t tbl);
        void handle_boundary_condition(std::string _caller, toml::node_t tbl);
        void handle_elmer(std::string _caller, toml::node_t tbl);
};

namespace toml {
    typedef std::pair<std::string, void (FixFea::*)(std::string, toml::node_t)> dict_item_t;
    typedef std::vector<dict_item_t> dict_t;
}

#endif