#ifndef OUT_H
#define OUT_H

#include <string>
#include <iostream>
#include <tuple>
#include <sys/stat.h>

#include "toml.hpp"

#define FLERR __FILE__,__LINE__


// Structure which would store the metadata
struct stat sb;


void error(std::string file, int line, std::string msg) {
    std::cout << "ERROR: (" << file << ":" << line << ") " << msg << "\n";
    exit(1);
}



namespace toml {
    typedef toml::v3::node_view<toml::v3::node> node_t;
    typedef toml::node_type var_type_t;

    const toml::v3::node_type string_t =  toml::node_type::string;
    const toml::v3::node_type integer_t =  toml::node_type::integer;
    const toml::v3::node_type none_t = toml::node_type::none;
    const toml::v3::node_type double_t =  toml::node_type::floating_point;

    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val);

    // errors instead of setting to default
    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var);
}



class Base {
    protected:
        std::string name;
        toml::table contents;

    public:        
        void write_contents(std::string& buffer) {
            for (auto it : this->contents) {
                std::cout << it.first << "\n";
            }
        }
        
        void set_contents(toml::node_t __tbl) {
            if (__tbl.type() != toml::node_type::table) 
                error(FLERR, name + " must be given a table");

            // converting to a table
            toml::table* _tbl = __tbl.as<toml::table>();
            this->contents = (*_tbl);
        }

        void set_contents(toml::table _tbl) {
            this->contents = _tbl;
        }

        std::size_t length() {
            return this->contents.size();
        }
};


class Constants : Base {
    Constants() {
        this->name = "Constants";
    }
};


class Simulation : Base {
    Simulation() {
        this->name = "Simulation";
    }
};


class Elmer {
    public:
        Elmer() {};

        Simulation simulation;
        Constants  constants;

    public:
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

    private:
        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

        void print(std::string _msg) { std::cout << _msg << "\n"; }

        /**
         * executes a table
        */
        template<typename dict>
        void FixFea::eval_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        /**
         * executes a table that wrapped as a node_t type
        */
        template<typename dict>
        void FixFea::eval_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

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

typedef std::pair<std::string, void (FixFea::*)(std::string, toml::node_t)> dict_item_t;
typedef std::vector<dict_item_t> dict_t;

#endif
