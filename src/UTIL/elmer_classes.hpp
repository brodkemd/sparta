#ifndef ELMER_CLASSES_H
#define ELMER_CLASSES_H

#include <algorithm>

#include "toml.hpp"

namespace elmer {
    util::string_t makeVarName(util::string_t _var_name) {
        std::replace(_var_name.begin(), _var_name.end(), '_', ' ');
        return _var_name;
    }

    /**
     * Base class that all elmer sections inherit from
    */
    class Section {
        private:
            
            util::string_t _name, _sep;

            // tab in the file
            util::string_t _tab = "  ";

            // ends the section
            util::string_t _end = "End";
            util::int_t _id;
            util::bool_t _include_count;

        public:
            toml::OrderedDict_t contents;
            Section() {}
            Section(util::string_t name, util::int_t id = toml::noInt, util::string_t sep = " = ", util::bool_t include_count = true) {
                this->_name = name;
                this->_sep = sep;
                this->_id = id;
                this->contents = toml::OrderedDict_t();
                this->_include_count = include_count;
            }

            void join(util::oFile& _buf) {
                if (this->_id == toml::noInt)
                    _buf << this->_name << "\n";
                else
                    _buf << this->_name << " " << _id << "\n";
                
                toml::Item_t keys = this->contents.getKeys();
                toml::Item_t vals = this->contents.getValues();
                for (util::int_t i = 0; i < keys.length(); i++) {
                    _buf << (_tab + makeVarName(keys[i].toString()));
                    if (vals[i].getType() == toml::LIST) {
                        if (this->_include_count)
                            _buf << "(" << std::to_string(vals[i].length()) << ")";
                        _buf << _sep << vals[i].toString(" ", true, 1000);
                    } else {
                        _buf << _sep << vals[i].toString(" ", true);
                    }
                    _buf << "\n";
                }
                // this->contents.toFile(_buf, _tab, this->_sep);
                _buf << this->_end << "\n";
            }

            util::int_t getId() {
                return this->_id;
            }

            toml::Item_t getItem(toml::Item_t key) {
                return this->contents.getItem(key);
            }

            void setItem(toml::Item_t key, toml::Item_t val) {
                this->contents.setItem(key, val);
            }
            void clear() { this->contents.clear(); }
            util::int_t length() { return this->contents.length(); }
    };
        // private:
        //     /**
        //      * these functions convert the inputted var into a string
        //     */
        //     void _varToString(util::string_t& _buf, util::string_t _var) {
        //          _buf = "\""+_var+"\"";
        //     }

        //     void _varToString(util::string_t& _buf, double_t _var) {
        //         _buf = util::dtos(_var);
        //     }

        //     void _varToString(util::string_t& _buf, util::int_t _var) {
        //         _buf = std::to_string(_var);
        //     }

        //     void _varToString(util::string_t& _buf, util::bool_t _var) {
        //         if (_var) _buf = "True";
        //         else _buf = "False";
        //     }

        // public:
        //     // set function, this is overriden
        //     virtual void set(toml::handler& _h) { return; }

        //     /**
        //      * these functions are outward facing and they build elmer file
        //      * lines from the inputs
        //     */ 
        //     void varToString(util::oFile& _buf, util::string_t _name, util::string_t _var) {
        //         if (_var == toml::noString) return;
        //         _buf << (_tab + _name + sep + "\""+_var+"\"" + "\n");
        //     }

        //     void varToString(util::oFile& _buf, util::string_t _name, double_t _var) {
        //         if (_var == toml::noDouble) return;
        //         _buf << (_tab + _name + sep) << _var << "\n";
        //     }

        //     void varToString(util::oFile& _buf, util::string_t _name, util::int_t _var) {
        //         if (_var == toml::noInt) return;
        //         _buf << (_tab + _name + sep) << _var << "\n";
        //     }

        //     void varToString(util::oFile& _buf, util::string_t _name, util::bool_t _var) {
        //         if (_var == toml::noBool) return;
        //         if (_var) _buf << (_tab + _name + sep + "True" + "\n");
        //         else      _buf << (_tab + _name + sep + "False" + "\n");
        //     }

        //     // for vectors
        //     template<typename T>
        //     void varToString(util::oFile& _buf, util::string_t _name, std::vector<T> _var, util::bool_t include_count = true, util::int_t wrap_every = 1000) {

        //         if (_var.size() == 0) return;

        //         util::string_t _temp;
        //         if (include_count)
        //             _buf << (_tab + _name + "(") << _var.size() << (")" + sep);
        //         else 
        //             _buf << (_tab + _name + sep);

        //         util::int_t count = 1;
        //         for (std::size_t i = 0; i < (_var.size() - 1); i++){
        //             _varToString(_temp, _var[i]);
        //             if (count >= wrap_every) {
        //                 _buf << "\n";
        //                 count = 0;
        //             }
        //             _buf<<(_temp + " ");
        //             count++;
        //         }
        //         _varToString(_temp, _var[_var.size() - 1]);
        //         _buf << (_temp + "\n");
        //     }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer header section
    // */
    // class Header : public Base {
    //     public:
    //         Header() : Base() {
    //             this->name = "Header";
    //             sep = " ";
    //         }

    //         std::vector<util::string_t> Mesh_DB;
    //         util::string_t Include_Path, Results_Directory;

    //         void set(toml::handler& _h) {
    //             _h.getAtPath(Mesh_DB,           "elmer.header.Mesh_DB", true);
    //             _h.getAtPath(Include_Path,      "elmer.header.Include_Path", false);
    //             _h.getAtPath(Results_Directory, "elmer.header.Results_Directory", true);
    //         }

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             // false tells it to not put array length in array definition in file
    //             varToString(_buf, "Mesh DB",           Mesh_DB, false);
    //             varToString(_buf, "Include Path",      Include_Path);
    //             varToString(_buf, "Results Directory", Results_Directory);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer constants section
    // */
    // class Constants : public Base {
    //     public:
    //         Constants() : Base() {
    //             this->name = "Constants";
    //             sep = " = ";
    //         }

    //         std::vector<double_t> Gravity;
    //         double_t Stefan_Boltzmann, Permittivity_of_Vacuum, Permeability_of_Vacuum, Boltzmann_Constant, Unit_Charge;


    //         void set(toml::handler& _h) {
    //             _h.getAtPath(Gravity,                 "elmer.constants.Gravity", true);
    //             _h.getAtPath(Stefan_Boltzmann,        "elmer.constants.Stefan_Boltzmann", true);
    //             _h.getAtPath(Permittivity_of_Vacuum,  "elmer.constants.Permittivity_of_Vacuum", true);
    //             _h.getAtPath(Permeability_of_Vacuum,  "elmer.constants.Permeability_of_Vacuum", true);
    //             _h.getAtPath(Boltzmann_Constant,      "elmer.constants.Boltzmann_Constant", true);
    //             _h.getAtPath(Unit_Charge,             "elmer.constants.Unit_Charge", true);
    //         }


    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Gravity",                Gravity);
    //             varToString(_buf, "Stefan Boltzmann",       Stefan_Boltzmann);
    //             varToString(_buf, "Permittivity of Vacuum", Permittivity_of_Vacuum);
    //             varToString(_buf, "Permeability of Vacuum", Permeability_of_Vacuum);
    //             varToString(_buf, "Boltzmann Constant",     Boltzmann_Constant);
    //             varToString(_buf, "Unit Charge",            Unit_Charge);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer simulation section
    // */
    // class Simulation : public Base {
    //     public:
    //         Simulation() : Base() {
    //             this->name = "Simulation";
    //             sep = " = ";
    //         }

    //         void set(toml::handler& _h) {
    //             _h.getDictAtPath(this->contents, "elmer.simulation");
    //         }


    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             this->contents.toFile(_buf, "  ", " = ", );
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer thermal_solver section
    // */
    // class Thermal_Solver : public Base {
    //     public:
    //         Thermal_Solver() : Base() {
    //             this->name = "solver";
    //             sep = " = ";
    //             id = 1;
    //             this->name += (" " + std::to_string(id));
    //         }

    //         std::vector<util::string_t> Procedure;
    //         util::string_t Equation, Variable, Exec_Solver, Linear_System_Solver, Linear_System_Iterative_Method, Linear_System_Preconditioning;
    //         util::int_t Linear_System_Residual_Output, Nonlinear_System_Max_Iterations, Nonlinear_System_Newton_After_Iterations, Linear_System_Max_Iterations, BiCGstabl_polynomial_degree, Nonlinear_System_Relaxation_Factor, Linear_System_Precondition_Recompute;
    //         util::bool_t Stabilize, Optimize_Bandwidth, Linear_System_Abort_Not_Converged, Calculate_Principal, Calculate_Stresses;
    //         double_t Steady_State_Convergence_Tolerance, Nonlinear_System_Convergence_Tolerance, Nonlinear_System_Newton_After_Tolerance, Linear_System_Convergence_Tolerance, Linear_System_ILUT_Tolerance;


    //         void set(toml::handler& _h) {
    //             _h.getAtPath(Equation,                                 "elmer.thermal_solver.Equation", false);
    //             _h.getAtPath(Procedure,                                "elmer.thermal_solver.Procedure", false);
    //             _h.getAtPath(Variable,                                 "elmer.thermal_solver.Variable", false);
    //             _h.getAtPath(Exec_Solver,                              "elmer.thermal_solver.Exec_Solver", false);
    //             _h.getAtPath(Stabilize,                                "elmer.thermal_solver.Stabilize", false);
    //             _h.getAtPath(Optimize_Bandwidth,                       "elmer.thermal_solver.Optimize_Bandwidth", false);
    //             _h.getAtPath(Steady_State_Convergence_Tolerance,       "elmer.thermal_solver.Steady_State_Convergence_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Convergence_Tolerance,   "elmer.thermal_solver.Nonlinear_System_Convergence_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Max_Iterations,          "elmer.thermal_solver.Nonlinear_System_Max_Iterations", false);
    //             _h.getAtPath(Nonlinear_System_Newton_After_Iterations, "elmer.thermal_solver.Nonlinear_System_Newton_After_Iterations", false);
    //             _h.getAtPath(Nonlinear_System_Newton_After_Tolerance,  "elmer.thermal_solver.Nonlinear_System_Newton_After_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Relaxation_Factor,       "elmer.thermal_solver.Nonlinear_System_Relaxation_Factor", false);
    //             _h.getAtPath(Linear_System_Solver,                     "elmer.thermal_solver.Linear_System_Solver", false);
    //             _h.getAtPath(Linear_System_Iterative_Method,           "elmer.thermal_solver.Linear_System_Iterative_Method", false);
    //             _h.getAtPath(Linear_System_Max_Iterations,             "elmer.thermal_solver.Linear_System_Max_Iterations", false);
    //             _h.getAtPath(Linear_System_Convergence_Tolerance,      "elmer.thermal_solver.Linear_System_Convergence_Tolerance", false);
    //             _h.getAtPath(BiCGstabl_polynomial_degree,              "elmer.thermal_solver.BiCGstabl_polynomial_degree", false);
    //             _h.getAtPath(Linear_System_Preconditioning,            "elmer.thermal_solver.Linear_System_Preconditioning", false);
    //             _h.getAtPath(Linear_System_ILUT_Tolerance,             "elmer.thermal_solver.Linear_System_ILUT_Tolerance", false);
    //             _h.getAtPath(Linear_System_Abort_Not_Converged,        "elmer.thermal_solver.Linear_System_Abort_Not_Converged", false);
    //             _h.getAtPath(Linear_System_Residual_Output,            "elmer.thermal_solver.Linear_System_Residual_Output", false);
    //         }


    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Equation",                                  Equation);
    //             varToString(_buf, "Procedure",                                 Procedure, false);
    //             varToString(_buf, "Variable",                                  Variable);
    //             varToString(_buf, "Exec Solver",                               Exec_Solver);
    //             varToString(_buf, "Stabilize",                                 Stabilize);
    //             varToString(_buf, "Optimize Bandwidth",                        Optimize_Bandwidth);
    //             varToString(_buf, "Steady State Convergence Tolerance",        Steady_State_Convergence_Tolerance);
    //             varToString(_buf, "Nonlinear System Convergence Tolerance",    Nonlinear_System_Convergence_Tolerance);
    //             varToString(_buf, "Nonlinear System Max Iterations",           Nonlinear_System_Max_Iterations);
    //             varToString(_buf, "Nonlinear System Newton After Iterations",  Nonlinear_System_Newton_After_Iterations);
    //             varToString(_buf, "Nonlinear System Newton After Tolerance",   Nonlinear_System_Newton_After_Tolerance);
    //             varToString(_buf, "Nonlinear System Relaxation Factor",        Nonlinear_System_Relaxation_Factor);
    //             varToString(_buf, "Linear System Solver",                      Linear_System_Solver);
    //             varToString(_buf, "Linear System Iterative Method",            Linear_System_Iterative_Method);
    //             varToString(_buf, "Linear System Max Iterations",              Linear_System_Max_Iterations);
    //             varToString(_buf, "Linear System Convergence Tolerance",       Linear_System_Convergence_Tolerance);
    //             varToString(_buf, "BiCGstabl polynomial degree",               BiCGstabl_polynomial_degree);
    //             varToString(_buf, "Linear System Preconditioning",             Linear_System_Preconditioning);
    //             varToString(_buf, "Linear System ILUT Tolerance",              Linear_System_ILUT_Tolerance);
    //             varToString(_buf, "Linear System Abort Not Converged",         Linear_System_Abort_Not_Converged);
    //             varToString(_buf, "Linear System Residual Output",             Linear_System_Residual_Output);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer thermal_solver section
    // */
    // class Elastic_Solver : public Base {
    //     public:
    //         Elastic_Solver() : Base() {
    //             this->name = "solver";
    //             sep = " = ";
    //             id = 2;
    //             this->name += (" " + std::to_string(id));
    //         }

    //         std::vector<util::string_t> Procedure;
    //         util::string_t Equation, Variable, Exec_Solver, Linear_System_Solver, Linear_System_Iterative_Method, Linear_System_Preconditioning;
    //         util::int_t Linear_System_Residual_Output, Nonlinear_System_Max_Iterations, Nonlinear_System_Newton_After_Iterations, Linear_System_Max_Iterations, BiCGstabl_polynomial_degree, Nonlinear_System_Relaxation_Factor, Linear_System_Precondition_Recompute;
    //         util::bool_t Stabilize, Optimize_Bandwidth, Linear_System_Abort_Not_Converged, Calculate_Principal, Calculate_Stresses;
    //         double_t Steady_State_Convergence_Tolerance, Nonlinear_System_Convergence_Tolerance, Nonlinear_System_Newton_After_Tolerance, Linear_System_Convergence_Tolerance, Linear_System_ILUT_Tolerance;


    //         void set(toml::handler& _h) {
    //             _h.getAtPath(Equation,                                 "elmer.elastic_solver.Equation", false);
    //             _h.getAtPath(Procedure,                                "elmer.elastic_solver.Procedure", false);
    //             _h.getAtPath(Calculate_Principal,                      "elmer.elastic_solver.Calculate_Principal", false);
    //             _h.getAtPath(Calculate_Stresses,                       "elmer.elastic_solver.Calculate_Stresses", false);
    //             _h.getAtPath(Variable,                                 "elmer.elastic_solver.Variable", false);
    //             _h.getAtPath(Exec_Solver,                              "elmer.elastic_solver.Exec_Solver", false);
    //             _h.getAtPath(Stabilize,                                "elmer.elastic_solver.Stabilize", false);
    //             _h.getAtPath(Optimize_Bandwidth,                       "elmer.elastic_solver.Optimize_Bandwidth", false);
    //             _h.getAtPath(Steady_State_Convergence_Tolerance,       "elmer.elastic_solver.Steady_State_Convergence_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Convergence_Tolerance,   "elmer.elastic_solver.Nonlinear_System_Convergence_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Max_Iterations,          "elmer.elastic_solver.Nonlinear_System_Max_Iterations", false);
    //             _h.getAtPath(Nonlinear_System_Newton_After_Iterations, "elmer.elastic_solver.Nonlinear_System_Newton_After_Iterations", false);
    //             _h.getAtPath(Nonlinear_System_Newton_After_Tolerance,  "elmer.elastic_solver.Nonlinear_System_Newton_After_Tolerance", false);
    //             _h.getAtPath(Nonlinear_System_Relaxation_Factor,       "elmer.elastic_solver.Nonlinear_System_Relaxation_Factor", false);
    //             _h.getAtPath(Linear_System_Solver,                     "elmer.elastic_solver.Linear_System_Solver", false);
    //             _h.getAtPath(Linear_System_Iterative_Method,           "elmer.elastic_solver.Linear_System_Iterative_Method", false);
    //             _h.getAtPath(Linear_System_Max_Iterations,             "elmer.elastic_solver.Linear_System_Max_Iterations", false);
    //             _h.getAtPath(Linear_System_Convergence_Tolerance,      "elmer.elastic_solver.Linear_System_Convergence_Tolerance", false);
    //             _h.getAtPath(BiCGstabl_polynomial_degree,              "elmer.elastic_solver.BiCGstabl_polynomial_degree", false);
    //             _h.getAtPath(Linear_System_Preconditioning,            "elmer.elastic_solver.Linear_System_Preconditioning", false);
    //             _h.getAtPath(Linear_System_ILUT_Tolerance,             "elmer.elastic_solver.Linear_System_ILUT_Tolerance", false);
    //             _h.getAtPath(Linear_System_Abort_Not_Converged,        "elmer.elastic_solver.Linear_System_Abort_Not_Converged", false);
    //             _h.getAtPath(Linear_System_Residual_Output,            "elmer.elastic_solver.Linear_System_Residual_Output", false);
    //             _h.getAtPath(Linear_System_Precondition_Recompute,     "elmer.elastic_solver.Linear_System_Precondition_Recompute", false);
    //         }


    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Equation",                                  Equation);
    //             varToString(_buf, "Procedure",                                 Procedure, false);
    //             varToString(_buf, "Variable",                                  Variable);
    //             varToString(_buf, "Exec Solver",                               Exec_Solver);
    //             varToString(_buf, "Stabilize",                                 Stabilize);
    //             varToString(_buf, "Optimize Bandwidth",                        Optimize_Bandwidth);
    //             varToString(_buf, "Steady State Convergence Tolerance",        Steady_State_Convergence_Tolerance);
    //             varToString(_buf, "Nonlinear System Convergence Tolerance",    Nonlinear_System_Convergence_Tolerance);
    //             varToString(_buf, "Nonlinear System Max Iterations",           Nonlinear_System_Max_Iterations);
    //             varToString(_buf, "Nonlinear System Newton After Iterations",  Nonlinear_System_Newton_After_Iterations);
    //             varToString(_buf, "Nonlinear System Newton After Tolerance",   Nonlinear_System_Newton_After_Tolerance);
    //             varToString(_buf, "Nonlinear System Relaxation Factor",        Nonlinear_System_Relaxation_Factor);
    //             varToString(_buf, "Linear System Solver",                      Linear_System_Solver);
    //             varToString(_buf, "Linear System Iterative Method",            Linear_System_Iterative_Method);
    //             varToString(_buf, "Linear System Max Iterations",              Linear_System_Max_Iterations);
    //             varToString(_buf, "Linear System Convergence Tolerance",       Linear_System_Convergence_Tolerance);
    //             varToString(_buf, "BiCGstabl polynomial degree",               BiCGstabl_polynomial_degree);
    //             varToString(_buf, "Linear System Preconditioning",             Linear_System_Preconditioning);
    //             varToString(_buf, "Linear System ILUT Tolerance",              Linear_System_ILUT_Tolerance);
    //             varToString(_buf, "Linear System Abort Not Converged",         Linear_System_Abort_Not_Converged);
    //             varToString(_buf, "Linear System Residual Output",             Linear_System_Residual_Output);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer equation section
    // */
    // class Equation : public Base {
    //     public:
    //         Equation() : Base() {
    //             this->name = "Equation";
    //             sep = " = ";
    //             id = 1;
    //             this->name += (" " + std::to_string(id));
    //         }

    //         std::vector<util::int_t> Active_Solvers;

    //         // void set(toml::handler& _h) {
    //         //     _h.getAtPath(Active_Solvers, "elmer.equation.Active_Solvers", true);
    //         // }

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Name",           this->name);
    //             varToString(_buf, "Active Solvers", Active_Solvers);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer material section
    // */
    // class Material : public Base {
    //     public:
    //         Material() : Base() {
    //             this->name = "Material";
    //             sep = " = ";
    //             id = 1;
    //             this->name += (" " + std::to_string(id));
    //         }

    //         double_t Poisson_ratio, Heat_Capacity, Density, Youngs_modulus, Heat_expansion_Coefficient, Sound_speed, Heat_Conductivity;

    //         void set(toml::handler& _h) {
    //             _h.getAtPath(Poisson_ratio,               "elmer.material.Poisson_ratio", true);
    //             _h.getAtPath(Heat_Capacity,               "elmer.material.Heat_Capacity", true);
    //             _h.getAtPath(Density,                     "elmer.material.Density", true);
    //             _h.getAtPath(Youngs_modulus,              "elmer.material.Youngs_modulus", true);
    //             _h.getAtPath(Heat_expansion_Coefficient,  "elmer.material.Heat_expansion_Coefficient", true);
    //             _h.getAtPath(Sound_speed,                 "elmer.material.Sound_speed", true);
    //             _h.getAtPath(Heat_Conductivity,           "elmer.material.Heat_Conductivity", true);
    //         }

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Name",                       this->name);
    //             varToString(_buf, "Poisson ratio",              Poisson_ratio);
    //             varToString(_buf, "Heat Capacity",              Heat_Capacity);
    //             varToString(_buf, "Density",                    Density);
    //             varToString(_buf, "Youngs modulus",             Youngs_modulus);
    //             varToString(_buf, "Heat expansion Coefficient", Heat_expansion_Coefficient);
    //             varToString(_buf, "Sound speed",                Sound_speed);
    //             varToString(_buf, "Heat Conductivity",          Heat_Conductivity);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer body section
    // */
    // class Body : public Base {
    //     public:
    //         Body(util::int_t _id) : Base() {
    //             this->name = "Body " + std::to_string(_id);
    //             sep = " = ";
    //             id = _id;
    //             Equation = Material = 1;
    //             Initial_condition = Body_Force = toml::noInt;
    //         }
    //         std::vector<util::int_t> Target_Bodies;
    //         util::int_t Initial_condition, Equation, Material, Body_Force;

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Target Bodies",      this->Target_Bodies, true);
    //             // varToString(_buf, "Name",               this->name);
    //             varToString(_buf, "Equation",           Equation);
    //             varToString(_buf, "Material",           Material);
    //             varToString(_buf, "Initial condition",  Initial_condition);
    //             varToString(_buf, "Body Force",         Body_Force);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // class Body_Force : public Base {
    //     public:
    //         Body_Force(util::int_t _id) {
    //             this->name = "Body Force " + std::to_string(_id);
    //             sep = " = ";
    //             id = _id;
    //             Stress_Bodyforce_1 = toml::noDouble;
    //             Stress_Bodyforce_2 = toml::noDouble;
    //             Stress_Bodyforce_3 = toml::noDouble;
    //             Density = toml::noDouble;
    //         }
    //         double_t Stress_Bodyforce_1, Stress_Bodyforce_2, Stress_Bodyforce_3, Density;

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Name",               this->name);
                
    //             // custom lines
    //             if (Stress_Bodyforce_1 != toml::noDouble) {
    //                 _buf << _tab << "Stress Bodyforce 1" << sep << "$ " << Stress_Bodyforce_1;
    //                 if (Density != toml::noDouble)
    //                     _buf << "*" << Density << "\n";
    //                 else
    //                     _buf << "*1.0" << "\n";
    //             }
    //             if (Stress_Bodyforce_2 != toml::noDouble) {
    //                 _buf << _tab << "Stress Bodyforce 2" << sep << "$ " << Stress_Bodyforce_2;
    //                 if (Density != toml::noDouble)
    //                     _buf << "*" << Density << "\n";
    //                 else
    //                     _buf << "*1.0" << "\n";
    //             }
    //             if (Stress_Bodyforce_3 != toml::noDouble) {
    //                 _buf << _tab << "Stress Bodyforce 3" << sep << "$ " << Stress_Bodyforce_3;
    //                 if (Density != toml::noDouble)
    //                     _buf << "*" << Density << "\n";
    //                 else
    //                     _buf << "*1.0" << "\n";
    //             }

    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer initial condition section
    // */
    // class Initial_Condition : public Base {
    //     public:
    //         Initial_Condition(util::int_t _id) : Base() {
    //             this->name = "Initial Condition " + std::to_string(_id);
    //             sep = " = ";
    //             id = _id;
    //             Temperature = toml::noDouble;
    //         }

    //         double_t Temperature, Velocity_1, Velocity_2, Velocity_3;

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Name",        this->name);
    //             varToString(_buf, "Temperature", Temperature);
    //             varToString(_buf, "Velocity 1",    Velocity_1);
    //             varToString(_buf, "Velocity 2",    Velocity_2);
    //             varToString(_buf, "Velocity 2",    Velocity_3);
    //             _buf<<(_end + "\n");
    //         }
    // };

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * handles elmer boundary condition section
    // */
    // class Boundary_Condition : public Base {
    //     public:
    //         Boundary_Condition(util::int_t _id) : Base() {
    //             this->name = "Boundary Condition " + std::to_string(_id);
    //             sep = " = ";
    //             id = _id;
    //             Heat_Flux_BC = toml::noBool;
    //             Heat_Flux    = toml::noDouble;
    //         }

    //         std::vector<util::int_t> Target_Boundaries;
    //         util::bool_t Heat_Flux_BC;
    //         double_t Heat_Flux;
    //         std::vector<double_t> stress_6_vector;

    //         void join(util::oFile& _buf) {
    //             _buf << (name + "\n");
    //             varToString(_buf, "Target Boundaries",  Target_Boundaries, true);
    //             varToString(_buf, "Name",               this->name);
    //             varToString(_buf, "Heat Flux BC",       Heat_Flux_BC);
    //             varToString(_buf, "Heat Flux",          Heat_Flux);
    //             varToString(_buf, "Stress",             stress_6_vector, true);
    //             _buf<<(_end + "\n");
    //         }
    // };
} // namespace elmer 


#endif