#ifndef ELMER_H
#define ELMER_H

#include <fstream>
#include <array>
#include "toml.hpp"

namespace elmer {
    const unsigned int boundary_size = 4; // number of elements in a boundary element including id
    const std::string indicator      = "Perm:"; // indicator for start of new timestep data section in data file

    // error command
    // NOTE: all functions that using this command must be wrapped in a try-catch statement that catches strings
    void error(std::string _msg) { throw _msg; }

    // used in the EXEC function below, represents data returned from a system command
    struct CommandResult {
        std::string output;
        int exitstatus;
    };


    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return commandResult containing STDOUT (not stderr) output & exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    static CommandResult EXEC(const std::string &command) {
        int exitcode = 0;
        std::array<char, 1048576> buffer {};
        std::string result;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result += std::string(buffer.data(), bytesread);
            }
        } catch (...) {
            pclose(pipe);
            error("unhandled exception occured");
        }
        exitcode = WEXITSTATUS(pclose(pipe));
        return CommandResult{result, exitcode};
    }


    void ltrim(std::string& _s) {
        int j = 0;
        for (j; j < _s.length(); j++) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(j, _s.length() - j);
    }


    void rtrim(std::string& _s) {
        int j = _s.length()-1;
        for (j; j >= 0; j--) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(0,  j+1);
    }


    void trim(std::string& _s) { rtrim(_s); ltrim(_s); }


    void split(std::string& _s, std::vector<std::string>& _v) {
        _v.clear();
        int last = 0, j = 0;
        for (j; j <= _s.length(); j++) {
            if ((_s[j] == ' ' || j == _s.length())) {
                if (j == last) { last++; continue; }
                _v.push_back(_s.substr(last, j-last));
                last = j+1;
            }
        }
    }


    /**
     * erases the file at the path inputted
    */ 
    void eraseFile(std::string filename) {
        std::ofstream ofs;
        ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
        if (!(ofs.is_open())) error(filename + " did not open");
        ofs.close();
    }

    /**
     * writes the inputted string to file at the path inputted
    */ 
    void writeFile(std::string filename, std::string& lines) {
        std::ofstream out(filename);
        if (!(out.is_open())) error(filename + " did not open");
        out << lines;
        out.close();
    }

    /**
     * counts the number of lines in the file at the path inputted
    */ 
    int count_lines_in_file(std::string _file) {
        int number_of_lines = 0;
        std::string line;
        std::ifstream myfile(_file);

        while (std::getline(myfile, line))
            ++number_of_lines;
        return number_of_lines;
    }

    /**
     * reads the contents of the file at the path inputted into the inputted
     * vector, this function is designed to read the data file
    */ 
    void readDataFile(std::string filename, std::vector<std::string>& lines) {
        std::ifstream myfile(filename);
        if (!(myfile.is_open())) error(filename + " did not open");

        std::string line;
        while (std::getline(myfile, line)) {
            if (line.length() >= indicator.length()) {
                // resets the vector if it detects a new indicator, this allows it to get only the new data
                if (line.substr(0, indicator.length()) == indicator){
                    lines.clear();
                    continue;
                }
            }
            lines.push_back(line + "\n");
        }
        myfile.close();
    }

    /**
     * reads the contents of the file at the path inputted into the inputted
     * vector
    */ 
    void readFile(std::string fileName, std::vector<std::string>& lines) {
        std::ifstream myfile(fileName);
        std::string line;
        while (std::getline(myfile, line)) {
            lines.push_back(line + "\n");
        }
        myfile.close();
    }

    /**
     * Gets the latest data for each node from the inputted file
    */
    void getLatestNodeData(std::string filename, std::vector<double>& data) {
        std::vector<std::string> lines;

        data.clear();

        // reads the data file
        readDataFile(filename, lines);

        for (std::string it : lines) {
            trim(it);
            // if the line only has one thing in it
            if (it.find(" ") == std::string::npos) {
                // convert to double and add to list
                data.push_back(std::stod(it));
            }
        }
    }

    /**
     * Loads the boundary data from the inputted file into the vector
    */
    void getBoundaryData(std::string filename, std::vector<std::array<int, boundary_size>>& data) {
        std::vector<std::string> lines, _split;
        std::array<int, boundary_size> arr;
        data.clear();

        readFile(filename, lines);

        int count = 1;
        for (std::string it : lines) {
            trim(it);
            if (it.length() == 0) continue;

            // splitting the line at spaces
            split(it, _split);

            if (_split[4] != (std::string)"303") 
                error("element is not a triangle in boundary file at line: " + std::to_string(count));

            // getting rid of stuff I do not need
            _split.erase(_split.begin()+1, _split.begin() + 5);

            // catching errors
            if (_split.size() != boundary_size)
                error("too many boundary elements to assign at line: " + std::to_string(count));

            // adding the data
            for (int i = 0; i < boundary_size; i++)
                arr[i] = std::stoi(_split[i]);

            // adding the data to the class vector
            data.push_back(arr);
            count++;
        }
    }

    /**
     * Base class that all elmer sections inherit from
    */
    class Base {
        protected:
            std::string name, sep;

            // tab in the file
            std::string _tab = "  ";

            // ends the section
            std::string _end = "End";
            int id;

            std::ostringstream double_converter;

            Base() {
                // setting up the double converter
                double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
            }

        private:
            /**
             * these functions convert the inputted var into a string
            */
            void _varToString(std::string& _buf, std::string _var) {
                 _buf = "\""+_var+"\"";
            }

            void _varToString(std::string& _buf, double _var) {
                double_converter.str("");
                double_converter << _var;
                _buf = double_converter.str();
            }

            void _varToString(std::string& _buf, int _var) {
                _buf = std::to_string(_var);
            }

            void _varToString(std::string& _buf, bool _var) {
                if (_var)
                    _buf = "True";
                else
                    _buf = "False";
            }

        public:
            // set function, this is overriden
            virtual void set(toml::handler& _h) { return; }

            /**
             * these functions are outward facing and they build elmer file
             * lines from the inputs
            */ 
            void varToString(std::string& _buf, std::string _name, std::string _var) {
                _buf = (_tab + _name + sep + "\""+_var+"\"" + "\n");
            }

            void varToString(std::string& _buf, std::string _name, double _var) {
                double_converter.str("");
                double_converter << _var;
                _buf = (_tab + _name + sep + double_converter.str() + "\n");
            }

            void varToString(std::string& _buf, std::string _name, int _var) {
                _buf = (_tab + _name + sep + std::to_string(_var) + "\n");
            }

            void varToString(std::string& _buf, std::string _name, bool _var) {
                if (_var)
                    _buf = (_tab + _name + sep + "True" + "\n");
                else
                    _buf = (_tab + _name + sep + "False" + "\n");
            }

            // for vectors
            template<typename T>
            void varToString(std::string& _buf, std::string _name, std::vector<T> _var, bool include_count = true) {

                if (_var.size() == 0) return;

                std::string _temp;
                if (include_count)
                    _buf = (_tab + _name + "(" + std::to_string(_var.size()) + ")" + sep);
                else 
                    _buf = (_tab + _name + sep);


                for (int i = 0; i < (_var.size() - 1); i++){
                    _varToString(_temp, _var[i]);
                    _buf+=(_temp + " ");
                }
                _varToString(_temp, _var[_var.size() - 1]);
                _buf += (_temp + "\n");
            }
    };

    /**
     * handles elmer header section
    */
    class Header : public Base {
        public:
            Header() : Base() {
                this->name = "Header";
                sep = " ";
            }

            std::vector<std::string> Mesh_DB;
            std::string Include_Path, Results_Directory;

            void set(toml::handler& _h) {
                _h.get_at_path(Mesh_DB, "elmer.header.Mesh_DB");
                _h.get_at_path(Include_Path, "elmer.header.Include_Path");
                _h.get_at_path(Results_Directory, "elmer.header.Results_Directory");
            }

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Mesh DB", Mesh_DB, false); _buf+=_temp;
                varToString(_temp, "Include Path", Include_Path); _buf+=_temp;
                varToString(_temp, "Results Directory", Results_Directory); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer constants section
    */
    class Constants : public Base {
        public:
            Constants() : Base() {
                this->name = "Constants";
                sep = " = ";
            }

            std::vector<double> Gravity;
            double Stefan_Boltzmann, Permittivity_of_Vacuum, Permeability_of_Vacuum, Boltzmann_Constant, Unit_Charge;


            void set(toml::handler& _h) {
                _h.get_at_path(Gravity, "elmer.constants.Gravity");
                _h.get_at_path(Stefan_Boltzmann, "elmer.constants.Stefan_Boltzmann");
                _h.get_at_path(Permittivity_of_Vacuum, "elmer.constants.Permittivity_of_Vacuum");
                _h.get_at_path(Permeability_of_Vacuum, "elmer.constants.Permeability_of_Vacuum");
                _h.get_at_path(Boltzmann_Constant, "elmer.constants.Boltzmann_Constant");
                _h.get_at_path(Unit_Charge, "elmer.constants.Unit_Charge");
            }


            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Gravity", Gravity); _buf+=_temp;
                varToString(_temp, "Stefan Boltzmann", Stefan_Boltzmann); _buf+=_temp;
                varToString(_temp, "Permittivity of Vacuum", Permittivity_of_Vacuum); _buf+=_temp;
                varToString(_temp, "Permeability of Vacuum", Permeability_of_Vacuum); _buf+=_temp;
                varToString(_temp, "Boltzmann Constant", Boltzmann_Constant); _buf+=_temp;
                varToString(_temp, "Unit Charge", Unit_Charge); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer simulation section
    */
    class Simulation : public Base {
        public:
            Simulation() : Base() {
                this->name = "Simulation";
                sep = " = ";
            }

            int Max_Output_Level, Steady_State_Max_Iterations, BDF_Order;
            std::string Coordinate_System, Simulation_Type, Timestepping_Method, Solver_Input_File, Output_File,Post_File;
            std::vector<int> Coordinate_Mapping, Output_Intervals, Timestep_intervals;
            std::vector<double> Timestep_Sizes;
            

            void set(toml::handler& _h) {
                _h.get_at_path(Max_Output_Level, "elmer.simulation.Max_Output_Level");
                _h.get_at_path(Coordinate_System, "elmer.simulation.Coordinate_System");
                _h.get_at_path(Coordinate_Mapping, "elmer.simulation.Coordinate_Mapping");
                _h.get_at_path(Simulation_Type, "elmer.simulation.Simulation_Type");
                _h.get_at_path(Steady_State_Max_Iterations, "elmer.simulation.Steady_State_Max_Iterations");
                _h.get_at_path(Output_Intervals, "elmer.simulation.Output_Intervals");
                _h.get_at_path(Timestep_intervals, "elmer.simulation.Timestep_intervals");
                _h.get_at_path(Timestep_Sizes, "elmer.simulation.Timestep_Sizes");
                _h.get_at_path(Timestepping_Method, "elmer.simulation.Timestepping_Method");
                _h.get_at_path(BDF_Order, "elmer.simulation.BDF_Order");
                _h.get_at_path(Solver_Input_File, "elmer.simulation.Solver_Input_File");
                _h.get_at_path(Output_File, "elmer.simulation.Output_File");
                _h.get_at_path(Post_File, "elmer.simulation.Post_File");
            }


            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Max Output Level", Max_Output_Level); _buf+=_temp;
                varToString(_temp, "Coordinate System", Coordinate_System); _buf+=_temp;
                varToString(_temp, "Coordinate Mapping", Coordinate_Mapping); _buf+=_temp;
                varToString(_temp, "Simulation Type", Simulation_Type); _buf+=_temp;
                varToString(_temp, "Steady State Max Iterations", Steady_State_Max_Iterations); _buf+=_temp;
                varToString(_temp, "Output Intervals", Output_Intervals); _buf+=_temp;
                varToString(_temp, "Timestep intervals", Timestep_intervals); _buf+=_temp;
                varToString(_temp, "Timestep Sizes", Timestep_Sizes); _buf+=_temp;
                varToString(_temp, "Timestepping Method", Timestepping_Method); _buf+=_temp;
                varToString(_temp, "BDF Order", BDF_Order); _buf+=_temp;
                varToString(_temp, "Solver Input File", Solver_Input_File); _buf+=_temp;
                varToString(_temp, "Output File", Output_File); _buf+=_temp;
                varToString(_temp, "Post File", Post_File); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer solver section
    */
    class Solver : public Base {
        public:
            Solver() : Base() {
                this->name = "Solver";
                sep = " = ";
                id = 1;
                this->name += (" " + std::to_string(id));
            }

            std::vector<std::string> Procedure;
            std::string Equation, Variable, Exec_Solver, Linear_System_Solver, Linear_System_Iterative_Method, Linear_System_Preconditioning;
            int Linear_System_Residual_Output, Nonlinear_System_Max_Iterations, Nonlinear_System_Newton_After_Iterations, Linear_System_Max_Iterations, BiCGstabl_polynomial_degree, Nonlinear_System_Relaxation_Factor;
            bool Stabilize, Optimize_Bandwidth, Linear_System_Abort_Not_Converged;
            double Steady_State_Convergence_Tolerance, Nonlinear_System_Convergence_Tolerance, Nonlinear_System_Newton_After_Tolerance, Linear_System_Convergence_Tolerance, Linear_System_ILUT_Tolerance;


            void set(toml::handler& _h) {
                _h.get_at_path(Equation, "elmer.solver.Equation");
                _h.get_at_path(Procedure, "elmer.solver.Procedure");
                _h.get_at_path(Variable, "elmer.solver.Variable");
                _h.get_at_path(Exec_Solver, "elmer.solver.Exec_Solver");
                _h.get_at_path(Stabilize, "elmer.solver.Stabilize");
                _h.get_at_path(Optimize_Bandwidth, "elmer.solver.Optimize_Bandwidth");
                _h.get_at_path(Steady_State_Convergence_Tolerance, "elmer.solver.Steady_State_Convergence_Tolerance");
                _h.get_at_path(Nonlinear_System_Convergence_Tolerance, "elmer.solver.Nonlinear_System_Convergence_Tolerance");
                _h.get_at_path(Nonlinear_System_Max_Iterations, "elmer.solver.Nonlinear_System_Max_Iterations");
                _h.get_at_path(Nonlinear_System_Newton_After_Iterations, "elmer.solver.Nonlinear_System_Newton_After_Iterations");
                _h.get_at_path(Nonlinear_System_Newton_After_Tolerance, "elmer.solver.Nonlinear_System_Newton_After_Tolerance");
                _h.get_at_path(Nonlinear_System_Relaxation_Factor, "elmer.solver.Nonlinear_System_Relaxation_Factor");
                _h.get_at_path(Linear_System_Solver, "elmer.solver.Linear_System_Solver");
                _h.get_at_path(Linear_System_Iterative_Method, "elmer.solver.Linear_System_Iterative_Method");
                _h.get_at_path(Linear_System_Max_Iterations, "elmer.solver.Linear_System_Max_Iterations");
                _h.get_at_path(Linear_System_Convergence_Tolerance, "elmer.solver.Linear_System_Convergence_Tolerance");
                _h.get_at_path(BiCGstabl_polynomial_degree, "elmer.solver.BiCGstabl_polynomial_degree");
                _h.get_at_path(Linear_System_Preconditioning, "elmer.solver.Linear_System_Preconditioning");
                _h.get_at_path(Linear_System_ILUT_Tolerance, "elmer.solver.Linear_System_ILUT_Tolerance");
                _h.get_at_path(Linear_System_Abort_Not_Converged, "elmer.solver.Linear_System_Abort_Not_Converged");
                _h.get_at_path(Linear_System_Residual_Output, "elmer.solver.Linear_System_Residual_Output");
            }


            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Equation", Equation); _buf+=_temp;
                varToString(_temp, "Procedure", Procedure); _buf+=_temp;
                varToString(_temp, "Variable", Variable); _buf+=_temp;
                varToString(_temp, "Exec Solver", Exec_Solver); _buf+=_temp;
                varToString(_temp, "Stabilize", Stabilize); _buf+=_temp;
                varToString(_temp, "Optimize Bandwidth", Optimize_Bandwidth); _buf+=_temp;
                varToString(_temp, "Steady State Convergence Tolerance", Steady_State_Convergence_Tolerance); _buf+=_temp;
                varToString(_temp, "Nonlinear System Convergence Tolerance", Nonlinear_System_Convergence_Tolerance); _buf+=_temp;
                varToString(_temp, "Nonlinear System Max Iterations", Nonlinear_System_Max_Iterations); _buf+=_temp;
                varToString(_temp, "Nonlinear System Newton After Iterations", Nonlinear_System_Newton_After_Iterations); _buf+=_temp;
                varToString(_temp, "Nonlinear System Newton After Tolerance", Nonlinear_System_Newton_After_Tolerance); _buf+=_temp;
                varToString(_temp, "Nonlinear System Relaxation Factor", Nonlinear_System_Relaxation_Factor); _buf+=_temp;
                varToString(_temp, "Linear System Solver", Linear_System_Solver); _buf+=_temp;
                varToString(_temp, "Linear System Iterative Method", Linear_System_Iterative_Method); _buf+=_temp;
                varToString(_temp, "Linear System Max Iterations", Linear_System_Max_Iterations); _buf+=_temp;
                varToString(_temp, "Linear System Convergence Tolerance", Linear_System_Convergence_Tolerance); _buf+=_temp;
                varToString(_temp, "BiCGstabl polynomial degree", BiCGstabl_polynomial_degree); _buf+=_temp;
                varToString(_temp, "Linear System Preconditioning", Linear_System_Preconditioning); _buf+=_temp;
                varToString(_temp, "Linear System ILUT Tolerance", Linear_System_ILUT_Tolerance); _buf+=_temp;
                varToString(_temp, "Linear System Abort Not Converged", Linear_System_Abort_Not_Converged); _buf+=_temp;
                varToString(_temp, "Linear System Residual Output", Linear_System_Residual_Output); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer equation section
    */
    class Equation : public Base {
        public:
            Equation() : Base() {
                this->name = "Equation";
                sep = " = ";
                id = 1;
                this->name += (" " + std::to_string(id));
            }

            std::vector<int> Active_Solvers;

            void set(toml::handler& _h) {
                _h.get_at_path(Active_Solvers, "elmer.equation.Active_Solvers");
            }

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Name", this->name); _buf+=_temp;
                varToString(_temp, "Active Solvers", Active_Solvers); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer material section
    */
    class Material : public Base {
        public:
            Material() : Base() {
                this->name = "Material";
                sep = " = ";
                id = 1;
                this->name += (" " + std::to_string(id));
            }

            double Poisson_ratio, Heat_Capacity, Density, Youngs_modulus, Heat_expansion_Coefficient, Sound_speed, Heat_Conductivity;

            void set(toml::handler& _h) {
                _h.get_at_path(Poisson_ratio, "elmer.material.Poisson_ratio");
                _h.get_at_path(Heat_Capacity, "elmer.material.Heat_Capacity");
                _h.get_at_path(Density, "elmer.material.Density");
                _h.get_at_path(Youngs_modulus, "elmer.material.Youngs_modulus");
                _h.get_at_path(Heat_expansion_Coefficient, "elmer.material.Heat_expansion_Coefficient");
                _h.get_at_path(Sound_speed, "elmer.material.Sound_speed");
                _h.get_at_path(Heat_Conductivity, "elmer.material.Heat_Conductivity");
            }

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Name", this->name); _buf+=_temp;
                varToString(_temp, "Poisson ratio", Poisson_ratio); _buf+=_temp;
                varToString(_temp, "Heat Capacity", Heat_Capacity); _buf+=_temp;
                varToString(_temp, "Density", Density); _buf+=_temp;
                varToString(_temp, "Youngs modulus", Youngs_modulus); _buf+=_temp;
                varToString(_temp, "Heat expansion Coefficient", Heat_expansion_Coefficient); _buf+=_temp;
                varToString(_temp, "Sound speed", Sound_speed); _buf+=_temp;
                varToString(_temp, "Heat Conductivity", Heat_Conductivity); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer body section
    */
    class Body : public Base {
        public:
            Body(int _id) : Base() {
                this->name = "Body " + std::to_string(_id);
                sep = " = ";
                id = _id;
                Equation = 1; Material = 1;
            }
            std::vector<int> Target_Bodies;
            int Initial_condition, Equation, Material;

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Target Bodies", this->Target_Bodies, true); _buf+=_temp;
                varToString(_temp, "Name", this->name); _buf+=_temp;
                varToString(_temp, "Equation", Equation); _buf+=_temp;
                varToString(_temp, "Material", Material); _buf+=_temp;
                varToString(_temp, "Initial condition", Initial_condition); _buf+=_temp;
                _buf+=(_end + "\n");
            }
    };

    /**
     * handles elmer initial condition section
    */
    class Initial_Condition : public Base {
        public:
            Initial_Condition(int _id) : Base() {
                this->name = "Initial Condition " + std::to_string(_id);
                sep = " = ";
                id = _id;
            }

            double Temperature;

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Name", this->name); _buf+=_temp;
                varToString(_temp, "Temperature", Temperature); _buf+=_temp;
                _buf+=(_end + "\n");
            }

    };

    /**
     * handles elmer boundary condition section
    */
    class Boundary_Condition : public Base {
        public:
            Boundary_Condition(int _id) : Base() {
                this->name = "Boundary Condition " + std::to_string(_id);
                sep = " = ";
                id = _id;
            }

            std::vector<int> Target_Boundaries;
            bool Heat_Flux_BC = true;
            double Heat_Flux;

            void join(std::string& _buf) {
                std::string _temp; _buf = (name + "\n");
                varToString(_temp, "Target Boundaries", Target_Boundaries); _buf+=_temp;
                varToString(_temp, "Name", this->name); _buf+=_temp;
                varToString(_temp, "Heat Flux BC", Heat_Flux_BC); _buf+=_temp;
                varToString(_temp, "Heat Flux", Heat_Flux); _buf+=_temp;
                _buf+=(_end + "\n");
            }

    };

    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            const std::string sep = "\n\n";

        public:
            Elmer() {
                this->header     = Header();
                this->simulation = Simulation();
                this->constants  = Constants();
                this->equation   = Equation();
                this->material   = Material();
                this->solver     = Solver();
            }


            void run() {
                // running the command
                // must have " 2>&1" at end to pipe stderr to stdout
                CommandResult command_result;
                try {
                    command_result = EXEC(exe+" "+this->sif+" 2>&1");
                } catch (std::exception& e) {
                    error(e.what());
                }
                // if the command did not succeed
                if (command_result.exitstatus)
                    error(command_result.output);
                
            }


            void small_join(std::string& _buf) {
                std::string _temp; _buf.clear();
                header.join(_temp); _buf+=_temp;
                simulation.join(_temp); _buf+=_temp;
                constants.join(_temp); _buf+=_temp;
                solver.join(_temp); _buf+=_temp;
                equation.join(_temp); _buf+=_temp;
                material.join(_temp); _buf+=_temp;
            }


            void join_conditions(std::string& _buf) {
                std::string _temp; _buf.clear();

                for (int i = 0; i < this->bodys.size(); i++){
                    this->bodys[i]->join(_temp); _buf+=_temp;
                }

                for (int i = 0; i < this->initial_conditions.size(); i++){
                    this->initial_conditions[i]->join(_temp); _buf+=_temp;
                }

                for (int i = 0; i < this->boundary_conditions.size(); i++){
                    this->boundary_conditions[i]->join(_temp); _buf+=_temp;
                }
            }


            void join(std::string& _buf) {
                std::string _temp; _buf.clear();
                small_join(_temp); _buf+=_temp;
                join_conditions(_temp); _buf+=_temp;

            }


            void loadNodeTemperatureData() {
                // get the temperatures per node from the file and putting them into the elmer variable, elmer needs to know this info
                getLatestNodeData(this->temperature_data_file, this->node_temperature_data);
            }


            void loadBoundaryData() {
                // loading the boundary data into the local variable
                getBoundaryData(this->meshDBstem + ".boundary", this->boundary_data);
            }


            void createInitCondsFromData() {
                // used later
                double avg;
                int _size;
                int size = this->node_temperature_data.size();

                this->initial_conditions.clear();
                this->bodys.clear();
                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (int i = 0; i < this->element_data.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    _size = this->element_data[i].size();
                    for (int j = 0; j < _size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        if (this->element_data[i][j]-1 >= size)
                            error("index out of bounds in creating initial conditions");
                        avg += (this->node_temperature_data[this->element_data[i][j]-1]/((double)_size));
                    }

                    // surface element
                    Body* _body = new Body(i+1);
                    _body->Initial_condition = i+1;
                    _body->Target_Bodies = {i+1};
                    Initial_Condition* _ic = new Initial_Condition(i+1);
                    _ic->Temperature = avg;

                    this->initial_conditions.push_back(_ic);
                    this->bodys.push_back(_body);
                }
            }


            void setupTemperatureDataFile() {
                // getting the number of nodes
                int count = elmer::count_lines_in_file(this->meshDBstem + ".nodes");

                // setting up object to convert doubles
                std::ostringstream double_converter;
                double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << this->base_temp;

                // getting string version of elmer.base_temp
                std::string base_temp_str = double_converter.str();

                // writing base temperature to file for each node
                std::ofstream output(this->temperature_data_file);
                for (int i = 0; i < count; i++) output << base_temp_str << "\n";
                output.close();
            }


            void makeSif() {
                std::string _buf;
                this->join(_buf);
                elmer::writeFile(this->sif, _buf);
            }


            void loadElements() {
                std::vector<std::string> _temp_data, lines, _split;
                std::vector<int> arr;
                this->element_data.clear();

                readFile(this->meshDBstem + ".elements", _temp_data);

                this->element_data.resize(_temp_data.size());

                for (std::string it : _temp_data) {
                    arr.clear();

                    trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    split(it, _split);

                    // getting rid of stuff I do not need
                    _split.erase(_split.begin() + 1, _split.begin() + 3);

                    // adding the data
                    for (int i = 1; i < _split.size(); i++)
                        arr.push_back(std::stoi(_split[i]));

                    // adding the data to the class vector
                    this->element_data[std::stoi(_split[0]) - 1] = arr;
                }
            }


            void set(toml::handler& _h) {
                // setting needed variables
                _h.get_at_path(this->meshDBstem,            "elmer.meshDBstem");
                _h.get_at_path(this->exe,                   "elmer.exe");
                _h.get_at_path(this->sif,                   "elmer.sif");
                _h.get_at_path(this->base_temp,             "elmer.base_temp");
                _h.get_at_path(this->temperature_data_file, "elmer.temperature_data_file");

                // each elmer section sets its own variables, so passing it the data structure
                this->header.set(_h);
                this->simulation.set(_h);
                this->constants.set(_h);
                this->equation.set(_h);
                this->material.set(_h);
                this->solver.set(_h);
            }


            std::vector<std::vector<int>> element_data;
            std::vector<double> node_temperature_data;
            std::vector<std::array<int, elmer::boundary_size>> boundary_data;
            std::string name, exe, sif, meshDBstem, temperature_data_file;
            double base_temp;

            Header     header;
            Simulation simulation;
            Constants  constants;
            Solver     solver;
            Equation   equation;
            Material   material;
            std::vector<Body*>               bodys;
            std::vector<Initial_Condition*>  initial_conditions;
            std::vector<Boundary_Condition*> boundary_conditions;
    };
}


#endif