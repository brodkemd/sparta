#ifndef ELMER_H
#define ELMER_H

#include "elmer_classes.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    static const char SEP = '\\';
#else
    static const char SEP = '/';
#endif

namespace elmer {
    const unsigned int boundary_size = 4; // number of elements in a boundary element including id
    const unsigned int dimension     = 3;

    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            const std::string sep = "\n\n";
            std::ostringstream double_converter;

        public:
            std::vector<std::vector<int>> element_data;
            std::vector<double> node_temperature_data, node_delta_x_data, node_delta_y_data, node_delta_z_data;
            std::vector<std::array<int,    elmer::boundary_size>> boundary_data;
            std::vector<std::array<double, elmer::dimension>>     node_data;
            std::string name, exe, sif, meshDB, node_data_file, dump_directory, boundary_file, element_file, node_file, node_temperature_file_ext, node_position_file_ext;
            double base_temp;

            Header          header;
            Simulation      simulation;
            Constants       constants;
            Thermal_Solver  thermal_solver;
            Elastic_Solver  elastic_solver;
            Equation        equation;
            Material        material;
            std::vector<Body*>               bodys;
            std::vector<Initial_Condition*>  initial_conditions;
            std::vector<Boundary_Condition*> boundary_conditions;

            Elmer() {
                this->header         = Header();
                this->simulation     = Simulation();
                this->constants      = Constants();
                this->equation       = Equation();
                this->material       = Material();
                this->thermal_solver = Thermal_Solver();
                this->elastic_solver = Elastic_Solver();

                // setting up the double converter
                double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10+2);
            }

            void set(toml::handler& _h) {
                // setting needed variables
                _h.get_at_path(this->meshDB,                "elmer.meshDB",     true);
                _h.get_at_path(this->exe,                   "elmer.exe", true);
                _h.get_at_path(this->sif,                   "elmer.sif", true);
                _h.get_at_path(this->base_temp,             "elmer.base_temp", true);
                _h.get_at_path(this->node_data_file,        "elmer.node_data_file", true);
                _h.get_at_path(this->dump_directory,        "elmer.dump_directory", true);
                _h.get_at_path(this->node_temperature_file_ext, "elmer.node_temperature_file_ext", true);
                _h.get_at_path(this->node_position_file_ext,    "elmer.node_position_file_ext", true);

                // each elmer section sets its own variables, so passing it the data structure
                this->header.set(_h);
                this->simulation.set(_h);
                this->constants.set(_h);
                this->equation.set(_h);
                this->material.set(_h);
                this->thermal_solver.set(_h);
                this->elastic_solver.set(_h);

                this->checks();
            }

            void checks() {
                struct stat sb;
                // making sure the elmer exe 
                if (!(stat(this->exe.c_str(),  &sb) == 0))
                    util::error("Illegal fix fea command, exe path does not exist");
                
                // making sure the mesh database is complete
                std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    if (!(stat((this->meshDB + SEP + "mesh." + exts[i]).c_str(),  &sb) == 0))
                        util::error("Illegal fix fea command, mesh database incomplete, " + (this->meshDB + SEP + "mesh." + exts[i]) + " does not exist");
                }

                // making sure base temperature is valid
                if (this->base_temp <= (double)0)
                    util::error("base temperature must be greater than 0");

                if (this->dump_directory.back() == SEP)
                    this->dump_directory = this->dump_directory.substr(0, this->dump_directory.length()-1);
            }

            /* ---------------------------------------------------------------------- */

            void setupDumpDirectory() {
                // making sure the mesh database is complete
                std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    util::copyFile(this->meshDB + "/mesh." + exts[i], this->dump_directory + "/mesh." + exts[i]);
                }

                this->boundary_file = this->dump_directory + SEP + "mesh.boundary";
                this->element_file = this->dump_directory + SEP + "mesh.elements";
                this->node_file = this->dump_directory + SEP + "mesh.nodes";
            }

            
            void setupNodeDataFileAndVectors() {
                // getting the number of nodes
                int count = util::count_lines_in_file(this->node_data_file);

                // setting up object to convert doubles
                this->double_converter.str("");
                this->double_converter << this->base_temp;

                // getting string version of elmer.base_temp
                std::string str = this->double_converter.str();

                // writing base temperature to file for each node
                std::ofstream output(this->node_data_file);
                if (!(output.is_open()))
                    util::error("could not open node data file");
                
                output << "Time:      1\n";
                output << "temperature\n";
                output << "Perm:      \n";
                for (int i = 0; i < count; i++) {
                    output << "   " << str << "\n";
                    this->node_temperature_data.push_back(this->base_temp);
                }

                output << "displacement 1\nPerm:    \n";
                this->double_converter.str("");
                this->double_converter << (double)0;
                str = double_converter.str();
                for (int i = 0; i < count; i++) {
                    output << "   " << str << "\n";
                    this->node_delta_x_data.push_back((double)0);
                }

                output << "displacement 2\nPerm:    \n";
                for (int i = 0; i < count; i++) {
                    output << "   " << str << "\n";
                    this->node_delta_y_data.push_back((double)0);
                }

                output << "displacement 3\nPerm:    \n";
                for (int i = 0; i < count; i++) {
                    output << "   " << str << "\n";
                    this->node_delta_z_data.push_back((double)0);
                }
                output.close();
                this->double_converter.str("");
            }


            void setup() {
                this->setupDumpDirectory();
                this->setupNodeDataFileAndVectors();
                this->loadBoundaryData();
                this->loadElements();
            }


            void run() {
                // running the command
                // must have " 2>&1" at end to pipe stderr to stdout
                util::CommandResult command_result;
                try {
                    command_result = util::EXEC(exe+" "+this->sif+" 2>&1");
                } catch (std::exception& e) {
                    util::error(e.what());
                }
                // if the command did not succeed
                if (command_result.exitstatus)
                    util::error(command_result.output);
                
            }


            void loadBoundaryData() {
                std::vector<std::string> lines, _split;
                std::array<int, boundary_size> arr;
                this->boundary_data.clear();

                util::readFile(this->boundary_file, lines);

                int count = 1;
                for (std::string it : lines) {
                    util::trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    util::split(it, _split, ' ');

                    if (_split[4] != (std::string)"303") 
                        util::error("element is not a triangle in boundary file at line: " + std::to_string(count));

                    // getting rid of stuff I do not need
                    _split.erase(_split.begin()+1, _split.begin() + 5);

                    // catching errors
                    if (_split.size() != boundary_size)
                        util::error("too many boundary elements to assign at line: " + std::to_string(count));

                    // adding the data
                    for (unsigned int i = 0; i < boundary_size; i++)
                        arr[i] = std::stoi(_split[i]);

                    // adding the data to the class vector
                    this->boundary_data.push_back(arr);
                    count++;
                }
            }


            void loadNodeData() {
                this->node_temperature_data.clear();
                this->node_delta_x_data.clear();
                this->node_delta_y_data.clear();
                this->node_delta_z_data.clear();

                long unsigned int i, start;
                std::vector<double>* v;
                std::vector<std::string> lines, split_line;
                std::string cur_key, line;
                
                util::readFile(this->node_data_file, lines);
                start = 0;
                for (i = 0; i < lines.size(); i++) {
                    util::split(lines[i], split_line, ' ');
                    if (split_line[0] == (std::string)"Time:")
                        start = i;
                }

                for (i = start+1; i < lines.size(); i++) {
                    line = lines[i];
                    util::trim(line);
                    util::split(line, split_line, ' ');
                    if (lines[i][0] == ' ' && split_line.size() == 1)
                        v->push_back(std::stod(line));
                    else if ("Perm:" == split_line[0]) {
                        cur_key = lines[i-1];
                        util::trim(cur_key);
                        if (cur_key == (std::string)"temperature")
                            v = &this->node_temperature_data;
                        else if (cur_key == (std::string)"displacement 1")
                            v = &this->node_delta_x_data;
                        else if (cur_key == (std::string)"displacement 2")
                            v = &this->node_delta_y_data;
                        else if (cur_key == (std::string)"displacement 3")
                            v = &this->node_delta_z_data;
                        else
                            util::error("got unknown key in data file: " + cur_key);
                    }
                }
            }


            void loadElements() {
                std::vector<std::string> _temp_data, _split;
                std::vector<int> arr;
                this->element_data.clear();

                util::readFile(this->element_file, _temp_data);

                this->element_data.resize(_temp_data.size());

                for (std::string it : _temp_data) {
                    arr.clear();

                    util::trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    util::split(it, _split, ' ');

                    // getting rid of stuff I do not need
                    _split.erase(_split.begin() + 1, _split.begin() + 3);

                    // adding the data
                    for (long unsigned int i = 1; i < _split.size(); i++)
                        arr.push_back(std::stoi(_split[i]));

                    // adding the data to the class vector
                    this->element_data[std::stoi(_split[0]) - 1] = arr;
                }
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
                for (long unsigned int i = 0; i < this->element_data.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    _size = this->element_data[i].size();
                    for (int j = 0; j < _size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        if (this->element_data[i][j]-1 >= size)
                            util::error("index out of bounds in creating initial conditions");
                        avg += (this->node_temperature_data[this->element_data[i][j]-1]/((double)_size));
                    }

                    // surface element
                    Body* _body = new Body(i+1);
                    _body->Initial_condition = i+1;
                    _body->Target_Bodies = {(int)i+1};
                    Initial_Condition* _ic = new Initial_Condition(i+1);
                    _ic->Temperature = avg;

                    this->initial_conditions.push_back(_ic);
                    this->bodys.push_back(_body);
                }
            }


            void makeSif() {
                std::ofstream _buf(this->sif);
                if (!(_buf.is_open())) util::error(this->sif + " did not open");
                this->join(_buf);                
                _buf.close();
            }


            // void loadNodeDataFromPostFile() {
            //     std::vector<std::string> data, v;
            //     std::array<double, dimension> arr;
            //     std::string lines;
            //     this->node_data.clear();
    
            //     std::ifstream f(this->node_data_file);
                
            //     if (f.is_open()) {
            //         std::ostringstream ss;
            //         ss << f.rdbuf();
            //         lines = ss.str();
            //     }

            //     f.close();

            //     int start = lines.find('\n', lines.find('#')+1)+1;
            //     int end   = lines.find('#', lines.find('#')+1)-1;
            //     lines = lines.substr(start,  end - start);

            //     split(lines, data, '\n');

            //     long unsigned int i, j;
            //     for (i = 0; i < data.size(); i++) {
            //         // std::cout << data[i] << "\n";
            //         trim(data[i]); split(data[i], v, ' ');
            //         if (v.size() != dimension)
            //             error("node element is not correct size at element: " + std::to_string(v.size()) + ", has size: " + std::to_string(v.size()));

            //         for (j = 0; j < dimension; j++)
            //             arr[j] = std::stod(v[j]);

            //         this->node_data.push_back(arr);
            //     }
            // }


            void small_join(std::ofstream& _buf) {
                header.join(_buf);
                simulation.join(_buf);
                constants.join(_buf);
                thermal_solver.join(_buf);
                elastic_solver.join(_buf);
                equation.join(_buf);
                material.join(_buf);
            }


            void join_conditions(std::ofstream& _buf) {
                for (long unsigned int i = 0; i < this->bodys.size(); i++)
                    this->bodys[i]->join(_buf);

                for (long unsigned int i = 0; i < this->initial_conditions.size(); i++)
                    this->initial_conditions[i]->join(_buf);

                for (long unsigned int i = 0; i < this->boundary_conditions.size(); i++)
                    this->boundary_conditions[i]->join(_buf);
            }


            void join(std::ofstream& _buf) {
                small_join(_buf);
                join_conditions(_buf);
            }


            void dump(long _timestep) {
                this->dumpNodeTemperatures(_timestep);
                this->dumpNodePositions(_timestep);
            }


            void dumpNodeTemperatures(long _timestep) {
                std::string filename = this->dump_directory + std::to_string(_timestep) + "." + this->node_temperature_file_ext;
                std::ofstream out(filename);
                if (!(out.is_open())) util::error(filename + " did not open");

                for (double it : this->node_temperature_data) {
                    this->double_converter.str("");
                    this->double_converter << it;
                    out << this->double_converter.str() << "\n";
                }
                out.close();
            }


            void dumpNodePositions(long _timestep) {
                std::string filename = this->dump_directory + std::to_string(_timestep) + "." + this->node_position_file_ext;
                std::ofstream out(filename);
                if (!(out.is_open())) util::error(filename + " did not open");

                if (this->node_delta_x_data.size() != this->node_delta_y_data.size() || this->node_delta_y_data.size() != this->node_delta_z_data.size() || this->node_delta_x_data.size() != this->node_delta_z_data.size())
                    util::error("position delta vectors do not have the same size");

                for (long unsigned int i = 0; i < this->node_delta_x_data.size(); i++) {
                    this->double_converter.str("");
                    this->double_converter << this->node_data[i][0] << " " << this->node_data[i][1] << " " << this->node_data[i][2];
                    out << this->double_converter.str() << "\n";
                }
                out.close();

            }

    };
}


#endif