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
    const long unsigned int node_element_size = 5;

    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            const std::string sep = "\n\n";
            std::ostringstream double_converter;

        public:
            std::vector<double> node_temperature_data, node_delta_x_data, node_delta_y_data, node_delta_z_data;
            std::vector<std::array<int, elmer::boundary_size>> boundaries;
            //                    temperature    x       y       z
            std::vector<std::tuple<double,     double, double, double>>  nodes;
            std::vector<std::vector<int>> elements;
            std::string name, exe, sif, meshDB, node_data_file, dump_directory, boundary_file, element_file, node_file, node_temperature_file_ext, node_position_file_ext;
            double base_temp;

            Header          header;
            Simulation      simulation;
            Constants       constants;
            Thermal_Solver  thermal_solver;
            Elastic_Solver  elastic_solver;
            Equation        equation;
            Material        material;
            std::vector<Body*>               bodies;
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
                _h.get_at_path(this->meshDB,                    "elmer.meshDB",                     true);
                _h.get_at_path(this->exe,                       "elmer.exe",                        true);
                _h.get_at_path(this->sif,                       "elmer.sif",                        true);
                _h.get_at_path(this->base_temp,                 "elmer.base_temp",                  true);
                _h.get_at_path(this->node_data_file,            "elmer.node_data_file",             true);
                _h.get_at_path(this->dump_directory,            "elmer.dump_directory",             true);
                _h.get_at_path(this->node_temperature_file_ext, "elmer.node_temperature_file_ext",  true);
                _h.get_at_path(this->node_position_file_ext,    "elmer.node_position_file_ext",     true);

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
                    util::copyFile(this->meshDB + SEP + "mesh." + exts[i], this->dump_directory + SEP + "mesh." + exts[i]);
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
                this->loadBoundary();
                this->loadElements();
            }


            void run() {
                // making the sif file
                std::ofstream _buf(this->sif);
                if (!(_buf.is_open())) util::error(this->sif + " did not open");
                this->join(_buf);                
                _buf.close();

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
                
                this->loadVarNodeData();

            }


            void getNodeAtIndex(int _index, int _boundary_index, double(&_point)[3]) {
                for (int j = 0; j < elmer::dimension; j++)
                    _point[j] = std::get<1>(this->nodes[this->boundaries[_index][_boundary_index]-1])[j];
            }


            void createConditionsFrom(double*& _qw, double*& _fx, double*& _fy, double*& _fz, double*& _shx, double*& _shy, double*& _shz, int _ntimesteps, double _timestep_size, int _length) {
                int i, _size, size;
                double avg;

                // clearing the vectors in elmer to make room for the new data
                this->bodies.clear();
                this->initial_conditions.clear();
                this->boundary_conditions.clear();

                this->simulation.Timestep_Sizes     = { _timestep_size };
                // modifying elmer variables to match sparta vairables
                this->simulation.Timestep_intervals = { _ntimesteps };
                this->simulation.Output_Intervals   = { _ntimesteps };

                /***
                 * boundary conditions
                 */
                for (i = 0; i < _length; i++) {
                    elmer::Boundary_Condition* bc = new elmer::Boundary_Condition(i+1);
                    bc->Heat_Flux = _qw[i]/_ntimesteps; // the average heat flux
                    this->boundary_conditions.push_back(bc);
                }

                /***
                 * initial conditions
                 */
                size = this->node_temperature_data.size();

                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (long unsigned int i = 0; i < this->elements.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    _size = this->elements[i].size();
                    for (int j = 0; j < _size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        if (this->elements[i][j]-1 >= size)
                            util::error("index out of bounds in creating initial conditions");
                        avg += (this->node_temperature_data[this->elements[i][j]-1]/((double)_size));
                    }

                    // surface element
                    Body* _body = new Body(i+1);
                    _body->Initial_condition = i+1;
                    _body->Target_Bodies = {(int)i+1};
                    Initial_Condition* _ic = new Initial_Condition(i+1);
                    _ic->Temperature = avg;

                    this->initial_conditions.push_back(_ic);
                    this->bodies.push_back(_body);
                }
            }


            void averageNodeTemperaturesInto(double*& _temperatures, int _length) {
                // making sure everything is consistent
                if ((long unsigned int)_length != this->boundaries.size())
                    util::error("boundary data does not match required size, required size " + std::to_string(_length) + ", size " + std::to_string(this->boundaries.size()));

                // used later
                double avg;
                
                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (long unsigned int i = 0; i < this->boundaries.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    for (unsigned int j = 1; j < boundary_size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        avg += this->node_temperature_data[this->boundaries[i][j]-1];
                    }
                    // computing the average by dividing the sum by the number of points and setting the surface element
                    _temperatures[this->boundaries[i][0]-1] = avg/(boundary_size - 1);
                }
            }

            void loadVarNodeData() {
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


            void loadNodeData() {
                std::vector<std::string> lines, split_line;
                std::tuple<int, std::array<double, 3>> data_item;

                util::readFile(node_file, lines);
                int count = 0;
                for (std::string& it : lines) {
                    count++;
                    util::trim(it);
                    util::split(it, split_line, ' ');
                    if (split_line.size() != node_element_size)
                        util::error("node element not correct size at line " + std::to_string(count) + ", should have size " + std::to_string(node_element_size));

                    if (std::stoi(split_line[0]) != count)
                        util::error("detected unordered node element at line " + std::to_string(count));

                    std::get<0>(data_item) = std::stoi(split_line[1]);
                    std::get<1>(data_item) = { std::stod(split_line[2]), std::stod(split_line[3]), std::stod(split_line[4]) };
                    this->nodes.push_back(data_item);
                }
            }


            void loadBoundary() {
                std::vector<std::string> lines, _split;
                std::array<int, boundary_size> arr;
                this->boundaries.clear();

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
                    this->boundaries.push_back(arr);
                    count++;
                }
            }


            void loadElements() {
                std::vector<std::string> _temp_data, _split;
                std::vector<int> arr;
                this->elements.clear();

                util::readFile(this->element_file, _temp_data);

                this->elements.resize(_temp_data.size());

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
                    this->elements[std::stoi(_split[0]) - 1] = arr;
                }
            }


            void join(std::ofstream& _buf) {
                header.join(_buf);
                simulation.join(_buf);
                constants.join(_buf);
                thermal_solver.join(_buf);
                elastic_solver.join(_buf);
                equation.join(_buf);
                material.join(_buf);

                for (long unsigned int i = 0; i < this->bodies.size(); i++)
                    this->bodies[i]->join(_buf);

                for (long unsigned int i = 0; i < this->initial_conditions.size(); i++)
                    this->initial_conditions[i]->join(_buf);

                for (long unsigned int i = 0; i < this->boundary_conditions.size(); i++)
                    this->boundary_conditions[i]->join(_buf);
            }


            void dump(long _timestep) {
                this->dumpNodeTemperatures(_timestep);
                this->dumpNodePositions(_timestep);
            }


            void dumpNodeTemperatures(long _timestep) {
                std::string filename = this->dump_directory + SEP + std::to_string(_timestep) + "." + this->node_temperature_file_ext;
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
                std::string filename = this->dump_directory + SEP + std::to_string(_timestep) + "." + this->node_position_file_ext;
                std::ofstream out(filename);
                if (!(out.is_open())) util::error(filename + " did not open");

                if (this->node_delta_x_data.size() != this->node_delta_y_data.size() || 
                    this->node_delta_y_data.size() != this->node_delta_z_data.size() || 
                    this->node_delta_x_data.size() != this->node_delta_z_data.size())
                    util::error("position delta vectors do not have the same size");

                this->double_converter.str("");
                for (long unsigned int i = 0; i < this->node_delta_x_data.size(); i++) {    
                    this->double_converter << std::get<1>(this->nodes[i])[0] << " " << std::get<1>(this->nodes[i])[1] << " " << std::get<1>(this->nodes[i])[2];
                    out << this->double_converter.str() << "\n";
                    this->double_converter.str("");
                }
                out.close();

            }

    };
}


#endif