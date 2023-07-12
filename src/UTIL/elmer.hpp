#ifndef ELMER_H
#define ELMER_H

#include "elmer_classes.hpp"

/* ---------------------------------------------------------------------- */

// defining the os path separator based on the detected os
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    static const char SEP = '\\';
#else
    static const char SEP = '/';
#endif

/* ---------------------------------------------------------------------- */

namespace elmer {
    const unsigned int boundary_size = 3; // number of elements in a boundary element
    const unsigned int dimension     = 3; // the dimension
    const std::size_t node_element_size = 5;

    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            const std::string sep = "\n\n";

        public:
            /* ---------------------------------------------------------------------- */

            std::vector<double> node_temperature_data;
            std::vector<std::array<int, elmer::boundary_size>> boundaries;
            // array format: x, y, z
            std::vector<std::array<double, 3>>  nodes;
            std::vector<std::vector<int>> elements;
            std::string name, exe, sif, meshDB, node_data_file, simulation_directory, boundary_file, element_file, node_file, node_temperature_file_ext, node_position_file_ext, print_intensity;
            double base_temp;
            bool gravity_on;

            Header          header;
            Simulation      simulation;
            Constants       constants;
            Thermal_Solver  thermal_solver;
            Elastic_Solver  elastic_solver;
            Equation        equation;
            Material        material;
            std::vector<Body_Force*>         body_forces;
            std::vector<Body*>               bodies;
            std::vector<Initial_Condition*>  initial_conditions;
            std::vector<Boundary_Condition*> boundary_conditions;

            /* ---------------------------------------------------------------------- */

            Elmer() {
                this->header         = Header();
                this->simulation     = Simulation();
                this->constants      = Constants();
                this->equation       = Equation();
                this->material       = Material();
                this->thermal_solver = Thermal_Solver();
                this->elastic_solver = Elastic_Solver();
            }

            /* ---------------------------------------------------------------------- */

            void set(toml::handler& _h) {
                ULOG("setting up elmer");
                // setting needed variables
                _h.getAtPath(this->meshDB,                    "elmer.meshDB",                     true);
                _h.getAtPath(this->exe,                       "elmer.exe",                        true);
                _h.getAtPath(this->sif,                       "elmer.sif",                        true);
                _h.getAtPath(this->base_temp,                 "elmer.base_temp",                  true);
                _h.getAtPath(this->simulation_directory,      "simulation_directory",             true);
                _h.getAtPath(this->node_temperature_file_ext, "elmer.node_temperature_file_ext",  true);
                _h.getAtPath(this->node_position_file_ext,    "elmer.node_position_file_ext",     true);
                _h.getAtPath(this->gravity_on,                "elmer.gravity_on",                 true);
                _h.getAtPath(this->print_intensity,           "elmer.print_intensity",            false);

                // each elmer section sets its own variables, so passing it the data structure
                // this->header.set(_h);
                this->simulation.set(_h);
                this->constants.set(_h);
                this->equation.set(_h);
                this->material.set(_h);
                this->thermal_solver.set(_h);
                this->elastic_solver.set(_h);

                this->checks();
            }

            /* ---------------------------------------------------------------------- */

            void checks() {
                ULOG("checking parameters");

                struct stat sb;
                // making sure the elmer exe 
                if (!(stat(this->exe.c_str(),  &sb) == 0))
                    UERR("Illegal fix fea command, exe path does not exist");
                
                // making sure the mesh database is complete
                std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    if (!(stat((this->meshDB + SEP + "mesh." + exts[i]).c_str(),  &sb) == 0))
                        UERR("Illegal fix fea command, mesh database incomplete, " + (this->meshDB + SEP + "mesh." + exts[i]) + " does not exist");
                }

                // making sure base temperature is valid
                if (this->base_temp <= 0.0)
                    UERR("base temperature must be greater than 0");

                if (this->simulation_directory.back() == SEP)
                    this->simulation_directory = this->simulation_directory.substr(0, this->simulation_directory.length()-1);
                
                if (this->print_intensity == toml::noString) this->print_intensity = "none";
                if (this->print_intensity != "none" && this->print_intensity != "limited" && this->print_intensity != "verbose")
                    UERR("invalid print intensity for elmer");
                
                // setting vars from provided values
                this->node_data_file = this->simulation_directory + SEP + this->simulation.Output_File;
                this->simulation.Output_File = this->node_data_file;
                
                this->header.Mesh_DB = { this->simulation_directory, "." };
                this->header.Include_Path = "";
                this->header.Results_Directory = this->simulation_directory;

                if (this->gravity_on) {
                    elmer::Body_Force* bf = new elmer::Body_Force(1);
                    if (this->constants.Gravity[0] != 0.0) {
                        bf->Stress_Bodyforce_1 = this->constants.Gravity[0]*this->constants.Gravity[(int)this->constants.Gravity.size()-1];
                        ULOG("gravity acts (at least in part) in the x-direction");
                    }
                    if (this->constants.Gravity[1] != 0.0) {
                        bf->Stress_Bodyforce_2 = this->constants.Gravity[1]*this->constants.Gravity[(int)this->constants.Gravity.size()-1];
                        ULOG("gravity acts (at least in part) in the y-direction");
                    }
                    if (this->constants.Gravity[2] != 0.0) {
                        bf->Stress_Bodyforce_3 = this->constants.Gravity[2]*this->constants.Gravity[(int)this->constants.Gravity.size()-1];
                        ULOG("gravity acts (at least in part) in the z-direction");
                    }
                    bf->Density = this->material.Density;
                    this->body_forces.push_back(bf);
                }
            }

            /* ---------------------------------------------------------------------- */

            void setupDumpDirectory() {
                ULOG("setting up dump directory");
                // making sure the mesh database is complete
                std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    util::copyFile(this->meshDB + SEP + "mesh." + exts[i], this->simulation_directory + SEP + "mesh." + exts[i]);
                }

                // setting files
                this->boundary_file = this->simulation_directory + SEP + "mesh.boundary";
                this->element_file = this->simulation_directory + SEP + "mesh.elements";
                this->node_file = this->simulation_directory + SEP + "mesh.nodes";
            }

            /* ---------------------------------------------------------------------- */
            
            void setupVectors() {
                ULOG("setting up vectors");
                // getting the number of nodes
                int count = util::countLinesInFile(this->node_file);

                // setting temperature vector to all base temperatures
                std::string str = util::dtos(this->base_temp);

                // writing base temperature to file for each node
                for (int i = 0; i < count; i++)
                    this->node_temperature_data.push_back(this->base_temp);
            }

            /* ---------------------------------------------------------------------- */

            void setup() {
                this->setupDumpDirectory();
                this->setupVectors();
                this->loadBoundaries();
                this->loadElements();
                this->loadNodes();
            }

            /* ---------------------------------------------------------------------- */

            void makeSif() {
                util::oFile _buf(this->sif);
                this->join(_buf);
            }

            /* ---------------------------------------------------------------------- */

            void run() {
                // making the sif file
                this->makeSif();

                ULOG("running");
                // running the command
                // must have " 2>&1" at end to pipe stderr to stdout
                int exit_status;
                try {
                    if (this->print_intensity == "none")
                        exit_status = util::quietExec(exe+" "+this->sif+" 2>&1");
                    else if (this->print_intensity == "limited")
                        exit_status = util::limitedExec(exe+" "+this->sif+" 2>&1", "MAIN: Time:", ":");
                    else
                        exit_status = util::verboseExec(exe+" "+this->sif+" 2>&1");
                } catch (std::exception& e) {
                    UERR(e.what());
                }
                // if the command did not succeed
                if (exit_status)
                    UERR("check elmer output");
                // else
                //     util::printToFile(command_result.output, 0);
                
                this->loadNodeData();
                this->updateNodeFile();

            }

            /* ---------------------------------------------------------------------- */

            void getNodePointAtIndex(int _index, int _boundary_index, double(&_point)[3]) {
                for (std::size_t j = 0; j < elmer::dimension; j++)
                    _point[j] = this->nodes[this->boundaries[_index][_boundary_index]-1][j];
            }

            /* ---------------------------------------------------------------------- */

            void createConditionsFrom(double*& _qw, double*& _px, double*& _py, double*& _pz, double*& _shx, double*& _shy, double*& _shz, int _ntimesteps, double _timestep_size, int _length) {
                int i, _size, size, _body_force;
                double avg;

                ULOG("creating conditions");

                if (this->gravity_on)
                    _body_force = 1;
                else
                    _body_force = toml::noInt;

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
                    elmer::Boundary_Condition* bc;
                    bc = new elmer::Boundary_Condition(2*i+1);
                    bc->Heat_Flux_BC = true;
                    bc->Heat_Flux = _qw[i]/_ntimesteps; // the average heat flux
                    this->boundary_conditions.push_back(bc);
                    
                    bc = new elmer::Boundary_Condition(2*(i+1));
                    // setting the stress tensor tensor components, averaging them
                    bc->stress_6_vector = {_px[i]/_ntimesteps, _py[i]/_ntimesteps, _pz[i]/_ntimesteps, _shx[i]/_ntimesteps, _shy[i]/_ntimesteps, _shz[i]/_ntimesteps};
                    this->boundary_conditions.push_back(bc);
                }

                /***
                 * initial conditions
                 */
                size = this->node_temperature_data.size();

                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (std::size_t i = 0; i < this->elements.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    _size = this->elements[i].size();
                    for (int j = 0; j < _size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        if (this->elements[i][j]-1 >= size)
                            UERR("index out of bounds in creating initial conditions");
                        avg += (this->node_temperature_data[this->elements[i][j]-1]/((double)_size));
                    }

                    // surface element
                    Body* _body = new Body(i+1);
                    _body->Initial_condition = i+1;
                    _body->Target_Bodies = {(int)i+1};
                    _body->Body_Force = _body_force;
                    Initial_Condition* _ic = new Initial_Condition(i+1);
                    _ic->Temperature = avg;

                    this->initial_conditions.push_back(_ic);
                    this->bodies.push_back(_body);
                }
            }

            /* ---------------------------------------------------------------------- */

            void averageNodeTemperaturesInto(double*& _temperatures, int _length) {
                // making sure everything is consistent
                if ((std::size_t)_length != this->boundaries.size())
                    UERR("boundary data does not match required size, required size " + std::to_string(_length) + ", size " + std::to_string(this->boundaries.size()));
                
                ULOG("averaging the per node temperatures into surf");

                // used later
                double avg;
                
                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (std::size_t i = 0; i < this->boundaries.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    avg = 0;
                    for (unsigned int j = 0; j < boundary_size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        avg += this->node_temperature_data[this->boundaries[i][j]-1];
                    }
                    // computing the average by dividing the sum by the number of points and setting the surface element
                    _temperatures[i] = avg/boundary_size;
                }
            }

            /* ---------------------------------------------------------------------- */

            void loadNodeData() {
                ULOG("loading node data");
                std::size_t i, start;
                std::vector<std::string> lines, split_line;
                std::string cur_key, line;
                util::readFile(this->node_data_file, lines);
                start = 0;
                for (i = 0; i < lines.size(); i++) {
                    util::split(lines[i], split_line, ' ');
                    if (split_line[0] == (std::string)"Time:") start = i;
                }
                std::vector<std::array<double, 3>> temp_nodes = this->nodes;

                int cur_index = 0;
                for (i = start+1; i < lines.size(); i++) {
                    line = lines[i];
                    util::trim(line);
                    util::split(line, split_line, ' ');
                    if (lines[i][0] == ' ' && split_line.size() == 1) {
                        if (cur_key == (std::string)"temperature") {
                            this->node_temperature_data[cur_index] = std::stod(line);
                        } else if (cur_key == (std::string)"displacement 1") {
                            this->nodes[cur_index][0] += std::stod(line);
                        } else if (cur_key == (std::string)"displacement 2") {
                            this->nodes[cur_index][1] += std::stod(line);
                        } else if (cur_key == (std::string)"displacement 3") {
                            this->nodes[cur_index][2] += std::stod(line);
                        } else
                            UERR("got unknown key in data file: " + cur_key);
                        cur_index++;
                    } else if ("Perm:" == split_line[0]) {
                        if (cur_index > 0) {
                            if ((std::size_t)cur_index != this->nodes.size()) UERR("did not get all indicies while loading data");
                        }
                        cur_index = 0;
                        cur_key = lines[i-1];
                        util::trim(cur_key);
                    }
                }

                // int j;
                // for (i = 0; i < temp_nodes.size(); i++) {
                //     for (j = 0; j < 3; j++)
                //         if (std::abs(nodes[i][j] - temp_nodes[i][j]) > 0.0) {
                //             ULOG("got difference");
                //             break;
                //         }
                // }
            }

            /* ---------------------------------------------------------------------- */

            void updateNodeFile() {
                ULOG("updating node file");
                util::oFile out(this->node_file);

                for (std::size_t i = 0; i < this->nodes.size(); i++)
                    out << i+1 << " " << -1 << " " << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }

            /* ---------------------------------------------------------------------- */
            
            /**
             * makes sure file from elmer object data
             * returns: the name of the surf file
            */
            std::string makeSpartaSurf() {
                ULOG("making surface file for sparta");
                std::size_t i;
                unsigned int j;
                std::string surf_file = this->simulation_directory + SEP + "mesh.surf";
                
                util::oFile f(surf_file);
                f << "# Surface element file written by SPARTA fea interface\n\n" << this->nodes.size() << " points\n" << this->boundaries.size() << " triangles\n\nPoints\n\n";

                for (i = 0; i < this->nodes.size(); i++) {
                    f << i+1;
                    for (j = 0; j < elmer::dimension; j++)
                        f << " " << this->nodes[i][j];
                    f << "\n";
                }

                f << "\nTriangles\n\n";
                for (i = 0; i < this->boundaries.size(); i++) {
                    f << i+1;
                    for (j = 0; j < elmer::boundary_size; j++)
                        f << " " << this->boundaries[i][j];
                    f << "\n";
                }

                return surf_file;
            }

            /* ---------------------------------------------------------------------- */

            void loadNodes() {
                ULOG("loading nodes");
                std::vector<std::string> lines, split_line;
                std::array<double, 3> data_item;

                util::readFile(node_file, lines);
                int count = 0;
                for (std::string& it : lines) {
                    count++;
                    util::trim(it);
                    util::split(it, split_line, ' ');
                    if (split_line.size() != node_element_size)
                        UERR("node element not correct size at line " + std::to_string(count) + ", should have size " + std::to_string(node_element_size));

                    if (std::stoi(split_line[0]) != count)
                        UERR("detected unordered node element at line " + std::to_string(count));

                    if (std::stoi(split_line[1]) != -1)
                        UERR("detected a value that is not -1 in the second column of the node file, this is not yet support");
                    

                    data_item = { std::stod(split_line[2]), std::stod(split_line[3]), std::stod(split_line[4]) };
                    this->nodes.push_back(data_item);
                }
            }

            /* ---------------------------------------------------------------------- */

            void loadBoundaries() {
                ULOG("loading boundaries");
                std::vector<std::string> lines, _split;
                std::array<int, boundary_size> arr;
                this->boundaries.clear();

                util::readFile(this->boundary_file, lines);

                int last_id = 0;
                int cur_id;
                int count = 1;
                for (std::string it : lines) {
                    util::trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    util::split(it, _split, ' ');

                    if (_split[4] != (std::string)"303") 
                        UERR("element is not a triangle in boundary file at line: " + std::to_string(count));

                    cur_id = std::stoi(_split[0]);
                    if ((cur_id - last_id) != 1)
                        UERR("got unordered boundary elements at line: " + std::to_string(count));
                    last_id = cur_id;

                    // getting rid of stuff I do not need
                    _split.erase(_split.begin(), _split.begin() + 5);

                    // catching errors
                    // if (_split.size() != boundary_size)
                    //     UERR("too many boundary elements to assign at line: " + std::to_string(count));

                    // adding the data
                    for (unsigned int i = 0; i < boundary_size; i++)
                        arr[i] = std::stoi(_split[i]);

                    // adding the data to the class vector
                    this->boundaries.push_back(arr);
                    count++;
                }
            }

            /* ---------------------------------------------------------------------- */

            void loadElements() {
                ULOG("loading elements");
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
                    for (std::size_t i = 1; i < _split.size(); i++)
                        arr.push_back(std::stoi(_split[i]));

                    // adding the data to the class vector
                    this->elements[std::stoi(_split[0]) - 1] = arr;
                }
            }

            /* ---------------------------------------------------------------------- */

            void join(util::oFile& _buf) {
                header.join(_buf);
                simulation.join(_buf);
                constants.join(_buf);
                thermal_solver.join(_buf);
                elastic_solver.join(_buf);
                equation.join(_buf);
                material.join(_buf);

                for (std::size_t i = 0; i < this->bodies.size(); i++)
                    this->bodies[i]->join(_buf);

                for (std::size_t i = 0; i < this->body_forces.size(); i++)
                    this->body_forces[i]->join(_buf);

                for (std::size_t i = 0; i < this->initial_conditions.size(); i++)
                    this->initial_conditions[i]->join(_buf);

                for (std::size_t i = 0; i < this->boundary_conditions.size(); i++)
                    this->boundary_conditions[i]->join(_buf);
            }

            /* ---------------------------------------------------------------------- */

            void dumpBefore(long _timestep) {
                this->dumpNodePositionsBefore(_timestep);
                this->dumpNodeTemperaturesBefore(_timestep);
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeTemperaturesBefore(long _timestep) {
                ULOG("dumping node temperatures before");
                std::string filename = this->simulation_directory + SEP + std::to_string(_timestep) + ".before_" + this->node_temperature_file_ext;
                util::oFile out(filename);

                for (double& it : this->node_temperature_data)
                    out << it << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodePositionsBefore(long _timestep) {
                ULOG("dumping node positions before");
                std::string filename = this->simulation_directory + SEP + std::to_string(_timestep) + ".before_" + this->node_position_file_ext;
                util::oFile out(filename);

                for (std::size_t i = 0; i < this->nodes.size(); i++) 
                    out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dump(long _timestep) {
                this->dumpNodeTemperatures(_timestep);
                this->dumpNodePositions(_timestep);
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeTemperatures(long _timestep) {
                ULOG("dumping node temperatures");
                std::string filename = this->simulation_directory + SEP + std::to_string(_timestep) + "." + this->node_temperature_file_ext;
                util::oFile out(filename);

                for (double& it : this->node_temperature_data)
                    out << it << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodePositions(long _timestep) {
                ULOG("dumping node positions");
                std::string filename = this->simulation_directory + SEP + std::to_string(_timestep) + "." + this->node_position_file_ext;
                util::oFile out(filename);

                for (std::size_t i = 0; i < this->nodes.size(); i++) 
                    out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }
    };
}

#endif