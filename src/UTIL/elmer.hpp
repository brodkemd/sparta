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

// dummy class to allow the use of a pointer later
namespace SPARTA_NS { class FixFea; }

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
            SPARTA_NS::FixFea* sparta;
            std::vector<util::string_t> run_these;
            std::vector<util::double_t> node_temperature_data;
            std::vector<std::array<util::int_t, elmer::boundary_size>> boundaries;
            // array format: x, y, z
            std::vector<std::array<util::double_t, elmer::dimension>>  nodes, node_velocity_data;
            std::vector<std::pair<util::int_t, std::vector<util::int_t>>> elements;
            toml::Item_t name, exe, sif, meshDB, node_data_file, simulation_directory, boundary_file, element_file, node_file, node_temperature_file_ext, node_position_file_ext, node_velocity_file_ext, print_intensity, gravity_on, base_temp, energy_threshold, pressure_threshold, shear_threshold;

            /* ---------------------------------------------------------------------- */

            Section header, simulation, constants;
            std::vector<Section*> solvers, materials, equations, body_forces, bodies, initial_conditions, boundary_conditions;

            /* ---------------------------------------------------------------------- */

        public:
            Elmer(SPARTA_NS::FixFea* _sparta) {
                // referencing the calling sparta instance
                this->sparta         = &*_sparta;

                // initing class variables
                this->header          = Section("Header", toml::noInt, " ", false);
                this->simulation      = Section("Simulation");
                this->constants       = Section("Constants");
                this->equations.clear(); this->materials.clear();
                this->body_forces.clear();  this->bodies.clear();
                this->solvers.clear();
                this->initial_conditions.clear();
                this->boundary_conditions.clear();
            }

            /* ---------------------------------------------------------------------- */

            void set(toml::handler& _h) {
                ULOG("setting elmer class parameters");
                // setting needed variables
                _h.getAtPath(this->meshDB,                    "elmer.meshDB",                     toml::STRING);
                _h.getAtPath(this->exe,                       "elmer.exe",                        toml::STRING);
                _h.getAtPath(this->sif,                       "elmer.sif",                        toml::STRING);
                _h.getAtPath(this->base_temp,                 "elmer.base_temp",                  toml::DOUBLE);
                _h.getAtPath(this->simulation_directory,      "simulation_directory",             toml::STRING);
                _h.getAtPath(this->node_temperature_file_ext, "elmer.node_temperature_file_ext",  toml::STRING);
                _h.getAtPath(this->node_position_file_ext,    "elmer.node_position_file_ext",     toml::STRING);
                _h.getAtPath(this->node_velocity_file_ext,    "elmer.node_velocity_file_ext",     toml::STRING);
                _h.getAtPath(this->gravity_on,                "elmer.gravity_on",                 toml::BOOL);
                _h.getAtPath(this->energy_threshold,          "sparta.energy_threshold",          toml::DOUBLE);
                _h.getAtPath(this->pressure_threshold,        "sparta.pressure_threshold",        toml::DOUBLE);
                _h.getAtPath(this->shear_threshold,           "sparta.shear_threshold",           toml::DOUBLE);
                _h.getAtPath(this->print_intensity,           "elmer.print_intensity",            toml::STRING);

                // each elmer section sets its own variables, so passing it the data structure
                // this->header.set(_h);
                _h.getDictAtPath(this->simulation.contents, "elmer.simulation");
                _h.getDictAtPath(this->constants.contents, "elmer.constants");

                std::vector<toml::Item_t> keys;

                _h.getDictKeysAtPath(keys, "elmer.material");
                for (auto key : keys) {
                    Section* material = new Section("Material", std::stoi(key.toString()));
                    _h.getDictAtPath(material->contents, "elmer.material." + key.toString());
                    this->materials.push_back(material);
                }
                ULOG("# of materials loaded: " + std::to_string(this->materials.size()));
                
                _h.getDictKeysAtPath(keys, "elmer.solver");
                for (auto key : keys) {
                    Section* solver = new Section("Solver", std::stoi(key.toString()));
                    _h.getDictAtPath(solver->contents, "elmer.solver." + key.toString());
                    this->solvers.push_back(solver);
                }
                ULOG("# of solvers loaded: " + std::to_string(this->solvers.size()));


                this->checks();
            }

            /* ---------------------------------------------------------------------- */
            
            /**
             * makes sure file from elmer object data
             * returns: the name of the surf file
            */
            util::string_t makeSpartaSurf() {
                ULOG("making surface file for sparta");
                std::size_t i;
                unsigned int j;
                util::string_t surf_file = this->simulation_directory.toString() + SEP + "mesh.surf";
                
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

            void getNodePointAtIndex(util::int_t _index, util::int_t _boundary_index, double(&_point)[3]) {
                for (std::size_t j = 0; j < elmer::dimension; j++)
                    _point[j] = this->nodes[this->boundaries[_index][_boundary_index]-1][j];
            }

            /* ---------------------------------------------------------------------- */

            void run() {
                this->createConditions();
                ULOG("# of boundary conditions: " + std::to_string(this->boundary_conditions.size()));
                ULOG("# of initial conditions and bodies: " + std::to_string(this->initial_conditions.size()));

                // setting the timestep size
                if (!(this->simulation.hasKey("Timestep_Sizes"))) {
                    ULOG("did not provide \"timestep size\", using Sparta's");
                    this->simulation.setItem("Timestep_Sizes", (std::vector<util::double_t>){ (util::double_t)this->sparta->update->dt });
                }

                // setting number of timesteps to run
                if (!(this->simulation.hasKey("Timestep_intervals"))) {
                    ULOG("did not provide \"number of timesteps\", using Sparta's");
                    this->simulation.setItem("Timestep_intervals", (std::vector<util::int_t>){ (util::int_t)this->sparta->run_every });
                }
                
                if (!(this->simulation.hasKey("Output_Intervals"))) {
                    ULOG("did not provide \"output interval\", using Sparta's");
                    this->simulation.setItem("Output_Intervals", (std::vector<util::int_t>){ (util::int_t)this->sparta->run_every });
                }

                // making the sif file
                this->makeSif();

                ULOG("running");
                
                // running the command
                // must have " 2>&1" at end to pipe stderr to stdout
                int exit_status = 0;
                try {
                    if (this->print_intensity == "none")
                        exit_status = util::quietExec(this->exe.toString() + " " + this->sif.toString() + " 2>&1");
                    else if (this->print_intensity == "limited")
                        exit_status = util::limitedExec(this->exe.toString() + " " + this->sif.toString() + " 2>&1", "MAIN: Time:", ":");
                    else
                        exit_status = util::verboseExec(this->exe.toString() + " " + this->sif.toString() + " 2>&1");
                } catch (std::exception& e) {
                    UERR(e.what());
                }
                // if the command did not succeed
                if (exit_status)
                    UERR("check elmer output");
                
                this->loadNodeData();
                this->updateNodeFile();
            }


            bool shouldRun() {
                util::int_t i; util::double_t sum;

                ULOG("checking if solver should be run");
                this->run_these.clear();

                // summing all of the vars and checking if they are above the threshold
                // returns true if above the threshold
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->qw[i]);
                if (sum > this->energy_threshold.toDouble())   { this->run_these.push_back("qw"); }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->px[i]);
                if (sum > this->pressure_threshold.toDouble()) { this->run_these.push_back("px");     }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->py[i]);
                if (sum > this->pressure_threshold.toDouble()) { this->run_these.push_back("py");     }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->pz[i]);
                if (sum > this->pressure_threshold.toDouble()) { this->run_these.push_back("pz");     }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->shx[i]);
                if (sum > this->shear_threshold.toDouble())    { this->run_these.push_back("shx");    }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->shy[i]);
                if (sum > this->shear_threshold.toDouble())    { this->run_these.push_back("shy");    }
                sum = 0;
                for (i = 0; i < this->sparta->nsurf; i++) sum+=std::abs(this->sparta->shz[i]);
                if (sum > this->shear_threshold.toDouble())    { this->run_these.push_back("shz");    }
                
                return this->run_these.size() > (std::size_t)0;
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

            void dumpBefore() { 
                this->dumpNodePositionsBefore();
                this->dumpNodeTemperaturesBefore(); 
                this->dumpNodeVelocitiesBefore();
            }
            
            /* ---------------------------------------------------------------------- */

            void dump() {
                this->dumpNodeTemperatures();
                this->dumpNodePositions();
                this->dumpNodeVelocities();
            }

        /* ---------------------------------------------------------------------- */

        protected:
            void checks() {
                ULOG("checking elmer class parameters");

                // making sure the dimension is correct
                if (this->sparta->domain->dimension != elmer::dimension)
                    UERR("Invalid dimension detected, must be " + std::to_string(elmer::dimension) + " dimensional");

                struct stat sb;
                // making sure the elmer exe 
                if (!(stat(this->exe.toString().c_str(),  &sb) == 0))
                    UERR("Illegal fix fea command, exe path does not exist");
                
                // making sure the mesh database is complete
                util::string_t exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    if (!(stat((this->meshDB.toString() + SEP + "mesh." + exts[i]).c_str(),  &sb) == 0))
                        UERR("Illegal fix fea command, mesh database incomplete, " + (this->meshDB.toString() + SEP + "mesh." + exts[i]) + " does not exist");
                }

                // making sure base temperature is valid
                if (this->base_temp.toDouble() <= 0.0)
                    UERR("base temperature must be greater than 0");

                if (this->simulation_directory.toString().back() == SEP)
                    this->simulation_directory = this->simulation_directory.toString().substr(0, this->simulation_directory.length()-1);
                
                // setting vars from provided values
                this->node_data_file = this->simulation_directory.toString() + SEP + this->simulation.getItem("Output_File").toString();
                this->simulation.setItem("Output_File", this->node_data_file);
                ULOG("Will load node data from: " + this->node_data_file.toString());

                toml::Item_t _temp_solvers = std::vector<toml::Item_t>({});
                for (auto& it : this->solvers) {
                    _temp_solvers.append(it->getId());
                }

                for (auto& it : this->equations) {
                    it->setItem("Active_Solvers",  _temp_solvers); // id for each solver of the class
                }
                
                this->header.setItem("Mesh_DB", std::vector<util::string_t>({ this->simulation_directory.toString(), "." }));
                this->header.setItem("Results_Directory", this->simulation_directory);

                if (this->gravity_on.toBool()) {
                    util::int_t count = 1;
                    for (auto it : this->materials) {
                        elmer::Section* bf = new elmer::Section("Body Force", it->getId());
                        ULOG("gravity acts as: " + this->constants.getItem("Gravity").toString());
                        if (this->constants.getItem("Gravity")[0].toDouble() != 0.0) {
                            (*bf).setItem("Stress_Bodyforce_1", this->constants.getItem("Gravity")[0].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                            //ULOG("gravity acts (at least in part) in the x-direction");
                        }
                        if (this->constants.getItem("Gravity")[1].toDouble() != 0.0) {
                            (*bf).setItem("Stress_Bodyforce_2", this->constants.getItem("Gravity")[1].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                            //ULOG("gravity acts (at least in part) in the y-direction");
                        }
                        if (this->constants.getItem("Gravity")[2].toDouble() != 0.0) {
                            (*bf).setItem("Stress_Bodyforce_3", this->constants.getItem("Gravity")[2].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                            //ULOG("gravity acts (at least in part) in the z-direction");
                        }
                        this->body_forces.push_back(bf);
                        count++;
                    }
                }
                // UERR("done");
            }

            /* ---------------------------------------------------------------------- */

            

            void setupDumpDirectory() {
                ULOG("setting up dump directory");
                // making sure the mesh database is complete
                util::string_t exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
                for (int i = 0; i < 4; i++) {
                    util::copyFile(this->meshDB.toString() + SEP + "mesh." + exts[i], this->simulation_directory.toString() + SEP + "mesh." + exts[i]);
                }

                // setting files
                this->boundary_file = this->simulation_directory.toString() + SEP + "mesh.boundary";
                this->element_file  = this->simulation_directory.toString() + SEP + "mesh.elements";
                this->node_file     = this->simulation_directory.toString() + SEP + "mesh.nodes";
            }

            /* ---------------------------------------------------------------------- */
            
            void setupVectors() {
                ULOG("setting up vectors");
                // getting the number of nodes
                long count = util::countLinesInFile(this->node_file.toString());

                std::array<util::double_t, elmer::dimension> temp;
                temp.fill((util::double_t)0);

                // writing base temperature to file for each node
                for (long i = 0; i < count; i++) {
                    this->node_temperature_data.push_back(this->base_temp.toDouble());
                    this->node_velocity_data.push_back(temp);
                }
            }

            /* ---------------------------------------------------------------------- */

            void makeSif() { util::oFile _buf(this->sif.toString()); this->join(_buf); }

            /* ---------------------------------------------------------------------- */

            void createConditions() {
                std::vector<std::vector<util::int_t>> ids;
                std::vector<util::int_t> _temp_vec, indicies;
                std::size_t _size, size;
                util::int_t index, j, k, _body_force;
                const util::int_t numDataPoints = 7;
                ids.clear(); _temp_vec.clear(); indicies.clear();

                util::int_t bc_count = 1;
                util::string_t _temp;
                for (util::string_t& it : this->run_these) _temp+=(" " + it);

                ULOG("variables over threshold:" + _temp);

                // clearing the vectors in elmer to make room for the new data
                this->bodies.clear();
                this->initial_conditions.clear();
                this->boundary_conditions.clear();

                /***
                 * boundary conditions
                */
                // this part is for the stresses
                // first determines which stress should be accounted for
                index = 0;
                if (util::find(run_these, std::string("qw"))  != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("px"))  != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("py"))  != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("pz"))  != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("shx")) != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("shy")) != util::npos)
                    indicies.push_back(index++);
                if (util::find(run_these, std::string("shz")) != util::npos)
                    indicies.push_back(index++);

                // if stress should be accounted for
                if (indicies.size() > 0) {
                    // get the stress values from sparta
                    double* arr[numDataPoints] = { &*this->sparta->qw, &*this->sparta->px, &*this->sparta->py, &*this->sparta->pz, &*this->sparta->shx, &*this->sparta->shy, &*this->sparta->shz };
                    // ULOG("adding stress boundary conditions");

                    std::vector<std::array<util::double_t, numDataPoints>> data;
                    std::array<util::double_t, numDataPoints> temp_data;
                    data.clear();
                    
                    // detecting unique values and saving them to the stress vector and
                    // their corresponding ids to the ids vector, this condenses the
                    // number of boundary conditions (removes boundary conditions with
                    // the same values)
                    for (util::int_t i = 0; i < this->sparta->nsurf; i++) {
                        temp_data.fill(0.0);
                        for (util::int_t& it : indicies) {
                            temp_data[it] = arr[it][i]/((util::double_t)this->sparta->run_every);
                        }

                        index = util::find(data, temp_data);
                        if (index == util::npos) {
                            _temp_vec.clear();
                            _temp_vec = {i+1};
                            ids.push_back(_temp_vec);
                            data.push_back(temp_data);
                        } else {
                            ids[index].push_back(i+1);
                        }
                    }

                    // adding the necessary boundary conditions
                    std::vector<util::double_t> _temp_vec;
                    for (std::size_t i = 0; i < data.size(); i++) {
                        elmer::Section* bc = new elmer::Section("Boundary Condition" ,bc_count++);
                        (*bc).setItem("Target_Boundaries", ids[i]);
                        (*bc).setItem("Heat_Flux", data[i][0]);
                        _temp_vec.clear();
                        for (j = 1; j < numDataPoints; j++)
                            _temp_vec.push_back(data[i][j]);
                        (*bc).setItem("Stress", _temp_vec);
                
                        this->boundary_conditions.push_back(bc);
                    }
                }
                ids.clear(); _temp_vec.clear();


                /***
                 * initial conditions
                */
                size = this->node_temperature_data.size();
    
                if (this->gravity_on.toBool()) _body_force = 1;
                else                  _body_force = toml::noInt;

                const util::int_t arr_size = elmer::dimension+2;
                std::vector<std::array<util::double_t, arr_size>> init_conds;
                std::array<util::double_t, arr_size> init_cond;
                init_conds.clear();

                // averages values for the nodes of a surface element and sets this average to the 
                // temperature of the surface element
                for (std::size_t i = 0; i < this->elements.size(); i++) {
                    // computes the average temperature of the nodes that make up the surface element
                    // this value is used to set the surface element temperature
                    init_cond.fill(0.0);
                    _size = this->elements[i].second.size();
                    init_cond[0] = this->elements[i].first;
                    for (j = 0; j < (signed)_size; j++) {
                        // gets the data point corresponding to node id and adds it to the rolling sum
                        if ((unsigned)this->elements[i].second[j]-1 >= size)
                            UERR("index out of bounds in creating initial conditions");
                        
                        init_cond[1] += (this->node_temperature_data[this->elements[i].second[j]-1]/((util::double_t)_size));
                        for (k = 0; k < elmer::dimension; k++)
                            init_cond[k+2] += (this->node_velocity_data[this->elements[i].second[j]-1][k]/((util::double_t)_size));
                    }
                    index = util::find(init_conds, init_cond);
                    if (index == util::npos) {
                        _temp_vec.clear();
                        _temp_vec = {(long)i+1};
                        ids.push_back(_temp_vec);
                        init_conds.push_back(init_cond);
                    } else {
                        ids[index].push_back(i+1);
                    }
                }
                
                for (std::size_t i = 0; i < init_conds.size(); i++) {
                    Section* _body = new Section("Body", i+1);
                    (*_body).setItem("Initial_Condition", (util::int_t)(i+1));
                    (*_body).setItem("Body_Force", _body_force);
                    (*_body).setItem("Material", (util::int_t)init_conds[i][0]);
                    (*_body).setItem("Target_Bodies", ids[i]);
                    Section* _ic = new Section("Initial Condition", i+1);
                    (*_ic).setItem("Temperature", init_conds[i][1]);
                    (*_ic).setItem("Velocity_1", init_conds[i][2]);
                    (*_ic).setItem("Velocity_2", init_conds[i][3]);
                    (*_ic).setItem("Velocity_3", init_conds[i][4]);

                    this->initial_conditions.push_back(_ic);
                    this->bodies.push_back(_body);
                }
            }

            /* ---------------------------------------------------------------------- */

            void loadNodeData() {
                ULOG("loading node data");
                std::size_t i, start;
                std::vector<util::string_t> lines, split_line;
                util::string_t cur_key, line;
                util::readFile(this->node_data_file.toString(), lines);
                start = 0;
                for (i = 0; i < lines.size(); i++) {
                    util::split(lines[i], split_line, ' ');
                    if (split_line[0] == (util::string_t)"Time:") start = i;
                }

                util::int_t cur_index = 0;
                for (i = start+1; i < lines.size(); i++) {
                    line = lines[i];
                    util::trim(line);
                    util::split(line, split_line, ' ');
                    if (lines[i][0] == ' ' && split_line.size() == 1) {
                        if (cur_key == (util::string_t)"temperature") {
                            this->node_temperature_data[cur_index] = std::stod(line);
                        } else if (cur_key == (util::string_t)"displacement 1") {
                            this->nodes[cur_index][0] += std::stod(line);
                        } else if (cur_key == (util::string_t)"displacement 2") {
                            this->nodes[cur_index][1] += std::stod(line);
                        } else if (cur_key == (util::string_t)"displacement 3") {
                            this->nodes[cur_index][2] += std::stod(line);
                        } else if (cur_key == (util::string_t)"velocity 1") {
                            this->node_velocity_data[cur_index][0] = std::stod(line);
                        } else if (cur_key == (util::string_t)"velocity 2") {
                            this->node_velocity_data[cur_index][1] = std::stod(line);
                        } else if (cur_key == (util::string_t)"velocity 3") {
                            this->node_velocity_data[cur_index][2] = std::stod(line);
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
            }

            /* ---------------------------------------------------------------------- */

            void updateNodeFile() {
                ULOG("updating node file");
                util::oFile out(this->node_file.toString());

                for (std::size_t i = 0; i < this->nodes.size(); i++)
                    out << i+1 << " " << -1 << " " << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void loadNodes() {
                ULOG("loading nodes");
                std::vector<util::string_t> lines, split_line;
                std::array<util::double_t, 3> data_item;

                util::readFile(node_file.toString(), lines);
                util::int_t count = 0;
                for (util::string_t& it : lines) {
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
                std::vector<util::string_t> lines, _split;
                std::array<util::int_t, boundary_size> arr;
                this->boundaries.clear();

                util::readFile(this->boundary_file.toString(), lines);
                util::oFile out(this->boundary_file.toString());

                util::int_t last_id = 0;
                util::int_t cur_id;
                util::int_t count = 1;
                for (util::string_t it : lines) {
                    util::trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    util::split(it, _split, ' ');

                    if (_split[4] != (util::string_t)"303") 
                        UERR("element is not a triangle in boundary file at line: " + std::to_string(count));

                    cur_id = std::stoi(_split[0]);
                    if ((cur_id - last_id) != 1)
                        UERR("got unordered boundary elements at line: " + std::to_string(count));
                    last_id = cur_id;

                    _split[1] = _split[0];
                    for (util::string_t it : _split) out << it + " ";
                    out << "\n";

                    // getting rid of stuff I do not need
                    _split.erase(_split.begin(), _split.begin() + 5);

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
                std::vector<util::string_t> _temp_data, _split;
                std::pair<util::int_t, std::vector<util::int_t>> arr;
                this->elements.clear();

                util::readFile(this->element_file.toString(), _temp_data);
                util::oFile out(this->element_file.toString());
                //this->elements.resize(_temp_data.size());

                util::int_t last_id = 0;
                util::int_t cur_id;
                util::int_t count = 1;
                for (util::string_t it : _temp_data) {
                    arr.second.clear();

                    util::trim(it);
                    if (it.length() == 0) continue;

                    // splitting the line at spaces
                    util::split(it, _split, ' ');

                    // getting rid of stuff I do not need
                    
                    arr.first = std::stoi(_split[1]);

                    cur_id = std::stoi(_split[0]);
                    if ((cur_id - last_id) != 1)
                        UERR("got unordered boundary elements at line: " + std::to_string(count));
                    
                    last_id = cur_id;

                    _split[1] = _split[0];
                    for (util::string_t it : _split) out << it + " ";
                    out << "\n";

                    _split.erase(_split.begin(), _split.begin() + 3);
                    
                    // adding the data
                    for (std::size_t i = 0; i < _split.size(); i++)
                        arr.second.push_back(std::stoi(_split[i]));

                    // adding the data to the class vector
                    this->elements.push_back(arr);
                }
            }

            /* ---------------------------------------------------------------------- */

            void join(util::oFile& _buf) {
                header.join(_buf);
                simulation.join(_buf);
                constants.join(_buf);

                std::size_t i;
                for (i = 0; i < this->solvers.size(); i++)
                    this->solvers[i]->join(_buf);

                for (i = 0; i < this->equations.size(); i++)
                    this->equations[i]->join(_buf);
                    
                for (i = 0; i < this->materials.size(); i++)
                    this->materials[i]->join(_buf);

                for (i = 0; i < this->bodies.size(); i++)
                    this->bodies[i]->join(_buf);

                for (i = 0; i < this->body_forces.size(); i++)
                    this->body_forces[i]->join(_buf);

                for (i = 0; i < this->initial_conditions.size(); i++)
                    this->initial_conditions[i]->join(_buf);

                for (i = 0; i < this->boundary_conditions.size(); i++)
                    this->boundary_conditions[i]->join(_buf);
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeTemperaturesBefore() {
                ULOG("dumping node temperatures before");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_temperature_file_ext.toString();
                util::oFile out(filename);

                for (util::double_t& it : this->node_temperature_data)
                    out << it << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeVelocitiesBefore() {
                ULOG("dumping node velocities before");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_velocity_file_ext.toString();
                util::oFile out(filename);

                for (auto& it : this->node_velocity_data)
                    out << it[0] << " " << it[1] << " " << it[2] << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodePositionsBefore() {
                ULOG("dumping node positions before");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_position_file_ext.toString();
                util::oFile out(filename);

                for (std::size_t i = 0; i < this->nodes.size(); i++) 
                    out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeTemperatures() {
                ULOG("dumping node temperatures");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_temperature_file_ext.toString();
                util::oFile out(filename);

                for (util::double_t& it : this->node_temperature_data)
                    out << it << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodePositions() {
                ULOG("dumping node positions");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_position_file_ext.toString();
                util::oFile out(filename);

                for (std::size_t i = 0; i < this->nodes.size(); i++) 
                    out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
            }

            /* ---------------------------------------------------------------------- */

            void dumpNodeVelocities() {
                ULOG("dumping node velocities");
                util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_velocity_file_ext.toString();
                util::oFile out(filename);

                for (auto& it : this->node_velocity_data)
                    out << it[0] << " " << it[1] << " " << it[2] << "\n";
            }
    };
}

#endif