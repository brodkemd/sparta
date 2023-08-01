#ifndef ELMER_H
#define ELMER_H

#include "elmer_classes.hpp"
#include "shell_server_config.h"

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
    /* ---------------------------------------------------------------------- 
     * 
     * 
     * Declarations
     * 
     * 
     * ---------------------------------------------------------------------- 
    */

    const unsigned int boundary_size = 3; // number of elements in a boundary element
    const unsigned int dimension     = 3; // the dimension
    const std::size_t node_element_size = 5;

    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            SPARTA_NS::FixFea* sparta;
            std::vector<util::double_t> node_temperature_data;
            std::vector<std::array<util::int_t, elmer::boundary_size>> boundaries;
            std::vector<util::bool_t> boundary_node_moved;
            
            // array format: x, y, z
            std::vector<std::array<util::double_t, elmer::dimension>> nodes, node_velocity_data, node_displacements;
            std::vector<std::pair<util::int_t, std::vector<util::int_t>>> elements;
            toml::Item_t exe, sif, meshDB, node_data_file, simulation_directory, shell_server_file, boundary_file, element_file, header_file, node_file, node_temperature_file_ext, node_position_file_ext, node_velocity_file_ext, print_intensity, gravity_on, base_temp;

            /* ---------------------------------------------------------------------- */

            Section header, simulation, constants;
            std::vector<Section*> solvers, materials, equations, body_forces, bodies, initial_conditions, boundary_conditions;

            /* ---------------------------------------------------------------------- */

        public:
            Elmer(SPARTA_NS::FixFea* _sparta);
            ~Elmer();
            void set(toml::handler& _h);
            /**
             * makes sure file from elmer object data
             * returns: the name of the surf file
            */
            util::string_t makeSpartaSurf();
            void averageNodeTemperaturesInto(double*& _temperatures, int _length);
            void getNodePointAtIndex(util::int_t _index, util::int_t _boundary_index, double(&_point)[elmer::dimension]);
            // void checkUpdatedAllBoundaryNodes();
            void createInitialConditions();
            void createBoundaryConditions();
            void run();
            bool shouldRun() { return true; }
            void setup();
            void dumpBefore();
            void dump();

        protected:
            void setupDumpDirectory();
            void setupVectors();
            void makeServerFile()
            void makeSif();
            void loadNodeData();
            void updateNodes();
            void updateNodeFile();
            void loadNodes();
            void loadBoundaries();
            void loadElements();
            void join(util::oFile& _buf);
            void dumpNodeTemperaturesBefore();
            void dumpNodeVelocitiesBefore();
            void dumpNodePositionsBefore();
            void dumpNodeTemperatures();
            void dumpNodePositions();
            void dumpNodeVelocities();
    };

    /* ---------------------------------------------------------------------- *
     * 
     * 
     * Definitions
     * 
     * 
     * ---------------------------------------------------------------------- *
    */

    Elmer::Elmer(SPARTA_NS::FixFea* _sparta) {
        // referencing the calling sparta instance
        this->sparta         = &*_sparta;

        // initing class variables
        this->header     = Section("Header", toml::noInt, " ", false);
        this->simulation = Section("Simulation");
        this->constants  = Section("Constants");

        this->equations.clear();    this->materials.clear();
        this->body_forces.clear();  this->bodies.clear();
        this->solvers.clear();      this->initial_conditions.clear();
        this->boundary_conditions.clear();

        this->makeServerFile();
        this->waitForServerStart();

    }

    /* ---------------------------------------------------------------------- */

    Elmer::~Elmer() {
        // this->makeExitFile();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::set(toml::handler& _h) {
        ULOG("setting elmer class parameters");
        // setting needed variables
        _h.getAtPath(this->meshDB,                    "elmer.meshDB",                     toml::STRING);
        _h.getAtPath(this->exe,                       "elmer.exe",                        toml::STRING);
        _h.getAtPath(this->sif,                       "elmer.sif",                        toml::STRING);
        _h.getAtPath(this->base_temp,                 "elmer.base_temp",                  toml::DOUBLE);
        _h.getAtPath(this->shell_server_file,         "shell_server_file",                toml::STRING);
        _h.getAtPath(this->simulation_directory,      "simulation_directory",             toml::STRING);
        _h.getAtPath(this->node_temperature_file_ext, "elmer.node_temperature_file_ext",  toml::STRING);
        _h.getAtPath(this->node_position_file_ext,    "elmer.node_position_file_ext",     toml::STRING);
        _h.getAtPath(this->node_velocity_file_ext,    "elmer.node_velocity_file_ext",     toml::STRING);
        _h.getAtPath(this->gravity_on,                "elmer.gravity_on",                 toml::BOOL);
        _h.getAtPath(this->print_intensity,           "elmer.print_intensity",            toml::STRING);

        // each elmer section sets its own variables, so passing it the data structure
        _h.getDictAtPath(this->simulation, "elmer.simulation");
        _h.getDictAtPath(this->constants, "elmer.constants");

        std::vector<toml::Item_t> keys;

        // getting the materials
        _h.getDictKeysAtPath(keys, "elmer.material");
        for (auto key : keys) {
            Section* material = new Section("Material", std::stoi(key.toString()));
            _h.getDictAtPath(*material, "elmer.material." + key.toString());
            this->materials.push_back(material);
        }
        ULOG("# of materials loaded: " + std::to_string(this->materials.size()));
        
        // getting the solvers
        _h.getDictKeysAtPath(keys, "elmer.solver");
        for (auto key : keys) {
            Section* solver = new Section("Solver", std::stoi(key.toString()));
            _h.getDictAtPath(*solver, "elmer.solver." + key.toString());
            this->solvers.push_back(solver);
        }
        ULOG("# of solvers loaded: " + std::to_string(this->solvers.size()));

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
        this->node_data_file = this->simulation.getItem("Output_File").toString();
        // this->simulation.setItem("Output_File", this->node_data_file);
        ULOG("Will load node data from: " + this->node_data_file.toString());

        toml::Item_t _temp_solvers = std::vector<toml::Item_t>({});
        for (auto& it : this->solvers)
            _temp_solvers.append(it->getId());

        Section* eqn = new Section("Equation", 1);
        eqn->setItem("Active_Solvers",  _temp_solvers); // id for each solver of the class
        this->equations.push_back(eqn);

        // setting the header
        this->header.setItem("Mesh_DB", std::vector<util::string_t>({ this->simulation_directory.toString(), "." }));
        this->header.setItem("Results_Directory", this->simulation_directory);

        // creating body force if gravity is on
        if (this->gravity_on.toBool()) {
            for (auto it : this->materials) {
                elmer::Section* bf = new elmer::Section("Body Force", it->getId());
                ULOG("gravity acts as: " + this->constants.getItem("Gravity").toString());
                if (this->constants.getItem("Gravity")[0].toDouble() != 0.0) {
                    (*bf).setItem("Stress_Bodyforce_1", this->constants.getItem("Gravity")[0].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                }
                if (this->constants.getItem("Gravity")[1].toDouble() != 0.0) {
                    (*bf).setItem("Stress_Bodyforce_2", this->constants.getItem("Gravity")[1].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                }
                if (this->constants.getItem("Gravity")[2].toDouble() != 0.0) {
                    (*bf).setItem("Stress_Bodyforce_3", this->constants.getItem("Gravity")[2].toDouble()*this->constants.getItem("Gravity")[3].toDouble()*(*it).getItem("Density").toDouble());
                }
                this->body_forces.push_back(bf);
            }
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::makeServerFile() {
        util::oFile server_file(this->shell_server_file.toString());
        server_file << SHELL_SERVER_STRING;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * makes sure file from elmer object data
     * returns: the name of the surf file
    */
    util::string_t Elmer::makeSpartaSurf() {
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

    void Elmer::averageNodeTemperaturesInto(double*& _temperatures, int _length) {
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
            for (unsigned int j = 0; j < elmer::boundary_size; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg += this->node_temperature_data[this->boundaries[i][j]-1];
            }
            // computing the average by dividing the sum by the number of points and setting the surface element
            _temperatures[i] = avg/elmer::boundary_size;
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::getNodePointAtIndex(util::int_t _index, util::int_t _boundary_index, double(&_point)[elmer::dimension]) {
        for (std::size_t j = 0; j < elmer::dimension; j++)
            _point[j] = this->nodes[this->boundaries[_index][_boundary_index]-1][j];
        //this->boundary_node_moved[_boundary_index] = true;
    }

    /* ---------------------------------------------------------------------- */

    // void checkUpdatedAllBoundaryNodes() {
    //     for (std::size_t i = 0; i < this->boundary_node_moved.size(); i++) {
    //         if (!(this->boundary_node_moved[i]))
    //             UERR("Boundary node did not get updated at index: " + std::to_string(i+1));
    //     }
    //     ULOG("Updated all nodes");
    //     boundary_node_moved.resize(this->boundaries.size(), false);
    // }

    /* ---------------------------------------------------------------------- */

    void Elmer::createInitialConditions() {
        ULOG("Creating Initial Conditions");
        this->bodies.clear();
        this->initial_conditions.clear();

        std::size_t _size;
        util::int_t _body_force, j, k, index;
        
        const std::size_t numDataPoints = elmer::dimension+2;
        std::vector<std::array<util::double_t, numDataPoints>> data; data.clear();
        std::array<util::double_t, numDataPoints> temp_data;
        std::vector<std::vector<util::int_t>> ids;  ids.clear();
        
        if (this->gravity_on.toBool()) _body_force = 1;
        else                           _body_force = toml::noInt;

        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        for (std::size_t i = 0; i < this->elements.size(); i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            _size = this->elements[i].second.size();
            temp_data.fill(0.0);
            temp_data[0] = this->elements[i].first;
            for (j = 0; j < (signed)_size; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                // if ((unsigned)this->elements[i].second[j]-1 >= size)
                //     UERR("index out of bounds in creating initial conditions");
                
                temp_data[1] += (this->node_temperature_data[this->elements[i].second[j]-1]/((util::double_t)_size));
                for (k = 0; k < elmer::dimension; k++)
                    temp_data[k+2] += (this->node_velocity_data[this->elements[i].second[j]-1][k]/((util::double_t)_size));
            }
            index = util::find(data, temp_data);
            if (index == util::npos) {
                ids.push_back((std::vector<util::int_t>){(util::int_t)i+1});
                data.push_back(temp_data);
            } else { ids[index].push_back(i+1); }
        }

        for (std::size_t i = 0; i < data.size(); i++) {
            Section* _body = new Section("Body", i+1);
            (*_body).setItem("Initial_Condition", (util::int_t)(i+1));
            (*_body).setItem("Body_Force",        _body_force);
            (*_body).setItem("Equation",          (util::int_t)1);
            (*_body).setItem("Material",          (util::int_t)data[i][0]);
            (*_body).setItem("Target_Bodies",     ids[i]);

            Section* _ic = new Section("Initial Condition", i+1);
            (*_ic).setItem("Temperature", data[i][1]);
            (*_ic).setItem("Velocity_1",  data[i][2]);
            (*_ic).setItem("Velocity_2",  data[i][3]);
            (*_ic).setItem("Velocity_3",  data[i][4]);

            this->initial_conditions.push_back(_ic);
            this->bodies.push_back(_body);
        }
        ULOG("# of initial conditions and bodies: " + std::to_string(this->initial_conditions.size()));
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::createBoundaryConditions() {
        ULOG("Creating Boundary Conditions");
        this->boundary_conditions.clear();

        util::int_t index;
        util::int_t bc_count = 1;
        const util::int_t numDataPoints = 7;
        // get the stress values from sparta
        double* arr[numDataPoints] = {
            &*this->sparta->qw,
            &*this->sparta->px,  &*this->sparta->py,  &*this->sparta->pz,
            &*this->sparta->shx, &*this->sparta->shy, &*this->sparta->shz
        };

        std::vector<std::array<util::double_t, numDataPoints>> data; data.clear();
        std::array<util::double_t, numDataPoints> temp_data;
        std::vector<util::double_t> _temp_vec; _temp_vec.clear();
        std::vector<std::vector<util::int_t>> ids;  ids.clear();
        
        
        // detecting unique values and saving them to the stress vector and
        // their corresponding ids to the ids vector, this condenses the
        // number of boundary conditions (removes boundary conditions with
        // the same values)
        util::int_t i, j;
        for (i = 0; i < this->sparta->nsurf; i++) {
            for (j = 0; j < numDataPoints; j++) temp_data[j] = arr[j][i];

            index = util::find(data, temp_data);
            if (index == util::npos) {
                ids.push_back((std::vector<util::int_t>){i+1});
                data.push_back(temp_data);
            } else { ids[index].push_back(i+1); }
        }

        // adding the necessary boundary conditions
        for (std::size_t i = 0; i < data.size(); i++) {
            elmer::Section* bc = new elmer::Section("Boundary Condition", bc_count++);
            (*bc).setItem("Target_Boundaries", ids[i]);
            (*bc).setItem("Heat_Flux", data[i][0]);
            _temp_vec.clear();
            for (j = 1; j < numDataPoints; j++)
                _temp_vec.push_back(data[i][j]);
            (*bc).setItem("Stress", _temp_vec);

            this->boundary_conditions.push_back(bc);
        }
        ULOG("# of boundary conditions: " + std::to_string(this->boundary_conditions.size()));
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::run() {
        // setting the timestep size
        bool has_timestep_sizes = true;
        if (!(this->simulation.hasKey("Timestep_Sizes"))) {
            has_timestep_sizes = false;
            ULOG("did not provide \"timestep size\", using Sparta's");
            this->simulation.setItem("Timestep_Sizes", (std::vector<util::double_t>){ (util::double_t)this->sparta->update->dt });
        } else
            ULOG("Using provided \"timestep size\"");

        // setting number of timesteps to run
        bool has_timestep_intervals = true;
        if (!(this->simulation.hasKey("Timestep_intervals"))) {
            has_timestep_intervals = false;
            ULOG("did not provide \"number of timesteps\", using Sparta's");
            this->simulation.setItem("Timestep_intervals", (std::vector<util::int_t>){ (util::int_t)this->sparta->run_every });
        } else
            ULOG("Using provided \"number of timesteps\"");
        
        bool has_output_intervals = true;
        if (!(this->simulation.hasKey("Output_Intervals"))) {
            has_output_intervals = false;
            ULOG("did not provide \"output interval\", using Sparta's");
            this->simulation.setItem("Output_Intervals", (std::vector<util::int_t>){ (util::int_t)this->sparta->run_every });
        } else
            ULOG("Using provided \"output interval\"");

        
        // making the sif file
        ULOG("Making sif file for elmer");
        this->makeSif();

        ULOG("running");

        // running the command
        // must have " 2>&1" at end to pipe stderr to stdout
        int exit_status = 0;
        try {
            // this->makeRunFile();
            // this->waitForDoneFile();
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
        
        if (!(has_timestep_sizes))
            this->simulation.removeKey("Timestep_Sizes");
        
        if (!(has_timestep_intervals))
            this->simulation.removeKey("Timestep_intervals");
        
        if (!(has_output_intervals))
            this->simulation.removeKey("Output_Intervals");
        
        this->loadNodeData();
        //this->checkLoadedNodeData();
        this->updateNodes();
        this->updateNodeFile();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::setup() {
        this->setupDumpDirectory();
        this->setupVectors();
        this->loadBoundaries();
        this->loadElements();
        this->loadNodes();

        boundary_node_moved.resize(this->boundaries.size(), false);
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpBefore() { 
        this->dumpNodePositionsBefore();
        this->dumpNodeTemperaturesBefore(); 
        this->dumpNodeVelocitiesBefore();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dump() {
        this->dumpNodeTemperatures();
        this->dumpNodePositions();
        this->dumpNodeVelocities();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::setupDumpDirectory() {
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
        this->header_file   = this->simulation_directory.toString() + SEP + "mesh.header";
    }

    /* ---------------------------------------------------------------------- */
    
    void Elmer::setupVectors() {
        // getting the number of nodes
        util::iFile h = util::iFile(this->header_file.toString());
        std::vector<util::string_t> lines, split;
        h.getLines(lines);
        util::split(lines[0], split, ' ');
        util::int_t count = std::stoi(split[0]);
        ULOG("Creating vectors with length: " + std::to_string(count));

        std::array<util::double_t, elmer::dimension> arr;
        arr.fill(0.0);

        // writing base temperature to file for each node
        for (util::int_t i = 0; i < count; i++) {
            this->node_temperature_data.push_back(this->base_temp.toDouble());
            this->node_velocity_data.push_back(arr);
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::makeSif() {
        util::oFile _buf(this->sif.toString());
        _buf << "! File created: " << util::getTime() << "\n";
        _buf << "! File made by process: " << util::_me << "\n\n";
        _buf << "Check Keywords warn\n";
        this->join(_buf);
    }

    /* ---------------------------------------------------------------------- */

    // void Elmer::checkLoadedNodeData() {
    //     util::int_t j;
    //     for (std::size_t i = 0; i < this->node_displacements.size(); i++) {
    //         for (j = 0; j < elmer::dimension; j++) {
    //             if (this->node_displacements[i][j] != 0.0) {
    //                 UERR("Did not load displacement values for node with id: " + std::to_string(i+1));
    //             }
    //         }
    //     }
    //     for (std::size_t i = 0; i < this->node_temperature_data.size(); i++) {
    //         if (node_temperature_data[i] <= 0.0) {
    //             UERR("Did not load temperature correct temperature value for node with id: " + std::to_string(i+1));
    //         }
    //     }
    //     ULOG("Node Parameters Loaded Successfully");
    // }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadNodeData() {
        ULOG("loading node data");
        std::size_t i, start;
        std::vector<util::string_t> lines, split_line;
        util::string_t cur_key, line;
        util::int_t cur_index, count, index;
        

        ULOG("Loading data from: " + this->node_data_file.toString());
        util::readFile(this->node_data_file.toString(), lines);
        toml::Dict_t<util::int_t> permutation_table;
        start = 0;
        const util::int_t list_length = (util::int_t)this->node_temperature_data.size();
        std::vector<util::bool_t> check_list;
        check_list.resize(list_length, true);

        // gets latest data and permutation table
        for (i = 0; i < lines.size(); i++) {
            util::trim(lines[i]);
            util::split(lines[i], split_line, ' ');
            if (split_line[0] == (util::string_t)"Time:") start = i;
            if (split_line[0] == (util::string_t)"Perm:") {
                if (split_line[1] == "use") {
                    ULOG("Got empty permutation table at line: " + std::to_string(i));
                    continue;
                } else {

                    util::split(lines[i-2], split_line, ' ');
                    if (split_line[0] != (util::string_t)"Time:")
                        UERR("got another permutation table from data file at line: " + std::to_string(i+1) + ", this is not currently support");
                    
                    ULOG("Loading permutation table at line (will overwrite existing): " + std::to_string(i+1));
                    permutation_table.clear();

                    for (i = i + 1; i < lines.size(); i++) {
                        util::trim(lines[i]);
                        util::split(lines[i], split_line, ' ');
                        if (split_line.size() < (std::size_t)2) {
                            i--;
                            break;
                        }
                        permutation_table.setKey(std::stoi(split_line[0]), std::stoi(split_line[1]));
                    }
                }
            }
        }

        cur_index = 0;
        check_list.resize(list_length, true);
        for (i = start+1; i < lines.size(); i++) {
            line = lines[i];
            util::trim(line);
            util::split(line, split_line, ' ');
            if ((std::isdigit(lines[i][0]) || lines[i][0] == '-') && split_line.size() == 1) {
                index = permutation_table.getKey(cur_index+1) - 1;
                if        (cur_key == (util::string_t)"temperature") {
                    this->node_temperature_data[index] = std::stod(line);
                } else if (cur_key == (util::string_t)"displacement 1") {
                    this->node_displacements[index][0] += std::stod(line);
                } else if (cur_key == (util::string_t)"displacement 2") {
                    this->node_displacements[index][1] += std::stod(line);
                } else if (cur_key == (util::string_t)"displacement 3") {
                    this->node_displacements[index][2] += std::stod(line);
                } else if (cur_key == (util::string_t)"velocity 1") {
                    this->node_velocity_data[index][0] = std::stod(line);
                } else if (cur_key == (util::string_t)"velocity 2") {
                    this->node_velocity_data[index][1] = std::stod(line);
                } else if (cur_key == (util::string_t)"velocity 3") {
                    this->node_velocity_data[index][2] = std::stod(line);
                } else
                    UERR("got unknown key in data file with name: " + cur_key);
                check_list[index] = true;
                cur_index++;
            } else if ((util::string_t)"Perm:" == split_line[0]) {
                if (cur_index > 0) {
                    if ((std::size_t)cur_index != this->nodes.size())
                        UERR("did not get all indicies while loading \"" + cur_key + "\" from data");
                    
                    count = 1;
                    for (auto it : check_list) {
                        if (!(it)) {
                            UERR("Did not get value at index \"" + std::to_string(count) + "\" while loading \"" + cur_key + "\"");
                        }
                        count++;
                    }
                    ULOG("Parameter \"" + cur_key + "\" completed checks");
                }

                check_list.resize(list_length, false);
                cur_index = 0;
                cur_key = lines[i-1];
                util::trim(cur_key);
                ULOG("Loading data from node data file for: " + cur_key);
            }
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::updateNodes() {
        ULOG("Updating nodes");
        util::int_t i, j;
        util::double_t before, after;
        util::bool_t got_no_difference;
        std::vector<util::bool_t> checks;
        checks.resize(this->nodes.size(), false);

        for (i = 0; i < (signed)this->nodes.size(); i++) {
            got_no_difference = false;
            for (j = 0; j < elmer::dimension; j++) {
                before = this->nodes[i][j];
                after  = before + this->node_displacements[i][j];
                if (before == after && this->node_displacements[i][j] != 0.0) {
                    got_no_difference = true;
                    break;
                }
            }

            checks[i] = true;

            if (got_no_difference) {
                ULOG("did not update node point: " + std::to_string(i+1) + ", remembering displacement");
            } else {
                for (j = 0; j < elmer::dimension; j++) {
                    this->nodes[i][j] += this->node_displacements[i][j];
                    this->node_displacements[i][j] = 0.0;
                }
            }
        }

        for (i = 0; i < (signed)checks.size(); i++) {
            if (!(checks[i])) {
                UERR("Did not load node at index: " + std::to_string(i+1));
            }
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::updateNodeFile() {
        ULOG("updating node file");
        util::oFile out(this->node_file.toString());

        for (std::size_t i = 0; i < this->nodes.size(); i++)
            out << i+1 << " " << -1 << " " << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadNodes() {
        ULOG("loading nodes");
        std::vector<util::string_t> lines, split_line;
        std::array<util::double_t, 3> data_item;
        
        this->nodes.clear();
        this->node_displacements.clear();
        ULOG("Loading raw data");
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
            data_item.fill(0.0);
            this->node_displacements.push_back(data_item);
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadBoundaries() {
        ULOG("loading boundaries");
        std::vector<util::string_t> lines, _split;
        std::array<util::int_t, boundary_size> arr;
        this->boundaries.clear();

        ULOG("Loading raw data");
        util::readFile(this->boundary_file.toString(), lines);
        ULOG("Itemizing boundary");
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

    void Elmer::loadElements() {
        ULOG("loading elements");
        std::vector<util::string_t> _temp_data, _split;
        std::pair<util::int_t, std::vector<util::int_t>> arr;
        this->elements.clear();

        ULOG("Loading raw data");
        util::readFile(this->element_file.toString(), _temp_data);
        ULOG("Itemizing elements");
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

    void Elmer::join(util::oFile& _buf) {
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

    void Elmer::dumpNodeTemperaturesBefore() {
        ULOG("dumping node temperatures before");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_temperature_file_ext.toString();
        util::oFile out(filename);

        for (util::double_t& it : this->node_temperature_data)
            out << it << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeVelocitiesBefore() {
        ULOG("dumping node velocities before");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_velocity_file_ext.toString();
        util::oFile out(filename);

        for (auto& it : this->node_velocity_data)
            out << it[0] << " " << it[1] << " " << it[2] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodePositionsBefore() {
        ULOG("dumping node positions before");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_position_file_ext.toString();
        util::oFile out(filename);

        for (std::size_t i = 0; i < this->nodes.size(); i++) 
            out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeTemperatures() {
        ULOG("dumping node temperatures");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_temperature_file_ext.toString();
        util::oFile out(filename);

        for (util::double_t& it : this->node_temperature_data)
            out << it << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodePositions() {
        ULOG("dumping node positions");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_position_file_ext.toString();
        util::oFile out(filename);

        for (std::size_t i = 0; i < this->nodes.size(); i++) 
            out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeVelocities() {
        ULOG("dumping node velocities");
        util::string_t filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + "." + this->node_velocity_file_ext.toString();
        util::oFile out(filename);

        for (auto& it : this->node_velocity_data)
            out << it[0] << " " << it[1] << " " << it[2] << "\n";
    }
}

#endif