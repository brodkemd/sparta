#ifndef ELMER_H
#define ELMER_H

#include "server.hpp"
#include "mpi.h"
#include <cstdlib>
#include <unordered_map>

#define ELMER_LINE_START "  "
#define ELMER_ARRAY_SEP " "
#define ELMER_SECTION_END "End"

#define ELMER_DIMENSION 3

#define BOUNDARY_ELEMENT_START_READ 1 // index for a line in the boundary file to start reading data from (including it)
#define BOUNDARY_ELEMENT_LINE_SIZE 8
#define BOUNDARY_ELEMENT_ARRAY_SIZE BOUNDARY_ELEMENT_LINE_SIZE-BOUNDARY_ELEMENT_START_READ
#define BOUNDARY_ELEMENT_ARRAY_NODE_START 4 // is an index starting from second column

// elements are not restricted to a size
#define ELEMENT_START_READ 1 // index for a line in the boundary file to start reading data from (including it)
#define ELEMENT_ARRAY_NODE_START 2 // is an index starting from second column

// difference of these two should equal ELMER_DIMENSION
#define NODE_START_READ 2 // index for a line in the boundary file to start reading data from (including it)
#define NODE_LINE_SIZE 5


/*
Notes:
    - On the elements file
        - the first column is removed when read
    - On the boundary file
        - the first column is removed when read
    - On the nodes file
        - the first and second column are remove when read

*/



/* ---------------------------------------------------------------------- */

// defining the os path separator based on the detected os
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    const char SEP = '\\';
#else
    const char SEP = '/';
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

    // const unsigned int boundary_size    = 7; // number of elements in a boundary element
    // const unsigned int dimension        = 3; // the dimension
    // const std::size_t node_element_size = 5;


    
    /**
     * Base class that is used to construct the sections in a sif file
    */
    class Section {
        private:
            std::string _name, _sep;

            // ends the section
            long _id;
            bool _include_count;
            std::vector<std::array<std::string, 2>> _content_pairs;

        public:            
            Section() {}
            Section(std::string name, long id = util::NO_INT, std::string sep = " = ", bool include_count = true) {
                this->_name = name;
                this->_sep = sep;
                this->_id = id;
                // this->contents = toml::OrderedDict_t();
                this->_include_count = include_count;
            }

            void join(util::oFile& _buf) {
                if (this->_id == util::NO_INT)
                    _buf << this->_name << "\n";
                else
                    _buf << this->_name << " " << _id << "\n";
                
 
                for (auto it : this->_content_pairs)
                    _buf << ELMER_LINE_START << it[0] << this->_sep << it[1] << "\n";

                _buf << ELMER_SECTION_END << "\n";
            }

            void addEquality(std::string var, long val) {
                this->_content_pairs.push_back({var, std::to_string(val)});
            }

            // void addEquality(std::string var, unsigned long val) {
            //     this->_content_pairs.push_back({var, std::to_string(val)});
            // }

            void addEquality(std::string var, double val) {
                this->_content_pairs.push_back({var, util::dtos(val)});
            }

            void addEquality(std::string var, std::vector<long>& val, bool include_count = true) {
                std::array<std::string, 2> buf;
                if (include_count) {
                    buf[0] = var + "(" + std::to_string(val.size()) + ")";
                }

                buf[1] = "";

                if (val.size()) {
                    buf[1] = std::to_string(val[0]);
                    for (std::size_t i = 1; i < val.size(); i++){
                        buf[1] += " " + std::to_string(val[i]);
                    }
                } else {
                    UERR("can not interpret vector with no length");
                }
                this->_content_pairs.push_back(buf);
            }

            void addEquality(std::string var, std::vector<double>& val, bool include_count = true) {
                std::array<std::string, 2> buf;
                if (include_count) {
                    buf[0] = var + "(" + std::to_string(val.size()) + ")";
                } else 
                    buf[0] = var;

                buf[1] = "";

                if (val.size()) {
                    buf[1] = util::dtos(val[0]);
                    for (std::size_t i = 1; i < val.size(); i++){
                        buf[1] += " " + util::dtos(val[i]);
                    }
                } else {
                    UERR("can not interpret vector with no length");
                }
                this->_content_pairs.push_back(buf);
            }

            void addEquality(std::string var, double*& arr, long len, bool include_count = true) {
                std::array<std::string, 2> buf;
                if (include_count) {
                    buf[0] = var + "(" + std::to_string(len) + ")";
                } else 
                    buf[0] = var;

                buf[1] = "";

                if (len) {
                    buf[1] = util::dtos(arr[0]);
                    for (long i = 1; i < len; i++){
                        buf[1] += " " + util::dtos(arr[i]);
                    }
                } else {
                    UERR("can not interpret array with no length");
                }
                this->_content_pairs.push_back(buf);
            }

            // long getId() { return this->_id; }
    };


    /**
     * the main elmer class, this handles almost everything to do with elmer
    */
    class Elmer {
        private:
            SPARTA_NS::FixFea* sparta;
            python::handler*   python;

            MPI_Comm world;
            elmer::Server server;

            long **elements, **boundary, *element_lengths, me, body_force_id, equation_id, material_id, num_nodes, num_elements, num_boundary;//, num_user_boundary_conditions, start_boundary_condition_count_at;

            bool *allow_heat_flux, *allow_forces, *clamped;

            double **nodes, **node_velocities, **node_displacements,  *node_temperatures, base_temp;
            
            char *sif, *meshDB, *node_data_file, *simulation_directory, *boundary_file, *element_file, *header_file, *node_file, *node_temperature_file_ext, *node_position_file_ext, *node_velocity_file_ext, *sif_file_config_str;
            
            bool gravity_on, has_been_initialized = false;

            /* ---------------------------------------------------------------------- */

            PyObject *header, *simulation, *equation, *body_force, *elmer, *material;//, *user_specified_boundary_conditions;
            std::vector<Section*> bodies, initial_conditions, boundary_conditions;

            /* ---------------------------------------------------------------------- */

        public:
            Elmer(SPARTA_NS::FixFea* _sparta, long _me, MPI_Comm _world, python::handler* _h);
            ~Elmer();

            /**
             * makes sure file from elmer object data
             * returns: the name of the surf file
            */
            void makeSpartaSurf(char* fname);
            void averageNodeTemperaturesInto(double*& _temperatures, long _length);
            void getNodePointAtIndex(long _index, long _boundary_index, double(&_point)[ELMER_DIMENSION]);
            void createInitialConditions();
            void createBoundaryConditions();
            void run();
            bool shouldRun() { return true; }
            // void dumpBefore();
            void dump();

        protected:
            void handleUserBoundaryConditions();
            // void setupDumpDirectory();
            void setupVectors();
            // void makeSif();
            void loadNodeData();
            void updateNodes();
            void updateNodeFile();
            void loadNodes();
            void loadBoundaries();
            void loadElements();
            // void join(util::oFile& _buf);
            // void dumpNodeTemperaturesBefore();
            // void dumpNodeVelocitiesBefore();
            // void dumpNodePositionsBefore();
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

    Elmer::Elmer(SPARTA_NS::FixFea* _sparta, long _me, MPI_Comm _world, python::handler* _h) {
        ULOG("setting up elmer class");
        // referencing the calling sparta instance
        this->sparta     = &*_sparta;
        this->me         = _me;
        this->world      = _world;
        this->python     = &*_h;
        this->server     = elmer::Server(_h);

        ULOG("setting up elmer class");
        
        // setting needed variables
        this->elmer = this->python->loadObjectWithSetupFromMain("elmer");
        python::loadAttrFromObjectAndConvert(this->elmer, "meshDB",                     this->meshDB);
        python::loadAttrFromObjectAndConvert(this->elmer, "sif",                        this->sif);
        python::loadAttrFromObjectAndConvert(this->elmer, "base_temperature",           this->base_temp);
        python::loadAttrFromObjectAndConvert(this->elmer, "node_temperature_file_ext",  this->node_temperature_file_ext);
        python::loadAttrFromObjectAndConvert(this->elmer, "node_position_file_ext",     this->node_position_file_ext);
        python::loadAttrFromObjectAndConvert(this->elmer, "node_velocity_file_ext",     this->node_velocity_file_ext);
        python::loadAttrFromObjectAndConvert(this->elmer, "gravity_on",                 this->gravity_on);
        python::loadAttrFromObjectAndConvert(this->elmer, "_header_file",               this->header_file);
        python::loadAttrFromObjectAndConvert(this->elmer, "_node_file",                 this->node_file);
        python::loadAttrFromObjectAndConvert(this->elmer, "_element_file",              this->element_file);
        python::loadAttrFromObjectAndConvert(this->elmer, "_boundary_file",             this->boundary_file);
        // python::loadAttrFromObjectAndConvert(this->elmer, "_start_boundary_condition_count_at", this->start_boundary_condition_count_at);

        // getting runtime information
        PyObject* runtime = this->python->loadObjectWithSetupFromMain("runtime");
        python::loadAttrFromObjectAndConvert(runtime, "simulation_directory", this->simulation_directory);

        // each elmer section sets its own variables, so passing it the data structure
        this->simulation                         = this->python->loadObjectWithSetupFromObject(this->elmer, "simulation");
        this->equation                           = this->python->loadObjectWithSetupFromObject(this->elmer, "equation");
        this->material                           = this->python->loadObjectWithSetupFromObject(this->elmer, "material");
        this->body_force                         = this->python->loadObjectWithSetupFromObject(this->elmer, "body_force");
        

        // getting ids from attributes, used to set conditions
        python::loadAttrFromObjectAndConvert(this->equation, "_id", this->equation_id);
        python::loadAttrFromObjectAndConvert(this->material, "_id", this->material_id);

        if (this->gravity_on)
            python::loadAttrFromObjectAndConvert(this->body_force, "_id", this->body_force_id);
        else
            this->body_force_id = util::NO_INT;

        // setting vars from provided values
        python::loadAttrFromObjectAndConvert(this->simulation, "Output_File", this->node_data_file);

        python::convertObjectToString(this->elmer, this->sif_file_config_str);
        long length = strlen(this->sif_file_config_str);
        if (length < 3)
            UERR("Str method for elmer class returned a string with no length");
        
        ULOG("Will load node data from: " + std::string(this->node_data_file));

        // creates all of the vectors
        this->setupVectors();
        this->handleUserBoundaryConditions();

        // setting up server
        this->server.makeServerFile();
        this->server.start();

    }

    /* ---------------------------------------------------------------------- */

    Elmer::~Elmer() {
        delete [] nodes;
        delete [] elements;
        delete [] element_lengths;
        delete [] boundary;       
        delete [] node_velocities;
        delete [] node_temperatures;
        delete [] node_displacements;
        delete [] allow_forces;
        delete [] allow_heat_flux;
        delete [] clamped;
        this->server.end();
    }

    /* ---------------------------------------------------------------------- */

    /**
     * makes sure file from elmer object data
     * returns: the name of the surf file
    */
    void Elmer::makeSpartaSurf(char* path) {
        ULOG("making surface file for sparta");
        // char* path;
        // int result = asprintf(&path, "%s%c%s", this->simulation_directory, SEP, fname);
        // if (result == -1) { UERR("asprintf failed for path formatting"); }
        
        util::oFile f(path);
        f << "# Surface element file written by SPARTA fea interface\n\n" << this->num_nodes << " points\n" << this->num_boundary << " triangles\n\nPoints\n\n";

        long i, j;
        for (i = 0; i < this->num_nodes; i++) {
            f << i+1;
            for (j = 0; j < ELMER_DIMENSION; j++)
                f << " " << this->nodes[i][j];
            f << "\n";
        }

        f << "\nTriangles\n\n";
        for (i = 0; i < this->num_boundary; i++) {
            f << i+1;
            for (j = BOUNDARY_ELEMENT_ARRAY_NODE_START; j < BOUNDARY_ELEMENT_ARRAY_SIZE; j++)
                f << " " << this->boundary[i][j];
            f << "\n";
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::averageNodeTemperaturesInto(double*& _temperatures, long _length) {
        // making sure everything is consistent
        if (_length != this->num_boundary)
            UERR("boundary data does not match required size, required size " + std::to_string(_length) + ", size " + std::to_string(this->num_boundary));
        
        ULOG("averaging the per node temperatures into surf");

        // used later
        double avg;
        long i, j;
        
        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        for (i = 0; i < this->num_boundary; i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            avg = 0;
            // gets the data point corresponding to node id and adds it to the rolling sum
            for (j = BOUNDARY_ELEMENT_ARRAY_NODE_START; j < BOUNDARY_ELEMENT_ARRAY_SIZE; j++)
                avg += this->node_temperatures[this->boundary[i][j]-1];

            // computing the average by dividing the sum by the number of points and setting the surface element, the denominator should be 3 in the current implementation
            _temperatures[i] = avg/(BOUNDARY_ELEMENT_ARRAY_SIZE-BOUNDARY_ELEMENT_ARRAY_NODE_START); 
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::getNodePointAtIndex(long _index, long _boundary_index, double(&_point)[ELMER_DIMENSION]) {
        for (std::size_t j = 0; j < ELMER_DIMENSION; j++)
            _point[j] = this->nodes[this->boundary[_index][_boundary_index]-1][j];
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::createInitialConditions() {
        ULOG("Creating Initial Conditions");
        this->bodies.clear();
        this->initial_conditions.clear();

        long i, j, k;
        long index;

        std::vector<std::array<double,  ELMER_DIMENSION+1>> data; data.clear();
        std::array<double, ELMER_DIMENSION+1> avg;
        std::vector<std::vector<long>> indices; indices.clear();

        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        ULOG("Setting surface temperatures");
        for (long i = 0; i < this->num_elements; i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            avg.fill(0.0);           
            for (j = ELEMENT_ARRAY_NODE_START; j < this->element_lengths[i]; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg[0] += this->node_temperatures[this->elements[i][j]-1];
                for (k = 0; k < ELMER_DIMENSION; k++)
                    avg[k+1] += this->node_velocities[this->elements[i][j]-1][k];
            }

            for (j = 0; j < ELMER_DIMENSION; j++)
                avg[j] /= (this->element_lengths[i] - ELEMENT_ARRAY_NODE_START);

            index = util::find(data, avg);
            if (index == util::npos) {
                indices.push_back((std::vector<long>){i});
                data.push_back(avg);
            } else
                indices[index].push_back(i);
        }

        util::oFile out(this->element_file);
        long count = 1;
        for (i = 0; i < (long)indices.size(); i++) {
            for (j = 0; j < (long)indices[i].size(); j++) {
                out << count << " " << i + 1;
                for (k = 1; k < this->element_lengths[indices[i][j]]; k++)
                    out << " " << this->elements[indices[i][j]][k];
                out << "\n";
                count++;
            }
        }
        out.close();

        ULOG("Formatting sections");
        for (i = 0; i < (long)data.size(); i++) {
            Section* ic = new Section("Initial Condition", i+1);
            ic->addEquality("Temperature", data[i][0]);
            ic->addEquality("Velocity 1", data[i][1]);
            ic->addEquality("Velocity 2", data[i][2]);
            ic->addEquality("Velocity 3", data[i][3]);

            Section* body = new Section("Body", i+1);
            
            body->addEquality("Initial Condition", i+1);
            
            if (body_force_id != util::NO_INT)
                body->addEquality("Body Force", this->body_force_id);

            body->addEquality("Equation", this->equation_id);
            body->addEquality("Material", this->material_id);
            body->addEquality("Target Bodies", i+1);

            this->initial_conditions.push_back(ic);
            this->bodies.push_back(body);
        }
        ULOG("# of initial conditions and bodies: " + std::to_string(this->initial_conditions.size()));
        
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::createBoundaryConditions() {
        ULOG("Creating Boundary Conditions");

        long len, i, j, k, count = 1;
        double* temp_arr = new double[7];

        

        // std::vector<long> user_specified_boundary_condition_ids(size);
        // std::vector<long> user_specified_boundary_condition_surface_group_ids;
        // std::unordered_map<long, long> user_specified_boundary_condition_surface_group_map;
        // std::vector<bool> user_specified_boundary_condition_allow_heat(size);


        std::vector<std::string> hashes(this->num_boundary); // num_boundary = sparta->nsurf

        ULOG("Calculating Hashes");
        for (i = 0; i < this->num_boundary; i++) {
            len = 0;
            if (allow_heat_flux[i]) {
                temp_arr[len] = this->sparta->qw[i];
                len += 1;
            }

            if (allow_forces[i]) {
                temp_arr[len]   = this->sparta->px[i];
                temp_arr[len+1] = this->sparta->py[i];
                temp_arr[len+2] = this->sparta->pz[i];
                temp_arr[len+3] = this->sparta->shx[i];
                temp_arr[len+4] = this->sparta->shy[i];
                temp_arr[len+5] = this->sparta->shz[i];
                len+=6;
            } else if (clamped[i]) {
                // adding two so that no other combination of add values can be the same
                // i.e. so that for each combination of bools, a unique len is produced,
                // allowing for a guaranteed unique hash
                temp_arr[len]   = 0.0;
                temp_arr[len+1] = 0.0;
                len+=2;
            }

            hashes[i] = util::doubleArrayToHash(temp_arr, len);
        }

        // detecting unique values and saving them to the stress vector and
        // their corresponding indices to the indices vector, this condenses the
        // number of boundary conditions (removes boundary conditions with
        // the same values)
        std::unordered_map<std::string, std::vector<long>> data_dict;
        ULOG("Setting up Boundary Map");

        // setting all the hashes in the map, this essentially filters
        // for unique values (of surface parameters)
        for (i = 0; i < this->num_boundary; i++) {
            data_dict[hashes[i]] = std::vector<long>({});
        }

        // this step adds the indicies that correspond to each unique 
        // value (of surface parameters)
        for (i = 0; i < this->num_boundary; i++) {
            data_dict[hashes[i]].push_back(i);
        }

        // updating the boundary file, grouping all faces that have the same value into 
        // the same boundary groups
        ULOG("Generating Boundary Conditions");
        util::oFile out(this->boundary_file);
        
        // Elmer starts counting at 1 (not 0), so set these to 1
        i = 1; count = 1;
        ULOG("Boundary Conditions start with id:" + std::to_string(i));

        // clearing for good measure
        this->boundary_conditions.clear();

        // iterate over the map of unique surface parameters
        for (const auto& [key, value] : data_dict) {
            // writing to the boundary file, count is the index of the surface element
            // in the file and i is the boundary group. The rest is added from the 
            // originally loaded boundary file
            for (j = 0; j < (long)value.size(); j++) {
                out << count << " " << i;
                for (k = 1; k < BOUNDARY_ELEMENT_ARRAY_SIZE; k++)
                    out << " " << this->boundary[value[j]][k];
                out << "\n";
                count++;
            }
            
            // the rest of the loop is straight forward, create a new boundary condition
            // and add it to the
            elmer::Section* bc = new elmer::Section("Boundary Condition", i);
            
            bc->addEquality("Target Boundaries", i);
            // all data points at any index in the values vector has the same
            // surface parameters (by construction), so just pick the first one
            if (allow_heat_flux[value[0]])
                bc->addEquality("Heat Flux", this->sparta->qw[value[0]]);

            if (allow_forces[value[0]]) {
                temp_arr[0] = this->sparta->px[value[0]];
                temp_arr[1] = this->sparta->py[value[0]];
                temp_arr[2] = this->sparta->pz[value[0]];
                temp_arr[3] = this->sparta->shx[value[0]];
                temp_arr[4] = this->sparta->shy[value[0]];
                temp_arr[5] = this->sparta->shz[value[0]];
                bc->addEquality("Stress", temp_arr, 6, true);
            } else if (clamped[value[0]]) {              
                bc->addEquality("Displacement 1", (long)0);
                bc->addEquality("Displacement 2", (long)0);
                bc->addEquality("Displacement 3", (long)0);
            }
            
            this->boundary_conditions.push_back(bc);
            i++;
        }
        delete [] temp_arr;

        out.close();

        ULOG("# of boundary conditions: " + std::to_string(this->boundary_conditions.size()));
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::run() {
        if (!has_been_initialized) {

            PyObject* var;
            var = python::loadAttrFromObject(this->elmer, "couple_ratio", python::PYFLOAT, true);
            if (Py_IsNone(var)) {
                var = python::loadAttrFromObject(this->simulation, "Timestep_Sizes", python::PYFLOAT, true);
                if (Py_IsNone(var)) {
                    ULOG("did not provide \"timestep size\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Timestep_Sizes", this->sparta->update->dt);
                } else
                    ULOG("Using provided \"timestep size\"");

                // setting number of timesteps to run
                var = python::loadAttrFromObject(this->simulation, "Timestep_intervals", python::PYINT, true);
                if (Py_IsNone(var)) {
                    ULOG("did not provide \"number of timesteps\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Timestep_intervals", (long)this->sparta->run_every);
                } else
                    ULOG("Using provided \"number of timesteps\"");
                

                var = python::loadAttrFromObject(this->simulation, "Output_Intervals", python::PYINT, true);
                if (Py_IsNone(var)) {
                    ULOG("did not provide \"output interval\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Output_Intervals", (long)this->sparta->run_every);
                } else
                    ULOG("Using provided \"output interval\"");

            } else {
                ULOG("Creating Elmer Time parameters from couple_ratio");
                double r;
                python::loadAttrFromObjectAndConvert(this->elmer, "couple_ratio", r);

                long num_elmer_timesteps = this->sparta->run_every/r;
                double elmer_dt = this->sparta->update->dt * (double)this->sparta->run_every/(double)num_elmer_timesteps;

                python::setAttrOfObject(this->simulation, "Timestep_intervals", num_elmer_timesteps);
                python::setAttrOfObject(this->simulation, "Output_Intervals",   num_elmer_timesteps);
                python::setAttrOfObject(this->simulation, "Timestep_Sizes",     elmer_dt);
            }
            // telling the code to no longer run this section of code
            this->has_been_initialized = true;

            // need to redo this because I set parameters, redoing this will add them to the string
            python::convertObjectToString(this->elmer, this->sif_file_config_str);
            long length = strlen(this->sif_file_config_str);
            if (length < 3)
                UERR("Str method for elmer class returned a string with no length");
        }

        // making the sif file
        ULOG("Making sif file for elmer");
        util::oFile _buf(this->sif);
        // _buf << "! File created: " << util::getTime() << "\n";
        _buf << "! File made by process: " << util::_me << "\n\n";
        _buf << "Check Keywords warn\n";
    
        _buf << this->sif_file_config_str;

        std::size_t i;

        for (i = 0; i < this->bodies.size(); i++)
            this->bodies[i]->join(_buf);

        for (i = 0; i < this->initial_conditions.size(); i++)
            this->initial_conditions[i]->join(_buf);

        for (i = 0; i < this->boundary_conditions.size(); i++)
            this->boundary_conditions[i]->join(_buf);

        _buf.close();
        ULOG("running");

        UERR("DONE");

        // running the command
        this->server.runCommand();
        this->server.waitForDoneFile();
        this->loadNodeData();
        this->updateNodes();
        this->updateNodeFile();
        
    }

    /* ---------------------------------------------------------------------- */

    

    /* ---------------------------------------------------------------------- */

    // void Elmer::dumpBefore() {
    //     this->dumpNodePositionsBefore();
    //     this->dumpNodeTemperaturesBefore(); 
    //     this->dumpNodeVelocitiesBefore();
    // }

    /* ---------------------------------------------------------------------- */

    void Elmer::dump() {
        this->dumpNodeTemperatures();
        this->dumpNodePositions();
        this->dumpNodeVelocities();
    }


    /*----------------------------------------------
    
    
    Private functions

    
    -----------------------------------------------*/


    void Elmer::handleUserBoundaryConditions() {
        
        // long temp_id, size_target_boundaries, target_boundary_id, i;
        long i, j, k, size_target_boundaries, target_boundary_group_id;
        
        PyObject *user_specified_boundary_conditions = python::loadAttrFromObject(this->elmer, "boundary_conditions", python::PYLIST);
        long size = (long)PyList_Size(user_specified_boundary_conditions);
        ULOG("Collecting User Specified Boundary Conditions: " + std::to_string(size));

        bool allow_heat_flux_at_boundary, allow_forces_at_boundary, clamped_at_boundary;
        PyObject *user_bc, *target_boundaries, *list_item;
        
        for (i = 0; i < size; i++) {
            // getting the boundary condition from the list, then getting its id
            user_bc = PyList_GetItem(user_specified_boundary_conditions, (Py_ssize_t)i);
            this->python->setupObject(user_bc);
            // python::loadAttrFromObjectAndConvert(user_bc, "_id",              temp_id);
            python::loadAttrFromObjectAndConvert(user_bc, "_allow_heat_flux", allow_heat_flux_at_boundary);
            python::loadAttrFromObjectAndConvert(user_bc, "_allow_forces",    allow_forces_at_boundary);
            python::loadAttrFromObjectAndConvert(user_bc, "_clamped",         clamped_at_boundary);
            
            // gathering the surface ids in the vector, so I know which ones to skip
            target_boundaries = python::loadAttrFromObject(user_bc, "Target_Boundaries", python::PYLIST);
            
            size_target_boundaries = (long)PyList_Size(target_boundaries);
            
            for (j = 0; j < size_target_boundaries; j++) {
                list_item = PyList_GetItem(target_boundaries, (Py_ssize_t)j);
                python::convertObjToLong(list_item, target_boundary_group_id);
                ULOG("Looking for:" + std::to_string(target_boundary_group_id));
                for (k = 0; k < this->num_boundary; k++) {
                    if (this->boundary[k][0] == target_boundary_group_id) {
                        // UERR("in loop");
                        this->allow_heat_flux[k] = allow_heat_flux_at_boundary;
                        this->allow_forces[k]    = allow_forces_at_boundary;
                        this->clamped[k]         = clamped_at_boundary;
                    }
                }
            }
        }
        // ULOG("DONE looking");
    }




    /* ---------------------------------------------------------------------- */

    // void Elmer::setupDumpDirectory() {
    //     ULOG("setting up dump directory");
    //     int result;
    //     char *from, *to;
    //     // making sure the mesh database is complete
    //     char* exts[4] = {(char*)"boundary", (char*)"nodes", (char*)"header", (char*)"elements"}; // list of component file extensions
    //     for (int i = 0; i < 4; i++) {
    //         result = asprintf(&from, "%s%c%s%s", this->meshDB, SEP, "mesh.", exts[i]);
    //         if (result == -1) { UERR("asprintf failed for \"from\" path"); }

    //         result = asprintf(&to, "%s%c%s%s", this->simulation_directory, SEP, "mesh.", exts[i]);
    //         if (result == -1) { UERR("asprintf failed for \"to\" path"); }
            
    //         util::copyFile(from, to);
    //     }

    //     // setting files
    //     result = asprintf(&this->boundary_file, "%s%c%s", this->simulation_directory, SEP, "mesh.boundary");
    //     if (result == -1) { UERR("asprintf failed for boundary_file path"); }
        
    //     result = asprintf(&this->element_file, "%s%c%s", this->simulation_directory, SEP, "mesh.elements");
    //     if (result == -1) { UERR("asprintf failed for element_file path"); }
        
    //     result = asprintf(&this->node_file, "%s%c%s", this->simulation_directory, SEP, "mesh.nodes");
    //     if (result == -1) { UERR("asprintf failed for node_file path"); }
        
    //     result = asprintf(&this->header_file, "%s%c%s", this->simulation_directory, SEP, "mesh.header");
    //     if (result == -1) { UERR("asprintf failed for header_file path"); }
    //     // this->boundary_file = this->simulation_directory.toString() + SEP + "mesh.boundary";
    //     // this->element_file  = this->simulation_directory.toString() + SEP + "mesh.elements";
    //     // this->node_file     = this->simulation_directory.toString() + SEP + "mesh.nodes";
    //     // this->header_file   = this->simulation_directory.toString() + SEP + "mesh.header";
    // }

    

    /* ---------------------------------------------------------------------- */

    // void Elmer::makeSif() {
        
    // }

    /* ---------------------------------------------------------------------- */

    // void Elmer::checkLoadedNodeData() {
    //     long j;
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
        std::vector<std::string> lines, split_line;
        std::string cur_key, line;
        long cur_index, count, index, i, start = 0;

        ULOG("Loading data from: " + std::string(this->node_data_file));
        util::readFile(this->node_data_file, lines);
        std::vector<std::array<long, 2>> permutation_table;

        std::vector<bool> check_list;
        check_list.resize(this->num_nodes, true);

        // gets latest data and permutation table
        for (i = 0; i < (long)lines.size(); i++) {
            util::trim(lines[i]);
            util::split(lines[i], split_line, ' ');
            if (split_line[0] == (std::string)"Time:") start = i;
            if (split_line[0] == (std::string)"Perm:") {
                if (split_line[1] == "use") {
                    ULOG("Got empty permutation table at line: " + std::to_string(i+1));
                    continue;
                } else {
                    util::split(lines[i-2], split_line, ' ');
                    if (split_line[0] != (std::string)"Time:")
                        UERR("got another permutation table from data file at line: " + std::to_string(i+1) + ", this is not currently support");
                    
                    ULOG("Loading permutation table at line (will overwrite existing): " + std::to_string(i+1));
                    permutation_table.clear();

                    for (i = i + 1; i < (long)lines.size(); i++) {
                        util::trim(lines[i]);
                        util::split(lines[i], split_line, ' ');
                        if (split_line.size() < (std::size_t)2) {
                            i--;
                            break;
                        }
                        permutation_table.push_back({std::stoi(split_line[0]), std::stoi(split_line[1])});
                    }
                }
            }
        }

        // loading information based on the permutation table
        cur_index = 0;
        check_list.resize(this->num_nodes, true);
        for (i = start+1; i < (long)lines.size(); i++) {
            line = lines[i];
            util::trim(line);
            util::split(line, split_line, ' ');
            if ((std::isdigit(lines[i][0]) || lines[i][0] == '-') && split_line.size() == 1) {
                for (std::size_t j = 0; j < permutation_table.size(); j++) {
                    if (permutation_table[j][0] == cur_index+1)
                        index = permutation_table[j][1] - 1;
                }

                if        (cur_key == (std::string)"temperature") {
                    this->node_temperatures[index] = std::stod(line);
                } else if (cur_key == (std::string)"displacement 1") {
                    this->node_displacements[index][0] += std::stod(line);
                } else if (cur_key == (std::string)"displacement 2") {
                    this->node_displacements[index][1] += std::stod(line);
                } else if (cur_key == (std::string)"displacement 3") {
                    this->node_displacements[index][2] += std::stod(line);
                } else if (cur_key == (std::string)"velocity 1") {
                    this->node_velocities[index][0] = std::stod(line);
                } else if (cur_key == (std::string)"velocity 2") {
                    this->node_velocities[index][1] = std::stod(line);
                } else if (cur_key == (std::string)"velocity 3") {
                    this->node_velocities[index][2] = std::stod(line);
                } else
                    UERR("got unknown key in data file with name: " + cur_key);
                check_list[index] = true;
                cur_index++;
            } else if ((std::string)"Perm:" == split_line[0]) {
                if (cur_index > 0) {
                    if (cur_index != this->num_nodes)
                        UERR("did not get all indicies while loading \"" + cur_key + "\" from data");
                    
                    count = 1;
                    for (auto it : check_list) {
                        if (!(it))
                            UERR("Did not get value at index \"" + std::to_string(count) + "\" while loading \"" + cur_key + "\"");
                        count++;
                    }
                    ULOG("Parameter \"" + cur_key + "\" completed checks");
                }

                check_list.resize(this->num_nodes, false);
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
        long i, j;
        double before, after;
        bool checks[this->num_nodes], got_no_difference;

        for (i = 0; i < (signed)this->num_nodes; i++) {
            checks[i] = false;

            // this block of code checks if a displacement is non-zero but
            // falls below below double precision when added to another value
            // so it does not affect the other value
            got_no_difference = false;
            for (j = 0; j < ELMER_DIMENSION; j++) {
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
                for (j = 0; j < ELMER_DIMENSION; j++) {
                    this->nodes[i][j] += this->node_displacements[i][j];
                    if (nodes[i][j] < this->sparta->domain->boxlo[j] || nodes[i][j] > this->sparta->domain->boxhi[j]) {
                        this->dump();
                        UERR("Detected node outside of bounds at index: " + std::to_string(i) + ", dumped data");
                    }
                    this->node_displacements[i][j] = 0.0;
                }
            }
        }

        for (i = 0; i < this->num_nodes; i++) {
            if (!(checks[i]))
                UERR("Did not load node at index: " + std::to_string(i+1));
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::updateNodeFile() {
        ULOG("updating node file");
        util::oFile out(this->node_file);

        for (long i = 0; i < this->num_nodes; i++)
            out << i+1 << " " << -1 << " " << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    }

    /* ---------------------------------------------------------------------- */
    
    void Elmer::setupVectors() {
        // getting the number of nodes
        util::iFile h = util::iFile(this->header_file);
        std::vector<std::string> lines, split;
        h.getLines(lines);
        h.close();

        // splitting lines[0] at spaces, saving result to split
        util::split(lines[0], split, ' ');
        
        this->num_nodes    = std::stol(split[0]);
        this->num_elements = std::stol(split[1]);
        this->num_boundary = std::stol(split[2]);

        ULOG("Creating node vector with length: " + std::to_string(this->num_nodes));
        this->nodes    = new double*[this->num_nodes];

        ULOG("Creating element vector with length: " + std::to_string(this->num_elements));
        this->elements = new long*[this->num_elements];
        this->element_lengths = new long[this->num_elements];

        ULOG("Creating boundary vector with length: " + std::to_string(this->num_boundary));
        this->boundary = new long*[this->num_boundary];

        // writing base temperature to file for each node
        ULOG("Creating nodal data vectors with length: " + std::to_string(this->num_nodes));
        this->node_temperatures  = new double[this->num_nodes];
        this->node_velocities    = new double*[this->num_nodes];
        this->node_displacements = new double*[this->num_nodes];

        long j;
        for (long i = 0; i < num_nodes; i++) {
            this->node_temperatures[i]  = this->base_temp;
            this->node_velocities[i]    = new double[ELMER_DIMENSION];
            this->node_displacements[i] = new double[ELMER_DIMENSION];
            // ensuring all are set to 0
            for (j = 0; j < ELMER_DIMENSION; j++) {
                this->node_velocities[i][j]    = 0.0;
                this->node_displacements[i][j] = 0.0;
            }
        }

        // user inputs arrays
        this->allow_forces    = new bool[this->num_boundary];
        this->allow_heat_flux = new bool[this->num_boundary];
        this->clamped         = new bool[this->num_boundary];

        for (long i = 0; i < this->num_boundary; i++) {
            this->allow_forces[i]    = true;
            this->allow_heat_flux[i] = true;
            this->clamped[i]         = false;
        }

        // loading the stuff from the elmer mesh database
        this->loadBoundaries();
        this->loadElements();
        this->loadNodes();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadNodes() {
        ULOG("Loading nodes");
        std::vector<std::string> lines, split_line;        
        util::readFile(this->node_file, lines);

        if ((long)lines.size() != this->num_nodes)
            UERR("Number of nodes in file is not equal to number of nodes in header");

        long j, count = 0;
        for (long i = 0; i < this->num_nodes; i++) {
            count++;
            util::trim(lines[i]);
            util::split(lines[i], split_line, ' ');
            if (split_line.size() != NODE_LINE_SIZE)
                UERR("node element not correct size at line " + std::to_string(count) + ", should have size " + std::to_string(NODE_LINE_SIZE));

            if (std::stol(split_line[0]) != count)
                UERR("detected unordered node element at line " + std::to_string(count));

            if (std::stoi(split_line[1]) != -1)
                UERR("detected a value that is not -1 in the second column of the node file, this is not yet support");

            this->nodes[i] = new double[ELMER_DIMENSION];

            for (j = NODE_START_READ; j < NODE_LINE_SIZE; j++)
                this->nodes[i][j-NODE_START_READ] = std::stod(split_line[j]);
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadBoundaries() {
        ULOG("Loading boundaries");
        std::vector<std::string> lines, _split;
        util::readFile(this->boundary_file, lines);
        
        if ((long)lines.size() != this->num_boundary)
            UERR("Number of boundary elements in file is not equal to number of boundary elements in header");

        long j, count = 0;
        for (long i = 0; i < this->num_boundary; i++) {
            count++;
            util::trim(lines[i]);
            util::split(lines[i], _split, ' ');

            if (_split[4] != (std::string)"303") 
                UERR("element is not a triangle in boundary file at line: " + std::to_string(count));

            if (std::stol(_split[0]) != count)
                UERR("got unordered boundary elements at line: " + std::to_string(count));

            if (_split.size() != BOUNDARY_ELEMENT_LINE_SIZE)
                UERR("Caught boundary element with incorrect number of entries");

            this->boundary[i] = new long[BOUNDARY_ELEMENT_ARRAY_SIZE];

            // adding the data
            for (j = BOUNDARY_ELEMENT_START_READ; j < BOUNDARY_ELEMENT_LINE_SIZE; j++)
                this->boundary[i][j-BOUNDARY_ELEMENT_START_READ] = std::stol(_split[j]);
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadElements() {
        ULOG("Loading elements");
        std::vector<std::string> lines, _split;
        util::readFile(this->element_file, lines);

        if ((long)lines.size() != this->num_elements) {
            UERR("Number of elements in file is not equal to number of elements in header");
        }

        long j, count = 0;
        for (long i = 0; i < this->num_elements; i++) {
            count++;
            util::trim(lines[i]);
            util::split(lines[i], _split, ' ');

            if (std::stol(_split[0]) != count)
                UERR("got unordered boundary elements at line: " + std::to_string(count));

            this->element_lengths[i] = (long)(_split.size() - ELEMENT_START_READ);
            this->elements[i]        = new long[this->element_lengths[i]];

            for (j = ELEMENT_START_READ; j < this->element_lengths[i]+ELEMENT_START_READ; j++)
                this->elements[i][j-ELEMENT_START_READ] = std::stol(_split[j]);
        }
    }

    /* ---------------------------------------------------------------------- */

    // void Elmer::dumpNodeTemperaturesBefore() {
    //     ULOG("dumping node temperatures before");
    //     char* filename;
    //     int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".before_", this->node_temperature_file_ext);
    //     if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }
    //     //std::string filename = std::string(this->simulation_directory) + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + node_temperature_file_ext;
    //     util::oFile out(filename);

    //     for (double& it : this->node_temperature_data)
    //         out << it << "\n";
    // }

    // /* ---------------------------------------------------------------------- */

    // void Elmer::dumpNodeVelocitiesBefore() {
    //     ULOG("dumping node velocities before");
    //     char* filename;
    //     int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".before_", this->node_velocity_file_ext);
    //     if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }
    //     //std::string filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_velocity_file_ext.toString();
    //     util::oFile out(filename);

    //     for (auto& it : this->node_velocity_data)
    //         out << it[0] << " " << it[1] << " " << it[2] << "\n";
    // }

    // /* ---------------------------------------------------------------------- */

    // void Elmer::dumpNodePositionsBefore() {
    //     ULOG("dumping node positions before");
    //     char* filename;
    //     int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".before_", this->node_position_file_ext);
    //     if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }
    //     // std::string filename = this->simulation_directory.toString() + SEP + std::to_string(this->sparta->update->ntimestep) + ".before_" + this->node_position_file_ext.toString();
    //     util::oFile out(filename);

    //     for (std::size_t i = 0; i < this->nodes.size(); i++) 
    //         out << this->nodes[i][0] << " " << this->nodes[i][1] << " " << this->nodes[i][2] << "\n";
    // }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeTemperatures() {
        ULOG("dumping node temperatures");
        char* filename;
        int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".", this->node_temperature_file_ext);
        if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }

        util::oFile out(filename);

        for (long i = 0; i < this->num_nodes; i++)
            out << this->node_temperatures[i] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodePositions() {
        ULOG("dumping node positions");
        char* filename;
        int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".", this->node_position_file_ext);
        if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }

        util::oFile out(filename);

        int j;
        for (long i = 0; i < this->num_nodes; i++) {
            out << this->nodes[i][0];
            for (j = 1; j < ELMER_DIMENSION; j++)
                out << " " << this->nodes[i][j];
            out << "\n";
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeVelocities() {
        ULOG("dumping node velocities");
        char* filename;
        int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, this->sparta->update->ntimestep, ".", this->node_velocity_file_ext);
        if (result == -1) { UERR("asprintf failed for dumpNodeVelocities"); }

        util::oFile out(filename);
        
        int j;
        for (long i = 0; i < this->num_nodes; i++) {
            out << this->node_velocities[i][0];
            for (j = 1; j < ELMER_DIMENSION; j++)
                out << " " << this->node_velocities[i][j];
            out << "\n";
        }
    }
}

#endif