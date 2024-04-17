#include "elmer.h"
#include <unordered_map>
#include "server.h"
#include "../fix_fea.h"

namespace elmer{
    Elmer::Elmer(SPARTA_NS::FixFea* _sparta, long _me, python::handler* _h) {
        // ULOG("setting up elmer class");
        // referencing the calling sparta instance
        this->sparta     = &*_sparta;
        this->me         = _me;
        this->python     = &*_h;
        this->server     = new elmer::Server(_h);
        
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
        python::loadAttrFromObjectAndConvert(this->elmer, "_elastic_solver_is_used",    this->elastic_solver_is_used);

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

        // creating a string from the elmer config
        python::convertObjectToString(this->elmer, this->sif_file_config_str);        
        long length = strlen(this->sif_file_config_str);
        if (length < 3)
            UERR("Str method for elmer class returned a string with no length");

        // creates all of the vectors
        this->setupVectors();
        this->handleUserBoundaryConditions();

        // setting up server
        this->server->makeServerFile();
        this->server->start();
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
        delete this->server;
    }

    /* ---------------------------------------------------------------------- */

    bool Elmer::shouldUpdateSurf() {
        return this->elastic_solver_is_used;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * makes sure file from elmer object data
     * returns: the name of the surf file
    */
    void Elmer::makeSpartaSurf(char* path) {
        ULOG("making surface file for sparta");
        
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

    void Elmer::makeSif() {
        if (!has_been_initialized) {
            PyObject* var;
            var = python::loadAttrFromObject(this->elmer, "couple_ratio", python::PYFLOAT, true);
            if (python::isNone(var)) {
                var = python::loadAttrFromObject(this->simulation, "Timestep_Sizes", python::PYFLOAT, true);
                if (python::isNone(var)) {
                    ULOG("did not provide \"timestep size\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Timestep_Sizes", *this->sparta->dt);
                } else
                    ULOG("Using provided \"timestep size\"");

                // setting number of timesteps to run
                var = python::loadAttrFromObject(this->simulation, "Timestep_intervals", python::PYINT, true);
                if (python::isNone(var)) {
                    ULOG("did not provide \"number of timesteps\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Timestep_intervals", (long)this->sparta->run_every);
                } else
                    ULOG("Using provided \"number of timesteps\"");
                

                var = python::loadAttrFromObject(this->simulation, "Output_Intervals", python::PYINT, true);
                if (python::isNone(var)) {
                    ULOG("did not provide \"output interval\", using Sparta's");
                    python::setAttrOfObject(this->simulation, "Output_Intervals", (long)this->sparta->run_every);
                } else
                    ULOG("Using provided \"output interval\"");

            } else {
                ULOG("Creating Elmer Time parameters from couple_ratio");
                double r;
                python::loadAttrFromObjectAndConvert(this->elmer, "couple_ratio", r);

                long num_elmer_timesteps = this->sparta->run_every/r;
                double elmer_dt = *this->sparta->dt * (double)this->sparta->run_every/(double)num_elmer_timesteps;

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

        // making conditions
        this->createInitialConditions();
        this->createBoundaryConditions();

        // making the sif file
        // ULOG("Making sif file for elmer");

        util::oFile _buf(this->sif);
        _buf << "! File created: " << util::getTime() << "\n";
        _buf << "! File made by process: " << util::_me << "\n\n";
        _buf << "Check Keywords warn\n";
    
        _buf << this->sif_file_config_str;

        std::size_t i;

        for (i = 0; i < this->initial_conditions.size(); i++)
            this->initial_conditions[i]->joinInto(_buf);

        for (i = 0; i < this->bodies.size(); i++)
            this->bodies[i]->joinInto(_buf);

        for (i = 0; i < this->boundary_conditions.size(); i++)
            this->boundary_conditions[i]->joinInto(_buf);

        _buf.close();

        char* filename; int result;
        result = asprintf(&filename, "%s%c%ld%s", this->simulation_directory, SEP, *this->sparta->timestep, ".sif");
        if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }

        // making a version of the sif file at the timestep
        util::copyFile(this->sif, std::string(filename));

        // copying element file to a timestep instance
        result = asprintf(&filename, "%s%c%ld%s", this->simulation_directory, SEP, *this->sparta->timestep, ".elements");
        if (result == -1) { UERR("asprintf failed for sif file"); }
        util::copyFile(this->element_file, filename);

        // copying node file to a timestep instance
        result = asprintf(&filename, "%s%c%ld%s", this->simulation_directory, SEP, *this->sparta->timestep, ".nodes");
        if (result == -1) { UERR("asprintf failed for sif file"); }
        util::copyFile(this->node_file, filename);

        result = asprintf(&filename, "%s%c%ld%s", this->simulation_directory, SEP, *this->sparta->timestep, ".boundary");
        if (result == -1) { UERR("asprintf failed for sif file"); }
        util::copyFile(this->boundary_file, filename);

    }

    /* ---------------------------------------------------------------------- */

    void Elmer::run() {
        this->makeSif();
        // running the command
        this->server->runCommand();
        this->server->waitForDoneFile();
        this->loadNodeData();
        if (this->shouldUpdateSurf()) {
            this->updateNodes();
            this->updateNodeFile();
        }
    }

    /* ---------------------------------------------------------------------- */

    // void Elmer::dumpBefore() {
    //     this->dumpNodePositionsBefore();
    //     this->dumpNodeTemperaturesBefore(); 
    //     this->dumpNodeVelocitiesBefore();
    // }

    /* ---------------------------------------------------------------------- */

    void Elmer::dump() {
        this->dumpNodeTemperatures();
        // this->dumpNodePositions();
        this->dumpNodeVelocities();
    }
}