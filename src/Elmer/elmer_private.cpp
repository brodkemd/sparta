#include "elmer.h"
#include "mpi.h"
#include <cstdlib>
#include <unordered_map>
#include "../fix_fea.h"

namespace elmer {
    void Elmer::setupVectors() {
        // getting the number of nodes
        std::vector<std::string> lines, split;
        util::iFile h = util::iFile(this->header_file);
        h.getLines(lines);
        h.close();

        // splitting lines[0] at spaces, saving result to split
        util::trim(lines[0]);
        util::splitStringAtWhiteSpace(lines[0], split);
        
        // getting the parameters from the mesh header file
        this->num_nodes    = std::stol(split[0]);
        this->num_elements = std::stol(split[1]);
        this->num_boundary = std::stol(split[2]);

        // ULOG("Creating node vector with length: " + std::to_string(this->num_nodes));
        this->nodes    = new double*[this->num_nodes];

        // ULOG("Creating element vector with length: " + std::to_string(this->num_elements));
        this->elements = new long*[this->num_elements];
        this->element_lengths = new long[this->num_elements];

        // ULOG("Creating boundary vector with length: " + std::to_string(this->num_boundary));
        this->boundary = new long*[this->num_boundary];

        // writing base temperature to file for each node
        // ULOG("Creating nodal data vectors with length: " + std::to_string(this->num_nodes));
        this->node_temperatures  = new double[this->num_nodes];
        this->node_velocities    = new double*[this->num_nodes];
        this->node_displacements = new double*[this->num_nodes];

        // sets each element in vectors
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

        // setting the default values
        for (long i = 0; i < this->num_boundary; i++) {
            this->allow_forces[i]    = false;
            this->allow_heat_flux[i] = true;
            this->clamped[i]         = false;
        }
        // overwrites is deformation is allowed
        if (this->shouldUpdateSurf()) {
            for (long i = 0; i < this->num_boundary; i++) {
                this->allow_forces[i]    = true;
            }
        }   
        
        // loading the stuff from the elmer mesh database
        this->loadBoundaries();
        this->loadElements();
        this->loadNodes();
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::handleUserBoundaryConditions() {
        
        // long temp_id, size_target_boundaries, target_boundary_id, i;
        long i, j, k, size_target_boundaries, target_boundary_group_id;
        
        // getting the boundary conditions list from the elmer class
        PyObject *user_specified_boundary_conditions = python::loadAttrFromObject(this->elmer, "boundary_conditions", python::PYLIST);

        // getting how long the list of boundary conditions from the elmer class is
        long size = (long)PyList_Size(user_specified_boundary_conditions);
        ULOG("Collecting User Specified Boundary Conditions: " + std::to_string(size));

        // declaring the variables
        bool allow_heat_flux_at_boundary, allow_forces_at_boundary, clamped_at_boundary;
        PyObject *user_bc, *target_boundaries, *list_item;
        
        // iterating over list of boundary conditions
        for (i = 0; i < size; i++) {
            // getting the boundary condition from the list, then getting its id
            user_bc = PyList_GetItem(user_specified_boundary_conditions, (Py_ssize_t)i);
            this->python->setupObject(user_bc);

            // getting the user specified parameters from the boundary condition
            python::loadAttrFromObjectAndConvert(user_bc, "_allow_heat_flux", allow_heat_flux_at_boundary);
            python::loadAttrFromObjectAndConvert(user_bc, "_allow_forces",    allow_forces_at_boundary);
            python::loadAttrFromObjectAndConvert(user_bc, "_clamped",         clamped_at_boundary);
            
            // gathering the surface ids in the vector, so I know which ones to skip
            target_boundaries = python::loadAttrFromObject(user_bc, "Target_Boundaries", python::PYLIST);
            
            // getting length of the list of boundary boundary conditions
            size_target_boundaries = (long)PyList_Size(target_boundaries);
            
            // iterating over list of target boundaries
            for (j = 0; j < size_target_boundaries; j++) {
                // getting the list item at index j
                list_item = PyList_GetItem(target_boundaries, (Py_ssize_t)j);

                // getting the boundary_group_id as a long
                python::convertObjToLong(list_item, target_boundary_group_id);

                // setting all indices in the user specified conditions vectors
                // for boundary elements with the same target_group_id
                for (k = 0; k < this->num_boundary; k++) {
                    // if the boundary element has the same id as the target
                    if (this->boundary[k][0] == target_boundary_group_id) {
                        // setting to the provided values from the user
                        this->allow_heat_flux[k] = allow_heat_flux_at_boundary;
                        this->allow_forces[k]    = allow_forces_at_boundary;
                        this->clamped[k]         = clamped_at_boundary;
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::createInitialConditions() {
        ULOG("Creating Initial Conditions");

        long i, j, k, count;
        std::vector<std::string> hashes(this->num_elements);
        double avg[ELMER_DIMENSION+1];

        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        // ULOG("Setting surface temperatures");
        for (long i = 0; i < this->num_elements; i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            for (j = 0; j < ELMER_DIMENSION+1; j++)
                avg[j] = 0;
            
            for (j = ELEMENT_ARRAY_NODE_START; j < this->element_lengths[i]; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg[0] += this->node_temperatures[this->elements[i][j]-1];
                for (k = 0; k < ELMER_DIMENSION; k++)
                    avg[k+1] += this->node_velocities[this->elements[i][j]-1][k];
            }

            // dividing by the number of data points (the average)
            for (j = 0; j < ELMER_DIMENSION+1; j++)
                avg[j] /= (this->element_lengths[i] - ELEMENT_ARRAY_NODE_START);
            
            hashes[i] = util::hashDoubleArray(avg, ELMER_DIMENSION+1);
        }

        std::unordered_map<std::string, std::vector<long>> data_dict;
        // ULOG("Setting up Inital Condition Map");

        // setting all the hashes in the map, this essentially filters
        // for unique values (of surface parameters)
        for (i = 0; i < this->num_elements; i++)
            data_dict[hashes[i]] = std::vector<long>({});

        // this step adds the indicies that correspond to each unique 
        // value (of surface parameters)
        for (i = 0; i < this->num_elements; i++)
            data_dict[hashes[i]].push_back(i);

        this->bodies.clear();
        this->initial_conditions.clear();

        // updating the element file
        util::oFile out(this->element_file);

        // Elmer starts counting at 1 (not 0), so set these to 1
        i = 1; count = 1;

        // iterate over the map of unique surface parameters
        for (const auto& [key, value] : data_dict) {
            // writing to the element file, count is the index of the surface element
            // in the file and i is the element group. The rest is added from the 
            // originally loaded element file
            for (j = 0; j < (long)value.size(); j++) {
                out << count << " " << i;
                for (k = 1; k < this->element_lengths[value[j]]; k++)
                    out << " " << this->elements[value[j]][k];
                out << "\n";
                count++;
            }

            // recompute values
            for (j = 0; j < ELMER_DIMENSION+1; j++)
                avg[j] = 0;
            
            for (j = ELEMENT_ARRAY_NODE_START; j < this->element_lengths[value[0]]; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg[0] += this->node_temperatures[this->elements[value[0]][j]-1];
                for (k = 0; k < ELMER_DIMENSION; k++)
                    avg[k+1] += this->node_velocities[this->elements[value[0]][j]-1][k];
            }

            for (j = 0; j < ELMER_DIMENSION+1; j++)
                avg[j] /= (this->element_lengths[value[0]] - ELEMENT_ARRAY_NODE_START);

            Section* ic = new Section("Initial Condition", i);
            ic->addEquality("Temperature", avg[0]);
            ic->addEquality("Velocity 1",  avg[1]);
            ic->addEquality("Velocity 2",  avg[2]);
            ic->addEquality("Velocity 3",  avg[3]);

            Section* body = new Section("Body", i);
            
            body->addEquality("Initial Condition", i, true);
            
            if (body_force_id != util::NO_INT)
                body->addEquality("Body Force", this->body_force_id);

            body->addEquality("Equation", this->equation_id);
            body->addEquality("Material", this->material_id);
            body->addEquality("Target Bodies", i, true);

            this->initial_conditions.push_back(ic);
            this->bodies.push_back(body);
            
            i++;
        }

        ULOG("# of initial conditions and bodies: " + std::to_string(this->initial_conditions.size()));
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::createBoundaryConditions() {
        ULOG("Creating Boundary Conditions");

        long len, i, j, k, count = 1;
        double* temp_arr = new double[7];
        std::vector<std::string> hashes(this->num_boundary); // num_boundary = sparta->nsurf

        // ULOG("Calculating Hashes");
        for (i = 0; i < this->num_boundary; i++) {
            len = 0;
            if (allow_heat_flux[i]) {
                temp_arr[len] = this->sparta->qw[i];
                len += 1;
            }

            if (allow_forces[i]) {
                // if (this->shouldUpdateSurf()) {
                // do not need this if here because this is guarded by python
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

            hashes[i] = util::hashDoubleArray(temp_arr, len);
        }

        // detecting unique values and saving them to the stress vector and
        // their corresponding indices to the indices vector, this condenses the
        // number of boundary conditions (removes boundary conditions with
        // the same values)
        std::unordered_map<std::string, std::vector<long>> data_dict;
        // ULOG("Setting up Boundary Map");

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
        // ULOG("Generating Boundary Conditions");
        util::oFile out(this->boundary_file);
        
        // Elmer starts counting at 1 (not 0), so set these to 1
        i = 1; count = 1;
        // ULOG("Boundary Conditions start with id:" + std::to_string(i));

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
            
            bc->addEquality("Target Boundaries", i, true);
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
                bc->addEquality("Displacement 1", (double)0);
                bc->addEquality("Displacement 2", (double)0);
                bc->addEquality("Displacement 3", (double)0);
            }
            
            this->boundary_conditions.push_back(bc);
            i++;
        }
        delete [] temp_arr;

        out.close();

        ULOG("# of boundary conditions: " + std::to_string(this->boundary_conditions.size()));
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadNodeData() {
        // ULOG("loading node data");
        std::string line, var;
        std::unordered_map<long, long> permutation_table;

        // ULOG("Initializing node data permutation table");
        for (long i = 0; i < num_nodes; i++)
            permutation_table[i] = 0;

        ULOG("Loading data from: " + std::string(node_data_file));
        std::ifstream file(node_data_file);

        if (file.is_open()) {
            long line_count = 1, counter = 0;
            bool loading_perm = false;

            std::vector<std::string> split_line;
            while (std::getline(file, line)) {
                // removes leading and trailing whitespace
                util::trim(line);

                // split line at whitespace
                util::splitStringAtWhiteSpace(line, split_line);

                // if the first character in a line is not numeric
                if (!isdigit(line[0]) && line[0] != '-') {
                    // switching over whatever the first word is
                    // handles if a new time step was found
                    if (split_line[0] == "Time:") {
                        // ULOG("Caught new data, replacing old with new");
                        ;
                    }
                    // handles if a new permutation table was found
                    else if (split_line[0] == "Perm:") {
                        // handles if it says "use previous"
                        if (split_line[1] == "use") {
                            // do not need to do anything
                            // ULOG("  Using previous permutation");
                            ;
                        }
                        // makes sure there is a value number of entries in the table
                        else if (std::stol(split_line[1]) == num_nodes) {
                            // telling the rest of the code that a permutation table is being loaded
                            // ULOG("  Loading new permutation table at line:" + std::to_string(line_count));
                            loading_perm = true;
                        }
                        // caught a bad table
                        else
                            UERR("invalid permutation table at line: " + std::to_string(line_count));
                    }
                    // if here then it is a variable name
                    else {
                        var = util::joinBySpaces(split_line);
                        // ULOG("  Loading from var: " + var);
                        counter = 0;
                    }
                }
                // if the first character in a line is numeric, then it is data, so loading it
                else {
                    // if the data is part of a permutation table
                    if (loading_perm) {
                        // getting the map, indexing from 0 that why there is -1 here
                        permutation_table[std::stol(split_line[0])-1] = std::stol(split_line[1])-1;
                        counter++;

                        // if the end of the table was reached, stop loading the table
                        if (counter >= num_nodes) {
                            loading_perm = false;
                            counter = 0;
                        }
                    }
                    // if the data is just data
                    else {
                        if (split_line.size() != 1)
                            UERR("Invalid line detected at line number:" + std::to_string(line_count));

                        // adding the data at the index in the vector that is currently being loaded
                        if        (var == (std::string)"temperature") {
                            this->node_temperatures[counter]      = std::stod(split_line[0]);
                        } 
                        else if (var == (std::string)"displacement 1") {
                            this->node_displacements[counter][0] += std::stod(split_line[0]); // rolling addition so that values that fall below precision are remembered
                        }
                        else if (var == (std::string)"displacement 2") {
                            this->node_displacements[counter][1] += std::stod(split_line[0]); // rolling addition so that values that fall below precision are remembered
                        }
                        else if (var == (std::string)"displacement 3") {
                            this->node_displacements[counter][2] += std::stod(split_line[0]); // rolling addition so that values that fall below precision are remembered
                        }
                        else if (var == (std::string)"velocity 1") {
                            this->node_velocities[counter][0]     = std::stod(split_line[0]);
                        }
                        else if (var == (std::string)"velocity 2") {
                            this->node_velocities[counter][1]     = std::stod(split_line[0]);
                        }
                        else if (var == (std::string)"velocity 3") {
                            this->node_velocities[counter][2]     = std::stod(split_line[0]);
                        }
                        else { 
                            UERR("got unknown key in data file with name: " + var);
                        }
                        counter++;
                    }
                }
                line_count++;
            }
        } else {
            UERR("Could not open the node data file");
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::updateNodes() {
        ULOG("Updating nodes");
        long i, j;
        double before, after;
        bool got_difference;

        for (i = 0; i < this->num_nodes; i++) {
            // this block of code checks if a displacement is non-zero but
            // falls below below double precision when added to another value
            // so it does not affect the other value
            got_difference = true;
            for (j = 0; j < ELMER_DIMENSION; j++) {
                before = this->nodes[i][j];
                after  = before + this->node_displacements[i][j];
                if (before == after && this->node_displacements[i][j] != 0.0) {
                    got_difference = false;
                    break;
                }
            }

            // if the displacement changes the value
            if (got_difference) {
                for (j = 0; j < ELMER_DIMENSION; j++) {
                    this->nodes[i][j] += this->node_displacements[i][j];
                    if (nodes[i][j] < *(this->sparta->boxlo + j) || nodes[i][j] > *(this->sparta->boxhi + j)) {
                        this->dump();
                        UERR("Detected node outside of bounds at index: " + std::to_string(i) + ", dumped data");
                    }
                    this->node_displacements[i][j] = 0.0;
                }
            } //else ULOG("did not update node point: " + std::to_string(i+1) + ", remembering displacement");
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
    
    void Elmer::handleNodeFileSplitLine(std::vector<std::string>& split, long line_number) {
        if (split.size() != NODE_LINE_SIZE)
            UERR("node element not correct size at line " + std::to_string(line_number) + ", should have size " + std::to_string(NODE_LINE_SIZE));

        if (std::stol(split[0]) != line_number)
            UERR("detected unordered node element at line " + std::to_string(line_number));

        if (std::stoi(split[1]) != -1)
            UERR("detected a value that is not -1 in the second column of the node file, this is not yet support");

        if (line_number > this->num_nodes)
            UERR("Too many entries in the node file");

        this->nodes[line_number-1] = new double[ELMER_DIMENSION];

        for (long j = NODE_START_READ; j < NODE_LINE_SIZE; j++)
            this->nodes[line_number-1][j-NODE_START_READ] = std::stod(split[j]);
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadNodes() {
        ULOG("Loading nodes");
        // iterates over the lines in the node file, trims whitespace, splits at whitespace
        // then pass resulting vector to handleNodeFileSplitLine
        ITERATE_OVER_NONEMPTY_LINES_FROM_FILE_AND_SPLIT(
            this->node_file,
            this->handleNodeFileSplitLine
        );
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::handleBoundaryFileSplitLine(std::vector<std::string>& split, long line_number) {
        if (split[4] != (std::string)"303") 
            UERR("element is not a triangle in boundary file at line: " + std::to_string(line_number));

        if (std::stol(split[0]) != line_number)
            UERR("got unordered boundary elements at line: " + std::to_string(line_number));

        if (split.size() != BOUNDARY_ELEMENT_LINE_SIZE)
            UERR("Caught boundary element with incorrect number of entries");

        if (line_number > this->num_boundary)
            UERR("Too many entries in the boundary file");

        this->boundary[line_number-1] = new long[BOUNDARY_ELEMENT_ARRAY_SIZE]; // line number starts at 1

        // adding the data
        for (int j = BOUNDARY_ELEMENT_START_READ; j < BOUNDARY_ELEMENT_LINE_SIZE; j++)
            this->boundary[line_number-1][j-BOUNDARY_ELEMENT_START_READ] = std::stol(split[j]);
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadBoundaries() {
        ULOG("Loading boundaries");
        // iterates over the lines in the boundary file, trims whitespace, splits at whitespace
        // then pass resulting vector to handleBoundaryFileSplitLine
        ITERATE_OVER_NONEMPTY_LINES_FROM_FILE_AND_SPLIT(
            this->boundary_file,
            this->handleBoundaryFileSplitLine
        );
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::handleElementFileSplitLine(std::vector<std::string>& split, long line_number) {
        if (std::stol(split[0]) != line_number)
            UERR("got unordered boundary elements at line: " + std::to_string(line_number));
        
        if (line_number > this->num_elements)
            UERR("Too many entries in the elements file");

        this->element_lengths[line_number-1] = (long)(split.size() - ELEMENT_START_READ);
        this->elements[line_number-1]        = new long[this->element_lengths[line_number-1]];

        for (long j = ELEMENT_START_READ; j < this->element_lengths[line_number-1]+ELEMENT_START_READ; j++)
            this->elements[line_number-1][j-ELEMENT_START_READ] = std::stol(split[j]);
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::loadElements() {
        ULOG("Loading elements");
        // iterates over the lines in the element file, trims whitespace, splits at whitespace
        // then pass resulting vector to handleElementFileSplitLine
        ITERATE_OVER_NONEMPTY_LINES_FROM_FILE_AND_SPLIT(
            this->element_file,
            this->handleElementFileSplitLine
        );
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeTemperatures() {
        ULOG("dumping node temperatures");
        char* filename;
        int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, *this->sparta->timestep, ".", this->node_temperature_file_ext);
        if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }

        util::oFile out(filename);

        for (long i = 0; i < this->num_nodes; i++)
            out << this->node_temperatures[i] << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::dumpNodeVelocities() {
        ULOG("dumping node velocities");
        char* filename;
        int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, *this->sparta->timestep, ".", this->node_velocity_file_ext);
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

    /* ---------------------------------------------------------------------- */

    // void Elmer::dumpNodePositions() {
    //     ULOG("dumping node positions");
    //     char* filename;
    //     int result = asprintf(&filename, "%s%c%ld%s%s", this->simulation_directory, SEP, *this->sparta->timestep, ".", this->node_position_file_ext);
    //     if (result == -1) { UERR("asprintf failed for dumpNodePositions"); }

    //     util::oFile out(filename);

    //     int j;
    //     for (long i = 0; i < this->num_nodes; i++) {
    //         out << this->nodes[i][0];
    //         for (j = 1; j < ELMER_DIMENSION; j++)
    //             out << " " << this->nodes[i][j];
    //         out << "\n";
    //     }
    // }

    /* ---------------------------------------------------------------------- */
    /*
    
     all of this is handled in python now

    */
    
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
}

