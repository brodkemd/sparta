#include "elmer.h"
#include <unordered_map>
#include "server.h"
#include "../fix_fea.h"
#include <string.h>

namespace elmer{
    Elmer::Elmer(SPARTA_NS::FixFea* _sparta, SPARTA_NS::index_t _me, python::handler* _h) {
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
        this->simulation = this->python->loadObjectWithSetupFromObject(this->elmer, "simulation");
        this->equation   = this->python->loadObjectWithSetupFromObject(this->elmer, "equation");
        this->material   = this->python->loadObjectWithSetupFromObject(this->elmer, "material");
        this->body_force = this->python->loadObjectWithSetupFromObject(this->elmer, "body_force");

        // getting ids from attributes, used to set conditions
        python::loadAttrFromObjectAndConvert(this->equation, "_id", this->equation_id);
        python::loadAttrFromObjectAndConvert(this->material, "_id", this->material_id);

        if (this->gravity_on)
            python::loadAttrFromObjectAndConvert(this->body_force, "_id", this->body_force_id);
        else
            this->body_force_id = util::NO_INDEX_T;

        // setting vars from provided values
        PyObject* temp = python::loadAttrFromObject(this->simulation, "Output_File", python::PYSTRING, true);
        if (python::isNone(temp)) {
            temp = python::loadAttrFromObject(this->simulation, "Post_File", python::PYSTRING, true);
            if (python::isNone(temp)) {
                UERR("Must specify output file or post file in elmer");
            }
        }
        python::convertObjectToString(temp, this->node_data_file);
            
        // python::loadAttrFromObjectAndConvert(this->simulation, "Output_File", this->node_data_file);

        // creating a string from the elmer config
        python::convertObjectToString(this->elmer, this->sif_file_config_str);        
        SPARTA_NS::index_t length = (SPARTA_NS::index_t)strlen(this->sif_file_config_str);
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
        delete [] cells;
        delete [] cell_types;
        delete [] offsets;
        delete [] elements_additional_info;
        delete [] boundary_additional_info;
        delete [] node_velocities;
        delete [] node_temperatures;
        delete [] node_displacements;
        delete [] allow_forces;
        delete [] allow_heat_flux;
        delete [] clamped;
        delete this->server;
    }

    /* ---------------------------------------------------------------------- */

    bool Elmer::shouldUpdateSurf() { return this->elastic_solver_is_used; }

    /* ---------------------------------------------------------------------- */

    /**
     * makes sure file from elmer object data
     * returns: the name of the surf file
    */
    void Elmer::makeSpartaSurf(char* path) {
        ULOG("Writing Sparta surf to file:" + std::string(path));
        util::oFile f(path);
        f << "# Surface element file written by SPARTA fea interface\n\n" << this->num_nodes << " points\n" << this->num_boundary << " triangles\n\nPoints\n\n";

        SPARTA_NS::index_t i, j;
        for (i = 0; i < this->num_nodes; i++) {
            f << i+1;
            for (j = 0; j < ELMER_DIMENSION; j++)
                f << " " << this->nodes[3*i+j];
            f << "\n";
        }

        // ULOG("Diff:" + std::to_string(offsets[this->num_elements+1]-offsets[this->num_elements]));

        f << "\nTriangles\n\n";
        for (i = 0; i < this->num_boundary; i++) {
            f << i+1;
            for (j = 0; j < BOUNDARY_ELEMENT_NODE_COUNT; j++)
                f << " " << this->cells[offsets[this->num_elements+i]+j];
            f << "\n";
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::averageNodeTemperaturesInto(double*& _temperatures, SPARTA_NS::index_t _length) {
        // making sure everything is consistent
        if (_length != this->num_boundary)
            UERR("boundary data does not match required size, required size " + std::to_string(_length) + ", size " + std::to_string(this->num_boundary));

        // used later
        double avg;
        SPARTA_NS::index_t i, j;

        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        for (i = 0; i < this->num_boundary; i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            avg = 0;
            // gets the data point corresponding to node id and adds it to the rolling sum
            for (j = 0; j < BOUNDARY_ELEMENT_NODE_COUNT; j++)
                avg += this->node_temperatures[this->cells[offsets[i+num_elements]+j]-1];

            // computing the average by dividing the sum by the number of points and setting the surface element, the denominator should be 3 in the current implementation
            _temperatures[i] = avg/(BOUNDARY_ELEMENT_ARRAY_SIZE-BOUNDARY_ELEMENT_ARRAY_NODE_START); 
        }
    }

    /* ---------------------------------------------------------------------- */

    void Elmer::getNodePointAtIndex(SPARTA_NS::index_t _index, SPARTA_NS::index_t _boundary_index, double(&_point)[ELMER_DIMENSION]) {
        for (std::size_t j = 0; j < ELMER_DIMENSION; j++)
            _point[j] = this->nodes[3*this->cells[_index+num_elements+_boundary_index]-1+j];
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

                long num_elmer_timesteps = this->sparta->run_every/((long)r);
                double elmer_dt = *this->sparta->dt * (double)this->sparta->run_every/(double)num_elmer_timesteps;

                python::setAttrOfObject(this->simulation, "Timestep_intervals", num_elmer_timesteps);
                python::setAttrOfObject(this->simulation, "Output_Intervals",   num_elmer_timesteps);
                python::setAttrOfObject(this->simulation, "Timestep_Sizes",     elmer_dt);
            }
            // telling the code to no longer run this section of code
            this->has_been_initialized = true;

            // need to redo this because I set parameters, redoing this will add them to the string
            python::convertObjectToString(this->elmer, this->sif_file_config_str);
            SPARTA_NS::index_t length = (SPARTA_NS::index_t)strlen(this->sif_file_config_str);
            if (length < 3)
                UERR("Str method for elmer class returned a string with no length");
        }

        // making conditions
        this->createInitialConditions();
        this->createBoundaryConditions();

        // making the sif file
        // ULOG("Making sif file for elmer");

        util::oFile _buf(this->sif);
        // _buf << "! File created: " << util::getTime() << "\n";
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
        if (this->shouldUpdateSurf())
            this->dumpNodeVelocities();
        return;
        // // return;
        // char* filename;
        // // setting everything up, not much so just run every time
        // // cVTU_init();

        // int result = asprintf(&filename, "%s%cFix_FEA_ID_%s_step_t%ld.vtu", this->simulation_directory, SEP, this->sparta->id, *this->sparta->timestep);
        // if (result == -1) { UERR("asprintf failed for vtu file"); }

        // double *points = new double[3*num_nodes];

        // SPARTA_NS::index_t i, j, num_connection_nodes, rolling_offset;
        // for (i = 0; i < num_nodes; i++) {
        //     memcpy(&(points[3*i]), nodes[i], 3);
        // }

        // SPARTA_NS::index_t *connections;
        // cVTU_cell_type_t *types = new cVTU_cell_type_t[num_elements];
        // num_connection_nodes = 0;
        // for (i = 0; i < num_elements; i++) {
        //     num_connection_nodes+=element_lengths[i] - ELEMENT_ARRAY_NODE_START;
        // }

        // connections = new SPARTA_NS::index_t[num_connection_nodes];

        // rolling_offset = 0;
        // for (i = 0; i < num_elements; i++) {
        //     types[i] = typeCodeToVTKEnumType.at(elements[i][1]);
        //     for (j = ELEMENT_ARRAY_NODE_START; j < element_lengths[i]; j++) {
        //         connections[rolling_offset + j] = elements[i][j];
        //     }
        //     rolling_offset+=element_lengths[i] - ELEMENT_ARRAY_NODE_START;
        // }


        // unsigned num_point_data = 1;
        // char** point_data_names = new char*[4];
        // int point_data_types[] = {cVTU_FLOAT64, cVTU_FLOAT64, cVTU_FLOAT64, cVTU_FLOAT64};
        // char* first_name = (char*)"temperature";
        // point_data_names[0] = new char[strlen(first_name)+1];
        // strcpy(point_data_names[0], first_name);

        // if (shouldUpdateSurf()) {
        //     num_point_data+=3;
        //     char* second_name = (char*)"displacement _";
        //     point_data_names[1] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[1], "displacement x");
        //     point_data_names[2] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[2], "displacement y");
        //     point_data_names[3] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[3], "displacement z");
        // }

        // void* point_data = malloc(num_nodes*num_point_data*sizeof(double));
        
        // memcpy(point_data, this->node_temperatures, num_nodes*sizeof(double));

        // for (i = 1; i < num_point_data; i++) {
        //     for (j = 0; j < num_nodes; j++) {
        //         memcpy(&(((char*)point_data)[(i*num_nodes+j)*sizeof(double)]), &this->node_velocities[j][i-1], sizeof(double));
        //     }
        // }

        // cVTU_writeUnstructuredGridAppended(filename, points, num_nodes, connections, num_elements, types, point_data, num_point_data, point_data_types, (const char**)point_data_names, NULL, 0, NULL, NULL);

        // free(point_data);

        // for (i = 0; i < num_point_data; i++) delete [] point_data_names[i];
        // delete [] point_data_names;
        // delete [] points;
        // delete [] connections;
        // delete [] types;

        // // opening file for binary writing
        // FILE *fp = fopen(filename, "w"); // Open in write mode
        // if (fp == NULL) {
        //     UERR("failed to open file dump file");
        // }

        // // used for a variety of things
        // cVTU_uint_t i, j, offset = 0;
        // char buf[100];
        // char* buf_pointer;
        // unsigned temp_unsigned;

        // // write the xml tag specifying the version
        // cVTU_XML(fp);

        // // temporary variable
        // buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
        // if (buf_pointer == NULL) {
        //     UERR("failed to get vtk type from type name the for the header");
        // }

        // struct cVTU_ArgType1 temp_args[4];
        // temp_args[0] = {       "type",                                            (char*)"UnstructuredGrid"};
        // temp_args[1] = {    "version",                                                         (char*)"0.1"};
        // temp_args[2] = {"header_type",                                                          buf_pointer};
        // temp_args[3] = { "byte_order", cVTU_isLittleEndian() ? cVTU_VTK_LITTLE_ENDIAN : cVTU_VTK_BIG_ENDIAN};

        // // writing the xml tag for the start of the vtk data
        // cVTU_startXMLSectionArgType1(fp, "VTKFile", temp_args, 4);

        // // starting the unstructured grid data
        // cVTU_startXMLSectionJustName(fp, "UnstructuredGrid");

        // struct cVTU_ArgType2 temp_args_2[2];
        // temp_args_2[0] = { "NumberOfCells", this->num_elements};
        // temp_args_2[1] = {"NumberOfPoints",    this->num_nodes};

        // // starts the xml tag that holds the geometry
        // cVTU_startXMLSectionArgType2(fp, "Piece", temp_args_2, 2);

        // // starting the points section
        // cVTU_startXMLSectionJustName(fp, "Points");

        // // temporary variable
        // buf_pointer = cVTU_getVTKTypeFromTypeName(double);
        // if (buf_pointer == NULL) {
        //     UERR("failed to get vtk type from type name the for the points");
        // }
        
        // temp_args[0] = {"NumberOfComponents",        (char*)"3"};
        // temp_args[1] = {            "format", (char*)"appended"};
        // temp_args[2] = {            "offset",               buf};
        // temp_args[3] = {              "type",       buf_pointer};

        // // creating the data array pointing to the point data in the appended data
        // sprintf(buf, cVTU_getInt_t_format(), offset);
        // cVTU_DataArray(fp, temp_args, 4);

        // // point_size is made from being 3D, using the number of points and the size
        // // of each element, shift up after so that the next data arrays know it
        // cVTU_uint_t points_size = 3*num_nodes*sizeof(double);
        // offset+=points_size+sizeof(cVTU_uint_t); // add header size because header starts data

        // // ending point section
        // cVTU_endXMLSection(fp, "Points");

        // // starting cells section
        // cVTU_startXMLSectionJustName(fp, "Cells");

        // // temporary variable
        // buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
        // if (buf_pointer == NULL) {
        //     UERR("failed to get vtk type from type name the for cell connectivity");
        // }

        // temp_args[0] = {  "Name", (char*)"connectivity"};
        // temp_args[1] = {"format",     (char*)"appended"};
        // temp_args[2] = {"offset",                   buf};
        // temp_args[3] = {  "type",           buf_pointer};

        // // creating the data array pointing to the connectivity data in the appended data
        // sprintf(buf, cVTU_getInt_t_format(), offset);
        // cVTU_DataArray(fp, temp_args, 4);

        // // connections_size is made using the number of connections in each element and the size
        // // of each element in the array, shift up after so that the next data arrays know it
        // cVTU_uint_t connections_size = 0;
        // for (i = 0; i < num_elements; i++) {
        //     temp_unsigned = cVTU_getElementTypeNodeCount(typeCodeToVTKEnumType.at(elements[i][1]));
        //     if (temp_unsigned == 0) {
        //         UERR("failed to get element node count");
        //     }
        //     connections_size+=temp_unsigned;
        // }
        // connections_size*=sizeof(cVTU_uint_t);
        // offset+=connections_size+sizeof(cVTU_uint_t); // add header size because header starts data

        // // temporary variable
        // buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
        // if (buf_pointer == NULL) {
        //     UERR("failed to get vtk type from typename for cell offsets");
        // }

        // temp_args[0] = {  "Name",  (char*)"offsets"};
        // temp_args[1] = {"format", (char*)"appended"};
        // temp_args[2] = {"offset",               buf};
        // temp_args[3] = {  "type",       buf_pointer};

        // // creating the data array pointing to the offsets data in the appended data
        // sprintf(buf, cVTU_getInt_t_format(), offset);
        // cVTU_DataArray(fp, temp_args, 4);

        // offset+=num_elements*sizeof(cVTU_uint_t)+sizeof(cVTU_uint_t); // add header size because header starts data

        // // temporary variable
        // buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_CELL_TYPE);
        // if (buf_pointer == NULL) {
        //     UERR("failed to get vtk type from typename for cell types");
        // }

        // temp_args[0] = {  "Name",    (char*)"types"};
        // temp_args[1] = {"format", (char*)"appended"};
        // temp_args[2] = {"offset",               buf};
        // temp_args[3] = {  "type",       buf_pointer};

        // // creating the data array pointing to the cell_types data in the appended data
        // sprintf(buf, cVTU_getInt_t_format(), offset);
        // cVTU_DataArray(fp, temp_args, 4);

        // offset+=num_elements*sizeof(cVTU_cell_type_t)+sizeof(cVTU_uint_t); // add header size because header starts data

        // // end cells section
        // cVTU_endXMLSection(fp, "Cells");

        // // starting the point data section
        // cVTU_startXMLSectionJustName(fp, "PointData");

        // unsigned num_point_data = 1;
        // char** point_data_names = new char*[4];
        // char* first_name = (char*)"temperature";
        // point_data_names[0] = new char[strlen(first_name)+1];
        // strcpy(point_data_names[0], first_name);

        // if (shouldUpdateSurf()) {
        //     num_point_data+=3;
        //     char* second_name = (char*)"displacement _";
        //     point_data_names[1] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[1], "displacement x");
        //     point_data_names[2] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[2], "displacement y");
        //     point_data_names[3] = new char[strlen(second_name)+1];
        //     strcpy(point_data_names[3], "displacement z");
        // }

        // // iterating over the point data
        // for (i = 0; i < num_point_data; i++) {
        //     // creating the data array pointing to the per point data in the appended data
        //     sprintf(buf, cVTU_getInt_t_format(), offset);

        //     // temporary variable
        //     buf_pointer = cVTU_getVTKTypeFromTypeName(double);
        //     if (buf_pointer == NULL) {
        //         UERR("failed to get vtk type from point data type");
        //     }

        //     temp_args[0] = {  "Name", (char*)point_data_names[i]};
        //     temp_args[1] = {"format",          (char*)"appended"};
        //     temp_args[2] = {"offset",                        buf};
        //     temp_args[3] = {  "type",                buf_pointer};

        //     // add data array
        //     cVTU_DataArray(fp, temp_args, 4);
        //     offset+=num_nodes*sizeof(double)+sizeof(cVTU_uint_t);
        // }

        // // end point data section
        // cVTU_endXMLSection(fp, "PointData");

        // // start cell data section
        // cVTU_startXMLSectionJustName(fp, "CellData");

        // // end cell data section
        // cVTU_endXMLSection(fp, "CellData");

        // // end geometry section
        // cVTU_endXMLSection(fp, "Piece");

        // // end unstructured grid section
        // cVTU_endXMLSection(fp, "UnstructuredGrid");

        // temp_args[0] = {"encoding", (char*)"raw"};

        // // starting the appended data with raw binary encoding
        // cVTU_startXMLSectionArgType1(fp, "AppendedData", temp_args, 1);
        // fprintf(fp, "_"); // format of file specifies data must start with "_"

        // // writing the header of the points size to the file
        // fwrite(&points_size, sizeof(cVTU_uint_t), 1, fp);

        // // writes points to the file
        // for (i = 0; i < num_nodes; i++)
        //     fwrite(nodes[i], sizeof(double), 3, fp);

        // // writing the header of the connections size to the file
        // fwrite(&connections_size, sizeof(cVTU_uint_t), 1, fp);

        // // writes connections to the file
        // for (i = 0; i < num_elements; i++)
        //     fwrite(elements[i], sizeof(cVTU_uint_t), element_lengths[i], fp);

        // cVTU_uint_t temp_header; // just holds info

        // // writing the size of "offsets" array as the header of offsets data
        // temp_header = num_elements*sizeof(cVTU_uint_t);
        // fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);

        // // writing the offsets to the file
        // cVTU_uint_t temp_offset = 0;
        // for (i = 0; i < num_elements; i++) {
        //     // incrementing offset count for cell
        //     temp_unsigned = cVTU_getElementTypeNodeCount(typeCodeToVTKEnumType.at(elements[i][1]));
        //     if (temp_unsigned == 0) {
        //        UERR("failed to get vtk element node code from cell type");
        //     }
        //     temp_offset += temp_unsigned;

        //     // write offset
        //     fwrite(&temp_offset, sizeof(cVTU_uint_t), 1, fp);
        // }

        // // write header of cell_types size
        // temp_header = num_elements*sizeof(cVTU_cell_type_t);        
        // fwrite((&temp_header), sizeof(cVTU_uint_t), 1, fp);

        // // writing the cell types to the file
        // for (i = 0; i < num_elements; i++)
        //     fwrite(&typeCodeToVTKEnumType.at(elements[i][1]), sizeof(cVTU_cell_type_t), 1, fp);
    
        // cVTU_uint_t rolling_offset = 0;

        // // getting size of this point data
        // temp_header = num_nodes*sizeof(double);
        
        // // write header
        // fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);
        
        // // write point data
        // fwrite(this->node_temperatures, sizeof(double), num_nodes, fp);

        // // update rolling offset with the length of the current data array
        // rolling_offset+=temp_header;
        // // iterating over the point data
        // if (shouldUpdateSurf()) {
        //     for (i = 0; i < 3; i++) {
        //         // write header
        //         fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);
                
        //         // write point data
        //         for (j = 0; j < num_nodes; j++)
        //             fwrite(&(node_velocities[j][i]), sizeof(double), 1, fp);

        //         // update rolling offset with the length of the current data array
        //         rolling_offset+=temp_header;
        //     }
        // }

        // // new line to end data
        // fprintf(fp, "\n");

        // // close out sections and the file
        // cVTU_endXMLSection(fp, "AppendedData");
        // cVTU_endXMLSection(fp, "VTKFile");
        // fclose(fp);
    }
}