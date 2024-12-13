#include "cVTU.h"

#include <stdio.h>
#include <string.h>

#define cVTU_VTK_BIG_ENDIAN (char*)"BigEndian"
#define cVTU_VTK_LITTLE_ENDIAN (char*)"LittleEndian"
#define cVTU_TYPE_ID_STR(type) #type
#define cVTU_getVTKTypeFromTypeName(type) _cVTU_getVTKTypeFromTypeName(cVTU_TYPE_ID_STR(type))
#define cVTU_SET_TYPE_FROM_SIZE(type, type_name) (sizeof(type) < 3) ? ((sizeof(type) == 2) ? (cVTU_##type_name##16) :  (cVTU_##type_name##8)) : ((sizeof(type) == 4) ? (cVTU_##type_name##32) : (cVTU_##type_name##64));

// set type enums based on their size
const int cVTU_INT     = cVTU_SET_TYPE_FROM_SIZE(int, INT);
const int cVTU_BIGINT  = cVTU_SET_TYPE_FROM_SIZE(SPARTA_NS::bigint, INT);
const int cVTU_DOUBLE  = (sizeof(double) == 8) ? cVTU_FLOAT64 : cVTU_FLOAT32;
const int cVTU_INDEX_T = cVTU_SET_TYPE_FROM_SIZE(SPARTA_NS::index_t, UINT);


// used to write int_t values to files, determined later
static char *cVTU_int_t_format;
static char *cVTU_error_msg = NULL;

// used for passing arg pairs to functions
struct cVTU_ArgType1 {
    const char* arg1;
    char*       arg2;
};

// used for passing arg pairs to functions
struct cVTU_ArgType2 {
    const char* arg1;
    cVTU_uint_t arg2;
};

// holds info on a cell type for vtk, its string name, its place in the enum,
// and how many nodes is has
struct cVTU_VTKCellType {
    const char  name[21];
    int        enum_name;
    unsigned   num_nodes;
};

// holds info on a data types used in vtk, its string name, its place in the enum,
// how many bytes it is, and its string equivalent name in C
struct cVTU_VTKType {
    const char     name[8];
    int          enum_name;
    unsigned          size;
    const char data_type[9];
};

// holds info on the valid types used in vtk
const struct cVTU_VTKType cVTU_VTKTypes[] = {
    {"Float32", cVTU_FLOAT32, 4,    "float"},
    {"Float64", cVTU_FLOAT64, 8,   "double"},
    {   "Int8",    cVTU_INT8, 1,   "int8_t"},
    {  "UInt8",   cVTU_UINT8, 1,  "uint8_t"},
    {  "Int16",   cVTU_INT16, 2,  "int16_t"},
    { "UInt16",  cVTU_UINT16, 2, "uint16_t"},
    {  "Int32",   cVTU_INT32, 4,  "int32_t"},
    { "UInt32",  cVTU_UINT32, 4, "uint32_t"},
    {  "Int64",   cVTU_INT64, 8,  "int64_t"},
    { "UInt64",  cVTU_UINT64, 8, "uint64_t"}
};

// holds info on the valid cell types that are permited, not all that
// vtk has to offer, but there are the one that I choose to consider
const struct cVTU_VTKCellType cVTU_VTKCellTypes[] = {
    {              "vertex",               cVTU_VERTEX,  1},
    {                "line",                 cVTU_LINE,  2},
    {            "triangle",             cVTU_TRIANGLE,  3},
    {                "quad",                 cVTU_QUAD,  4},
    {               "tetra",                cVTU_TETRA,  4},
    {          "hexahedron",           cVTU_HEXAHEDRON,  8},
    {               "wedge",                cVTU_WEDGE,  6},
    {             "pyramid",              cVTU_PYRAMID,  5},
    {      "quadratic_edge",       cVTU_QUADRATIC_EDGE,  3},
    {  "quadratic_triangle",   cVTU_QUADRATIC_TRIANGLE,  6},
    {      "quadratic_quad",       cVTU_QUADRATIC_QUAD,  8},
    {     "quadratic_tetra",      cVTU_QUADRATIC_TETRA, 10},
    {"quadratic_hexahedron", cVTU_QUADRATIC_HEXAHEDRON, 20}
};

const char* cVTU_getInt_t_format() {
    return cVTU_int_t_format;
}

// gets the vtk data type name from the name in C
char* _cVTU_getVTKTypeFromTypeName(const char* type_name) {
    for (unsigned i = 0; i < sizeof(cVTU_VTKTypes)/sizeof(struct cVTU_VTKType); i++) {
        if (strcmp(cVTU_VTKTypes[i].data_type, type_name) == 0) return (char*)cVTU_VTKTypes[i].name;
    }
    return NULL;
}

// gets the vtk data type name from the enum value
char* cVTU_getVTKTypeFromEnum(const int enum_name) {
    for (unsigned i = 0; i < sizeof(cVTU_VTKTypes)/sizeof(struct cVTU_VTKType); i++) {
        if (cVTU_VTKTypes[i].enum_name == enum_name) return (char*)cVTU_VTKTypes[i].name;
    }
    return NULL;
}

// gets the vtk data type size from the enum value
unsigned cVTU_getVTKTypeSizeFromEnum(const int enum_name) {
    for (unsigned i = 0; i < sizeof(cVTU_VTKTypes)/sizeof(struct cVTU_VTKType); i++) {
        if (cVTU_VTKTypes[i].enum_name == enum_name) return cVTU_VTKTypes[i].size;
    }
    return 0;
}

// gets the number of nodes making up a cell type from the cell type enum value
unsigned cVTU_getElementTypeNodeCount(const cVTU_cell_type_t element_type) {
    for (unsigned i = 0; i < sizeof(cVTU_VTKCellTypes)/sizeof(struct cVTU_VTKCellType); i++) {
        if (cVTU_VTKCellTypes[i].enum_name == element_type) return cVTU_VTKCellTypes[i].num_nodes;
    }
    return 0;
}

/*
* Sets up some attributes
*/
void cVTU_init() {
    // setup stuff for sparta
    switch (sizeof(cVTU_uint_t)) {
        case 1:
            cVTU_int_t_format = (char*)"%hhu";
            break;
        case 2:
            cVTU_int_t_format = (char*)"%hu";
            break;
        case 4:
            cVTU_int_t_format = (char*)"%u";
            break;
        case 8:
            cVTU_int_t_format = (char*)"%llu";
            break;
    }
}

/*
* Determines if the computer uses little endian format for its data
*/
unsigned cVTU_isLittleEndian() {
    int test_num = 1;
    char *ptr = (char*)&test_num;
    return (ptr[0] == 1); 
}

/*
* writes starting xml tag to the file
*/
void cVTU_XML(FILE *fp) { fprintf(fp, "<?xml version=\"1.0\"?>\n"); }

/*
* Writes an xml tag with the name to the file with the provided parameters with
* string, string format
*/
void cVTU_startXMLSectionJustName(FILE *fp, const char* name) { fprintf(fp, "<%s>\n", name); }

/*
* Writes an xml tag with the name to the file with the provided parameters with
* string, string format
*/
void cVTU_startXMLSectionArgType1(FILE *fp, const char* name, struct cVTU_ArgType1 *parameters, const cVTU_uint_t len_parameters) {
    fprintf(fp, "<%s", name);
    for (cVTU_uint_t i = 0; i < len_parameters; i++)
        fprintf(fp, " %s=\"%s\"", parameters[i].arg1, parameters[i].arg2);
    fprintf(fp, ">\n");
}

/*
* Writes an xml tag with the name to the file with the provided parameters with
* string, int_t format
*/
void cVTU_startXMLSectionArgType2(FILE *fp, const char* name, struct cVTU_ArgType2 *parameters, const cVTU_uint_t len_parameters) {
    struct cVTU_ArgType1 args[len_parameters];
    for (cVTU_uint_t i = 0; i < len_parameters; i++) {
        args[i].arg1 = parameters[i].arg1;
        args[i].arg2 = (char*)malloc(100*sizeof(char));
        sprintf(args[i].arg2, cVTU_int_t_format, parameters[i].arg2);
    }
    cVTU_startXMLSectionArgType1(fp, name, args, len_parameters);
    for (cVTU_uint_t i = 0; i < len_parameters; i++) 
        free(args[i].arg2);
}

/*
* Write end section with the given name to the file
*/
void cVTU_endXMLSection(FILE *fp, const char* name) { fprintf(fp, "</%s>\n", name); }

/*
* Writes an inline xml tag with the name to the file with the provided parameters 
* with string, string format
*/
void cVTU_inlineXMLSection(FILE *fp, const char* name, struct cVTU_ArgType1 *parameters, const cVTU_uint_t len_parameters) {
    fprintf(fp, "<%s", name);
    for (cVTU_uint_t i = 0; i < len_parameters; i++)
        fprintf(fp, " %s=\"%s\"", parameters[i].arg1, parameters[i].arg2);
    fprintf(fp, "/>\n");
}

/*
* Writes a <DataArray .../> tag with the name to the file with the provided
* parameters with string, string format
*/
void cVTU_DataArray(FILE *fp, struct cVTU_ArgType1 *parameters, const cVTU_uint_t len_parameters) {
    cVTU_inlineXMLSection(fp, "DataArray", parameters, len_parameters);
}

/*
* Writes a <PDataArray .../> tag with the name to the file with the provided
* parameters with string, string format
*/
void cVTU_PDataArray(FILE *fp, struct cVTU_ArgType1 *parameters, const cVTU_uint_t len_parameters) {
    cVTU_inlineXMLSection(fp, "PDataArray", parameters, len_parameters);
}

/**
 * returns error message if error is encountered
*/
char* cVTU_checkForErrors() { return cVTU_error_msg; }

/**
 Writes an unstructured grid to the file with the given info
 Inputs:
 - filename         : ends in ".vtu", file to write to
 - points           : double array with size 3*num_points*sizeof(double), 
                      i.e. list all points with components in one array
 - num_points       : number of points, not the length of the points array,
                      just the number of 3 component points, i.e.
                      num_points = sizeof(points)/3
 - connections      : connection_t array that contains the connectivity of each cell
                      to the points in the points array. Just a list of 
                      connection_t's (unsigned long's)
 - num_cells        : number of cells in the connections array.
 - cell_types       : cell_type_t array containing the types of each cell in connections
                      array, i.e. array containing QUAD, LINE, ....
 - point_data       : array with length num_point_data*num_points, each sub array of
                      length num_points is interpreted as data, specify the type of 
                      this data in the point_data_types array
 - num_point_data   : number of data sets that correspond to all the points, it
                      is also the number of arrays of length num_points within 
                      point_data
 - point_data_types : array of ints specifying the type of each sub-array
                      within the point_data array, options are DOUBLE, ...
 - point_data_names : array of c-strings with length num_point_data containing 
                      the names of the point_data sets
 - cell_data        : array with length num_cell_data*num_cells, each sub array of
                      length num_cells is interpreted as data, specify the type of 
                      this data in the cell_data_types array
 - num_cell_data    : number of data sets that correspond to all the cells, it
                      is also the number of arrays of length num_cells within 
                      cell_data
 - cell_data_types  : array of ints specifying the type of each sub-array
                      within the cell_data array, options are DOUBLE, ...
 - cell_data_names  : array of c-strings with length num_cell_data containing 
                      the names of the cell_data sets
*/
void cVTU_writeUnstructuredGridAppended(const char* filename, 
        const double   *points, const cVTU_uint_t     num_points, const cVTU_uint_t    *connections, const cVTU_uint_t        num_cells, const cVTU_cell_type_t *cell_types, 
        const void* point_data, const cVTU_uint_t num_point_data, const int *point_data_types, const char** point_data_names, 
        const void*  cell_data, const cVTU_uint_t  num_cell_data, const int  *cell_data_types, const char**  cell_data_names) 
    {
    // setting everything up, not much so just run every time
    cVTU_init();

    // opening file for binary writing
    FILE *fp = fopen(filename, "w"); // Open in write mode
    if (fp == NULL) {
        cVTU_error_msg = (char*)"failed to open file";
        return;
    }

    // used for a variety of things
    cVTU_uint_t i, offset = 0;
    char buf[100];
    char* buf_pointer;
    unsigned temp_unsigned;

    // write the xml tag specifying the version
    cVTU_XML(fp);

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from type name the for the header";
        return;
    }

    struct cVTU_ArgType1 temp_args[4];
    temp_args[0] = {       "type",                                            (char*)"UnstructuredGrid"};
    temp_args[1] = {    "version",                                                         (char*)"0.1"};
    temp_args[2] = {"header_type",                                                          buf_pointer};
    temp_args[3] = { "byte_order", cVTU_isLittleEndian() ? cVTU_VTK_LITTLE_ENDIAN : cVTU_VTK_BIG_ENDIAN};

    // writing the xml tag for the start of the vtk data
    cVTU_startXMLSectionArgType1(fp, "VTKFile", temp_args, 4);

    // starting the unstructured grid data
    cVTU_startXMLSectionJustName(fp, "UnstructuredGrid");

    struct cVTU_ArgType2 temp_args_2[2];
    temp_args_2[0] = { "NumberOfCells",  num_cells};
    temp_args_2[1] = {"NumberOfPoints", num_points};

    // starts the xml tag that holds the geometry
    cVTU_startXMLSectionArgType2(fp, "Piece", temp_args_2, 2);

    // starting the points section
    cVTU_startXMLSectionJustName(fp, "Points");

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(double);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from type name the for the points";
        return;
    }
    
    temp_args[0] = {"NumberOfComponents",        (char*)"3"};
    temp_args[1] = {            "format", (char*)"appended"};
    temp_args[2] = {            "offset",               buf};
    temp_args[3] = {              "type",       buf_pointer};

    // creating the data array pointing to the point data in the appended data
    sprintf(buf, cVTU_int_t_format, offset);
    cVTU_DataArray(fp, temp_args, 4);

    // point_size is made from being 3D, using the number of points and the size
    // of each element, shift up after so that the next data arrays know it
    cVTU_uint_t points_size = 3*num_points*sizeof(double);
    offset+=points_size+sizeof(cVTU_uint_t); // add header size because header starts data

    // ending point section
    cVTU_endXMLSection(fp, "Points");

    // starting cells section
    cVTU_startXMLSectionJustName(fp, "Cells");

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from type name the for cell connectivity";
        return;
    }

    temp_args[0] = {  "Name", (char*)"connectivity"};
    temp_args[1] = {"format",     (char*)"appended"};
    temp_args[2] = {"offset",                   buf};
    temp_args[3] = {  "type",           buf_pointer};

    // creating the data array pointing to the connectivity data in the appended data
    sprintf(buf, cVTU_int_t_format, offset);
    cVTU_DataArray(fp, temp_args, 4);

    // connections_size is made using the number of connections in each element and the size
    // of each element in the array, shift up after so that the next data arrays know it
    cVTU_uint_t connections_size = 0;
    for (i = 0; i < num_cells; i++) {
        temp_unsigned = cVTU_getElementTypeNodeCount(cell_types[i]);
        if (temp_unsigned == 0) {
            cVTU_error_msg = (char*)"failed to get element node count";
            return;
        }
        connections_size+=temp_unsigned;
    }
    connections_size*=sizeof(cVTU_uint_t);
    offset+=connections_size+sizeof(cVTU_uint_t); // add header size because header starts data

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from typename for cell offsets";
        return;
    }

    temp_args[0] = {  "Name",  (char*)"offsets"};
    temp_args[1] = {"format", (char*)"appended"};
    temp_args[2] = {"offset",               buf};
    temp_args[3] = {  "type",       buf_pointer};

    // creating the data array pointing to the offsets data in the appended data
    sprintf(buf, cVTU_int_t_format, offset);
    cVTU_DataArray(fp, temp_args, 4);

    offset+=num_cells*sizeof(cVTU_uint_t)+sizeof(cVTU_uint_t); // add header size because header starts data

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_CELL_TYPE);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from typename for cell types";
        return;
    }

    temp_args[0] = {  "Name",    (char*)"types"};
    temp_args[1] = {"format", (char*)"appended"};
    temp_args[2] = {"offset",               buf};
    temp_args[3] = {  "type",       buf_pointer};

    // creating the data array pointing to the cell_types data in the appended data
    sprintf(buf, cVTU_int_t_format, offset);
    cVTU_DataArray(fp, temp_args, 4);

    offset+=num_cells*sizeof(cVTU_cell_type_t)+sizeof(cVTU_uint_t); // add header size because header starts data

    // end cells section
    cVTU_endXMLSection(fp, "Cells");

    // starting the point data section
    cVTU_startXMLSectionJustName(fp, "PointData");

    // if point data is provided
    if (point_data != NULL && num_point_data > 0) {
        // checks
        if (point_data_names == NULL) {
            cVTU_error_msg = (char*)"point_data_names can not be null";
            return;
        }
        if (point_data_types == NULL) {
            cVTU_error_msg = (char*)"point_data_types can not be null";
            return;
        }

        // iterating over the point data
        for (i = 0; i < num_point_data; i++) {
            // creating the data array pointing to the per point data in the appended data
            sprintf(buf, cVTU_int_t_format, offset);

            // temporary variable
            buf_pointer = cVTU_getVTKTypeFromEnum(point_data_types[i]);
            if (buf_pointer == NULL) {
                cVTU_error_msg = (char*)"failed to get vtk type from point data type";
                return;
            }

            temp_args[0] = {  "Name", (char*)point_data_names[i]};
            temp_args[1] = {"format",          (char*)"appended"};
            temp_args[2] = {"offset",                        buf};
            temp_args[3] = {  "type",                buf_pointer};

            // add data array
            cVTU_DataArray(fp, temp_args, 4);
            // add header size because header starts data in array
            temp_unsigned = cVTU_getVTKTypeSizeFromEnum(point_data_types[i]);
            if (temp_unsigned == 0) {
                cVTU_error_msg = (char*)"failed to get vtk type size from point data type";
                return;
            }
            offset+=num_points*temp_unsigned+sizeof(cVTU_uint_t);
        }
    }

    // end point data section
    cVTU_endXMLSection(fp, "PointData");

    // start cell data section
    cVTU_startXMLSectionJustName(fp, "CellData");

    // if cell data is provided
    if (cell_data != NULL && num_cell_data > 0) {
        // checks
        if (cell_data_names == NULL)  {
            cVTU_error_msg = (char*)"cell_data_names can not be null";
            return;
        }
        if (cell_data_types == NULL) {
            cVTU_error_msg = (char*)"cell_data_types can not be null";
            return;
        }

        // creating the data array pointing to the per cell data in the appended data
        for (i = 0; i < num_cell_data; i++) {
            // creating the data array pointing to the per cell data in the appended data
            sprintf(buf, cVTU_int_t_format, offset);

            // temporary variable
            buf_pointer = cVTU_getVTKTypeFromEnum(cell_data_types[i]);
            if (buf_pointer == NULL) {
                cVTU_error_msg = (char*)"failed to get vtk type from cell data type";
                return;
            }

            temp_args[0] = {  "Name", (char*)cell_data_names[i]};
            temp_args[1] = {"format",         (char*)"appended"};
            temp_args[2] = {"offset",                       buf};
            temp_args[3] = {  "type",               buf_pointer};

            // add data array
            cVTU_DataArray(fp, temp_args, 4);

            // add header size because header starts data in array
            temp_unsigned = cVTU_getVTKTypeSizeFromEnum(cell_data_types[i]);
            if (temp_unsigned == 0) {
                cVTU_error_msg = (char*)"failed to get vtk type size from cell data type";
                return;
            }
            offset+=num_cells*temp_unsigned+sizeof(cVTU_uint_t);
        }
    }

    // end cell data section
    cVTU_endXMLSection(fp, "CellData");

    // end geometry section
    cVTU_endXMLSection(fp, "Piece");

    // end unstructured grid section
    cVTU_endXMLSection(fp, "UnstructuredGrid");

    temp_args[0] = {"encoding", (char*)"raw"};

    // starting the appended data with raw binary encoding
    cVTU_startXMLSectionArgType1(fp, "AppendedData", temp_args, 1);
    fprintf(fp, "_"); // format of file specifies data must start with "_"

    // writing the header of the points size to the file
    fwrite(&points_size, sizeof(cVTU_uint_t), 1, fp);

    // writes points to the file
    fwrite(points, sizeof(double), points_size/sizeof(double), fp);

    // writing the header of the connections size to the file
    fwrite(&connections_size, sizeof(cVTU_uint_t), 1, fp);

    // writes connections to the file
    fwrite(connections, sizeof(cVTU_uint_t), connections_size/sizeof(cVTU_uint_t), fp);

    cVTU_uint_t temp_header; // just holds info

    // writing the size of "offsets" array as the header of offsets data
    temp_header = num_cells*sizeof(cVTU_uint_t);
    fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);

    // writing the offsets to the file
    cVTU_uint_t temp_offset = 0;
    for (i = 0; i < num_cells; i++) {
        // incrementing offset count for cell
        temp_unsigned = cVTU_getElementTypeNodeCount(cell_types[i]);
        if (temp_unsigned == 0) {
            cVTU_error_msg = (char*)"failed to get vtk element node code from cell type";
            return;
        }
        temp_offset += temp_unsigned;

        // write offset
        fwrite(&temp_offset, sizeof(cVTU_uint_t), 1, fp);
    }

    // write header of cell_types size
    temp_header = num_cells*sizeof(cVTU_cell_type_t);        
    fwrite((&temp_header), sizeof(cVTU_uint_t), 1, fp);

    // writing the cell types to the file
    fwrite(cell_types, sizeof(cVTU_cell_type_t), num_cells, fp);

    // if point data is provided
    if (point_data != NULL && num_point_data > 0) {
        cVTU_uint_t rolling_offset = 0;
        // iterating over the point data
        for (i = 0; i < num_point_data; i++) {
            // getting size of this point data
            temp_unsigned = cVTU_getVTKTypeSizeFromEnum(point_data_types[i]);
            if (temp_unsigned == 0) {
                cVTU_error_msg = (char*)"failed to get vtk type size from point data type";
                return;
            }
            temp_header = num_points*temp_unsigned;
            
            // write header
            fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);
            
            // write point data
            fwrite(&(((char*)point_data)[rolling_offset]), temp_unsigned, num_points, fp);

            // update rolling offset with the length of the current data array
            rolling_offset+=temp_header;
        }
    }

    // if cell data is provided
    if (cell_data != NULL && num_cell_data > 0) {
        cVTU_uint_t rolling_offset = 0;
        // creating the data array pointing to the per cell data in the appended data
        for (i = 0; i < num_cell_data; i++) {
            // getting size of this cell data
            temp_unsigned = cVTU_getVTKTypeSizeFromEnum(cell_data_types[i]);
            if (temp_unsigned == 0) {
                cVTU_error_msg = (char*)"failed to get vtk type size from cell data type";
                return;
            }
            temp_header = num_cells*temp_unsigned;
            
            // write header
            fwrite(&temp_header, sizeof(cVTU_uint_t), 1, fp);
            
            // write cell data
            fwrite(&(((char*)cell_data)[rolling_offset]), temp_unsigned, num_cells, fp);
            
            // update rolling offset with the length of the current data array
            rolling_offset+=temp_header;
        }
    }

    // new line to end data
    fprintf(fp, "\n");

    // close out sections and the file
    cVTU_endXMLSection(fp, "AppendedData");
    cVTU_endXMLSection(fp, "VTKFile");
    fclose(fp);
}

/**
 Writes an unstructured grid to the file with the given info
 Inputs:
 - filename         : ends in ".pvtu", file to write to

 - num_point_data   : number of data sets that correspond to all the points, it
                      is also the number of arrays of length num_points within 
                      point_data
 - point_data_types : array of ints specifying the type of each sub-array
                      within the point_data array, options are DOUBLE, ...
 - point_data_names : array of c-strings with length num_point_data containing 
                      the names of the point_data sets
 - num_cell_data    : number of data sets that correspond to all the cells, it
                      is also the number of arrays of length num_cells within 
                      cell_data
 - cell_data_types  : array of ints specifying the type of each sub-array
                      within the cell_data array, options are DOUBLE, ...
 - cell_data_names  : array of c-strings with length num_cell_data containing 
                      the names of the cell_data sets
*/
void cVTU_writeParallelUnstructuredGridHeadFile(
                    const char* filename, const char **data_filenames, const cVTU_uint_t num_data_files,
        const cVTU_uint_t num_point_data, const int *point_data_types,    const char** point_data_names, 
        const cVTU_uint_t  num_cell_data, const int  *cell_data_types,    const char**  cell_data_names)
    {
    // setting everything up, not much so just run every time
    cVTU_init();

    // opening file for binary writing
    FILE *fp = fopen(filename, "w"); // Open in write mode
    if (fp == NULL) {
        cVTU_error_msg = (char*)"failed to open file";
        return;
    }

    // used for a variety of things
    cVTU_uint_t i = 0;
    char* buf_pointer;

    // write the xml tag specifying the version
    cVTU_XML(fp);

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(cVTU_UINT_TYPE);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from type name the for the header";
        return;
    }

    struct cVTU_ArgType1   temp_args[3];

    temp_args[0] = {       "type",                                           (char*)"PUnstructuredGrid"};
    temp_args[1] = {    "version",                                                         (char*)"0.1"};
    temp_args[2] = { "byte_order", cVTU_isLittleEndian() ? cVTU_VTK_LITTLE_ENDIAN : cVTU_VTK_BIG_ENDIAN};

    // writing the xml tag for the start of the vtk data
    cVTU_startXMLSectionArgType1(fp, "VTKFile", temp_args, 3);

    temp_args[0] = {"GhostLevel", (char*)"0"};

    // starting the unstructured grid data
    cVTU_startXMLSectionArgType1(fp, "PUnstructuredGrid", temp_args, 1);

    // starting the points section
    cVTU_startXMLSectionJustName(fp, "PPoints");

    // temporary variable
    buf_pointer = cVTU_getVTKTypeFromTypeName(double);
    if (buf_pointer == NULL) {
        cVTU_error_msg = (char*)"failed to get vtk type from type name the for the points";
        return;
    }

    temp_args[0] = {"NumberOfComponents",        (char*)"3"};
    temp_args[1] = {              "type",       buf_pointer};
    
    // creating the data array pointing to the point data in the appended data
    cVTU_PDataArray(fp, temp_args, 2);

    // ending point section
    cVTU_endXMLSection(fp, "PPoints");

    // write each file with the connectivity info and other data as a piece
    for (i = 0; i < num_data_files; i++) {
        temp_args[0] = {"Source", (char*)data_filenames[i]};
        // starting cells section
        cVTU_inlineXMLSection(fp, "Piece", temp_args, 1);
    }

    // starting the point data section
    cVTU_startXMLSectionJustName(fp, "PPointData");

    // if point data is provided
    if (point_data_names != NULL && num_point_data > 0) {
        // checks
        if (point_data_types == NULL) {
            cVTU_error_msg = (char*)"point_data_types can not be null";
            return;
        }

        // iterating over the point data
        for (i = 0; i < num_point_data; i++) {
            temp_args[0] = {  "Name", (char*)point_data_names[i]};
            temp_args[1] = {  "type",                buf_pointer};
            // add data array
            cVTU_PDataArray(fp, temp_args, 2);
        }
    }

    // end point data section
    cVTU_endXMLSection(fp, "PPointData");

    // start cell data section
    cVTU_startXMLSectionJustName(fp, "PCellData");

    // if cell data is provided
    if (cell_data_names != NULL && num_cell_data > 0) {
        // checks
        if (cell_data_types == NULL) {
            cVTU_error_msg = (char*)"cell_data_types can not be null";
            return;
        }

        // creating the data array pointing to the per cell data in the appended data
        for (i = 0; i < num_cell_data; i++) {
            // temporary variable
            buf_pointer = cVTU_getVTKTypeFromEnum(cell_data_types[i]);
            if (buf_pointer == NULL) {
                cVTU_error_msg = (char*)"failed to get vtk type from cell data type";
                return;
            }

            temp_args[0] = {  "Name", (char*)cell_data_names[i]};
            temp_args[1] = {  "type",               buf_pointer};

            // add data array
            cVTU_PDataArray(fp, temp_args, 2);
        }
    }

    // end cell data section
    cVTU_endXMLSection(fp, "PCellData");

    // end unstructured grid section
    cVTU_endXMLSection(fp, "PUnstructuredGrid");
    // close out sections and the file
    cVTU_endXMLSection(fp, "VTKFile");
    fclose(fp);
}


#undef cVTU_INT_TYPE
#undef cVTU_CELL_TYPE
#undef cVTU_VTK_BIG_ENDIAN
#undef cVTU_VTK_LITTLE_ENDIAN
#undef cVTU_TYPE_ID_STR
#undef cVTU_getVTKTypeFromTypeName

/* ----------------------------------------------------------------------

    end cVTU stuff

   ---------------------------------------------------------------------- */