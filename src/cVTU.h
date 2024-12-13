#ifndef cVTU_H
#define cVTU_H

#include <stdlib.h>
#include "spatype.h"

#define cVTU_UINT_TYPE SPA_INDEX_T
#define cVTU_CELL_TYPE uint8_t

typedef cVTU_UINT_TYPE      cVTU_uint_t;
typedef cVTU_CELL_TYPE cVTU_cell_type_t;

enum{
    cVTU_FLOAT32,
    cVTU_FLOAT64,
    cVTU_INT8,
    cVTU_UINT8,
    cVTU_INT16,
    cVTU_UINT16,
    cVTU_INT32,
    cVTU_UINT32,
    cVTU_INT64,
    cVTU_UINT64
};

extern const int cVTU_INT, cVTU_BIGINT, cVTU_DOUBLE, cVTU_INDEX_T;

enum{
    cVTU_VERTEX         = 1,
    cVTU_LINE           = 3,
    cVTU_TRIANGLE       = 5,
    cVTU_QUAD           = 9,
    cVTU_TETRA,
    cVTU_HEXAHEDRON     = 12,
    cVTU_WEDGE,
    cVTU_PYRAMID,
    cVTU_QUADRATIC_EDGE = 21,
    cVTU_QUADRATIC_TRIANGLE,
    cVTU_QUADRATIC_QUAD,
    cVTU_QUADRATIC_TETRA,
    cVTU_QUADRATIC_HEXAHEDRON
};


// gets the vtk data type name from the enum value
char* cVTU_getVTKTypeFromEnum(const int enum_name);


// gets the vtk data type size from the enum value
unsigned cVTU_getVTKTypeSizeFromEnum(const int enum_name);


// gets the number of nodes making up a cell type from the cell type enum value
unsigned cVTU_getElementTypeNodeCount(const cVTU_cell_type_t element_type);


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
        const double   *points, const cVTU_uint_t     num_points, const cVTU_uint_t *connections,   const cVTU_uint_t num_cells, const cVTU_cell_type_t *cell_types, 
        const void* point_data, const cVTU_uint_t num_point_data,    const int *point_data_types, const char** point_data_names, 
        const void*  cell_data, const cVTU_uint_t  num_cell_data,    const int  *cell_data_types, const char**  cell_data_names);

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
        const cVTU_uint_t  num_cell_data, const int  *cell_data_types,    const char**  cell_data_names);


/**
 * returns error message if error is encountered
*/
char* cVTU_checkForErrors();

#endif