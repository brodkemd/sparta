/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This file was contributed by Marek Brodke, a student at the University
   of Cincinnati
------------------------------------------------------------------------- */

#include "error.h"
#include "modify.h"
#include "compute.h"
#include "grid.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "memory.h"


#include <string.h>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>

#include "fix_paraview.h"
#include "cVTU.h"
#include "Elmer/util.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

// DO NOT CHANGE THIS
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32
#define _START_ARGS 6

#define SEP '/'

#define PRINT(...) if (comm->me == 0) { if (screen) { fprintf(screen, __VA_ARGS__); } if (logfile) { fprintf(logfile, __VA_ARGS__); } }
#define PRINT_PROC(...) if (screen) { fprintf(screen, __VA_ARGS__); } if (logfile) { fprintf(logfile, __VA_ARGS__); }
// #define PRINT_ALL_PROCS(num_procs, ...) MPI_Barrier(world); for (int __i = 0; __i < num_procs; __i++) { if (comm->me == __i) { PRINT_PROC(__VA_ARGS__); } }
#define ERROR(msg) error->all(FLERR,msg);

// do not change these
unsigned START_ARGS;
enum{GRID,SURF};
enum{_INT,_DOUBLE};


/**
 * NOTE: when setting char*, I used new char[SOME_VALUE+1], the 1 ensures the string is null terminated
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
 * command format:
 *        fix ID paraview group NEVERY PREFIX STYLE PER_STYLE_ARGS ARGS
 * style = grid or surf
 * COLLECT_INTO_ONE_FILE? = yes or no
*/
FixParaview::FixParaview(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    PRINT("Fix Paraview:\n");

    // initialize to the default starting point of the arguments
    START_ARGS = _START_ARGS;

    // checking if the number of args is correct
    if (narg < START_ARGS) { ERROR("Illegal fix paraview command, not enough arguments"); }

    PRINT("  id: %s\n", id);

    // how often to run
    nevery = atoi(arg[3]);

    // set the prefix as to where to save data to
    this->prefix = new char[strlen(arg[4])+1];
    strcpy(this->prefix, arg[4]);

    // trimming trailing / is it is on the path
    if (prefix[strlen(prefix)-1] == '/')
        prefix[strlen(prefix)-1] = '\0';

    // ensures data path exists
    struct stat info;
    if (stat(prefix, &info) != 0) {
        ERROR("Could not find directory at the provided prefix");
    } else {
         if (S_ISDIR(info.st_mode)) {
            PRINT("  Saving to the directory: %s\n", prefix);
        } else {
            ERROR("provided prefix is not a directory");
        }
    }

    // set style and get the groupbit for geometry to look at
    if (strcmp(arg[5], "grid") == 0) {
        style = GRID;

        // getting the grid cells we are looking
        int igroup = grid->find_group(arg[2]);
        if (igroup < 0) ERROR("fix paraview group ID does not exist");
        groupbit = grid->bitmask[igroup];

        // tell code to collect all data into one file if the user wants
        if (strcmp(arg[6], "yes") == 0) {
            collect_into_one_file = true;
        } else if (strcmp(arg[6], "no") == 0) {
            collect_into_one_file = false;
        } else {
            ERROR("invalid value for collect_into_on_file in fix paraview");
        }

        // has 1 optional arg, shift starting index for rest of args by 1
        START_ARGS+=1;

    } else if (strcmp(arg[5], "surf") == 0) {
        style = SURF;
        
        // get surface group we are looking for
        int igroup = surf->find_group(arg[2]);
        if (igroup < 0) ERROR("fix paraview group ID does not exist");
        groupbit = surf->bitmask[igroup];

        // force collection to main process for surface, many are so small it does not need the
        // same special treatment a grid does
        collect_into_one_file = true;

    } else {
        ERROR("invalid style for fix paraview");
    }

    // expand args, particularly, expanded at "*"
    int expand = 0; char **expanded_args;
    num_expanded_args = (index_t)input->expand_args(narg-START_ARGS, &arg[START_ARGS], 1, expanded_args);

    // if the arguments were expanded
    if (expanded_args != &arg[START_ARGS]) expand = 1;

    // allocate field arrays
    funcs      = new FnPtrPack[num_expanded_args];
    argindex   = new int[num_expanded_args];
    func_types = new int[num_expanded_args];
    data_names = new char*[num_expanded_args];
    ids        = new char*[num_expanded_args];

    // set to default value
    there_is_variable = false;

    // handles different styles,
    if (style == GRID) {
        handleGridArgs(num_expanded_args, expanded_args);
    } else {
        handleSurfArgs(num_expanded_args, expanded_args);
    }

    // if wildcard expansion occurred, free earg memory from expand_args()
    // wait to do this until after column string is created
    if (expand) {
        for (int i = 0; i < num_expanded_args; i++) 
            delete [] expanded_args[i];
        memory->sfree(expanded_args);
    }
}

/**
 * handles grid arguments
*/
void FixParaview::handleGridArgs(int argc, char** argv) {
    // loop over args
    for (int j = 0; j < argc; j++) {
        // remembering arg name as a data name, labels the data in the vtu file
        data_names[j] = new char[strlen(argv[j])+1];
        strcpy(data_names[j], argv[j]);

        // handle the different arguments
        if (strcmp(argv[j], "id") == 0) {
            funcs[j] = &FixParaview::handleGridId;
            if (sizeof(cellint) == sizeof(smallint)) {
                func_types[j] = cVTU_INT;
            } else {
                func_types[j] = cVTU_BIGINT;
            }
        } else if (strcmp(argv[j], "split") == 0) {
            funcs[j] = &FixParaview::handleGridSplit;
            func_types[j] = cVTU_INT;
        } else if (strcmp(argv[j], "proc") == 0) {
            funcs[j] = &FixParaview::handleGridProc;
            func_types[j] = cVTU_INT;
        } else if (strcmp(argv[j], "vol") == 0) {
            funcs[j] = &FixParaview::handleGridVol;
            func_types[j] = cVTU_DOUBLE;

        // compute value = c_ID
        // if no trailing [], then index = 0, else index = int between []
        } else if (strncmp(argv[j], "c_", 2) == 0) {
            // Adding pointer to function to handle compute values
            funcs[j]      = &FixParaview::handleGridCompute;
            func_types[j] = cVTU_DOUBLE;

            // copying up until [
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix, &argv[j][2]);
            char *ptr = strchr(suffix,'[');

            // get the int between [] if [ was found, otherwise, set to 0
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']') {
                    ERROR("Invalid attribute in fix paraview grid command");
                }
                argindex[j] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j] = 0;

            // checks
            n = modify->find_compute(suffix);
            if (n < 0) ERROR("Could not find fix paraview grid compute ID");
            if (modify->compute[n]->per_grid_flag == 0)
                ERROR("fix paraview compute does not compute per-grid info");
            if (argindex[j] == 0 && modify->compute[n]->size_per_grid_cols != 0)
                ERROR("fix paraview compute does not calculate per-grid vector");
            if (argindex[j] > 0 && modify->compute[n]->size_per_grid_cols == 0)
                ERROR("fix paraview compute does not calculate per-grid array");
            if (argindex[j] > 0 && argindex[j] > modify->compute[n]->size_per_grid_cols)
                ERROR("fix paraview compute array is accessed out-of-range");

            // remembering the id to call the compute later
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

        // fix value = f_ID
        // if no trailing [], then index = 0, else index = int between []
        // if index = 0 and compute stores array, expand to one value per column
        } else if (strncmp(argv[j],"f_",2) == 0) {
            // Adding pointer to function to handle fix values
            funcs[j]      = &FixParaview::handleGridFix;
            func_types[j] = cVTU_DOUBLE;

            // copying up until [
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix,&argv[j][2]);
            char *ptr = strchr(suffix,'[');

            // get the int between [] if [ was found, otherwise, set to 0
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']') {
                    ERROR("Invalid attribute in fix paraview command");
                }
                argindex[j] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j] = 0;

            // checks
            n = modify->find_fix(suffix);
            if (n < 0) ERROR("Could not find fix paraview fix ID");
            if (modify->fix[n]->per_grid_flag == 0)
                ERROR("fix paraview fix does not compute per-grid info");
            if (argindex[j] == 0 && modify->fix[n]->size_per_grid_cols != 0)
                ERROR("fix paraview fix does not calculate a per-grid vector");
            if (argindex[j] > 0 && modify->fix[n]->size_per_grid_cols == 0)
                ERROR("fix paraview fix does not calculate per-grid array");
            if (argindex[j] > 0 && argindex[j] > modify->fix[n]->size_per_grid_cols)
                ERROR("fix paraview fix array is accessed out-of-range");

            // remembering the id to call the fix later
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

        // variable value = v_name
        } else if (strncmp(argv[j],"v_",2) == 0) {
            // Adding pointer to function to handle variables
            funcs[j] = &FixParaview::handleGridVariable;
            func_types[j] = cVTU_DOUBLE;

            // copy variable name into suffix
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix, &argv[j][2]);

            // set to default value of 0
            argindex[j] = 0;

            // checks
            n = input->variable->find(suffix);
            if (n < 0) ERROR("Could not find fix paraview variable name");
            if (input->variable->grid_style(n) == 0)
                ERROR("fix paraview variable is not grid-style variable");

            // remember the id for later use to call the variable
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

            // tell rest of code that there is a variable
            there_is_variable = true;

        } else {
            ERROR("invalid argument passed to fix paraview");
        }
    }
}

/**
 * handles surf arguments
*/
void FixParaview::handleSurfArgs(int argc, char **argv) {
    // loop over args
    for (int j = 0; j < argc; j++) {
        // remembering arg name as a data name, labels the data in the vtu file
        data_names[j] = new char[strlen(argv[j])+1];
        strcpy(data_names[j], argv[j]);

        // handle the different arguments
        if (strcmp(argv[j], "id") == 0) {
            funcs[j] = &FixParaview::handleSurfId;
            if (sizeof(surfint) == sizeof(smallint)) {
                func_types[j] = cVTU_INT;
            } else {
                func_types[j] = cVTU_BIGINT;
            }
        } else if (strcmp(argv[j], "type") == 0) {
            funcs[j] = &FixParaview::handleSurfType;
            func_types[j] = cVTU_INT;
        } else if (strcmp(argv[j], "area") == 0) {
            funcs[j] = &FixParaview::handleSurfArea;
            func_types[j] = cVTU_DOUBLE;

        // custom surf vector or array
        // if no trailing [], then arg is set to 0, else arg is int between []
        } else if (strncmp(argv[j], "s_", 2) == 0) {
            // Adding pointer to function to handle custom variables
            funcs[j] = &FixParaview::handleSurfCustom;

            // copy custom variable name into suffix
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix,&argv[j][2]);
            char *ptr = strchr(suffix,'[');

            // get the int between [] if [ was found, otherwise, set to 0
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']')
                    ERROR("Invalid attribute in fix paraview command");
                argindex[j] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j] = 0;


            n = surf->find_custom(suffix);

            // checks
            if (n < 0) ERROR("Could not find fix paraview custom attribute");
            
            // converting from other enum type names to the enum type names in this file
            int temp_type = surf->etype[n];
            switch (temp_type)  {
                case _INT:
                    func_types[j] = cVTU_INT;
                    break;
                case _DOUBLE:
                    func_types[j] = cVTU_DOUBLE;
                    break;
                default:
                    ERROR("unsupported type for surf vector")
                    break;
            }

            // checks
            if (argindex[j] == 0 && surf->esize[n] > 0)
                ERROR("fix paraview custom attribute does not store per-surf vector");
            if (argindex[j] > 0 && surf->esize[n] == 0)
                ERROR("fix paraview custom attribute does not store per-surf array");
            if (argindex[j] > 0 && argindex[j] > surf->esize[n])
                ERROR("fix paraview custom attribute is accessed out-of-range");

            // remembering the id to call the custom surface variable later
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

        // compute value = c_ID
        // if no trailing [], then index = 0, else index = int between []
        } else if (strncmp(argv[j], "c_", 2) == 0) {
            // Adding pointer to function to handle compute values
            funcs[j] = &FixParaview::handleSurfCompute;
            func_types[j] = cVTU_DOUBLE;

            // copying up until [
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix,&argv[j][2]);
            char *ptr = strchr(suffix,'[');

            // get the int between [] if [ was found, otherwise, set to 0
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']') ERROR("Invalid attribute in fix paraview command");
                argindex[j] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j] = 0;

            // checks
            n = modify->find_compute(suffix);
            if (n < 0) ERROR("Could not find fix paraview compute ID");
            if (surf->implicit)
                ERROR("Cannot use fix paraview compute with implicit surfs");
            if (modify->compute[n]->per_surf_flag == 0)
                ERROR("fix paraview compute does not compute per-surf info");
            if (argindex[j]== 0 && modify->compute[n]->size_per_surf_cols != 0)
                ERROR("fix paraview compute does not compute per-surf vector");
            if (argindex[j] > 0 && modify->compute[n]->size_per_surf_cols == 0)
                ERROR(
                        "fix paraview compute does not calculate per-surf array");
            if (argindex[j] > 0 &&
                argindex[j] > modify->compute[n]->size_per_surf_cols)
                ERROR("fix paraview compute array is accessed out-of-range");

            // remembering the id to call the compute later
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

        // fix value = f_ID
        // if no trailing [], then index = 0, else index = int between []
        // if index = 0 and compute stores array, expand to one value per column
        } else if (strncmp(argv[j], "f_", 2) == 0) {
            // Adding pointer to function to handle fix values
            funcs[j] = &FixParaview::handleSurfFix;
            func_types[j] = cVTU_DOUBLE;

            // copying up until [
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix,&argv[j][2]);
            char *ptr = strchr(suffix,'[');
        
            // get the int between [] if [ was found, otherwise, set to 0
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']') ERROR("Invalid attribute in fix paraview command");
                argindex[j] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j] = 0;

            // checks
            n = modify->find_fix(suffix);
            if (n < 0) ERROR("Could not find fix paraview fix ID");
            if (surf->implicit)
                ERROR("Cannot use fix paraview fix with implicit surfs");
            if (modify->fix[n]->per_surf_flag == 0)
                ERROR("fix paraview fix does not compute per-surf info");
            if (argindex[j]== 0 && modify->fix[n]->size_per_surf_cols != 0)
                ERROR("fix paraview fix does not compute per-surf vector");
            if (argindex[j] > 0 && modify->fix[n]->size_per_surf_cols == 0)
                ERROR("fix paraview fix does not compute per-surf array");
            if (argindex[j] > 0 && argindex[j] > modify->fix[n]->size_per_surf_cols)
                ERROR("fix paraview fix array is accessed out-of-range");

            // remembering the id to call the fix later
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

        // variable value = v_name
        } else if (strncmp(argv[j],"v_",2) == 0) {
            // Adding pointer to function to handle variables
            funcs[j] = &FixParaview::handleSurfVariable;
            func_types[j] = cVTU_DOUBLE;

            // copy variable name into suffix
            int n = strlen(argv[j]);
            char *suffix = new char[n+1];
            strcpy(suffix,&argv[j][2]);

            // set to default value of 0
            argindex[j] = 0;

            // checks
            n = input->variable->find(suffix);
            if (n < 0) ERROR("Could not find fix paraview variable name");
            if (input->variable->surf_style(n) == 0)
                ERROR("fix paraview variable is not surf-style variable");

            // remember the id for later use to call the variable
            ids[j] = new char[strlen(suffix)+1];
            strcpy(ids[j], suffix);
            delete [] suffix;

            // tell rest of code that there is a variable
            there_is_variable = true;

        } else {
            ERROR( "invalid argument passed to fix paraview");
        }
    }
}

/* ---------------------------------------------------------------------- */

/**
 * cleaning up the class and freeing up memory 
 */
FixParaview::~FixParaview() {
    delete [] prefix;
    delete [] funcs;
    delete [] ids;
    delete [] argindex;
    delete [] func_types;
    for (int i = 0; i < num_expanded_args; i++) delete [] data_names[i];
    delete [] data_names;
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixParaview::setmask() { return 0 | END_OF_STEP; }

/* ---------------------------------------------------------------------- */

/**
 * allocates memory, loads the initial data, and performs some checks
 */
void FixParaview::init() {
    // if should dump at step 0, call end_of_step at init
    if (dump_at_0) end_of_step();
}

/* ---------------------------------------------------------------------- */

/*
 * Calculates base^exponent (in traditional math sense of ^)
*/
index_t power(index_t base, index_t exponent) {
    if (exponent == 0) return 1;
    if (exponent == 1) return base;

    index_t temp = power(base, exponent / 2);
    if (exponent % 2 == 0) return temp * temp;
    else return base * temp * temp;
}

/* ---------------------------------------------------------------------- */

/*
 * Communication algorithm, pairs up processes and builds a "tree"
 * of messages ending at the 0 process
*/
// num_points       = ACTUAL_NUMBER_OF_POINTS/3
// num_cells  = ACTUAL_NUMBER_OF_CONNECTIONS/num_elements_per_cell_num
void FixParaview::collectDataFromProcesses(index_t step, index_t num, double *&points, index_t &num_points, index_t *&connections, index_t &num_cells, index_t num_elements_per_cell_num, void *&data, index_t *&offsets, index_t &offsets_last_index) {
    // if the number is 1, then all data has been collected onto the main process
    // namely, the 0 process
    if (num == 1) return;

    index_t i, j, k, h;
   
    // testing for even or odd, used later
    int offset = num%2;

    // used a lot so calculating it here
    index_t power_val = power(2, step);
    
    // essentially, this loop performs the communication between pairs of processes
    for (i = 0;  i < num-offset; i+=2) {
        // if a process has this id, it sends the info,
        // this process id was choosen careful so that the proper 
        // pairs of process communicate and all of the data ends up on the 0 process
        // BE CAREFUL if you change it, very easy to royally mesh things up
        if (comm->me == (i+1)*power_val) {
            // sending points
            MPI_Send(&num_points,            1, MPI_SPARTA_INDEX_T, i*power_val, 0, world);
            MPI_Send(     points, 3*num_points,         MPI_DOUBLE, i*power_val, 0, world);
            
            // sending connections
            MPI_Send(&num_cells,                                    1, MPI_SPARTA_INDEX_T, i*power_val, 0, world);
            MPI_Send(connections, num_elements_per_cell_num*num_cells, MPI_SPARTA_INDEX_T, i*power_val, 0, world);

            // sending offsets
            MPI_Send(&offsets_last_index,                    1, MPI_SPARTA_INDEX_T, i*power_val, 0, world);
            MPI_Send(            offsets, offsets_last_index+1, MPI_SPARTA_INDEX_T, i*power_val, 0, world);

            // sending data
            MPI_Send(data, offsets[offsets_last_index], MPI_BYTE, i*power_val, 0, world);
        
        // if a process has this id, it receives the info and performs the calculations
        } else if (comm->me == i*power_val) {
            index_t buf_num_points, buf_num_connections, buf_offsets_last_index;

            // receiving points
            MPI_Recv(&buf_num_points, 1, MPI_SPARTA_INDEX_T, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // creating buffer to hold points while the original array is resized
            double *buf_points;
            buf_points = (double*)memory->smalloc(3*num_points*sizeof(double), "fixParaview::iter::buf_points");

            // saving points data to buffer
            memcpy(buf_points, points, 3*num_points*sizeof(double));

            // free then realloc to new size including the size of points from the other process
            memory->sfree(points);
            points = (double*)memory->smalloc(3*(num_points+buf_num_points)*sizeof(double), "fixParaview::iter::points");

            // moving back from buffer
            memcpy(points, buf_points, 3*num_points*sizeof(double));

            // freeing the buffer
            memory->sfree(buf_points);

            // receiving new points starting at the end of the original points array
            // essentially appending them to the original points array
            MPI_Recv(&(points[3*num_points]), 3*buf_num_points, MPI_DOUBLE, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // receiving connections
            MPI_Recv(&buf_num_connections, 1, MPI_SPARTA_INDEX_T, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // creating buffer to hold points while the original array is resized
            index_t *buf_connections;
            buf_connections = (index_t*)memory->smalloc(num_elements_per_cell_num*num_cells*sizeof(index_t), "fixParaview::iter::buf_connections");

            // saving connections data to buffer
            memcpy(buf_connections, connections, num_elements_per_cell_num*num_cells*sizeof(index_t));

            // free then realloc to new size including the size of points from the other process
            memory->sfree(connections);
            connections = (index_t*)memory->smalloc(num_elements_per_cell_num*(num_cells+buf_num_connections)*sizeof(index_t), "fixParaview::iter::connections");

            // move the data buffer back to the connections
            memcpy(connections, buf_connections, num_elements_per_cell_num*num_cells*sizeof(index_t));

            // freeing the buffer
            memory->sfree(buf_connections);

            // receiving new points starting at the end of the original connections array
            // essentially appending them to the original connections array
            MPI_Recv(&(connections[num_elements_per_cell_num*num_cells]), num_elements_per_cell_num*buf_num_connections, MPI_SPARTA_INDEX_T, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // starting from the end of the original data points
            for (j = num_cells; j < num_cells+buf_num_connections; j++) {
                // shifting the indicies up to account for new data
                // NOTE: this shift is 1/3 the actual shift, this is on purpose
                // later on there will be multiplication by 3 for assigning values
                for (k = 0; k < num_elements_per_cell_num; k++) connections[num_elements_per_cell_num*j+k]+=num_points;
            }

            // incrementing the sizes to account for the new data
            num_points     +=buf_num_points;
            num_cells+=buf_num_connections;

            // checks
            for (j = 0; j < num_cells; j++) {
                for (k = 0; k < num_elements_per_cell_num; k++) {
                    if (connections[num_elements_per_cell_num*j+k] >= num_points || connections[num_elements_per_cell_num*j+k] < 0) {
                        PRINT_PROC("%d : index = %ld\n", comm->me, connections[num_elements_per_cell_num*j+k]);
                        ERROR("invalid index found in connections, exceeds valid index bounds");
                    }
                }
            }
            for (j = num_cells-buf_num_connections; j < num_cells; j++) {
                for (k = 0; k < num_elements_per_cell_num; k++) {
                    if (connections[num_elements_per_cell_num*j+k] - (num_points-buf_num_points) < 0)
                        ERROR("invalid index found in connections, added index was not shifted of added data");
                }
            }

            // filtering the points
            std::unordered_map<std::string, index_t> points_map;
            std::vector<std::string> hashes(num_points);

            // hashing the points, remembering hashes to reduce cost and replace a permutation table
            for (j = 0; j < num_points; j++) {
                hashes[j] = util::hashDoubleArray(&points[3*j], 3);
                points_map[hashes[j]] = j;
            }

            // getting the number of unique points after the hashing
            num_points = (index_t)points_map.size();

            // initializing storage array for the new/unique points
            buf_points = (double*)memory->smalloc(3*num_points*sizeof(double), "fixParaview::iter::buf_points");

            // iterating over the map of unique points and adding the unique
            // points to the buffer, also setting the value in the map
            // to the new index in the buffer array
            j = 0;
            for (auto& [key, value] : points_map) {
                for (k = 0; k < 3; k++) buf_points[3*j+k] = points[3*value+k];
                value = j;
                j++;
            }

            // shifting the indicies of the connections from the indicies in points to the indicies
            // in the buffer
            for (j = 0; j < num_cells; j++) {
                for (k = 0; k < num_elements_per_cell_num; k++)
                    connections[num_elements_per_cell_num*j+k] = points_map[hashes[connections[num_elements_per_cell_num*j+k]]];
            }

            // free up memory
            points_map.clear(); hashes.clear();

            // resizing points to the same size as new_points
            memory->sfree(points);
            points = (double*)memory->smalloc(3*num_points*sizeof(double), "fixParaview::iter::points");

            // transferring buffer into points
            memcpy(points, buf_points, 3*num_points*sizeof(double));
            
            // cleaning up
            memory->sfree(buf_points);

            /*
                data handling from now on 
            */

            // receive last index from other process            
            MPI_Recv(&buf_offsets_last_index, 1, MPI_SPARTA_INDEX_T, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // creating buffer to hold offsets while the original array is resized
            if (buf_offsets_last_index != offsets_last_index)
                ERROR("process have different output parameters, this is not allowed");

            // initialize a buffer for the data
            void* buf_data;
            buf_data = memory->smalloc(offsets[offsets_last_index]*sizeof(char), "fixParaview::iter::buf_data"); // char has the size of a byte, thats why sizeof(char)

            // copy into buffer
            memcpy(buf_data, data, offsets[offsets_last_index]*sizeof(char));

            // free to resize later
            memory->sfree(data);

            // create buffer for offsets from other process
            index_t *buf_offsets;
            buf_offsets = (index_t*)memory->smalloc((offsets_last_index+1)*sizeof(index_t), "fixParaview::iter::buf_offsets");

            // receiving offsets from other process
            MPI_Recv(buf_offsets, offsets_last_index+1, MPI_SPARTA_INDEX_T, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);
            
            // realloc data
            data = memory->smalloc((offsets[offsets_last_index]+buf_offsets[offsets_last_index])*sizeof(char), "fixParaview::iter::data");

            // copying data from the buffer back into the original, this essentially treats every block of data as
            // consisting of two pieces, the original data is added to the first piece
            for (j = 0; j < offsets_last_index; j++)
                memcpy(&(((char*)data)[offsets[j]+buf_offsets[j]]), &(((char*)buf_data)[offsets[j]]), offsets[j+1]-offsets[j]);

            // resize buf for new data
            memory->sfree(buf_data);
            buf_data = memory->smalloc(buf_offsets[offsets_last_index]*sizeof(char), "fixParaview::iter::buf_data");

            // receiving new data into end of resized data
            MPI_Recv(buf_data, buf_offsets[offsets_last_index], MPI_BYTE, (i+1)*power_val, 0, world, MPI_STATUS_IGNORE);

            // copying data from the buffer back into the original, this essentially treats every block of data as
            // consisting of two pieces, the original data is added to the second piece
            for (j = 0; j < offsets_last_index; j++)
                memcpy(&(((char*)data)[offsets[j+1]+buf_offsets[j]]), &(((char*)buf_data)[buf_offsets[j]]), buf_offsets[j+1]-buf_offsets[j]);

            // free up no longer needed data
            memory->sfree(buf_data);

            // account for the new offsets
            for (j = 0; j < offsets_last_index+1; j++)
                offsets[j] += buf_offsets[j];
            
            // free up no longer needed data
            memory->sfree(buf_offsets);
        }
    }
    MPI_Barrier(world);

    // offset can be 0 or 1. offset being 1 is the weird case, there are
    // an odd number of processes, adding offset is critical
    num = num/2 + offset; 

    // performing the recursion
    collectDataFromProcesses(step+1, num, points, num_points, connections, num_cells, num_elements_per_cell_num, data, offsets, offsets_last_index);
}

void FixParaview::pointSwitch(index_t _index, index_t _element_index, double _arr[3]) {
    // some friendly macros to avoid typing too much
    #define XLO(_index) grid->cells[_index].lo[0];
    #define YLO(_index) grid->cells[_index].lo[1];
    #define ZLO(_index) grid->cells[_index].lo[2];
    #define XHI(_index) grid->cells[_index].hi[0];
    #define YHI(_index) grid->cells[_index].hi[1];
    #define ZHI(_index) grid->cells[_index].hi[2];
    #define LINE_V1X(_index) lines[_index].p1[0]
    #define LINE_V1Y(_index) lines[_index].p1[1]
    #define LINE_V2X(_index) lines[_index].p2[0]
    #define LINE_V2Y(_index) lines[_index].p2[1]
    #define TRI_V1X(_index) tris[_index].p1[0]
    #define TRI_V1Y(_index) tris[_index].p1[1]
    #define TRI_V1Z(_index) tris[_index].p1[2]
    #define TRI_V2X(_index) tris[_index].p2[0]
    #define TRI_V2Y(_index) tris[_index].p2[1]
    #define TRI_V2Z(_index) tris[_index].p2[2]
    #define TRI_V3X(_index) tris[_index].p3[0]
    #define TRI_V3Y(_index) tris[_index].p3[1]
    #define TRI_V3Z(_index) tris[_index].p3[2]
    
    // do not change the order below, they are added in this order
    // to correspond with how vtk orders its elements
    // this switch is big but it does not hurt performance because with modern
    // compilers switch statements are constant time (always take same amount
    // of time to evaluate), unlike if statements
    switch (_element_index) {
        case 0:
            _arr[0] = XLO(_index);
            _arr[1] = YLO(_index);
            _arr[2] = 0.0;
            break;
        case 1:
            _arr[0] = XLO(_index);
            _arr[1] = YHI(_index);
            _arr[2] = 0.0;
            break;
        case 2:
            _arr[0] = XHI(_index);
            _arr[1] = YHI(_index);
            _arr[2] = 0.0;
            break;
        case 3:
            _arr[0] = XHI(_index);
            _arr[1] = YLO(_index);
            _arr[2] = 0.0;
            break;
        case 4:
            _arr[0] = XLO(_index);
            _arr[1] = YLO(_index);
            _arr[2] = ZLO(_index);
            break;
        case 5:
            _arr[0] = XLO(_index);
            _arr[1] = YHI(_index);
            _arr[2] = ZLO(_index);
            break;
        case 6:
            _arr[0] = XLO(_index);
            _arr[1] = YHI(_index);
            _arr[2] = ZHI(_index);
            break;
        case 7:
            _arr[0] = XLO(_index);
            _arr[1] = YLO(_index);
            _arr[2] = ZHI(_index);
            break;
        case 8:
            _arr[0] = XHI(_index);
            _arr[1] = YLO(_index);
            _arr[2] = ZLO(_index);
            break;
        case 9:
            _arr[0] = XHI(_index);
            _arr[1] = YHI(_index);
            _arr[2] = ZLO(_index);
            break;
        case 10:
            _arr[0] = XHI(_index);
            _arr[1] = YHI(_index);
            _arr[2] = ZHI(_index);
            break;
        case 11:
            _arr[0] = XHI(_index);
            _arr[1] = YLO(_index);
            _arr[2] = ZHI(_index);
            break;
        case 12:
            _arr[0] = LINE_V1X(_index);
            _arr[1] = LINE_V1Y(_index);
            _arr[2] = 0.0;
            break;
        case 13:
            _arr[0] = LINE_V2X(_index);
            _arr[1] = LINE_V2Y(_index);
            _arr[2] = 0.0;
            break;
        case 14:
            _arr[0] = TRI_V1X(_index);
            _arr[1] = TRI_V1Y(_index);
            _arr[2] = TRI_V1Z(_index);
            break;
        case 15:
            _arr[0] = TRI_V2X(_index);
            _arr[1] = TRI_V2Y(_index);
            _arr[2] = TRI_V2Z(_index);
            break;
        case 16:
            _arr[0] = TRI_V3X(_index);
            _arr[1] = TRI_V3Y(_index);
            _arr[2] = TRI_V3Z(_index);
            break;
        default:
            ERROR("invalid value encountered in switch"); break;
    }

    // cleaning up defines
    #undef XLO
    #undef YLO
    #undef ZLO
    #undef XHI
    #undef YHI
    #undef ZHI
    #undef LINE_V1X
    #undef LINE_V1Y
    #undef LINE_V2X
    #undef LINE_V2Y
    #undef TRI_V1X
    #undef TRI_V1Y
    #undef TRI_V1Z
    #undef TRI_V2X
    #undef TRI_V2Y
    #undef TRI_V2Z
    #undef TRI_V3X
    #undef TRI_V3Y
    #undef TRI_V3Z
}

/**
 * end of step
 */
void FixParaview::end_of_step() {
    std::unordered_map<std::string, index_t> points_map;
    double point[3], *points;
    index_t i,j,k, remainder, index, *connections, *data_offsets, points_map_size;
    void* data; // holds per cell data

    // holds number of nodes for cell type, 3 for tri, 2 for line, 4 for 2d grid cell,
    // 8 for 3d grid cell
    unsigned num_nodes_per_cell;

    // offset for switch in pointSwitch
    unsigned cell_count_offset;

    // handles different styles of output
    if (style == GRID) {
        // handling grid, so set to the number of grid cells I own
        total_num_cells    = grid->nlocal;
        
        // 4 for 2d, 8 for 3d
        num_nodes_per_cell = 4*(domain->dimension - 1);

        // this is special for the switch in pointSwitch
        // 0 for 2d, 4 for 3d
        cell_count_offset  = 4*(domain->dimension - 2);
    } else {
        // handling surf, so set to the number of surf cells I own
        total_num_cells    = surf->nown;

        // 2 for 2d, 3 for 3d
        num_nodes_per_cell = domain->dimension;

        // this is special for the switch in pointSwitch
        // 12 for 2d, 14 for 3d
        cell_count_offset  = 2*(domain->dimension - 2)+12; // +12 shift account for grid stuff in switch structure
    }

    // number of cells owned by process
    num_cells = 0;

    // handles different styles of output
    if (style == GRID) {
        // getting cells that belong to the specified group
        for (i = 0; i < total_num_cells; i++) {
            if (!(grid->cinfo[i].mask & groupbit)) continue;
            // if (grid->cells[i].nsplit > 1) continue; // i do not know what this does
            num_cells++;
        }
    } else {
        // temporary index variable
        index_t m;

        // if 2d do lines
        if (domain->dimension == 2) {
            // handle different type of surf
            if (surf->distributed && !surf->implicit)
                lines = surf->mylines;
            else
                lines = surf->lines;
            
            // loop over cells, calculate the index and 
            // getting cells that belong to the specified group
            for (i = 0; i < total_num_cells; i++) {
                if (!surf->distributed)
                    m = (unsigned)comm->me + i*((unsigned)comm->nprocs);
                else
                    m = i;

                // checking to make sure it is in the allowed group
                // if so, increment count
                if (lines[m].mask & groupbit)
                    num_cells++;
            }
        // if 3d do tris
        } else {
            // handle different type of surf
            if (surf->distributed && !surf->implicit)
                tris = surf->mytris;
            else
                tris = surf->tris;
            
            // loop over cells, calculate the index and 
            // getting cells that belong to the specified group
            for (i = 0; i < total_num_cells; i++) {
                if (!surf->distributed)
                    m = (unsigned)comm->me + i*((unsigned)comm->nprocs);
                else
                    m = i;
                
                // checking to make sure it is in the allowed group
                // if so, increment count
                if (tris[m].mask & groupbit)
                    num_cells++;
            }
        }
    }

    // maps the normal indicies, i.e. 0,1,... to the indicies of the valid cells in grid->cells or lines/tris
    cell_index_map = (index_t*)memory->smalloc(num_cells*sizeof(index_t), "fixParaview/end_of_step/cell_index_map");

    // set for counter
    j = 0;

    // handles different styles of output
    if (style == GRID) {
        // getting cells that below to the specified group
        for (i = 0; i < total_num_cells; i++) {
            if (!(grid->cinfo[i].mask & groupbit)) continue;
            // if (grid->cells[i].nsplit > 1) continue;
            cell_index_map[j] = i;
            j++;
        }
    } else {
        // temporary index variable
        index_t m;

        // if 2d do lines
        if (domain->dimension == 2) {        
            // loop over cells, calculate the index and 
            // getting cells that belong to the specified group    
            for (i = 0; i < total_num_cells; i++) {
                if (!surf->distributed)
                    m = (unsigned)comm->me + i*((unsigned)comm->nprocs);
                else
                    m = i;
                
                // checking to make sure it is in the allowed group
                // if so, add it to the map
                if (!(lines[m].mask & groupbit)) continue;
                cell_index_map[j] = m;
                j++;
            }
        // if 3d do tris
        } else {
            // loop over cells, calculate the index and 
            // getting cells that belong to the specified group
            for (i = 0; i < total_num_cells; i++) {
                if (!surf->distributed)
                    m = (unsigned)comm->me + i*((unsigned)comm->nprocs);
                else
                    m = i;
                
                // checking to make sure it is in the allowed group
                // if so, add it to the map
                if (!(tris[m].mask & groupbit)) continue;
                cell_index_map[j] = m;
                j++;
            }
        }
    }

    // quick check
    if (j != num_cells) { ERROR("number of cells changed will generating the index map"); }

    // hashing all of the points of each cell, this removes duplicate points
    // first, iterating over all of the cells
    for (i = 0; i < num_cells; i++) {
        // looping over each vertex, pointSwitch yields the double[3] of the
        // point at the desired index, this index corresponds the vtk ordering
        // system
        for (j = 0; j < num_nodes_per_cell; j++) {
            pointSwitch(cell_index_map[i], j+cell_count_offset, point);
            // adding hashed point all with a unique long that indicates
            // both its index in the cells array and its location in the cell
            // to be used with pointSwitch
            points_map[util::hashDoubleArray(point, 3)] = num_nodes_per_cell*i+j;
        }
    }
    
    // allocating the number of unique points found in the hashing process
    num_points = (index_t)points_map.size();
    points = (double*)memory->smalloc(3*num_points*sizeof(double), "fixParaview/end_of_step/points");

    i = 0; // used as an iteration counter in the loop below

    // iterating over the filtered points, putting them into the points array,
    // and setting the value
    for (auto& [key, value] : points_map) {
        // used to figure which location in the box this point is,
        // the pointSwitch macro contains the ordering of the points
        remainder = value%num_nodes_per_cell;

        // index in the array of grid cells
        index = (value-remainder)/num_nodes_per_cell;

        // getting the point at the index with the location in cell
        // specified with the remainder
        pointSwitch(index, remainder+cell_count_offset, point);

        // assigning the point to the array
        memcpy(&(points[3*i]), point, 3*sizeof(double));

        // setting the value to the index of the start of the values in the points
        // array, to get all of the points, start here and get the next indicies as well
        value = i;
        i++;
    }
    
    // allocating the same number of connections as grid cells
    connections = (index_t*)memory->smalloc(num_nodes_per_cell*num_cells*sizeof(index_t), "fixParaview/end_of_step/connections");

    // iterating over the grid points
    for (i = 0; i < num_cells; i++) {

        // looping num_elements_per_cell_num times to get the num_elements_per_cell_num points of each cell
        for (j = 0; j < num_nodes_per_cell; j++) {
            // getting the point from the cell at the index, i, and location
            // in the cell, j, saving into point
            pointSwitch(cell_index_map[i], j+cell_count_offset, point);

            // getting the index (actually 1/3 of it) of this point in the points array
            // from the map using the points hash
            connections[num_nodes_per_cell*i+j] = points_map[util::hashDoubleArray(point, 3)];
        }
    }
    // getting the size for error checking later
    points_map_size = (index_t)points_map.size();
    
    // no longer needed, so free up the memory
    points_map.clear();

    // create offsets array for remembering locations in a data
    data_offsets = (index_t*)memory->smalloc((num_expanded_args+1)*sizeof(index_t), "fixParaview/end_of_step/offsets");
    data_offsets[0] = 0; // set first offset to zero because data starts at 0 index

    // setup offsets
    for (i = 0; i < num_expanded_args; i++) {
        // set next offset to the size of the data at the previous offset
        data_offsets[i+1] = num_cells*cVTU_getVTKTypeSizeFromEnum(func_types[i])+data_offsets[i];
    }
    
    // create data array
    data = memory->smalloc(data_offsets[num_expanded_args], "fixParaview/end_of_step/data");

    // create buffer for variables to use, if there are any
    if (there_is_variable)
        vbuf = (double*)memory->smalloc(total_num_cells*sizeof(double), "fixParaview/end_of_step/vbuf");
    
    // loop over functions, getting all the data
    for (i = 0; i < num_expanded_args; i++) {
        (this->*funcs[i])(i, data_offsets, data);
    }

    // free up variable buffer if it was made
    if (there_is_variable)
        memory->sfree(vbuf);
    
    MPI_Barrier(world);
    
    // if we should collect data from all processes into one large file
    if (collect_into_one_file) {
        // PRINT("Collecting Into One File\n");
        // this function is quite complex, it essentially sends the point and connection data
        // from each process to process 0 by incrementally collecting data from pairs of processes
        collectDataFromProcesses(0, comm->nprocs, points, num_points, connections, num_cells, num_nodes_per_cell, data, data_offsets, num_expanded_args);

        // some checks for my sake
        if (comm->me == 0 && comm->nprocs > 1) {
            if (num_points == points_map_size)
                ERROR("did not update points");
            if (num_cells  == total_num_cells)
                ERROR("did not update connections");
        }

        // main process writes the file
        if (comm->me == 0) {
            // formatting file name
            char* filename; int result;
            result = asprintf(&filename, "%s%cFixID%s_step_t%ld.vtu", prefix, SEP, id, update->ntimestep);
            if (result == -1) { ERROR("asprintf failed for vtu file"); }
            
            // create cell_type array
            cVTU_cell_type_t *cell_types = (cVTU_cell_type_t*)memory->smalloc(num_cells*sizeof(cVTU_cell_type_t), "fixParaview/end_of_step/cell_types");

            // temporary cell_type storage
            int per_cell_type;

            // handles different style of cell
            if (style == GRID) {
                // if 3d, hexahedron, if 2d, quadrilateral
                per_cell_type = (domain->dimension == 3) ? cVTU_HEXAHEDRON : cVTU_QUAD;
            } else {
                // if 3d, triangle, if 2d, line
                per_cell_type = (domain->dimension == 3) ? cVTU_TRIANGLE : cVTU_LINE;
            }

            // setting the cell types
            for (i = 0; i < num_cells; i++)
                cell_types[i] = per_cell_type;

            // writing file
            cVTU_writeUnstructuredGridAppended(
                filename, // file to write to
                points, num_points, connections, num_cells, cell_types, // defines the geometry
                NULL, 0, NULL, NULL, // point data is nothing
                data, num_expanded_args, func_types, (const char**)data_names // cell data
            );

            // detects if an error occurred in the vtu writer, and returns the error message
            char* msg = cVTU_checkForErrors();
            if (msg != NULL) {
                ERROR(msg);
            }

            // no longer needed
            memory->sfree(cell_types);
        }
    // if we are not collecting, each process writes its own file
    } else {
        // formatting file name
        char* filename; int result;

        // main process makes directory
        if (comm->me == 0) {
            // change directory name based on style
            result = asprintf(&filename, "%s%cFixID%s_step_t%ld", prefix, SEP, id, update->ntimestep);
            if (result == -1) { ERROR("asprintf failed for step directory"); }

            struct stat sb;
            // if the step data directory does not exist
            if (!(stat(filename, &sb) == 0 && S_ISDIR(sb.st_mode))) {
                // Create directory with permissions 0755 (drwxr-xr-x)
                result = mkdir(filename, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
                // check if directory was made
                if (result != 0) { ERROR("could not make step directory"); }
            }
        }

        // ensures directory is made first
        MPI_Barrier(world);

        // change filename based on style
        result = asprintf(&filename, "%s%cFixID%s_step_t%ld%crank_%d.vtu", prefix, SEP, id, update->ntimestep, SEP, comm->me);
        if (result == -1) { ERROR("asprintf failed for per rank vtu file"); }
        
        // create cell_type array
        cVTU_cell_type_t *cell_types = (cVTU_cell_type_t*)memory->smalloc(num_cells*sizeof(cVTU_cell_type_t), "fixParaview/end_of_step/cell_types");

        // temporary cell_type storage
        int per_cell_type;

        // handles different style of cell
        if (style == GRID) {
            // if 3d, hexahedron, if 2d, quadrilateral
            per_cell_type = (domain->dimension == 3) ? cVTU_HEXAHEDRON : cVTU_QUAD;
        } else {
            // if 3d, triangle, if 2d, line
            per_cell_type = (domain->dimension == 3) ? cVTU_TRIANGLE : cVTU_LINE;
        }

        // setting the cell types
        for (i = 0; i < num_cells; i++)
            cell_types[i] = per_cell_type;

        // writing file
        cVTU_writeUnstructuredGridAppended(
            filename, // file to write to
            points, num_points, connections, num_cells, cell_types, // defines the geometry
            NULL, 0, NULL, NULL, // point data is nothing
            data, num_expanded_args, func_types, (const char**)data_names // cell data
        );

        // detects if an error occurred in the vtu writer, and returns the error message
        char* msg = cVTU_checkForErrors();
        if (msg != NULL) {
            ERROR(msg);
        }

        // no longer needed
        memory->sfree(cell_types);

        // only the 0 process writes the parallel vtu head file
        if (comm->me == 0) {
            // generate the filenames based on each process
            char **process_file_names = new char*[comm->nprocs];
            for (i = 0; i < comm->nprocs; i++) {
                // change filename based on style
                result = asprintf(&filename, "FixID%s_step_t%ld%crank_%d.vtu", id, update->ntimestep, SEP, i);
                if (result == -1) { ERROR("asprintf failed for per rank vtu file"); }
                
                // copy file name into process filenames
                process_file_names[i] = new char[strlen(filename)+1];
                strcpy(process_file_names[i], filename);
            }

            // change filename based on style
            result = asprintf(&filename, "%s%cFixID%s_step_t%ld.pvtu", prefix, SEP, id, update->ntimestep);
            if (result == -1) { ERROR("asprintf failed for pvtu file"); }

            // write pvtu file
            cVTU_writeParallelUnstructuredGridHeadFile(
                filename, (const char**)process_file_names, comm->nprocs,
                0, NULL, NULL,
                num_expanded_args, func_types, (const char**)data_names
            );

            // detects if an error occurred in the vtu writer, and returns the error message
            char* msg = cVTU_checkForErrors();
            if (msg != NULL) { ERROR(msg); }

            // clean up
            for (i = 0; i < comm->nprocs; i++) delete [] process_file_names[i];
            delete [] process_file_names;
        }
    }
    MPI_Barrier(world);

    // free allocated variables
    memory->sfree(points);
    memory->sfree(connections);
    memory->sfree(data_offsets);
    memory->sfree(data);
    memory->sfree(cell_index_map);
}

// macro for later use
#define COPY_VAL_INTO_DATA memcpy(&((char*)data)[offsets[index_in_offsets]+i*sizeof(val)], &val, sizeof(val));

/**
 * Handles per grid cell id
*/
void FixParaview::handleGridId(   int index_in_offsets, index_t *offsets, void* data) {
    cellint val;
    // iterate over valid cells and get their id
    for (index_t i = 0; i < total_num_cells; i++) {
        val = grid->cells[cell_index_map[i]].id;
        COPY_VAL_INTO_DATA
    }
}

void FixParaview::handleGridSplit(int index_in_offsets, index_t *offsets, void* data) {
    int val;
    for (int i = 0; i < total_num_cells; i++) {
        val = -(grid->cells[cell_index_map[i]].nsplit) + 1;
        COPY_VAL_INTO_DATA
    }
}

void FixParaview::handleGridProc( int index_in_offsets, index_t *offsets, void* data) {
    int val;
    for (int i = 0; i < total_num_cells; i++) {
        val = grid->cells[cell_index_map[i]].proc;
        COPY_VAL_INTO_DATA

    }
}

void FixParaview::handleGridVol( int index_in_offsets, index_t *offsets, void* data) {
    double val;
    for (int i = 0; i < total_num_cells; i++) {
        val = grid->cinfo[cell_index_map[i]].volume;
        COPY_VAL_INTO_DATA
    }
}

void FixParaview::handleGridCompute( int index_in_offsets, index_t *offsets, void* data) {
    int index = argindex[index_in_offsets];

    // finding the compute, if it exists
    int icompute = modify->find_compute(ids[index_in_offsets]);
    if (icompute < 0) { ERROR("Could not find fix paraview compute ID"); }

    Compute *c = modify->compute[icompute];// compute[field2index[index_in_offsets]];

    if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
    }

    // if one of post_process flags is set,
    //   invoke post_process_grid() or invoke post_process_tally()
    // else extract from compute's vector_grid and array_grid directly
    // dump buf only stores values for grid cells with particles
    //   use cpart indices to extract needed subset

    if (c->post_process_grid_flag)
        c->post_process_grid(index,1,NULL,NULL,NULL,1);
    else if (c->post_process_isurf_grid_flag)
        c->post_process_isurf_grid();

    index_t i;
    double val;
    if (index == 0 || c->post_process_grid_flag) {
        double *vector = c->vector_grid;
        for (i = 0; i < num_cells; i++) {
            val = vector[cell_index_map[i]];
            COPY_VAL_INTO_DATA
        }
    } else {
        index--;
        double **array = c->array_grid;
        for (i = 0; i < num_cells; i++) {
            val = array[cell_index_map[i]][index];
            COPY_VAL_INTO_DATA
        }
    }
}

void FixParaview::handleGridFix(     int index_in_offsets, index_t *offsets, void* data) {
    int index = argindex[index_in_offsets];
    
    // finding the fix, if it exists
    int ifix = modify->find_fix(ids[index_in_offsets]);
    if (ifix < 0) { ERROR("Could not find fix paraview fix ID"); }

    // quick check for compatable dump times
    if (nevery % modify->fix[ifix]->per_grid_freq)
      ERROR("fix paraview and fix not computed at compatible times");

    // getting the fix
    Fix *f = modify->fix[ifix];

    index_t i;
    double val;
    if (index == 0) {
        for (i = 0; i < num_cells; i++) {
            val = f->vector_grid[cell_index_map[i]];
            COPY_VAL_INTO_DATA
        }
    } else {
        index--;
        for (i = 0; i < num_cells; i++) {
            val = f->array_grid[cell_index_map[i]][index];
            COPY_VAL_INTO_DATA
        }
    }
}

void FixParaview::handleGridVariable(int index_in_offsets, index_t *offsets, void* data) {
    // getting the variable
    int ivariable = input->variable->find(ids[index_in_offsets]);
    if (ivariable < 0) ERROR("Could not find fix paraview variable name");

    // get the variable data into vbuf 
    input->variable->compute_grid(ivariable,vbuf,1,0);

    double val;
    for (index_t i = 0; i < num_cells; i++) {
        val = vbuf[cell_index_map[i]];
        COPY_VAL_INTO_DATA
    }
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void FixParaview::handleSurfCompute(int index_in_offsets, index_t* offsets, void* data) {
    int index = argindex[index_in_offsets];

    // finding the compute, if it exists
    int icompute = modify->find_compute(ids[index_in_offsets]);
    if (icompute < 0) { ERROR("Could not find fix paraview compute ID"); }

    Compute *c = modify->compute[icompute];

    if (!(c->invoked_flag & INVOKED_PER_SURF)) {
        c->compute_per_surf();
        c->invoked_flag |= INVOKED_PER_SURF;
    }

    c->post_process_surf();

    index_t i;
    double val;
    if (index == 0 || c->post_process_grid_flag) {
        double *vector = c->vector_surf;
        for (i = 0; i < num_cells; i++) {
            val = vector[cell_index_map[i]];
            COPY_VAL_INTO_DATA
        }
    } else {
        index--;
        double **array = c->array_surf;
        for (i = 0; i < num_cells; i++) {
            val = array[cell_index_map[i]][index];
            COPY_VAL_INTO_DATA
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixParaview::handleSurfFix(int index_in_offsets, index_t* offsets, void* data) {
    int index = argindex[index_in_offsets];
    
    // finding the fix, if it exists
    int ifix = modify->find_fix(ids[index_in_offsets]);
    if (ifix < 0) { ERROR("Could not find fix paraview fix ID"); }

    // quick check for compatable dump times
    if (nevery % modify->fix[ifix]->per_surf_freq)
      ERROR("fix paraview and fix not computed at compatible times");

    // getting the fix
    Fix *f = modify->fix[ifix];

    index_t i;
    double val;
    if (index == 0) {
        for (i = 0; i < num_cells; i++) {
            val = f->vector_surf[cell_index_map[i]];
            COPY_VAL_INTO_DATA
        }
    } else {
        index--;
        for (i = 0; i < num_cells; i++) {
            val = f->array_surf[cell_index_map[i]][index];
            COPY_VAL_INTO_DATA
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixParaview::handleSurfVariable(int index_in_offsets, index_t* offsets, void* data) {
    // getting the variable
    int ivariable = input->variable->find(ids[index_in_offsets]);
    if (ivariable < 0) ERROR("Could not find fix paraview variable name");

    // get the variable data into vbuf 
    input->variable->compute_surf(ivariable,vbuf,1,0);

    double val;
    for (index_t i = 0; i < num_cells; i++) {
        val = vbuf[cell_index_map[i]];
        COPY_VAL_INTO_DATA
    }
}

/* ----------------------------------------------------------------------
   extraction of custom surf attribute
------------------------------------------------------------------------- */

void FixParaview::handleSurfCustom(int index_in_offsets, index_t* offsets, void* data) {
    int index = surf->find_custom(ids[index_in_offsets]);
    if (index < 0) ERROR("Could not find fix paraview custom attribute");

  // for now, custom data only allowed for explicit all
  // so custom data is nlocal in length, not nown
  // when enable distributed, commented out lines replace 2 previous ones
    index_t i;
    // int m;
    if (surf->etype[index] == cVTU_INT) {
        int val;
        if (surf->esize[index] == 0) {
            int *vector = surf->eivec[surf->ewhich[index]];
            for (int i = 0; i < num_cells; i++) {
                // m = comm->me + i*comm->nprocs;
                val = vector[cell_index_map[i]];
                COPY_VAL_INTO_DATA
                //buf[n] = vector[clocal[i]];
            }
        } else {
            int icol = argindex[index_in_offsets]-1;
            int **array = surf->eiarray[surf->ewhich[index]];
            for (i = 0; i < num_cells; i++) {
                val = array[cell_index_map[i]][icol];
                COPY_VAL_INTO_DATA
                //buf[n] = array[clocal[i]][icol];
            }
        }
    } else {
        double val;
        if (surf->esize[index] == 0) {
            double *vector = surf->edvec[surf->ewhich[index]];
            for (i = 0; i < num_cells; i++) {
                val = vector[cell_index_map[i]];
                //buf[n] = vector[clocal[i]];
                COPY_VAL_INTO_DATA
            }
        } else {
            int icol = argindex[index_in_offsets]-1;
            double **array = surf->edarray[surf->ewhich[index]];
            for (i = 0; i < num_cells; i++) {
                val = array[cell_index_map[i]][icol];
                COPY_VAL_INTO_DATA
                //buf[n] = array[clocal[i]][icol];
            }
        }
    }
}

/* ----------------------------------------------------------------------
   one method for every attribute fix paraview can output
   the surf property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void FixParaview::handleSurfId(int index_in_offsets, index_t* offsets, void* data) {
    index_t i; surfint val;
    // NOTE: surfint (bigint) won't fit in double in some cases
    if (domain->dimension == 2) {
        for (i = 0; i < num_cells; i++) {
            val = lines[cell_index_map[i]].id;
            COPY_VAL_INTO_DATA
        }
    } else {
        for (i = 0; i < num_cells; i++) {
            val = tris[cell_index_map[i]].id;
            COPY_VAL_INTO_DATA
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixParaview::handleSurfType(int index_in_offsets, index_t* offsets, void* data) {
    index_t i;
    int val;
    if (domain->dimension == 2) {
        for (i = 0; i < num_cells; i++) {
            val = lines[cell_index_map[i]].type;
            COPY_VAL_INTO_DATA
        }
    } else if (domain->dimension == 3) {
        for (i = 0; i < num_cells; i++) {
            val = tris[cell_index_map[i]].type;
            COPY_VAL_INTO_DATA
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixParaview::handleSurfArea(int index_in_offsets, index_t* offsets, void* data) {
    index_t i;
    double val;
    if (domain->dimension == 2) {
        if (domain->axisymmetric) {
            for (i = 0; i < num_cells; i++) {
                val = surf->axi_line_size(&lines[cell_index_map[i]]);
                COPY_VAL_INTO_DATA;
            }
        } else {
            for (i = 0; i < num_cells; i++) {
                val = surf->line_size(&lines[cell_index_map[i]]);
                COPY_VAL_INTO_DATA;
            }
        }
    } else if (domain->dimension == 3) {
        double tmp;
        for (i = 0; i < num_cells; i++) {
            val = surf->tri_size(&tris[cell_index_map[i]],tmp); // tmp does not do anything here
            COPY_VAL_INTO_DATA
        }
    }
}