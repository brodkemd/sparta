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

#include "fix_fea.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};  // several files

/* ---------------------------------------------------------------------- */

/**
 * inputs, format INDEX : DESCRIPTION or VALUE
 * 0  : id
 * 1  : fea
 * 2  : nevery, execute every this many time steps
 * 3  : elmer_exe, path to the executable for elmer
 * 4  : sif_file, path the sif file to run with elmer
 * 5  : surf_file, path or name of file to dump surface to
 * 6  : compute-id, id for compute 
 * 7  : group-id, group ID for which surface elements to perform calculation on
 * 8  : mix-ID, mixture ID for particles to perform calculation on
 * 9  : dump-id, 
 * 10 : dump-file, file to dump surf info to
 * 11 : select-ID, what surface elements to dump
 * 12 : dump compute format, this is c_[*]
 */
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    // checking if the number of args is correct
    if (narg > 13)
        error->all(FLERR,"Illegal fix fea command, too many inputs, format \"fix id fea EXE FILE\"");
    else if (narg < 13)
        error->all(FLERR,"Illegal fix fea command, too few inputs, format \"fix id fea EXE FILE\"");

    // making sure there is a surface to analyze
    if (!surf->exist)
        error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // getting when to run from the args to this fix command, nevery
    this->nevery = atoi(arg[2]);
    if (this->nevery <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");

    // getting the executable and file path from the args passed to this command
    char* exe_path  = arg[3];
    char* sif_path  = arg[4];
    char* surf_path = arg[5];

    // Structure which would store the metadata
    struct stat sb;

    // Calls the function with path as argument
    // If the file/directory exists at the path returns 0
    // If block executes if path exists
    if (!(stat(exe_path,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, exe path does not exist");
    if (!(stat(sif_path,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, sif path does not exist");

    // command style: compute id surf group-id mix-id args
    char* compute_args[COMPUTE_SURF_ARGS_SIZE] = {
               arg[6],
        (char*)"surf",
               arg[7],
               arg[8], 
        (char*)"etot"
    };

    // adding the surf path to the surf args
    this->surf_args[0] = surf_path;

    // dump id surf select-id nevery output_file id c_id[*]
    char* dump_args[DUMP_SURF_ARGS_SIZE] = {
               arg[9],
        (char*)"fea",
               arg[11],
        (char*)std::to_string(this->nevery).c_str(),
               arg[10],
        (char*)"id",
               arg[12]
    };

    // adding the needed compute 
    modify->add_compute(COMPUTE_SURF_ARGS_SIZE, compute_args);
    output->add_dump(DUMP_SURF_ARGS_SIZE, dump_args);

    // initializing the surface writing class
    this->writer = new WriteSurf(sparta);

    // must have " 2>&1" at end to pipe stderr to stdout
    std::string command = std::string(exe_path)+" "+std::string(sif_path)+" 2>&1";

    // dump arguments for the dump fea command
    char* dump_fea_modify_args[DUMP_FEA_MODIFY_ARGS_SIZE] = {
        (char*)"command",
        (char*)command.c_str()
    };

    // modifying the fea dump to add command to run
    output->dump[output->ndump - 1]->modify_params(DUMP_FEA_MODIFY_ARGS_SIZE, dump_fea_modify_args);
}

/* ---------------------------------------------------------------------- */

/**
 * deleting the writer on destruction
 */
FixFea::~FixFea() { delete this->writer; }

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP /* | START_OF_STEP */; }


/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() { this->writer->command(SURF_ARGS_SIZE, this->surf_args); }