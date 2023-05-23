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
#include "error.h"
#include "surf.h"
#include <bits/stdc++.h> 

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    // checking if the number of args is correct
    if (narg > 6)
        error->all(FLERR,"Illegal fix fea command, too many inputs, format \"fix id fea EXE FILE\"");
    else if (narg < 6)
        error->all(FLERR,"Illegal fix fea command, too few inputs, format \"fix id fea EXE FILE\"");

    // making sure there is a surface to analyze
    if (!surf->exist)
        error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // getting when to run from the args to this fix command, nevery
    nevery = atoi(arg[2]);
    if (nevery <= 0)
        error->all(FLERR,"Illegal fix fea command, nevery <= 0");

    // getting the executable and file path from the args passed to this command
    const char* exe_path =  arg[3];
    const char* sif_path =  arg[4];
    this->surf_path = arg[5];

    // Structure which would store the metadata
    struct stat sb;

    // Calls the function with path as argument
    // If the file/directory exists at the path returns 0
    // If block executes if path exists
    if (!(stat(exe_path,  &sb) == 0)) error->all(FLERR,"Illegal fix fea command, exe path does not exist");
    if (!(stat(sif_path,  &sb) == 0)) error->all(FLERR,"Illegal fix fea command, sif path does not exist");
    if (!(stat(this->surf_path, &sb) == 0)) error->all(FLERR,"Illegal fix fea command, surf path does not exist");

    // 
    this->writer = new WriteSurf(sparta);

    this->surf_id = strcat(arg[0], arg[1]);

    // must have " 2>&1" at end to pipe stderr to stdout
    this->command = std::string(exe_path) + " " + std::string(sif_path) + " 2>&1";
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->writer, this->surf_path, this->surf_id;
    this->command.clear();
    return;
}

/* ---------------------------------------------------------------------- */

int FixFea::setmask() {
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}


void FixFea::end_of_step() {
    error->message(FLERR, "Running FEA");

    char* cmds[3] = {strcat(,  ), (char*)'g', (char*)'g'};
    writer->command(3, cmds);

    CommandResult command_result = EXEC(this->command);

    if (command_result.exitstatus) {
        fprintf(logfile, command_result.output.c_str());
        error->all(FLERR, "fix fea failed, see sparta log file");
    }
}