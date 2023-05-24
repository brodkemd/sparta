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

// void genRandomID(int size, char* outstr) {
//     srand((unsigned)time(NULL) * getpid());
//     char alphanum[] = "0123456789_abcdefghijklmnopqrstuvwxyz";
//     for (int i = 0; i < size; i++) outstr[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
// }


/**
 * Execute system command and get STDOUT result.
 * Regular system() only gives back exit status, this gives back output as well.
 * @param command system command to execute
 * @return commandResult containing STDOUT (not stderr) output & exitstatus
 * of command. Empty if command failed (or has no output). If you want stderr,
 * use shell redirection (2&>1).
 */
static CommandResult EXEC(const std::string &command) {
    int exitcode = 0;
    std::array<char, 1048576> buffer {};
    std::string result;

    FILE *pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("popen() failed!");
    }
    try {
        std::size_t bytesread;
        while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
            result += std::string(buffer.data(), bytesread);
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    exitcode = WEXITSTATUS(pclose(pipe));
    return CommandResult{result, exitcode};
}

/* ---------------------------------------------------------------------- */

/***
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

    this->surf_args[0] = surf_path;
    // this->surf_dump_file = arg[10];
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

    char* dump_modify_args[DUMP_MODIFY_ARGS_SIZE] = {
               dump_args[0],
        (char*)"flush",
        (char*)"yes"
    };

    modify->add_compute(COMPUTE_SURF_ARGS_SIZE, compute_args);
    output->add_dump(DUMP_SURF_ARGS_SIZE, dump_args);
    // output->modify_dump(DUMP_MODIFY_ARGS_SIZE, dump_modify_args);

    // error->all(FLERR, this->surf_dump_id);
    // // modifying the latest dump so that the file is flushed every dump
    // for (int i = 0; i < output->ndump; i++) {
    //     error->message(FLERR, output->dump[i]->id);
    // }
    // output->dump[output->ndump-1]->flush_flag = 1;//modify_params(DUMP_MODIFY_ARGS_SIZE, dump_modify_args);

    this->writer = new WriteSurf(sparta); // initializing the surface writing class

    //this->surf_dump_id = arg[9];

    // must have " 2>&1" at end to pipe stderr to stdout
    this->command = std::string(exe_path) + " " + std::string(sif_path) + " 2>&1";

    // delete [] exe_path, sif_path, surf_path;
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->writer;
    this->command.clear();
    return;
}

/* ---------------------------------------------------------------------- */

int FixFea::setmask() {
    int mask = 0;
    mask |= END_OF_STEP;
    //mask |= START_OF_STEP;
    return mask;
}

// void FixFea::start_of_step() {
//     error->message(FLERR, "Running FEA");
//     return;
// }

void FixFea::end_of_step() {
    writer->command(SURF_ARGS_SIZE, this->surf_args);
    CommandResult command_result = EXEC(this->command);

    if (command_result.exitstatus) {
        fprintf(logfile, command_result.output.c_str());
        error->all(FLERR, "fix fea failed, see sparta log file");
    }
    // remove(this->dump_file);
}