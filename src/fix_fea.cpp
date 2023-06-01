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
#include "stdlib.h"
#include "string.h"
#include <array>

using namespace SPARTA_NS;

#define SB_SI 5.670374419e-8
#define SB_CGS 5.670374419e-5

enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};

struct CommandResult {
    std::string output;
    int exitstatus;
};

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

/**
 * inputs, format INDEX : DESCRIPTION or VALUE
 * 0  : id
 * 1  : fea
 * 2  : nevery, execute every this many time steps
 * 3  : elmer_exe, path to the executable for elmer
 * 4  : sif_file, path the sif file to run with elmer 
 * 5  : group-id, group ID for which surface elements to perform calculation on
 * 6  : mix-ID, mixture ID for particles to perform calculation on
 * 7  : surf-ID, group ID for which surface elements to consider
 * 8  : Tsurf_file, path to surface temperature file 
 * 9  : Tsurf_file_format,
            elmer = elmer data file format
            sparta = sparta surf data file format
 * 10 : emisurf, emissivity of the surface (unitless, 0 < emisurf <= 1)
 * 11 : customID, name of a custom per-surf variable to create
 */
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    // checking if the number of args is correct
    if (narg > 12)
        error->all(FLERR,"Illegal fix fea command, too many inputs");
    else if (narg < 12)
        error->all(FLERR,"Illegal fix fea command, too few inputs");

    // making sure there is a surface to analyze
    if (!surf->exist)
        error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // getting when to run from the args to this fix command, nevery
    this->nevery = atoi(arg[2]);
    if (this->nevery <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");

    // getting the executable and file path from the args passed to this command
    char* exe_path = arg[3];
    char* sif_path = arg[4];

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
        (char*)"CCCC",
        (char*)"surf",
               arg[5],
               arg[6],
        (char*)"etot"
    };

    // adding the needed compute
    modify->add_compute(COMPUTE_SURF_ARGS_SIZE, compute_args);

    /********  start temperature stuff  **********/

    if (surf->implicit)
        error->all(FLERR,"Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        error->all(FLERR,"Cannot use fix fea with distributed surfs");

    int igroup = surf->find_group(arg[7]);
    if (igroup < 0)
        error->all(FLERR,"Fix fea group ID does not exist");
    
    groupbit = surf->bitmask[igroup];

    source = COMPUTE;
    // int n = strlen(arg[4]);
    // id_qw = new char[n];
    // strcpy(id_qw,&arg[4][2]);

    // char *ptr = strchr(id_qw,'[');
    // if (ptr) {
    //     if (id_qw[strlen(id_qw)-1] != ']')
    //         error->all(FLERR,"Invalid source in fix fea command");
    //     qwindex = atoi(ptr+1);
    //     *ptr = '\0';
    // } else
    
    // setting to zero because the compute create above only has one compute value, namely etot
    qwindex = 0;

    // error checks
    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    icompute = modify->ncompute - 1;// modify->find_compute(id_qw);
    cqw = modify->compute[icompute];
    if (icompute < 0)
        error->all(FLERR,"Could not find fix fea compute ID");
    if (cqw->per_surf_flag == 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf info");
    if (qwindex == 0 && cqw->size_per_surf_cols > 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf vector");
    if (qwindex > 0 && cqw->size_per_surf_cols == 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf array");
    if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
        error->all(FLERR,"Fix fea compute array is accessed out-of-range");

    this->twall_file = arg[8];
    
    // Structure which would store the metadata
    struct stat sb;

    // if the twall file does not exist
    if (!(stat(twall_file,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, twall_file does not exist");

    emi = input->numeric(FLERR,arg[10]);
    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR,"Fix fea emissivity must be > 0.0 and <= 1");

    // parsing the file format
    if (strcmp(arg[9],"elmer") == 0) {
        this->file_format_flag = 0;
    } else if (strcmp(arg[9],"sparta") == 0) {
        this->file_format_flag = 1;
    } else {
        error->all(FLERR, "Fix fea temperature file format must be elmer or sparta");
    }

    int n = strlen(arg[11]) + 1;
    char *id_custom = new char[n];
    strcpy(id_custom,arg[11]);

    // create per-surf temperature vector

    tindex = surf->add_custom(id_custom,DOUBLE,0);
    delete [] id_custom;

    // prefactor and threshold in Stefan/Boltzmann equation
    // units of prefactor (SI) is K^4 / (watt - m^2)
    // same in 3d vs 2d, since SPARTA treats 2d cell volume as 1 m in z

    int dimension = domain->dimension;

    if (strcmp(update->unit_style,"si") == 0) {
        prefactor = 1.0 / (emi * SB_SI);
        threshold = 1.0e-6;
    } else if (strcmp(update->unit_style,"cgs") == 0) {
        prefactor = 1.0 / (emi * SB_CGS);
        threshold = 1.0e-3;
    }

    // trigger setup of list of owned surf elements belonging to surf group
    firstflag = 1;

    // initialize data structure
    tvector_me = NULL;

    // if the main process, read the file
    if (comm->me == 0) this->file_handler = true;

    // getting number of processes
    MPI_Comm_size(world, &nprocs);

    /********  end temperature stuff  **********/



    // // command style: compute id surf group-id mix-id args
    // char* fix_args[COMPUTE_SURF_ARGS_SIZE] = {
    //     (char*)"FFFF",
    //     (char*)"surf/temp/dynamic",
    //            arg[5],
    //     (char*)std::to_string(this->nevery).c_str(),
    //     (char*)"c_CCCC"
    // };

    // // dump id surf select-id nevery output_file id c_id[*]
    // char* dump_args[DUMP_SURF_ARGS_SIZE] = {
    //     (char*)"DDDD",
    //     (char*)"fea",
    //     (char*)"dummy",
    //     (char*)std::to_string(this->nevery).c_str(),
    //     (char*)"dummy.dump", // this dump command does not actually write anything
    //     (char*)"id",
    //     (char*)"c_CCCC"
    // };
    
    // output->add_dump(DUMP_SURF_ARGS_SIZE-1, dump_args);

    // initializing the surface writing class
    // this->writer = new WriteSurf(sparta);

    // must have " 2>&1" at end to pipe stderr to stdout
    this->command = std::string(exe_path)+" "+std::string(sif_path)+" 2>&1";


    // // dump arguments for the dump fea command
    // char* dump_fea_modify_args[DUMP_FEA_MODIFY_ARGS_SIZE] = {
    //     (char*)"command",
    //     (char*)command.c_str()
    // };

    std::cout << modify->compute[this->icompute]->id << "\n";


    // error->all(FLERR, std::to_string(this->compute_index).c_str());
    // modifying the fea dump to add command to run
    // output->dump[output->ndump - 1]->modify_params(DUMP_FEA_MODIFY_ARGS_SIZE, dump_fea_modify_args);

    // debug_msg();
}

/* ---------------------------------------------------------------------- */

/**
 * deleting the writer on destruction
 */
// FixFea::~FixFea() {} // { delete this->writer; }

FixFea::~FixFea() {
    delete [] id_qw;
    memory->destroy(tvector_me);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP | START_OF_STEP; }



void FixFea::init() {
    // one-time initialization of temperature for all surfs in custom vector
    if (!firstflag) return;
    firstflag = 0;

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;

    for (int i = 0; i < nlocal; i++)
        tvector[i] = twall[i];

    // allocate per-surf vector for explicit all surfs
    memory->create(tvector_me,nlocal,"fea:tvector_me");
}

void FixFea::start_of_step() {
}

/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() {
    int i,m,mask;
    double qw;

    int me = comm->me;
    // int nprocs = comm->nprocs;
    int dimension = domain->dimension;

    // access source compute or fix
    // set new temperature via Stefan-Boltzmann eq for nown surfs I own
    // use Twall if surf is not in surf group or eng flux is too small
    // compute/fix output is just my nown surfs, indexed by M
    // store in tvector_me = all nlocal surfs, indexed by I

    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;

    int nlocal = surf->nlocal;

    memset(tvector_me, 0, nlocal*sizeof(double));

    if (qwindex == 0) {
        double *vector;
        if (source == COMPUTE) {
            cqw->post_process_surf();
            vector = cqw->vector_surf;
        } else vector = fqw->vector_surf;

        m = 0;
        for (i = me; i < nlocal; i += nprocs) {
            if (dimension == 3) mask = tris[i].mask;
            else mask = lines[i].mask;
            if (!(mask & groupbit)) tvector_me[i] = twall[i];
            else {
                /////////////// remove this
                if (twall[i] == 0) { error->all(FLERR, "twall not set correctly"); }
                //////////////

                qw = vector[m];
                if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
                else tvector_me[i] = twall[i];
            }
            m++;
        }
    } else {
        double **array;
        if (source == COMPUTE) {
            cqw->post_process_surf();
            array = cqw->array_surf;
        } else array = fqw->array_surf;

        int icol = qwindex-1;

        m = 0;
        for (i = me; i < nlocal; i += nprocs) {
            if (dimension == 3) mask = tris[i].mask;
            else mask = lines[i].mask;
            if (!(mask & groupbit)) tvector_me[i] = twall[i];
            else {
                /////////////// remove this
                if (twall[i] == 0) { error->all(FLERR, "twall not set correctly"); }
                ///////////////
                qw = array[m][icol];
                if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
                else tvector_me[i] = twall[i];
            }
            m++;
        }
    }

    // Allreduce tvector_me with my owned surfs to tvector custom variable
    // so that all procs know new temperature of all surfs
    // NOTE: could possibly just Allreduce a vector size of surface group
    double *tvector = surf->edvec[surf->ewhich[tindex]];
    MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);

    MPI_Barrier(world);

    if (me == 0) {
        write_elmer_file();

        // if a command was provided
        if (this->command.length() > 0) {
            // running the command and blocks until it is completed
            
            CommandResult command_result = EXEC(this->command);
            

            // if the command did not succeed
            if (command_result.exitstatus) {
                // writing to the logfile and erroring out
                fprintf(logfile, command_result.output.c_str());
                error->all(FLERR, "fix fea failed, see sparta log file");
            }
            std::cout << "done running command\n";
        }
    }
    MPI_Barrier(world);
} // { this->writer->command(SURF_ARGS_SIZE, this->surf_args); }