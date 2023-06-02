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

#include "fix_fea.h"
#include "stdlib.h"
#include "string.h"

#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sys/stat.h>

using namespace SPARTA_NS;

#define SB_SI 5.670374419e-8
#define SB_CGS 5.670374419e-5

enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};

// used in the EXEC function below, represents data returned from a system command
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
 * 7  : Tsurf_file, path to surface temperature file 
 * 8  : Tsurf_file_format,
            elmer = elmer data file format
            sparta = sparta surf data file format
 * 9  : emisurf, emissivity of the surface (unitless, 0 < emisurf <= 1)
 * 10 : customID, name of a custom per-surf variable to create
 * 11 : meshDBstem, stem path to the elmer mesh database
 */
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    this->print("Setting up fix fea:", false);
    
    // checking if the number of args is correct
    if (narg > 12) error->all(FLERR,"Illegal fix fea command, too many inputs");
    else if (narg < 12) error->all(FLERR,"Illegal fix fea command, too few inputs");

    // making sure there is a surface to analyze
    if (!surf->exist) error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // getting when to run from the args to this fix command, nevery
    this->nevery = atoi(arg[2]);
    if (this->nevery <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");
    this->print("running every: " + std::to_string(this->nevery) + " steps");

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
    this->print("Elmer exe path: " + std::string(exe_path));

    if (!(stat(sif_path,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, sif path does not exist");
    this->print("sif file path: " + std::string(sif_path));

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
    icompute = modify->ncompute - 1;// modify->find_compute(id_qw);

    this->print("added surf compute with id: " + std::string(modify->compute[icompute]->id));

    /********  start temperature stuff  **********/

    if (surf->implicit)
        error->all(FLERR,"Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        error->all(FLERR,"Cannot use fix fea with distributed surfs");

    int igroup = surf->find_group(arg[5]);
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
    
    // setting to 1 because the compute create above only has one compute value, namely etot
    qwindex = 1;

    // error checks
    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    
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

    this->twall_file = arg[7];

    // if the twall file does not exist
    if (!(stat(twall_file,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, twall_file does not exist");
    this->print("Surf temperature file path: " + std::string(twall_file));

    emi = input->numeric(FLERR,arg[9]);
    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR,"Fix fea emissivity must be > 0.0 and <= 1");

    // parsing the file format
    if (strcmp(arg[8],"elmer") == 0) {
        this->file_format_flag = 0;
    } else if (strcmp(arg[8],"sparta") == 0) {
        this->file_format_flag = 1;
    } else {
        error->all(FLERR, "Fix fea temperature file format must be elmer or sparta");
    }
    this->print("Surf temperature file format: " + std::string(arg[8]));

    int n = strlen(arg[10]) + 1;
    char *id_custom = new char[n];
    strcpy(id_custom,arg[10]);

    // create per-surf temperature vector
    tindex = surf->add_custom(id_custom,DOUBLE,0);
    this->print("Created temperature variable: " + std::string(id_custom));
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
    
    // getting the mesh database path stem from the args
    this->file_stem = std::string(arg[11]);

    this->print("Checking mesh database at: " + file_stem);
    
    // making sure all of the component files are a part of the database
    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->file_stem + "." + exts[i]).c_str(),  &sb) == 0))
            error->all(FLERR,("Illegal fix fea command, mesh database incomplete, " + (this->file_stem + "." + exts[i]) + " does not exist").c_str());
    }


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

    // std::cout << modify->compute[this->icompute]->id << "\n";


    // error->all(FLERR, std::to_string(this->compute_index).c_str());
    // modifying the fea dump to add command to run
    // output->dump[output->ndump - 1]->modify_params(DUMP_FEA_MODIFY_ARGS_SIZE, dump_fea_modify_args);

    this->load_boundary();

    // debug_msg();
    error->all(FLERR, "done setting up fix fea");
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete [] id_qw;
    memory->destroy(tvector_me);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the start and end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP | START_OF_STEP; }

/* ---------------------------------------------------------------------- */

void FixFea::init() {
    // one-time initialization of temperature for all surfs in custom vector
    if (!firstflag) return;
    firstflag = 0;

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;

    for (int i = 0; i < nlocal; i++) tvector[i] = twall[i];

    // allocate per-surf vector for explicit all surfs
    memory->create(tvector_me,nlocal,"fea:tvector_me");
}

/* ---------------------------------------------------------------------- */

/**
 * Loading wall temperature info on start of step
*/
void FixFea::start_of_step() {
    int nlocal = surf->nlocal;

    // reading the info from the file
    memset(twall, 0, nlocal*sizeof(double));

    // waits till all of the processes get here
    MPI_Barrier(world);

    // if this process is supposed to read the files, reading the temperature file
    if (this->file_handler) {
        // loading data from the temperature file, this is usually spit out from elmer
        this->load_data();

        if (nlocal != this->boundary_data.size())
            error->all(FLERR, "boundary data does not match required size");

        // used later
        double avg;
        
        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        for (int i = 0; i < this->boundary_data.size(); i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            avg = 0;
            for (int j = 1; j < BOUNDARY_DATA_SIZE; j++) {
                // std::cout << this->data[std::stoi(boundary_data[i][j])-1] << " ";
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg += this->data[boundary_data[i][j]-1];
            }
            // std::cout << "\n";
            // computing the average by dividing the sum by the number of points and setting the
            // surface element
            this->twall[i] = avg/(BOUNDARY_DATA_SIZE - 1);
        }

        // checking to make sure all values where set
        for (int i = 0; i < nlocal; i++) {
            if (this->twall[i] == (double)0) error->all(FLERR, "wall temperature not set correctly");
        }

        // sending to other processes
        for (int i = 1; i < nprocs; i++) {
            MPI_Send(&twall, nlocal, MPI_DOUBLE, i, 1, world);
        }
    } else {
        // master process is 0, blocks until data is received
        MPI_Recv(&twall, nlocal, MPI_DOUBLE, 0, 1, world, MPI_STATUS_IGNORE);
    }

    // waits till all of the processes get here
    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

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
    if (this->file_handler) {
        // write_elmer_file();

        // if a command was provided
        if (this->command.length() > 0) {
            // running the command
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
}

/* ---------------------------------------------------------------------- */

/**
 * loads data from the elmer output file
*/
void FixFea::load_sif(std::string sif_path) {
    this->print("Loading sif format from: " + sif_path);

    // Read from the text file
    std::ifstream sif_file;

    // opening the node file
    sif_file.open(sif_path);

    // temporary vector to store the split strings
    std::vector<std::string> v;

    // clearing data
    this->sif_format.clear();

    if (sif_file.is_open()) {
        sif_file >> this->sif_format;
    } else {
        // catching if the file did not open
        error->all(FLERR, ((std::string)"sif file did not open, " + sif_path).c_str());
    }

    // Close the file
    sif_file.close();
}

/* ---------------------------------------------------------------------- */

/**
 * loads data from the elmer output file
*/
void FixFea::load_data() {
    this->print("Loading temperature data from: " + this->data_file);
    // std::cout << "loading data\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream data_file;

    // opening the node file
    data_file.open(this->data_file);

    // temporary vector to store the split strings
    std::vector<std::string> v;

    // clearing data
    this->data.clear();

    if (data_file.is_open()) {
        bool go = false;

        // reading the file
        while (data_file) {
            // getting the latest line
            std::getline(data_file, line);

            // trimming whitespaces off of the line
            boost::algorithm::trim(line);
            
            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));
            
            // if good to record data
            if (go) {
                if (v.size() == 1) {
                    // adding the data to the class list 
                    this->data.push_back(std::stod(v[0]));
                }
            } else go = (v[0] == (std::string)"Perm:");
        }
    } else {
        // catching if the file did not open
        error->all(FLERR, ((std::string)"data file did not open, " + this->data_file).c_str());
    }

    // Close the file
    data_file.close();
}

/* ---------------------------------------------------------------------- */

/**
 * Loads boundary element ids and nodes from boundary file
*/
void FixFea::load_boundary() {
    this->print("Loading boundary from: " + this->file_stem + ".boundary");

    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream boundary_file;

    // opening the boundary file
    boundary_file.open(this->file_stem + ".boundary");

    // temporary vector to store the split strings
    std::vector<std::string> v;
    std::array<double, 4> v_int;

    // clearing for good measure
    this->boundary_data.clear();

    if (boundary_file.is_open()) {
        unsigned int count = 1;
        // reading the file
        while (boundary_file) {
            v.clear();

            // getting the latest line
            std::getline(boundary_file, line);

            boost::algorithm::trim(line);
            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));

            if (v[4] != (std::string)"303") {
                error->all(FLERR, ((std::string)("element is not a triangle in boundary file at line " + std::to_string(count))).c_str());
            }

            // if the line contains a certain number of values, erase the first 5
            if (v.size() != BOUNDARY_LINE_SIZE) {
                v.erase(v.begin()+1, v.begin() + 5);
            }

            if (v.size() != BOUNDARY_DATA_SIZE) error->all(FLERR, "too many boundary elements to assign");

            for (int i = 0; i < BOUNDARY_DATA_SIZE; i++) {
                v_int[i] = std::stoi(v[i]);
            }

            // adding the data to the class vector
            this->boundary_data.push_back(v_int);

            // incrementing line counter
            count++;
        }
    } else error->all(FLERR, "boundary file did not open");

    // Close the file
    boundary_file.close();
}

/* ---------------------------------------------------------------------- */

/**
 * custom printing for this class
 */
void FixFea::print(std::string str, bool indent, std::string end) {
    std::string space;
    if (indent) space = "  ";
    else space = "";

    if (comm->me == 0) {
        if (screen)  fprintf(screen,(space + str + end).c_str());
        if (logfile) fprintf(logfile,(space + str + end).c_str());
    }
}