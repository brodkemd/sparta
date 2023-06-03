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


/* ----------------------------------------------------------------------

Elmer class methods

---------------------------------------------------------------------- */

/**
 * Constructor takes output file and the error instance used in the FixFea class.
 * The error instance is used to display error messages
 * Do not need to check validity of output file path because is already checked in
 * fixfea class
*/
FixFea::Elmer::Elmer(std::string output_file, Error*& _error) {
    this->sif_path = output_file;
    this->error = _error;
}

/**
 * Writes collected commands to the output file. Will write the optional input
 * "header" first if it is provided
*/
void FixFea::Elmer::write(std::string header) {
    // Read from the text file
    std::ofstream sif_file;

    // this->error->message(FLERR, ("Writing sif file: " + sif_path).c_str());

    // opening the node file
    sif_file.open(sif_path);

    // making sure it is open
    if (sif_file.is_open()) {
        // if the header is provided, write it
        if (header.length()) sif_file << header << "\n\n";

        // writing each command in the format required by elmer
        for (int i = 0; i < this->commands.size(); i++) {
            sif_file << this->commands[i][0] << "\n";
            for (int j = 1; j < this->commands[i].size()-1; j++) {
                sif_file << this->tab << this->commands[i][j] << "\n";
            }
            sif_file << this->commands[i][this->commands[i].size() - 1] << "\n\n";
        }
    } else {
        // catching if the file did not open
        error->all(FLERR, ((std::string)"sif file did not open in Elmer class, " + sif_path).c_str());
    }

    // Close the file
    sif_file.close();

    // clear no longer needed commands
    this->commands.clear();
}

/**
 * Adds a boundary condition
*/
void FixFea::Elmer::Boundary_Condition(int _n, std::vector<std::string> args) {
    this->_add_section("Boundary Condition", std::to_string(_n), args);
}

/**
 * Generic method to add a section, such as the "Boundary Condition" section above
*/
void FixFea::Elmer::_add_section(std::string _name, std::string _n, std::vector<std::string> args) {
    // clearing the vector for good measure
    this->v.clear();

    // if _n was provided, add it, _n is "" if it is not provided
    if (_n.length()) this->v.push_back(_name + " " + _n);
    else             this->v.push_back(_name);

    // adding the arguments to the temporary vector
    for (int i = 0; i < args.size(); i++) this->v.push_back(args[i]);
    
    // adding the required "End" keyword to the end of the temp vector
    this->v.push_back("End");

    // adding the temp vector to the larger vector, each temp vector represents a command section
    // like "Boundary Condition" or "Header"
    this->commands.push_back(v);
}


/* ----------------------------------------------------------------------
    
    ConfigPaser class methods

---------------------------------------------------------------------- */

/**
 * Constructor, sets config file and reads it in
*/
FixFea::ConfigParser::ConfigParser(std::string _file_name, Error*& _error, bool ignore_var_name_case) {
    // referencing error
    this->error = _error;

    // making sure the config path exists
    struct stat sb;
    if (!(stat(_file_name.c_str(),  &sb) == 0))
        this->error->all(FLERR, "fix fea config path does not exist");

    this->file_name = _file_name;
    
    // loading the config file
    this->read_file();
}

/**
 * overloading [] to access the class vector while keeping it private
*/
std::pair<std::string, std::string> &FixFea::ConfigParser::operator[](int index) {
    return this->contents[index];
}

/**
 * returns the size of the class vector
*/
std::size_t FixFea::ConfigParser::size() {
    return this->contents.size();
}

/**
 * reads the config file
*/
void FixFea::ConfigParser::read_file(bool ignore_var_name_case) {
    // opening the config file
    std::ifstream file;
    file.open(this->file_name);

    // making sure the file is open, throws error if it isn't
    if (file.is_open()) {
        // temporary variables used in reading data
        std::vector<std::string> lines;
        std::stringstream buffer;
        std::vector<std::string> v;
        std::pair<std::string, std::string> temp;

        // counts line numbers
        int count = 0;

        // reading the file
        buffer << file.rdbuf();

        // splitting the file contents at newlines
        boost::split(lines, buffer.str(), boost::is_any_of("\n"));
        
        // clearing vars for good measure
        buffer.clear();
        this->contents.clear();

        for (std::string line : lines) {
            // trimming whitespaces off of the line
            boost::algorithm::trim(line);

            // filters empty lines
            if (!(line.length()) || line[0] == '#') continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));

            // removing empty strings from vector
            for (unsigned int i = 0; i < v.size(); i++) {
                if ( v.at(i) == "" ) {
                    //remove element if empty string
                    v.erase(v.begin() + i);
                    i--;
                }
            }

            // making sure line is valid
            if (v.size() != 2)
                this->error->all(FLERR, ("fix fea config line " + std::to_string(count+1) + " invalid").c_str());
            
            // cleaning up the args
            boost::algorithm::trim(v[0]); boost::algorithm::trim(v[1]);
            
            // if case insensitive for the variable names
            if (ignore_var_name_case) boost::algorithm::to_lower(v[0]);
            
            // adding the arguments
            temp.first = v[0]; temp.second = v[1];

            // adding to class content vector
            this->contents.push_back(temp);
            count++;
        }
    } else this->error->all(FLERR, "fix fea config file did not open");
}


/* ----------------------------------------------------------------------
    
    FixFea class methods

 ---------------------------------------------------------------------- */


/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
 * to si in sif file used as the format
 * 
 * Command format
 * fix id fea config_file_path
 * 
 * Config file format:
 * arg1 value1
 * arg2 value2
 *  .     .
 *  .     .
 * 
 * Available args
 * nevery,      execute every this many time steps, int
 * exe,         path to the executable for elmer, string
 * sif,         path the sif file to run with elmer, string
 * groupID,     group ID for which surface elements to perform calculation on, string
 * mixID,       mixture ID for particles to perform calculation on, string
 * tsurf_file,  path to surface temperature file, string
 * emi,         emissivity of the surface (unitless, 0 < emisurf <= 1), double
 * customID,    name of a custom per-surf variable to create, string
 * meshDBstem,  stem path to the elmer mesh database, string
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    this->print("Setting up fix fea:", false);

    // checking if the number of args is correct
    if (narg > 3)
        error->all(FLERR,"Illegal fix fea command, too many inputs");
    else if (narg < 3)
        error->all(FLERR,"Illegal fix fea command, too few inputs");

    // making sure there is a surface to analyze
    if (!surf->exist)
        error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // making sure the surface is the correct type
    if (surf->implicit)
        error->all(FLERR,"Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        error->all(FLERR,"Cannot use fix fea with distributed surfs");

    // Structure which would store the metadata
    struct stat sb;

    // parsing args
    ConfigParser inst{std::string(arg[2]), this->error};
    std::pair<std::string, std::string> it;

    std::vector<std::string> required_args = {
        "nevery", "exe", "sif", "groupid", "mixid", "tsurf_file", "emi", "customid", "meshdbstem"
    };

    // going through the args
    for (int i = 0; i < inst.size(); i++) {
        it = inst[i];
        
        // how often to run this command
        if (it.first == "nevery") {
            // getting when to run from the args to this fix command, nevery
            this->nevery = std::stoi(it.second);
            if (this->nevery <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");
            this->print("running every: " + std::to_string(this->nevery) + " steps");

        // name of the custom variable to create
        } else if (it.first == "customid") {
            // getting the variable name to create
            int n = it.second.length() + 1;
            char *id_custom = new char[n];
            strcpy(id_custom, it.second.c_str());

            this->customID = std::string(id_custom);

            // create per-surf temperature vector
            tindex = surf->add_custom(id_custom,DOUBLE,0);
            this->print("Created temperature variable: " + this->customID);
            delete [] id_custom;

        // path to the elmer exe
        } else if (it.first == "exe") {
            this->exe_path = it.second;

            // Calls the function with path as argument
            // If the file/directory exists at the path returns 0
            // If block executes if path exists
            if (!(stat(this->exe_path.c_str(),  &sb) == 0))
                error->all(FLERR,"Illegal fix fea command, exe path does not exist");
            this->print("Elmer exe path: " + this->exe_path);

        // path to elmer sif file
        } else if (it.first == "sif") {
            this->sif_path = it.second;

            // making sure the sif path exists
            if (!(stat(this->sif_path.c_str(),  &sb) == 0))
                error->all(FLERR,"Illegal fix fea command, sif path does not exist");
            this->print("sif file path: " + this->sif_path);

            // making the elmer class
            this->elmer = new Elmer(this->sif_path, this->error);

        // path to elmer mesh database file
        } else if (it.first == "meshdbstem") {
            this->meshDBstem = it.second;

            // making sure all of the component files are a part of the database
            std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
            for (int i = 0; i < 4; i++) {
                if (!(stat((this->meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
                    error->all(FLERR,("Illegal fix fea command, mesh database incomplete, " + (this->meshDBstem + "." + exts[i]) + " does not exist").c_str());
            }
            this->print("Using mesh database at: " + this->meshDBstem);

        // path to surface temperature file
        } else if (it.first == "tsurf_file") {
            this->tsurf_file = it.second;

            // if the twall file does not exist
            if (!(stat(this->tsurf_file.c_str(),  &sb) == 0))
                error->all(FLERR,"Illegal fix fea command, tsurf_file does not exist");
            this->print("Surf temperature file path: " + this->tsurf_file);
        
        // surface group id
        } else if (it.first == "groupid") {
            this->groupID = it.second;

            // get the surface group
            int igroup = surf->find_group(this->groupID.c_str());
            if (igroup < 0)
                error->all(FLERR,"Fix fea group ID does not exist");
            
            groupbit = surf->bitmask[igroup];
            this->print("Surface group id: " + this->groupID);

        // emissivity of the surface
        } else if (it.first == "emi") {
            // getting the surface emissivity
            emi = input->numeric(FLERR, it.second.c_str());
            if (emi <= 0.0 || emi > 1.0)
                error->all(FLERR,"Fix fea emissivity must be > 0.0 and <= 1");

            this->print("surface Emissivity: " + std::to_string(this->emi));

        // gas mixture id
        } else if (it.first == "mixid") {
            this->mixID = it.second;
            this->print("Mixture id: " + this->groupID);

        } else error->all(FLERR, ("Invalid arg: " + it.first).c_str());


        for (int j = 0; j < required_args.size(); j++) {
            if (required_args[j] == it.first) {
                required_args.erase(required_args.begin() + j);
                break;
            }
        }
    }
    if (required_args.size()) error->all(FLERR, "not all args inputted");

    // command style: compute id surf group-id mix-id args
    char* compute_args[COMPUTE_SURF_ARGS_SIZE] = {
        (char*)"CCCC",
        (char*)"surf",
        (char*)this->groupID.c_str(),
        (char*)this->mixID.c_str(),
        (char*)"etot"
    };

    // adding the needed compute
    modify->add_compute(COMPUTE_SURF_ARGS_SIZE, compute_args);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    icompute = modify->ncompute - 1;// modify->find_compute(id_qw);
    this->print("added surf compute with id: " + std::string(modify->compute[icompute]->id));

    // char* fix_ave_args[4] = {
    //     (char*)"AAAA",
    //     (char*)"ave/surf",
    //            arg[5],
    //     (char*)"1"
    // }

    // modify->add_fix();

    

    /********  start temperature stuff  **********/

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

   

    // prefactor and threshold in Stefan/Boltzmann equation
    // units of prefactor (SI) is K^4 / (watt - m^2)
    // same in 3d vs 2d, since SPARTA treats 2d cell volume as 1 m in z
    // int dimension = domain->dimension;

    // setting variables based on unit system
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
    
    // command style: 
    char* surf_collide_args[4] = {
        (char*)"SSSS",
        (char*)"diffuse",
        (char*)("s_" + this->customID).c_str(),
        (char*)"0.5"
    };

    // adding surf collision model needed
    this->surf->add_collide(4, surf_collide_args);

    // command style: 
    char* surf_modify_args[3] = {
        (char*)"all",
        (char*)"collide",
        (char*)"SSSS"
    };

    // modifying params
    this->surf->modify_params(3, surf_modify_args);   

    // must have " 2>&1" at end to pipe stderr to stdout
    this->command = std::string(this->exe_path)+" "+std::string(this->sif_path)+" 2>&1";

    // loading the boundary data
    this->load_boundary();

    // getting the sif file format from the provided
    this->load_sif(this->sif_path);

    // temporary
    error->all(FLERR, "done setting up fix fea");
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete [] id_qw;
    delete this->elmer;
    memory->destroy(tvector_me);
    memory->destroy(twall);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the start and end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP | START_OF_STEP; }

/* ---------------------------------------------------------------------- */

void FixFea::init() {
    this->print("fix fea init");
    // one-time initialization of temperature for all surfs in custom vector
    if (!firstflag) return;
    firstflag = 0;

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;

    memory->create(twall,nlocal,"fea:twall");
    this->start_of_step();

    for (int i = 0; i < nlocal; i++) tvector[i] = twall[i];
    this->print("fix fea init creating memory");
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
                avg += this->data[this->boundary_data[i][j]-1];
            }
            // std::cout << "\n";
            // computing the average by dividing the sum by the number of points and setting the
            // surface element
            this->twall[this->boundary_data[i][0]-1] = avg/(BOUNDARY_DATA_SIZE - 1);
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
    // NOTE: all data is put into tvector
    double *tvector = surf->edvec[surf->ewhich[tindex]];
    MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);

    MPI_Barrier(world);
    if (this->file_handler) {
        this->print("Generating Boundary Conditions");
        for (int i = 0; i < nlocal; i++) {
            this->elmer->Boundary_Condition(i+1, {
                "Target Boundaries(1) = " + std::to_string(this->boundary_data[i][0]-1),
                "Temperature = " + std::to_string(tvector[this->boundary_data[i][0]-1])
            });
        }
        this->print("writing to elemer file");
        this->elmer->write(this->sif_format);
        this->print("done writing to elemer file");


        // running the command
        error->all(FLERR, "would have run command");
        CommandResult command_result = EXEC(this->command);
        
        // if the command did not succeed
        if (command_result.exitstatus) {
            // writing to the logfile and erroring out
            fprintf(logfile, command_result.output.c_str());
            error->all(FLERR, "fix fea failed, see sparta log file");
        }
        std::cout << "done running command\n";
        
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
    this->print("Loading temperature data from: " + this->tsurf_file);
    // std::cout << "loading data\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream data_file;

    // opening the node file
    data_file.open(this->tsurf_file);

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
        error->all(FLERR, ((std::string)"data file did not open, " + this->tsurf_file).c_str());
    }

    // Close the file
    data_file.close();
}

/* ---------------------------------------------------------------------- */

/**
 * Loads boundary element ids and nodes from boundary file
 */
void FixFea::load_boundary() {
    this->print("Loading boundary from: " + this->meshDBstem + ".boundary");

    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream boundary_file;

    // opening the boundary file
    boundary_file.open(this->meshDBstem + ".boundary");

    // temporary vector to store the split strings
    std::vector<std::string> v;
    std::array<int, 4> v_int;

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

            // getting rid of stuff I do not need
            v.erase(v.begin()+1, v.begin() + 5);

            // catching errors
            // this->print("Boundary size: " + std::to_string(v.size()));
            if (v.size() != BOUNDARY_DATA_SIZE) error->all(FLERR, "too many boundary elements to assign");

            // adding the data
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
        if (screen)  fprintf(screen, (space + str + end).c_str());
        if (logfile) fprintf(logfile,(space + str + end).c_str());
    }
}