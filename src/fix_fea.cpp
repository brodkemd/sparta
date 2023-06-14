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

#include "write_surf.h"
#include "surf.h"
#include "modify.h"
#include "output.h"
#include "dump.h"
#include "compute.h"
#include "comm.h"
#include "memory.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"

// #include "toml.h"

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



void FixFea::handle_both(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("emi", &FixFea::handle_emi),
        std::make_pair("tsurf_file", &FixFea::handle_tsurf_file)
    };

    this->run_table(_caller, "both", tbl, options);
}


void FixFea::handle_emi(std::string _caller, toml::node_t val) {
    toml::set(_caller, "emi", toml::double_t, val, this->emi);
    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");

    this->print("emi = " + std::to_string(this->emi));
}


void FixFea::handle_tsurf_file(std::string _caller, toml::node_t val) {
    toml::set(_caller, "tsurf_file", toml::string_t, val, this->tsurf_file);

    if (!(stat(this->tsurf_file.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, tsurf_file path does not exist");

    this->print("tsurf_file = " + this->tsurf_file);
}



void FixFea::handle_sparta(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("nevery", &FixFea::handle_nevery),
        std::make_pair("groupID", &FixFea::handle_groupID),
        std::make_pair("mixID", &FixFea::handle_mixID),
        std::make_pair("customID", &FixFea::handle_customID)
    };

    this->run_table(_caller, "sparta", tbl, options);
}


void FixFea::handle_nevery(std::string _caller, toml::node_t val) {
    toml::set(_caller, "nevery", toml::integer_t, val, this->run_every);
    if (this->run_every <= 0)
        error->all(FLERR,"Illegal fix fea command, nevery <= 0");
    this->print("nevery = " + std::to_string(this->run_every));
}


void FixFea::handle_groupID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "groupID", toml::string_t, val, this->groupID);
    // get the surface group
    int igroup = surf->find_group(this->groupID.c_str());
    if (igroup < 0)
        error->all(FLERR,"Fix fea group ID does not exist");
    
    groupbit = surf->bitmask[igroup];
    this->print("groupID = " + this->groupID);
}


void FixFea::handle_mixID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "mixID", toml::string_t, val, this->mixID);
    int imix = particle->find_mixture((char*)this->mixID.c_str());
    if (imix < 0)
        error->all(FLERR,"Compute thermal/grid mixture ID does not exist");

    ngroup = particle->mixture[imix]->ngroup;
    this->print("mixID = " + this->mixID);
}


void FixFea::handle_customID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "customID", toml::string_t, val, this->customID);
    this->tindex = surf->add_custom((char*)this->customID.c_str(), DOUBLE, 0);
    this->print("customID = " + this->customID);
}


void FixFea::handle_elmer(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("exe", &FixFea::handle_exe),
        std::make_pair("sif", &FixFea::handle_sif),
        std::make_pair("meshDBstem", &FixFea::handle_meshDBstem),
        std::make_pair("header", &FixFea::handle_header),
        std::make_pair("simulation", &FixFea::handle_simulation),
        std::make_pair("constants", &FixFea::handle_constants),
        std::make_pair("solver", &FixFea::handle_solver),
        std::make_pair("equation", &FixFea::handle_equation),
        std::make_pair("material", &FixFea::handle_material),
        std::make_pair("body", &FixFea::handle_body),
        std::make_pair("initial_condition", &FixFea::handle_initial_condition),
        std::make_pair("boundary_condition", &FixFea::handle_boundary_condition)
    };

    this->run_table(_caller, "elmer", tbl, options);
}


void FixFea::handle_exe(std::string _caller, toml::node_t val) {
    toml::set(_caller, "exe", toml::string_t, val, this->elmer.exe);

    if (!(stat(this->elmer.exe.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, exe path does not exist");
    
    this->print("exe = " + this->elmer.exe);
}

void FixFea::handle_sif(std::string _caller, toml::node_t val) {
    toml::set(_caller, "sif", toml::string_t, val, this->elmer.sif);
    
    if (!(stat(this->elmer.sif.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, sif path does not exist");
    
    this->print("sif = " + this->elmer.sif);
}

void FixFea::handle_meshDBstem(std::string _caller, toml::node_t val) {
    toml::set(_caller, "meshDBstem", toml::string_t, val, this->elmer.meshDBstem);

    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer.meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error->all(FLERR, ("Illegal fix fea command, mesh database incomplete, " + (this->elmer.meshDBstem + "." + exts[i]) + " does not exist").c_str());
    }
    this->print("meshDBstem = " + this->elmer.meshDBstem);
}


void FixFea::handle_header(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".header";
    toml::table_value_parser(_caller, __tbl, this->elmer.header.contents, "");
}


void FixFea::handle_simulation(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".simulation";
    table_value_parser(_caller, __tbl, this->elmer.simulation.contents);
}


void FixFea::handle_constants(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".constants";
    table_value_parser(_caller, __tbl, this->elmer.constants.contents);
}


void FixFea::handle_solver(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".solver";
    id_table_value_parser(_caller, tbl, this->elmer.solvers);
}


void FixFea::handle_equation(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".equation";
    id_table_value_parser(_caller, tbl, this->elmer.equations);
}


void FixFea::handle_material(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".material";
    id_table_value_parser(_caller, tbl, this->elmer.materials);
}


void FixFea::handle_body(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".body";
    id_table_value_parser(_caller, tbl, this->elmer.bodys);
}


void FixFea::handle_initial_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".initial_condition";
    id_table_value_parser(_caller, tbl, this->elmer.initial_conditions);
}


void FixFea::handle_boundary_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".boundary_condition";
    id_table_value_parser(_caller, tbl, this->elmer.boundary_conditions);
}




/*-------------------------------------------------


Utility functions


-------------------------------------------------*/

/**
 * executes a table
*/
template<typename dict>
void FixFea::run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options) {
    if (_tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    if (_caller.length())
        _caller = _caller + "." + _name;
    else
        _caller = _name;

    for (toml::dict_item_t it : _options) {
        toml::node_t val = _tbl[it.first];
        (this->*it.second)(_caller, val);
        _tbl.erase(it.first);
    }

    if (_tbl.size()) {
        std::string msg;
        if (_caller.length())
            msg = "Invalid args in section \"" + _caller + "\":\n";
        else
            msg = "Invalid section:\n";
        for (auto it : _tbl) {
            msg+=("-> " + std::string(it.first.str()) + "\n");
        }
        error(FLERR, msg);
    }
}


/**
 * executes a table that wrapped as a node_t type
*/
template<typename dict>
void FixFea::run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options) {

    if (__tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    // converting to a table
    toml::table* _tbl = __tbl.as<toml::table>();

    // calling the other definition because there is nothing else unique needed
    run_table(_caller, _name, (*_tbl), _options);
}

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
    
    FixFea class methods

 ---------------------------------------------------------------------- */

void FixFea::handle_emi(std::string _caller, toml::node_t val) {
    // std::cout << "--> handing: emi\n";
    toml::set(_caller, "emi", toml::double_t, val, this->emi);
    this->print("emi = " + std::to_string(this->emi));
}

void FixFea::handle_tsurf_file(std::string _caller, toml::node_t val) {
    // std::cout << "--> handing: tsurf_file\n";
    toml::set(_caller, "tsurf_file", toml::string_t, val, this->tsurf_file);

    struct stat sb;
    if (!(stat(this->tsurf_file.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, tsurf_file does not exist");
    this->print("tsurf_file = " + this->tsurf_file);
}

void FixFea::handle_both(std::string _caller, toml::node_t tbl) {
    // std::cout << "-> handing: both\n";
    toml::dict_t options = {
        std::make_pair("emi", &FixFea::handle_emi),
        std::make_pair("tsurf_file", &FixFea::handle_tsurf_file)
    };

    this->run_table(_caller, "both", tbl, options);
}


void FixFea::handle_sparta(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("nevery", &FixFea::handle_nevery),
        std::make_pair("groupID", &FixFea::handle_groupID)
    };

    this->run_table(_caller, "sparta", tbl, options);
}

void FixFea::handle_nevery(std::string _caller, toml::node_t val) {
    int64_t temp_nevery;
    toml::set(_caller, "nevery", toml::integer_t, val, temp_nevery);
    
    if (temp_nevery-INT_MIN <= (uint64_t)INT_MAX-INT_MIN)
        this->run_every = temp_nevery;
    else
        error->all(FLERR, "can not convert int64_t to int, the int64_t is too large");

    // setting so that this command runs continuously to compute an average
    this->nevery = 1;

    // error checking
    if (this->run_every <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");
    // this->print("running every: " + std::to_string(this->run_every) + " steps");
    this->print("nevery = " + std::to_string(this->run_every));
}

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
    this->print("Setting up fix fea:", 0);

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

    if (!(stat(arg[2],  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, toml config file does not exist");
    this->print("Loading config from: " + std::string(arg[2]));

    toml::table tbl;
    try {
        tbl = toml::parse_file(arg[2]);
    } catch (const toml::parse_error& err) {
        error->all(FLERR, ("Parsing failed:\n" + std::string(err.description())).c_str());
    }

    toml::dict_t base_options = {
        std::make_pair("both", &FixFea::handle_both)
    };

    try {
        // std::cout << toml::json_formatter(tbl) << "\n";
        this->run_table("", "", tbl, base_options);
    } catch(std::string e) {
        error->all(FLERR, e.c_str());
    }
    
    error->all(FLERR, "Done");


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
    if (strcmp(update->unit_style, "si") == 0) {
        prefactor = 1.0 / (emi * SB_SI);
        threshold = 1.0e-6;
    } else if (strcmp(update->unit_style, "cgs") == 0) {
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
        (char*)std::to_string(this->emi).c_str() // (char*)"0.5"
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

    // char* group_args[4] = {
    //     "GGGG",
    //     "grid",
    //     "intersect",
    //     this->groupID.c_str(),
    //     this->
    // }

    // surf->group()
    
    // temporary
    // error->all(FLERR, "done setting up fix fea");
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete [] id_qw;
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

    // must have " 2>&1" at end to pipe stderr to stdout
    this->command = std::string(this->exe_path)+" "+std::string(this->sif_path)+" 2>&1";

    // loading the boundary data
    this->load_boundary();

    // getting the sif file format from the provided
    this->load_sif(this->sif_path);

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;
    this->last_nlocal = nlocal;

    memory->create(twall,nlocal,"fea:twall");
    memory->create(this->qw_avg, nlocal, "fea:qw_avg");
    memset(this->qw_avg, 0, nlocal*sizeof(double));

    // this->start_of_step(true);

    for (int i = 0; i < nlocal; i++) tvector[i] = twall[i];
    // this->print("fix fea init creating memory");

    // allocate per-surf vector for explicit all surfs
    memory->create(tvector_me,nlocal,"fea:tvector_me");

    // error->all(FLERR, "")
}

/* ---------------------------------------------------------------------- */

bool FixFea::run_condition() {
    return update->ntimestep % this->run_every == 0;
}

/**
 * Loading wall temperature info on start of step
 */
void FixFea::start_of_step() {
    // only running if a valid time step
    // if ((update->ntimestep + 1) % this->run_every == 0 && !force_run) return;
    //this->print("Timestep: " + std::to_string(update->ntimestep) + ", bool: " + std::to_string(this->run_condition()), 4);
    if (!(this->run_condition())) return;

    // this->print("Running start of step: " + std::to_string(update->ntimestep-1), 4);

    // if (!force_run) this->print("Running start of step: " + std::to_string(update->ntimestep));
    // else            this->print("Running start of step");

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

        // clearing no longer needed data
        this->data.clear();

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
    if (this->last_nlocal != nlocal)
        error->all(FLERR, "detected surface change, this is not allowed with fix fea command");

    // if (qwindex == 0) {
    //     double *vector;
    //     if (source == COMPUTE) {
    //         cqw->post_process_surf();
    //         vector = cqw->vector_surf;
    //     } else vector = fqw->vector_surf;

    //     m = 0;
    //     for (i = me; i < nlocal; i += nprocs) {
    //         if (dimension == 3) mask = tris[i].mask;
    //         else mask = lines[i].mask;
    //         if (!(mask & groupbit)) tvector_me[i] = twall[i];
    //         else {
    //             qw = vector[m];
    //             if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
    //             else tvector_me[i] = twall[i];
    //         }
    //         m++;
    //     }
    // } else {
    double **array;
    //if (source == COMPUTE) {
    cqw->post_process_surf();
    array = cqw->array_surf;
    //} else array = fqw->array_surf;

    int icol = qwindex-1;

    // m = 0;
    // for (i = me; i < nlocal; i += nprocs) {
    //     if (dimension == 3) mask = tris[i].mask;
    //     else mask = lines[i].mask;

    //     if (!(mask & groupbit)) tvector_me[i] = twall[i];
    //     else {
    //         qw = array[m][icol];
    //         if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
    //         else tvector_me[i] = twall[i];
    //     }
    //     m++;
    // }

    m = 0;
    for (i = me; i < nlocal; i += nprocs) {
        if (dimension == 3) mask = tris[i].mask;
        else mask = lines[i].mask;

        if (mask & groupbit) {
            qw = array[m][icol];
            if (qw > threshold) this->qw_avg[i] += qw/this->run_every;
            // else tvector_me[i] = twall[i];
        }
        m++;
    }

    if (this->run_condition()) {
        // memset(tvector_me, 0, nlocal*sizeof(double));

        MPI_Barrier(world);
        if (this->file_handler) {
            this->elmer.boundary_conditions.resize(nlocal);
            for (i = me; i < nlocal; i += nprocs) {
                
            }
        }
        MPI_Barrier(world);


        // m = 0;
        // for (i = me; i < nlocal; i += nprocs) {
        //     if (dimension == 3) mask = tris[i].mask;
        //     else mask = lines[i].mask;
        //     if (!(mask & groupbit)) tvector_me[i] = twall[i];
        //     else {
        //         qw = qw_avg[i];
        //         if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
        //         else tvector_me[i] = twall[i];
        //     }
        //     m++;
        // }

        // Allreduce tvector_me with my owned surfs to tvector custom variable
        // so that all procs know new temperature of all surfs
        // NOTE: could possibly just Allreduce a vector size of surface group
        // NOTE: all data is put into tvector
        double *tvector = surf->edvec[surf->ewhich[tindex]];
        MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);

        MPI_Barrier(world);
        if (this->file_handler) {
            this->print("Generating Boundary Conditions", 4);

            this->elmer.boundary_conditions.clear();
            this->elmer.initial_conditions.clear();

            for (int i = 0; i < nlocal; i++) {
                this->elmer->Boundary_Condition(i+1, {
                    "Target Boundaries(1) = " + std::to_string(this->boundary_data[i][0]-1),
                    "Temperature = " + std::to_string(tvector[this->boundary_data[i][0]-1])
                });
            }
            this->print("writing to elemer file", 4);
            this->elmer->write(this->sif_format);
            // this->print("done writing to elemer file", 4);


            // // running the command
            // error->all(FLERR, "would have run command");
            // CommandResult command_result = EXEC(this->command);
            
            // // if the command did not succeed
            // if (command_result.exitstatus) {
            //     // writing to the logfile/screen and erroring out
            //     if (logfile) {
            //         fprintf(logfile, command_result.output.c_str());
            //         error->all(FLERR, "fix fea failed, see sparta log file");
            //     } else if (screen) {
            //         fprintf(screen, command_result.output.c_str());
            //         error->all(FLERR, "fea exited with the above error");
            //     } else error->all(FLERR, "no place to output command error to");
            // }
            
        }
        MPI_Barrier(world);
    }
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
    this->print("Loading temperature data from: " + this->tsurf_file, 4);
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
void FixFea::print(std::string str, int num_indent, std::string end) {
    std::string space = "";
    for (int i = 0; i < num_indent; i++)
        space += "  ";

    if (comm->me == 0) {
        if (screen)  fprintf(screen, (space + str + end).c_str());
        if (logfile) fprintf(logfile,(space + str + end).c_str());
    }
}