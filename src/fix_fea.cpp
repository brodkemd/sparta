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

#include "error.h"
#include "surf.h"
#include "modify.h"
#include "dump.h"
#include "compute.h"
#include "comm.h"
#include "memory.h"
#include "domain.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "update.h"

#include "UTIL/elmer.hpp"

using namespace SPARTA_NS;


enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};


#define START_TRY try {
#define END_TRY } catch (std::string _msg) { error->all(FLERR, _msg.c_str()); } catch (std::exception& e) { error->all(FLERR, e.what()); } catch (...) { error->all(FLERR, "unidentified error occurred"); }


/* ----------------------------------------------------------------------
    
    FixFea class methods

 ---------------------------------------------------------------------- */


/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    /** temporary Variables used during construction **/
    char** arr;
    int size;
    struct stat sb;
    std::string groupID, mixID, customID;
    std::vector<std::string> compute_args, surf_collide_args, surf_modify_args;

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

    // making sure provided config path exists
    if (!(stat(arg[2],  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, toml config file does not exist");
    this->print(("Loading config from: " + std::string(arg[2])).c_str());

    // making new elmer class
    this->elmer = new elmer::Elmer();

    // wrapping in a try catch defined above
    START_TRY

    // parsing the file and generating a data structure
    toml::handler s(arg[2]);

    // getting the val at the path (second input) and setting the variable (first input) to that variable
    s.get_at_path(this->emi,         "sparta.emi");
    s.get_at_path(this->run_every,   "sparta.run_every");
    s.get_at_path(this->nevery,      "sparta.nevery");
    s.get_at_path(this->qwindex,     "sparta.qwindex");
    s.get_at_path(this->threshold,   "sparta.threshold");
    s.get_at_path(groupID,           "sparta.groupID");
    s.get_at_path(mixID,             "sparta.mixID");
    s.get_at_path(customID,          "sparta.customID");
    s.get_at_path(compute_args,      "sparta.compute");
    s.get_at_path(surf_collide_args, "sparta.surf_collide");
    s.get_at_path(surf_modify_args,  "sparta.surf_modify");  
    
    // letting elmer handle its variable setting
    this->elmer->set(s);

    END_TRY

    // making sure emi is valid
    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");

    this->print(("emi = " + std::to_string(this->emi)).c_str());

    if (threshold < (double)0)
        error->all(FLERR, "Fix fea threshold must be greater than 0.0");
    this->print(("threshold = " + std::to_string(this->threshold)).c_str());


    // making sure base temperature is valid
    if (this->elmer->base_temp <= (double)0)
        error->all(FLERR, "base temperature must be greater than 0");
    this->print(("base temperature = " + std::to_string(this->elmer->base_temp)).c_str());

    // will make temperature data file, do not need to check if it exists
    // if (!(stat(this->temperature_data_file.c_str(),  &sb) == 0))
    //     error->all(FLERR, "Illegal fix fea command, temperature_data_file path does not exist");

    this->print(("temperature_data_file = " + this->elmer->temperature_data_file).c_str());

    // checking run_every
    if (this->run_every <= 0)
        error->all(FLERR,"Illegal fix fea command, run_every <= 0");

    this->print(("run_every = " + std::to_string(this->run_every)).c_str());

    // checking nevery
    if (this->nevery <= 0)
        error->all(FLERR,"Illegal fix fea command, nevery <= 0");
    if (this->run_every % this->nevery != 0)
        error->all(FLERR, "Illegal fix fea command, run_every must be a multiple of nevery");

    this->print(("nevery = " + std::to_string(this->nevery)).c_str());


    // get the surface group
    int igroup = surf->find_group(groupID.c_str());
    if (igroup < 0)
        error->all(FLERR,"Fix fea group ID does not exist");
    
    groupbit = surf->bitmask[igroup];
    this->print(("groupID = " + groupID).c_str());

    // getting the mixture
    int imix = particle->find_mixture((char*)mixID.c_str());
    if (imix < 0)
        error->all(FLERR,"Compute thermal/grid mixture ID does not exist");

    // ngroup = particle->mixture[imix]->ngroup;
    this->print(("mixID = " + mixID).c_str());
    
    // creating the custom variable for the surface temperature
    this->tindex = surf->add_custom((char*)customID.c_str(), DOUBLE, 0);
    this->print(("customID = " + customID).c_str());

    // making sure the elmer exe 
    if (!(stat(this->elmer->exe.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, exe path does not exist");
    
    this->print(("exe = " + this->elmer->exe).c_str());


    // if (!(stat(this->elmer->sif.c_str(),  &sb) == 0))
    //     error->all(FLERR, "Illegal fix fea command, sif path does not exist");
    
    // this->print("sif = " + this->elmer->sif);
    

    // making sure the mesh database is complete
    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer->meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error->all(FLERR, ("Illegal fix fea command, mesh database incomplete, " + (this->elmer->meshDBstem + "." + exts[i]) + " does not exist").c_str());
    }
    this->print(("meshDBstem = " + this->elmer->meshDBstem).c_str());    

    surf_collide_args[2] = "s_" + surf_collide_args[2];
    size = toml::vec_to_arr(surf_collide_args, arr);
    this->surf->add_collide(size, arr);

    // adding compute
    size = toml::vec_to_arr(compute_args, arr);
    modify->add_compute(size, arr);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;// modify->find_compute(id_qw);

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

    this->print(("added surf compute with id: " + std::string(modify->compute[icompute]->id)).c_str());

    // prefactor and threshold in Stefan/Boltzmann equation
    // units of prefactor (SI) is K^4 / (watt - m^2)
    // same in 3d vs 2d, since SPARTA treats 2d cell volume as 1 m in z
    this->dimension = domain->dimension;

    // setting variables based on unit system
    // if (strcmp(update->unit_style, "si") == 0) 
    //     threshold = 1.0e-6;
    // else if (strcmp(update->unit_style, "cgs") == 0) 
    //     threshold = 1.0e-3;

    // getting number of processes
    MPI_Comm_size(world, &nprocs);
    this->print(("Number of procs: " + std::to_string(nprocs)).c_str());

    // adding surf collision model needed
    // adding s_ to temperature variable, this is required
    

    // adding collision by modifying surf
    size = toml::vec_to_arr(surf_modify_args, arr);
    this->surf->modify_params(size, arr);
   
    // deleting no longer needed var
    delete [] arr;

    // modifying elmer variables to match sparta vairables
    this->elmer->simulation.Timestep_intervals = {this->run_every};
    this->elmer->simulation.Output_Intervals = {this->run_every};
    this->elmer->simulation.Output_File = this->elmer->temperature_data_file;

    // do not know what this does, but I think it is needed
    /// firstflag = 1;

    // initing var
    qw_avg_me = NULL;

    // this->print("done constructing");

}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->elmer;
    // memory->destroy(tvector_me);
    // memory->destroy(twall);
    memory->destroy(qw_avg_me);
    memory->destroy(qw_avg);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP; }// | START_OF_STEP; }

/* ---------------------------------------------------------------------- */

/**
 * Loads temperature data from file and sets needed variables
 * Must only run on process 0 
*/
void FixFea::load_temperatures() {
    // getting the number of surface elements
    int nlocal = surf->nlocal;

    this->print("loading temperatures");

    START_TRY
    // loading per node temperature data
    this->elmer->loadNodeTemperatureData();    
    END_TRY

    // the wall temperature variable
    double *tvector = surf->edvec[surf->ewhich[tindex]];

    // making sure everything is consistent
    if (nlocal != this->elmer->boundary_data.size())
        error->all(FLERR, ("boundary data does not match required size, required size " + std::to_string(nlocal) + ", size " + std::to_string(this->elmer->boundary_data.size())).c_str());

    // used later
    double avg;
    
    // averages values for the nodes of a surface element and sets this average to the 
    // temperature of the surface element
    for (int i = 0; i < this->elmer->boundary_data.size(); i++) {
        // computes the average temperature of the nodes that make up the surface element
        // this value is used to set the surface element temperature
        avg = 0;
        for (int j = 1; j < elmer::boundary_size; j++) {
            // gets the data point corresponding to node id and adds it to the rolling sum
            avg += this->elmer->node_temperature_data[this->elmer->boundary_data[i][j]-1];
        }
        // computing the average by dividing the sum by the number of points and setting the surface element
        tvector[this->elmer->boundary_data[i][0]-1] = avg/(elmer::boundary_size - 1);
    }

    // checking to make sure all values where set to something non zero
    for (int i = 0; i < nlocal; i++) {
        if (tvector[i] == (double)0)
            error->all(FLERR, "wall temperature not set correctly");
    }
}

void FixFea::init() {
    // number of surface elements
    int nlocal = surf->nlocal;
    this->last_nlocal = nlocal;


    memory->create(this->qw_avg, nlocal, "fea:qw_avg");
    memset(this->qw_avg, 0, nlocal*sizeof(double));

    memory->create(this->qw_avg_me, nlocal, "fea:qw_avg_me");
    memset(this->qw_avg_me, 0, nlocal*sizeof(double));

    // loading the boundary data
    if (comm->me == 0) {
        START_TRY
        this->elmer->setupTemperatureDataFile();
        this->elmer->loadBoundaryData();
        this->elmer->loadElements();
        this->load_temperatures();
        END_TRY
    }

    // waits for all processes to get here
    MPI_Barrier(world);

    // this->start_of_step(true);
    // for (int i = 0; i < nlocal; i++) tvector[i] = twall[i];
    // this->print("done initing");
}

/* ---------------------------------------------------------------------- */

/**
 * Condition to run elmer on
*/
bool FixFea::run_condition() {
    return update->ntimestep % this->run_every == 0;
}


/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() {
    int i,m,mask;
    double qw;
    modify->addstep_compute(update->ntimestep);
    // access source compute or fix
    // set new temperature via Stefan-Boltzmann eq for nown surfs I own
    // use Twall if surf is not in surf group or eng flux is too small
    // compute/fix output is just my nown surfs, indexed by M
    // store in tvector_me = all nlocal surfs, indexed by I
    //Surf::Line *lines = surf->lines;
    //Surf::Tri *tris = surf->tris;

    int nlocal = surf->nlocal;
    if (this->last_nlocal != nlocal)
        error->all(FLERR, "detected surface change, this is not allowed with fix fea command");

    double **array;
    cqw->post_process_surf();
    array = cqw->array_surf;

    int icol = qwindex-1;

    // this->print(std::to_string(particle->nlocal).c_str());

    m = 0;
    for (i = comm->me; i < nlocal; i += nprocs) {
        if (dimension == 3) mask = surf->tris[i].mask;
        else mask = surf->lines[i].mask;

        if (mask & groupbit) {
            this->qw_avg_me[i] = this->qw_avg_me[i] + array[m][icol];
        }
        m++;
    }

    if (this->run_condition()) {
        MPI_Allreduce(this->qw_avg_me, this->qw_avg, nlocal, MPI_DOUBLE, MPI_SUM, world);

        if (comm->me == 0) {
            int i;
            double sum = 0;
            for (i = 0; i < nlocal; i++) sum=sum+std::abs(qw_avg[i]);

            if (sum > this->threshold) {
                this->print(("Got index, " + std::to_string(i) + ", with value " + std::to_string(qw_avg[i])).c_str());
                this->print("Setting up Elmer");
                this->elmer->simulation.Timestep_Sizes = {update->dt};

                this->elmer->bodys.clear();
                this->elmer->initial_conditions.clear();
                this->elmer->boundary_conditions.clear();

                START_TRY

                for (int i = 0; i < nlocal; i++) {
                    elmer::Boundary_Condition* bc = new elmer::Boundary_Condition(i+1);
                    bc->Heat_Flux = this->qw_avg[i]/this->run_every;
                    this->elmer->boundary_conditions.push_back(bc);
                }

                this->elmer->createInitCondsFromData();
                this->elmer->makeSif();
                this->elmer->run();
                
                // done running so reloading temperatures
                this->load_temperatures();
                
                END_TRY

            } else {
                this->print("Skipping elmer, no flux detected");
            }

            START_TRY
            this->print("dumping node temperatures");
            this->elmer->dumpNodeTemperatures(update->ntimestep);
            END_TRY

        }
        // error->all(FLERR, "generated everything");
        for (int i = 0; i < nlocal; i++) this->qw_avg_me[i] = (double)0;
        MPI_Barrier(world);
    }
}

/* ---------------------------------------------------------------------- */


/**
 * custom printing for this class
*/
void FixFea::print(const char* str, int num_indent, const char* end) {
    if (comm->me == 0) {
        std::string space = "";
        for (int i = 0; i < num_indent; i++)
            space += "  ";

        if (screen)  fprintf(screen,  "%s%s%s", space.c_str(), str, end);
        if (logfile) fprintf(logfile, "%s%s%s", space.c_str(), str, end);
    }
}