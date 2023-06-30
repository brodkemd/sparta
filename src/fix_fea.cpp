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
// #include "move_surf.h"

#include "dump.h"
#include "compute.h"
#include "grid.h"
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
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files


#define INVOKED_PER_SURF 32
#define START_TRY try {
#define END_TRY } catch (std::string _msg) { error->all(FLERR, _msg.c_str()); } catch (std::exception& e) { error->all(FLERR, e.what()); } catch (...) { error->all(FLERR, "unidentified error occurred"); }


static std::ostringstream double_converter;

/* ----------------------------------------------------------------------
    
    FixFea class methods

 ---------------------------------------------------------------------- */

/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    /** temporary Variables used during construction **/
    char** arr;
    int size, qwindex;
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
    s.get_at_path(this->threshold,   "sparta.threshold");
    s.get_at_path(qwindex,           "sparta.qwindex");
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
    
    // making sure the mesh database is complete
    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer->meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error->all(FLERR, ("Illegal fix fea command, mesh database incomplete, " + (this->elmer->meshDBstem + "." + exts[i]) + " does not exist").c_str());
    }
    this->print(("meshDBstem = " + this->elmer->meshDBstem).c_str());    

    this->dimension = domain->dimension;
    if (this->dimension != 3)
        error->all(FLERR, "Can not use fix fea with 2d simulation");

    // adding surf collision model needed
    // adding s_ to temperature variable, this is required
    surf_collide_args[2] = "s_" + surf_collide_args[2];
    size = toml::vec_to_arr(surf_collide_args, arr);
    this->surf->add_collide(size, arr);

    // adding compute
    size = toml::vec_to_arr(compute_args, arr);
    modify->add_compute(size, arr);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;

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
    icol = qwindex-1;

    this->print(("added surf compute with id: " + std::string(modify->compute[icompute]->id)).c_str());

    // getting number of processes
    MPI_Comm_size(world, &nprocs);
    this->print(("Number of procs: " + std::to_string(nprocs)).c_str());
    

    // adding collision by modifying surf
    size = toml::vec_to_arr(surf_modify_args, arr);
    this->surf->modify_params(size, arr);
   
    // deleting no longer needed var
    delete [] arr;

    // modifying elmer variables to match sparta vairables
    this->elmer->simulation.Timestep_intervals = { this->run_every };
    this->elmer->simulation.Output_Intervals   = { this->run_every };
    this->elmer->simulation.Output_File        =   this->elmer->temperature_data_file;

    // initing var
    qw_avg_me = NULL;

    // telling the compute surf etot to run
    modify->addstep_compute_all(1);

    // setting up the double converter
    double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10+2);


    // setting moving surf

    // connecting the surface is a good idea to keep it water tight
    // so setting it to one
    this->connectflag = 1;

    // scalar_flag = 1;
    // global_freq = 1;
    // movesurf = new MoveSurf(sparta);
    // movesurf->mode = 1;
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->elmer;
    memory->destroy(qw_avg_me);
    memory->destroy(qw_avg);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP; }

/* ---------------------------------------------------------------------- */

/**
 * Loads temperature data from file and sets needed variables
 * Must only run on process 0 
*/
void FixFea::load_temperatures() {
    // getting the number of surface elements
    int nlocal = surf->nlocal;

    START_TRY
    // loading per node temperature data
    this->print("loading temperatures");
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

/* ---------------------------------------------------------------------- */

void FixFea::init() {
    // number of surface elements
    int nlocal = surf->nlocal;
    this->nsurf = nlocal;

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
}

/* ---------------------------------------------------------------------- */

/**
 * Condition to run elmer on
*/
bool FixFea::run_condition() { return update->ntimestep % this->run_every == 0; }

/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() {
    int i,m;
    double qw,*p1,*p2,*p3;
    
    // number of surface elements
    int nlocal = surf->nlocal;
    if (this->nsurf != nlocal)
        error->all(FLERR, "detected surface change, this is not allowed with fix fea command");

    // required to run the compute
    modify->clearstep_compute();
    if (!(cqw->invoked_flag & INVOKED_PER_SURF)) {
        cqw->compute_per_surf();
        cqw->invoked_flag |= INVOKED_PER_SURF;
    }
    cqw->post_process_surf();

    // getting the data from the compute
    double **array = cqw->array_surf;

    // adding the heat flux to the running average
    m = 0;
    for (i = comm->me; i < nsurf; i += nprocs) {
        if (surf->tris[i].mask & groupbit) this->qw_avg_me[i]+=array[m][icol];
        m++;
    }

    // if should run elmer
    if (this->run_condition()) {
        // summing all of the heat fluxes into one vector, namely qw_avg
        MPI_Allreduce(this->qw_avg_me, this->qw_avg, nsurf, MPI_DOUBLE, MPI_SUM, world);

        // only run on the main process
        if (comm->me == 0) {
            // computing the total heat flux to test if it is less than the threshold
            double sum = 0;
            for (i = 0; i < nsurf; i++) sum+=std::abs(qw_avg[i]);

            // if it less than the threshold, do not run elmer
            if (sum > this->threshold) {
                // messages
                double_converter.str(""); double_converter << sum;
                this->print(("Got sum, " + double_converter.str()).c_str());
                this->print("Setting up Elmer");
                
                // setting the elmer timestep to the sparta timestep
                this->elmer->simulation.Timestep_Sizes = { update->dt };

                // clearing the vectors in elmer to make room for the new data
                this->elmer->bodys.clear();
                this->elmer->initial_conditions.clear();
                this->elmer->boundary_conditions.clear();

                // wrapping in try to catch if exception is,
                // elmer throws exceptions if it encounters an error
                START_TRY
                this->print("creating boundary conditions");
                for (i = 0; i < nsurf; i++) {
                    elmer::Boundary_Condition* bc = new elmer::Boundary_Condition(i+1);
                    bc->Heat_Flux = this->qw_avg[i]/this->run_every; // the average heat flux
                    this->elmer->boundary_conditions.push_back(bc);
                }

                this->print("creating initial conditions");
                this->elmer->createInitCondsFromData();
                
                this->print("making sif");
                this->elmer->makeSif();
                
                this->print("running sif");
                this->elmer->run();
                
                // done running so reloading temperatures
                this->load_temperatures();

                /**
                * moving the surface
                */

                /* -------------------------------
                 loading new surface 
                   ------------------------------- */
                this->move_surf();
                
                
                END_TRY

            } else this->print("Skipping elmer, no flux detected");

            START_TRY
            this->print("dumping node temperatures");
            this->elmer->dumpNodeTemperatures(update->ntimestep);
            END_TRY

        }
        // resetting the heat flux array
        for (i = 0; i < nsurf; i++) this->qw_avg_me[i] = (double)0;
        MPI_Barrier(world);
    }
    // telling the compute to run on the next timestep
    modify->addstep_compute(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixFea::move_surf() {
    int i;
    // sort particles
    if (particle->exist) particle->sort();
    if (connectflag && groupbit != 1) this->connect_3d_pre();

    // this part does the moving
    // for (i = 0; i < nsurf; i++) {
    //     if (!(surf->tris[i].mask & groupbit)) continue;
    //     p1 = surf->tris[i].p1;
    //     p2 = surf->tris[i].p2;
    //     p3 = surf->tris[i].p3;
    // }
    
    if (connectflag && groupbit != 1) this->connect_3d_post();

    surf->compute_tri_normal(0);

    // check that all points are still inside simulation box
    surf->check_point_inside(0);

    // assign split cell particles to parent split cell
    // assign surfs to grid cells
    grid->unset_neighbors();
    grid->remove_ghosts();

    if (grid->nsplitlocal) {
        Grid::ChildCell *cells = grid->cells;
        int nglocal = grid->nlocal;
        for (int icell = 0; icell < nglocal; icell++)
        if (cells[icell].nsplit > 1)
            grid->combine_split_cell_particles(icell,1);
    }

    grid->clear_surf();
    grid->surf2grid(1,0);

    // re-setup owned and ghost cell info
    grid->setup_owned();
    grid->acquire_ghosts();
    grid->reset_neighbors();
    comm->reset_neighbors();

    // flag cells and corners as OUTSIDE or INSIDE
    grid->set_inout();
    grid->type_check(0);

    // remove particles as needed due to surface move
    // set ndeleted for scalar output
    if (particle->exist) ndeleted = this->remove_particles();

    // notify all classes that store per-grid data that grid may have changed
    grid->notify_changed();
}

/* ---------------------------------------------------------------------- */

void FixFea::connect_3d_pre() {
    // hash for corner points of moved triangles
    // key = corner point
    // value = global index (0 to 3*Ntri-1) of the point
    // NOTE: could prealloc hash to correct size here
    hash = new MyHash();

    // add moved points to hash
    double *p1,*p2,*p3;
    OnePoint3d key;

    Surf::Tri *tris = surf->tris;
    int nsurf = surf->nsurf;

    for (int i = 0; i < nsurf; i++) {
        if (!(tris[i].mask & groupbit)) continue;
        p1 = tris[i].p1;
        p2 = tris[i].p2;
        p3 = tris[i].p3;
        key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = p1[2];
        if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+0;
        key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = p2[2];
        if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+1;
        key.pt[0] = p3[0]; key.pt[1] = p3[1]; key.pt[2] = p3[2];
        if (hash->find(key) == hash->end()) (*hash)[key] = 3*i+2;
    }
}

/* ---------------------------------------------------------------------- */

void FixFea::connect_3d_post() {
    // check if non-moved points are in hash
    // if so, set their coords to matching point
    // set pselect for newly moved points so remove_particles() will work
    int m,value,j,jwhich;
    double *p[3],*q;
    OnePoint3d key;

    Surf::Tri *tris = surf->tris;
    int nsurf = surf->nsurf;

    for (int i = 0; i < nsurf; i++) {
        if (tris[i].mask & groupbit) continue;
        p[0] = tris[i].p1;
        p[1] = tris[i].p2;
        p[2] = tris[i].p3;

        for (m = 0; m < 3; m++) {
            key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = p[m][2];
            if (hash->find(key) != hash->end()) {
                value = (*hash)[key];
                j = value/3;
                jwhich = value % 3;
                if (jwhich == 0) q = tris[j].p1;
                else if (jwhich == 1) q = tris[j].p2;
                else q = tris[j].p3;
                p[m][0] = q[0];
                p[m][1] = q[1];
                p[m][2] = q[2];
                if (m == 0) pselect[3*i] = 1;
                else if (m == 1) pselect[3*i+1] = 1;
                else pselect[3*i+2] = 1;
            }
        }
    }

    // free the hash

    delete hash;
}

/* ---------------------------------------------------------------------- */

bigint FixFea::remove_particles() {
    int isurf,nsurf;
    surfint *csurfs;

    Grid::ChildCell *cells = grid->cells;
    Grid::ChildInfo *cinfo = grid->cinfo;
    int nglocal = grid->nlocal;
    int delflag = 0;

    for (int icell = 0; icell < nglocal; icell++) {
        // cell is inside surfs
        // remove particles in case it wasn't before
        if (cinfo[icell].type == INSIDE) {
            if (cinfo[icell].count) delflag = 1;
            particle->remove_all_from_cell(cinfo[icell].first);
            cinfo[icell].count = 0;
            cinfo[icell].first = -1;
            continue;
        }

        // cell has surfs or is split
        // if m < nsurf, loop over csurfs did not finish
        // which means cell contains a moved surf, so delete all its particles
        if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
            nsurf = cells[icell].nsurf;
            csurfs = cells[icell].csurfs;

            int m;

            for (m = 0; m < nsurf; m++) {
                isurf = csurfs[m];
                if (pselect[3*isurf]) break;
                if (pselect[3*isurf+1]) break;
                if (pselect[3*isurf+2]) break;
            }

            if (m < nsurf) {
                if (cinfo[icell].count) delflag = 1;
                particle->remove_all_from_cell(cinfo[icell].first);
                cinfo[icell].count = 0;
                cinfo[icell].first = -1;
            }
        }

        if (cells[icell].nsplit > 1)
            grid->assign_split_cell_particles(icell);
    }

    int nlocal_old = particle->nlocal;
    if (delflag) particle->compress_rebalance();
    bigint delta = nlocal_old - particle->nlocal;
    bigint ndeleted;
    MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    return ndeleted;
}

/* ---------------------------------------------------------------------- */

/**
 * custom printing for this class
*/
void FixFea::print(const char* str, int num_indent, const char* end) {
    // only prints on the main process
    if (comm->me == 0) {
        std::string space = "";
        for (int i = 0; i < num_indent; i++) space += "  ";
        if (screen)  fprintf(screen,  "%s%s%s", space.c_str(), str, end);
        if (logfile) fprintf(logfile, "%s%s%s", space.c_str(), str, end);
    }
}