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
    util::_screen  = &*screen;
    util::_logfile = &*logfile;
    util::_me = comm->me;

    util::print("Setting up fix fea:", 0);

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
    util::print("Loading config from: " + std::string(arg[2]));

    // making new fea class
    this->fea = new elmer::Elmer();

    // wrapping in a try catch defined above
    START_TRY

    // parsing the file and generating a data structure
    toml::handler s(arg[2]);

    // getting the val at the path (second input) and setting the variable (first input) to that variable
    s.get_at_path(this->run_every,          "sparta.run_every",         true);
    s.get_at_path(this->nevery,             "sparta.nevery",            true);
    s.get_at_path(this->energy_threshold,   "sparta.energy_threshold",  true);
    s.get_at_path(this->force_threshold,    "sparta.force_threshold",   true);
    s.get_at_path(this->shear_threshold,    "sparta.shear_threshold",   true);
    s.get_at_path(this->connectflag,        "sparta.connect",           true);
    s.get_at_path(groupID,                  "sparta.groupID",           true);
    s.get_at_path(mixID,                    "sparta.mixID",             true);
    s.get_at_path(customID,                 "sparta.customID",          true);
    s.get_at_path(compute_args,             "sparta.compute",           true);
    s.get_at_path(surf_collide_args,        "sparta.surf_collide",      true);
    s.get_at_path(surf_modify_args,         "sparta.surf_modify",       true);  

    // letting fea handle its variable setting
    this->fea->set(s);

    END_TRY
    
    this->tindex = surf->add_custom((char*)customID.c_str(),DOUBLE,0);
    // adding surf collision model needed
    // adding s_ to temperature variable, this is required
    size = toml::vec_to_arr(surf_collide_args, arr);
    this->surf->add_collide(size, arr);

    // adding compute
    size = toml::vec_to_arr(compute_args, arr);
    modify->add_compute(size, arr);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;

    // error checks
    this->cqw = modify->compute[icompute];
    if (icompute < 0)
        error->all(FLERR,"Could not find fix fea compute ID");
    if (cqw->per_surf_flag == 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf info");
    // if (qwindex == 0 && cqw->size_per_surf_cols > 0)
    //     error->all(FLERR,"Fix fea compute does not compute per-surf vector");
    // if (qwindex > 0 && cqw->size_per_surf_cols == 0)
    //     error->all(FLERR,"Fix fea compute does not compute per-surf array");
    // if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
    //     error->all(FLERR,"Fix fea compute array is accessed out-of-range");
    // this->icol = qwindex-1;

    energy_loc = util::find(compute_args, (std::string)"etot")-4;
    if (energy_loc < 0)
        error->all(FLERR, "etot is not provided as compute arg");
    util::print("Energy flux index = " + std::to_string(energy_loc));
    
    std::vector<std::string> opts = {"fx", "fy", "fz"};
    for (std::size_t i = 0; i < opts.size(); i++) {
        force_locs[i] = util::find(compute_args, opts[i])-4;
        if (force_locs[i] < 0)
            error->all(FLERR, (opts[i] + " is not provided as compute arg").c_str());
    }

    util::print("fx index = " + std::to_string(force_locs[0]));
    util::print("fy index = " + std::to_string(force_locs[1]));
    util::print("fz index = " + std::to_string(force_locs[2]));

    opts = {"shx", "shy", "shz"};
    for (std::size_t i = 0; i < opts.size(); i++) {
        shear_locs[i] = util::find(compute_args, opts[i])-4;
        if (shear_locs[i] < 0)
            error->all(FLERR, (opts[i] + " is not provided as compute arg").c_str());
    }
    
    util::print("shx index = " + std::to_string(shear_locs[0]));
    util::print("shy index = " + std::to_string(shear_locs[1]));
    util::print("shz index = " + std::to_string(shear_locs[2]));

    std::vector<int> _temp_indicies = {energy_loc, force_locs[0], force_locs[1], force_locs[2], shear_locs[0], shear_locs[1], shear_locs[2]};
    int max = util::max(_temp_indicies);
    if (max > 0 && max > cqw->size_per_surf_cols)
        error->all(FLERR,"Fix fea compute array is accessed out-of-range");

    util::print("added surf compute with id: " + std::string(modify->compute[icompute]->id));

    // getting number of processes
    MPI_Comm_size(world, &nprocs);
    util::print("nprocs = " + std::to_string(nprocs));

    // adding collision by modifying surf
    size = toml::vec_to_arr(surf_modify_args, arr);
    this->surf->modify_params(size, arr);

    // initing vars
    this->qw_avg = NULL;
    this->qw_avg_me = NULL;
    this->fx_avg = NULL;
    this->fx_avg_me = NULL;
    this->fy_avg = NULL;
    this->fy_avg_me = NULL;
    this->fz_avg = NULL;
    this->fz_avg_me = NULL;
    this->shx_avg = NULL;
    this->shx_avg_me = NULL;
    this->shy_avg = NULL;
    this->shy_avg_me = NULL;
    this->shz_avg = NULL;
    this->shz_avg_me = NULL;

    // telling the compute surf etot to run
    modify->addstep_compute_all(1);

    // initial output
    this->ndeleted = 0;

    // deleting no longer needed var
    delete [] arr;
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->fea;
    memory->destroy(this->pselect);
    memory->destroy(this->qw_avg_me);
    memory->destroy(this->qw_avg);
    memory->destroy(this->fx_avg_me);
    memory->destroy(this->fx_avg);
    memory->destroy(this->fy_avg_me);
    memory->destroy(this->fy_avg);
    memory->destroy(this->fz_avg_me);
    memory->destroy(this->fz_avg);
    memory->destroy(this->shx_avg_me);
    memory->destroy(this->shx_avg);
    memory->destroy(this->shy_avg_me);
    memory->destroy(this->shy_avg);
    memory->destroy(this->shz_avg_me);
    memory->destroy(this->shz_avg);
    surf->remove_custom(this->tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP; }

/* ---------------------------------------------------------------------- */

void FixFea::init() {
    // number of surface elements
    // int nlocal = surf->nlocal;
    this->nsurf = surf->nlocal;

    memory->create(this->qw_avg, this->nsurf, "fea:qw_avg");
    memset(this->qw_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->qw_avg_me, this->nsurf, "fea:qw_avg_me");
    memset(this->qw_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->fx_avg, this->nsurf, "fea:fx_avg");
    memset(this->fx_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->fx_avg_me, this->nsurf, "fea:fx_avg_me");
    memset(this->fx_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->fy_avg, this->nsurf, "fea:fy_avg");
    memset(this->fy_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->fy_avg_me, this->nsurf, "fea:fy_avg_me");
    memset(this->fy_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->fx_avg, this->nsurf, "fea:fx_avg");
    memset(this->fx_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->fz_avg_me, this->nsurf, "fea:fz_avg_me");
    memset(this->fz_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->shx_avg, this->nsurf, "fea:shx_avg");
    memset(this->shx_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->shx_avg_me, this->nsurf, "fea:shx_avg_me");
    memset(this->shx_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->shy_avg, this->nsurf, "fea:shy_avg");
    memset(this->shy_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->shz_avg_me, this->nsurf, "fea:shz_avg_me");
    memset(this->shz_avg_me, 0, this->nsurf*sizeof(double));
    
    memory->create(this->shz_avg, this->nsurf, "fea:shz_avg");
    memset(this->shz_avg, 0, this->nsurf*sizeof(double));

    memory->create(this->shz_avg_me, this->nsurf, "fea:shz_avg_me");
    memset(this->shz_avg_me, 0, this->nsurf*sizeof(double));

    memory->create(this->pselect,3*surf->nsurf,"fea:pselect");
    memset(this->pselect, 0, 3*this->nsurf*sizeof(int));

    util::print("done creating memory");

    // loading the boundary data
    if (comm->me == 0) {
        START_TRY
        this->fea->setup();
        this->update_temperatures();
        // no need to update surf here
        END_TRY
    }

    // waits for all processes to get here
    MPI_Barrier(world);
    util::print("done initing");
}

/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() {
    util::print("end of step");
    int i,m;
    
    // number of surface elements
    if (this->nsurf != surf->nlocal)
        error->all(FLERR, "detected surface change, this is not allowed with fix fea command");

    util::print("setting up and running compute");
    // required to run the compute
    modify->clearstep_compute();
    if (!(cqw->invoked_flag & INVOKED_PER_SURF)) {
        cqw->compute_per_surf();
        cqw->invoked_flag |= INVOKED_PER_SURF;
    }
    cqw->post_process_surf();

    // adding the heat flux to the running average
    util::print("averaging");
    m = 0;
    for (i = comm->me; i < nsurf; i += nprocs) {
        if (!(surf->tris[i].mask & groupbit))
        
        this->qw_avg_me[i] +=cqw->array_surf[m][energy_loc];
        util::print("energy");
        this->shy_avg_me[i]+=cqw->array_surf[m][shear_locs[1]];
        util::print("x force");
        this->fy_avg_me[i] +=cqw->array_surf[m][force_locs[1]];
        util::print("y force");
        this->fz_avg_me[i] +=cqw->array_surf[m][force_locs[2]];
        util::print("z force");
        this->shx_avg_me[i]+=cqw->array_surf[m][shear_locs[0]];
        util::print("shx force");
        this->shy_avg_me[i]+=cqw->array_surf[m][shear_locs[1]];
        util::print("shy force");
        this->shz_avg_me[i]+=cqw->array_surf[m][shear_locs[2]];
        util::print("shz force");
        m++;
    }
    util::print("done averaging");

    // if should run fea
    if (this->run_condition()) {
        util::print("running");
        util::print("mpi reducing vars");
        // summing all of the heat fluxes into one vector, namely qw_avg
        MPI_Allreduce(this->qw_avg_me,  this->qw_avg,  nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->fx_avg_me,  this->fx_avg,  nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->fy_avg_me,  this->fy_avg,  nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->fz_avg_me,  this->fz_avg,  nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shx_avg_me, this->shx_avg, nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shy_avg_me, this->shy_avg, nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shz_avg_me, this->shz_avg, nsurf, MPI_DOUBLE, MPI_SUM, world);
        

        // only run on the main process
        if (comm->me == 0) {
            std::string var_name;
            // if it less than the threshold, do not run fea
            if (this->checkVarSums(var_name)) {
                // messages
                util::print("got sum over threshold for variable: " + var_name);
                util::print("Setting up fea");

                // wrapping in try to catch if exception is,
                // fea throws exceptions if it encounters an error
                START_TRY

                this->fea->createConditionsFrom(
                    qw_avg, 
                    fx_avg, fy_avg, fz_avg, 
                    shx_avg, shy_avg, shz_avg, 
                    this->run_every, update->dt, this->nsurf
                );

                this->fea->run();

                /* -------------------------------
                 loading new surface and temperatures
                   ------------------------------- */
                this->update_temperatures();
                this->update_surf();

                END_TRY

            } else util::print("Skipping fea, no params detected");

            START_TRY
            util::print("dumping node temperatures and points");
            this->fea->dump(update->ntimestep);
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

double FixFea::compute_scalar() { return (double) this->ndeleted; }

/* ---------------------------------------------------------------------- */

bool FixFea::checkVarSums(std::string& _name) {
    int i;
    double sum;

    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(qw_avg[i]);
    if (sum > this->energy_threshold) {
        _name = "energy";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(fx_avg[i]);
    if (sum > this->force_threshold) {
        _name = "fx";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(fy_avg[i]);
    if (sum > this->force_threshold) {
        _name = "fy";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(fz_avg[i]);
    if (sum > this->force_threshold) {
        _name = "fz";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(shx_avg[i]);
    if (sum > this->shear_threshold) {
        _name = "shx";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(shy_avg[i]);
    if (sum > this->shear_threshold) {
        _name = "shy";
        return true;
    }
    sum = 0;
    for (i = 0; i < nsurf; i++) sum+=std::abs(shz_avg[i]);
    if (sum > this->shear_threshold) {
        _name = "shz";
        return true;
    }
    return false;
}

/**
 * Condition to run fea on
*/
bool FixFea::run_condition() { return update->ntimestep % this->run_every == 0; }

/* ---------------------------------------------------------------------- */


/**
 * Loads temperature data from file and sets needed variables
 * Must only run on process 0 
*/
void FixFea::update_temperatures() {
    // the wall temperature variable
    double *tvector = surf->edvec[surf->ewhich[tindex]];

    START_TRY
    // loading per node temperature data
    util::print("loading temperatures");
    this->fea->averageNodeTemperaturesInto(tvector, this->nsurf);
    END_TRY

    // checking to make sure all values where set to something non zero
    for (int i = 0; i < this->nsurf; i++) {
        if (tvector[i] == (double)0)
            error->all(FLERR, "wall temperature not set correctly");
    }
}

/* ---------------------------------------------------------------------- */

void FixFea::update_surf() {

    // if (this->fea->node_data.size() != (std::size_t)this->nsurf)
    //     error->all(FLERR, "detected a surface change while moving the surface");

    // sort particles
    if (particle->exist) particle->sort();
    if (this->connectflag && this->groupbit != 1) this->connect_3d_pre();

    /**** this part does the moving */
    unsigned int i;

    // resetting pselect
    for (i = 0; i < (unsigned)3*this->nsurf; i++) this->pselect[i] = 0;

    // moving points
    for (i = 0; i < (unsigned)this->nsurf; i++) {
        if (!(surf->tris[i].mask & this->groupbit)) continue;

        this->fea->getNodePointAtIndex(i, 1, surf->tris[i].p1);
        this->pselect[3*i] = 1; // saying the first point of the triangle moved

        this->fea->getNodePointAtIndex(i, 2, surf->tris[i].p2);
        this->pselect[3*i+1] = 1; // saying the second point of the triangle moved

        this->fea->getNodePointAtIndex(i, 3, surf->tris[i].p3);
        this->pselect[3*i+2] = 1; // saying the third point of the triangle moved
    }
    
    if (this->connectflag && this->groupbit != 1) this->connect_3d_post();

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
        for (int icell = 0; icell < nglocal; icell++) {
            if (cells[icell].nsplit > 1)
                grid->combine_split_cell_particles(icell,1);
        }
    }

    grid->clear_surf();
    grid->surf2grid(1,0);

    // checks
    surf->check_point_near_surf_3d();
    surf->check_watertight_3d();

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
    if (particle->exist) this->ndeleted = this->remove_particles();

    // notify all classes that store per-grid data that grid may have changed
    grid->notify_changed();
}

/* ---------------------------------------------------------------------- */

void FixFea::connect_3d_pre() {
    // hash for corner points of moved triangles
    // key = corner point
    // value = global index (0 to 3*Ntri-1) of the point
    // NOTE: could prealloc hash to correct size here
    this->hash = new MyHash();

    // add moved points to hash
    double *p1,*p2,*p3;
    OnePoint3d key;

    for (int i = 0; i < this->nsurf; i++) {
        if (!(surf->tris[i].mask & this->groupbit)) continue;
        p1 = surf->tris[i].p1;
        p2 = surf->tris[i].p2;
        p3 = surf->tris[i].p3;
        key.pt[0] = p1[0]; key.pt[1] = p1[1]; key.pt[2] = p1[2];
        if (this->hash->find(key) == this->hash->end()) (*this->hash)[key] = 3*i+0;
        key.pt[0] = p2[0]; key.pt[1] = p2[1]; key.pt[2] = p2[2];
        if (this->hash->find(key) == this->hash->end()) (*this->hash)[key] = 3*i+1;
        key.pt[0] = p3[0]; key.pt[1] = p3[1]; key.pt[2] = p3[2];
        if (this->hash->find(key) == this->hash->end()) (*this->hash)[key] = 3*i+2;
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

    // Surf::Tri *tris = surf->tris;

    for (int i = 0; i < nsurf; i++) {
        if (surf->tris[i].mask & groupbit) continue;
        p[0] = surf->tris[i].p1;
        p[1] = surf->tris[i].p2;
        p[2] = surf->tris[i].p3;

        for (m = 0; m < 3; m++) {
            key.pt[0] = p[m][0]; key.pt[1] = p[m][1]; key.pt[2] = p[m][2];
            if (hash->find(key) != hash->end()) {
                value = (*hash)[key];
                j = value/3;
                jwhich = value % 3;
                if (jwhich == 0) q = surf->tris[j].p1;
                else if (jwhich == 1) q = surf->tris[j].p2;
                else q = surf->tris[j].p3;
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
    int isurf, cell_nsurf;
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
        // if m < cell_nsurf, loop over csurfs did not finish
        // which means cell contains a moved surf, so delete all its particles
        if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
            cell_nsurf = cells[icell].nsurf;
            csurfs = cells[icell].csurfs;

            int m;

            for (m = 0; m < cell_nsurf; m++) {
                isurf = csurfs[m];
                if (pselect[3*isurf]) break;
                if (pselect[3*isurf+1]) break;
                if (pselect[3*isurf+2]) break;
            }

            if (m < cell_nsurf) {
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