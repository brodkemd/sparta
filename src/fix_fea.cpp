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
#include "input.h"
#include "read_surf.h"

#include "Elmer/elmer.hpp"

using namespace SPARTA_NS;


// DO NOT CHANGE THESE
enum{INT,DOUBLE};                     // several files
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP}; // several files

/* ---------------------------------------------------------------------- */

// DO NOT CHANGE THIS
#define INVOKED_PER_SURF 32

// THESE CAN CHANGE BUT BE VERY CAREFUL
#define SERR(_msg) error->all(FLERR, _msg)
#define START_TRY try {
#define END_TRY } catch (std::string _msg) { error->allNoFileAndLine(_msg.c_str()); } catch (std::exception& e) { error->all(FLERR, e.what()); } catch (...) { error->all(FLERR, "unidentified error occurred"); }

/* ---------------------------------------------------------------------- */

/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    /** temporary Variables used during construction **/
    char** arr;
    util::int_t size;
    struct stat sb;
    toml::Item_t customID, temp_run_every, temp_nevery, temp_connectflag;
    toml::Item_t compute_args, surf_collide_args, surf_modify_args;
    // making new fea class
    this->fea = new elmer::Elmer(this);

    util::_screen  = &*screen;
    util::_logfile = &*logfile;
    util::_me      = comm->me;

    ULOG("Setting up fix fea");

    // checking if the number of args is correct
    if (narg > 3)
        SERR("Illegal fix fea command, too many inputs");
    else if (narg < 3)
        SERR("Illegal fix fea command, too few inputs");

    // making sure provided config path exists
    if (!(stat(arg[2],  &sb) == 0))
        SERR("Illegal fix fea command, toml config file does not exist");
    ULOG("Loading config from: " + std::string(arg[2]));

    // wrapping in a try catch defined above
    START_TRY

    // parsing the file and generating a data structure
    toml::handler s(arg[2]);

    // getting the val at the path (second input) and setting the variable (first input) to that variable
    s.getAtPath(temp_run_every,           "sparta.run_every",         toml::INT);
    s.getAtPath(temp_nevery,              "sparta.nevery",            toml::INT);
    s.getAtPath(temp_connectflag,         "sparta.connect",           toml::BOOL);
    s.getAtPath(customID,                 "sparta.customID",          toml::STRING);
    s.getAtPath(compute_args,             "sparta.compute",           toml::LIST);
    s.getAtPath(surf_collide_args,        "sparta.surf_collide",      toml::LIST);
    s.getAtPath(surf_modify_args,         "sparta.surf_modify",       toml::LIST);  

    this->run_every   = temp_run_every.toInt();
    this->nevery      = temp_nevery.toInt();
    this->connectflag = temp_connectflag.toBool();

    // sets up elmer
    if (comm->me == 0) {
        // letting fea handle its variable setting
        this->fea->set(s);
        this->fea->setup();
        this->fea->dump();
    }

    END_TRY
    MPI_Barrier(world);

    // makes a surface file if none is provided
    if (!(surf->exist)) {
        this->loadSurf();
    } else
        UERR("can not use fix fea with existing surface, no guarantee it will match with the elmer body");

    if (!grid->exist)
        error->all(FLERR,"Cannot use fix fea before grid is defined");

    // making sure the surface is the correct type
    if (surf->implicit)
        SERR("Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        SERR("Cannot use fix fea with distributed surfs");
    
    // getting the index of the surface variable added
    this->tindex = surf->add_custom((char*)customID.toString().c_str(),DOUBLE,0);
    ULOG("Add surface variable with name: " + customID.toString());
    
    // adding surf collision model needed
    size = toml::listToCharArray(surf_collide_args, arr);
    this->surf->add_collide(size, arr);
    ULOG("added surface collision model");

    // adding compute
    size = toml::listToCharArray(compute_args, arr);
    modify->add_compute(size, arr);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;

    // error checks
    this->cqw = modify->compute[icompute];
    if (icompute < 0)
        SERR("Could not find fix fea compute ID");
    if (cqw->per_surf_flag == 0)
        SERR("Fix fea compute does not compute per-surf info");

    // getting the index in the compute of the energy value
    energy_loc = util::find(compute_args.toVector(), (toml::Item_t)"etot")-4;
    if (energy_loc == util::npos)
        SERR("etot is not provided as compute arg");
    // ULOG("Energy flux index = " + std::to_string(energy_loc));

    // getting the indicies for the pressures (stresses) from the compute
    std::vector<util::string_t> opts = {"px", "py", "pz"};
    for (std::size_t i = 0; i < opts.size(); i++) {
        force_locs[i] = util::find(compute_args.toVector(), (toml::Item_t)opts[i])-4;
        if (force_locs[i] == util::npos)
            SERR((opts[i] + " is not provided as compute arg").c_str());
    }

    // getting the indicies for the shear stresses from the compute
    opts = {"shx", "shy", "shz"};
    for (std::size_t i = 0; i < opts.size(); i++) {
        shear_locs[i] = util::find(compute_args.toVector(), (toml::Item_t)opts[i])-4;
        if (shear_locs[i] == util::npos)
            SERR((opts[i] + " is not provided as compute arg").c_str());
    }

    // checking to make sure the indicies just detected do not go above the bounds of the compute array
    std::vector<util::int_t> _temp_indicies = {
        energy_loc, force_locs[0], force_locs[1], force_locs[2], shear_locs[0], shear_locs[1], shear_locs[2]
    };

    util::int_t max = util::max(_temp_indicies);
    if (max > 0 && max > cqw->size_per_surf_cols)
        SERR("Fix fea compute array is accessed out-of-range");

    ULOG("added surf compute with id: " + std::string(modify->compute[icompute]->id));

    // getting number of processes
    MPI_Comm_size(world, &nprocs);
    ULOG("running on: " + std::to_string(nprocs) + " processes");

    // adding collision by modifying surf
    size = toml::listToCharArray(surf_modify_args, arr);
    this->surf->modify_params(size, arr);
    ULOG("modified surface");

    // initing vars
    this->qw     = NULL;
    this->qw_me  = NULL;
    this->px     = NULL;
    this->px_me  = NULL;
    this->py     = NULL;
    this->py_me  = NULL;
    this->pz     = NULL;
    this->pz_me  = NULL;
    this->shx    = NULL;
    this->shx_me = NULL;
    this->shy    = NULL;
    this->shy_me = NULL;
    this->shz    = NULL;
    this->shz_me = NULL;

    // telling the compute surf etot to run
    modify->addstep_compute_all(1);

    // deleting no longer needed var
    delete [] arr;

    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

/**
 * cleaning up the class and freeing up memory 
 */
FixFea::~FixFea() {
    // deleting the pointer
    delete this->fea;

    // freeing up memory
    memory->destroy(this->pselect);
    memory->destroy(this->qw_me);
    memory->destroy(this->qw);
    memory->destroy(this->px_me);
    memory->destroy(this->px);
    memory->destroy(this->py_me);
    memory->destroy(this->py);
    memory->destroy(this->pz_me);
    memory->destroy(this->pz);
    memory->destroy(this->shx_me);
    memory->destroy(this->shx);
    memory->destroy(this->shy_me);
    memory->destroy(this->shy);
    memory->destroy(this->shz_me);
    memory->destroy(this->shz);

    // removing the variable from the surface
    surf->remove_custom(this->tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the end of each timestep 
 */
int FixFea::setmask() { return 0 | END_OF_STEP; }

/* ---------------------------------------------------------------------- */

/**
 * allocates memory, loads the initial data, and performs some checks
 */
void FixFea::init() {
    // number of surface elements
    this->nsurf = surf->nlocal;

    // creating the memory
    memory->create(this->qw,        this->nsurf, "fea:qw");
    memory->create(this->qw_me,     this->nsurf, "fea:qw_me");
    memory->create(this->px,        this->nsurf, "fea:px");
    memory->create(this->px_me,     this->nsurf, "fea:px_me");
    memory->create(this->py,        this->nsurf, "fea:py");
    memory->create(this->py_me,     this->nsurf, "fea:py_me");
    memory->create(this->pz,        this->nsurf, "fea:px");
    memory->create(this->pz_me,     this->nsurf, "fea:pz_me");
    memory->create(this->shx,       this->nsurf, "fea:shx");
    memory->create(this->shx_me,    this->nsurf, "fea:shx_me");
    memory->create(this->shy,       this->nsurf, "fea:shy");
    memory->create(this->shy_me,    this->nsurf, "fea:shz_me");
    memory->create(this->shz,       this->nsurf, "fea:shz");
    memory->create(this->shz_me,    this->nsurf, "fea:shz_me");
    memory->create(this->pselect, 3*surf->nsurf, "fea:pselect"); // needs "surf->nsurf"

    // setting all of the entries to zero
    memset(this->qw,      0,   this->nsurf*sizeof(double));
    memset(this->qw_me,   0,   this->nsurf*sizeof(double));
    memset(this->px,      0,   this->nsurf*sizeof(double));
    memset(this->px_me,   0,   this->nsurf*sizeof(double));
    memset(this->py,      0,   this->nsurf*sizeof(double));
    memset(this->py_me,   0,   this->nsurf*sizeof(double));
    memset(this->pz,      0,   this->nsurf*sizeof(double));
    memset(this->pz_me,   0,   this->nsurf*sizeof(double));
    memset(this->shx,     0,   this->nsurf*sizeof(double));
    memset(this->shx_me,  0,   this->nsurf*sizeof(double));
    memset(this->shy,     0,   this->nsurf*sizeof(double));
    memset(this->shy_me,  0,   this->nsurf*sizeof(double));
    memset(this->shz,     0,   this->nsurf*sizeof(double));
    memset(this->shz_me,  0,   this->nsurf*sizeof(double));
    memset(this->pselect, 0, 3*this->nsurf*sizeof(int));

    // loading the boundary data
    if (comm->me == 0) {
        START_TRY

        // sets the temperatures of the surface elements
        this->updateTemperatures();
        
        // no need to update surf here
        END_TRY
    }

    // waits for all processes to get here
    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step gets values from the compute and sends them
 * to elmer
 */
void FixFea::end_of_step() {
    //ULOG("end of step");
    util::int_t i, m;
    
    // number of surface elements
    if (this->nsurf != this->surf->nlocal)
        SERR("detected surface change, this is not allowed with fix fea command");

    // required to run the compute
    this->modify->clearstep_compute();
    if (!(this->cqw->invoked_flag & INVOKED_PER_SURF)) {
        this->cqw->compute_per_surf();
        this->cqw->invoked_flag |= INVOKED_PER_SURF;
    }
    this->cqw->post_process_surf();

    // adding the heat flux and all of the stress to the running averages
    m = 0;
    for (i = comm->me; i < this->nsurf; i += this->nprocs) {
        if (!(surf->tris[i].mask & groupbit))
        this->qw_me[i] +=this->cqw->array_surf[m][this->energy_loc];
        this->px_me[i] +=this->cqw->array_surf[m][this->force_locs[0]];
        this->py_me[i] +=this->cqw->array_surf[m][this->force_locs[1]];
        this->pz_me[i] +=this->cqw->array_surf[m][this->force_locs[2]];
        this->shx_me[i]+=this->cqw->array_surf[m][this->shear_locs[0]];
        this->shy_me[i]+=this->cqw->array_surf[m][this->shear_locs[1]];
        this->shz_me[i]+=this->cqw->array_surf[m][this->shear_locs[2]];
        m++;
    }
    MPI_Barrier(world);

    // if should run fea
    if (this->runCondition()) {
        ULOG("got true run condition, running fea");
        double denominator = ((double)this->run_every)/((double)this->nevery);
        // ULOG("Dividing by: " + std::to_string(denominator));
        // averaging the values over time range
        m = 0;
        for (i = comm->me; i < this->nsurf; i += this->nprocs) {
            if (!(surf->tris[i].mask & groupbit))
            this->qw_me[i]  /= denominator;
            this->px_me[i]  /= denominator;
            this->py_me[i]  /= denominator;
            this->pz_me[i]  /= denominator;
            this->shx_me[i] /= denominator;
            this->shy_me[i] /= denominator;
            this->shz_me[i] /= denominator;
            m++;
        }

        // summing all of the per process arrays into shared arrays
        MPI_Allreduce(this->qw_me,  this->qw,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->px_me,  this->px,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->py_me,  this->py,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->pz_me,  this->pz,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shx_me, this->shx, this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shy_me, this->shy, this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(this->shz_me, this->shz, this->nsurf, MPI_DOUBLE, MPI_SUM, world);

        // only run on the main process
        if (comm->me == 0) {
            // wrapping in try to catch if exception is,
            // fea throws exceptions if it encounters an error
            START_TRY
            // this->fea->dumpBefore();
            // checks to see if elmer should be run
            if (this->fea->shouldRun()) {
                this->fea->createInitialConditions();
                this->fea->createBoundaryConditions();
                // runs the fea solver
                this->fea->run();

                // loading new temperatures and surface points from the fea result
                this->updateTemperatures();
                this->updateSurf();
            } else ULOG("Skipping fea, it did not detect sufficient values");

            // dumps the data obtained from the fea solver
            this->fea->dump();
            END_TRY

        }
        // resetting arrays to all zeros
        memset(this->qw_me,  0, this->nsurf*sizeof(double));
        memset(this->px_me,  0, this->nsurf*sizeof(double));
        memset(this->py_me,  0, this->nsurf*sizeof(double));
        memset(this->pz_me,  0, this->nsurf*sizeof(double));
        memset(this->shx_me, 0, this->nsurf*sizeof(double));
        memset(this->shy_me, 0, this->nsurf*sizeof(double));
        memset(this->shz_me, 0, this->nsurf*sizeof(double));
    }
    // telling the compute to run on the next timestep
    modify->addstep_compute(update->ntimestep+1);
    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

/**
 * Condition to run fea on
*/
bool FixFea::runCondition() { return update->ntimestep % this->run_every == 0; }

/* ---------------------------------------------------------------------- */

/**
 * Loads temperature data from file and sets needed variables
 * Must only run on process 0 and wrapped in try
*/
void FixFea::updateTemperatures() {
    // the wall temperature variable
    double *tvector = surf->edvec[surf->ewhich[tindex]];

    // loading per node temperature data
    ULOG("updating surface temperatures");

    // averages node temperatures on a per surface basis
    this->fea->averageNodeTemperaturesInto(tvector, this->nsurf);

    // checking to make sure all values where set to something non zero
    for (int i = 0; i < this->nsurf; i++) {
        if (tvector[i] == 0.0) SERR("wall temperature not set correctly");
    }
}

/* ---------------------------------------------------------------------- */

/**
 * Updates the surface in sparta using the data from fea
*/
void FixFea::updateSurf() {
    ULOG("updating surface positions");
    // sort particles
    if (particle->exist) particle->sort();

    // connects the surfs, makes it water tight
    if (this->connectflag && this->groupbit != 1) this->connect3dPre();

    /**** this part does the moving */
    unsigned int i;

    // resetting pselect
    memset(this->pselect, 0, 3*this->nsurf*sizeof(int));

    // moving points
    for (i = 0; i < (unsigned)this->nsurf; i++) {
        if (!(surf->tris[i].mask & this->groupbit)) continue;

        this->fea->getNodePointAtIndex(i, 0, surf->tris[i].p1);
        this->pselect[3*i] = 1; // saying the first point of the triangle moved

        this->fea->getNodePointAtIndex(i, 1, surf->tris[i].p2);
        this->pselect[3*i+1] = 1; // saying the second point of the triangle moved

        this->fea->getNodePointAtIndex(i, 2, surf->tris[i].p3);
        this->pselect[3*i+2] = 1; // saying the third point of the triangle moved
    }
    //this->fea->checkUpdatedAllBoundaryNodes();
    
    // connects the surfs, makes it water tight
    if (this->connectflag && this->groupbit != 1) this->connect3dPost();

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
    if (particle->exist) {
        bigint ndeleted = this->removeParticles();
        ULOG("number of particles deleted: " + std::to_string(ndeleted));
    }

    // notify all classes that store per-grid data that grid may have changed
    grid->notify_changed();
}

/* ---------------------------------------------------------------------- */

void FixFea::connect3dPre() {
    ULOG("connecting surface before move");
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

void FixFea::connect3dPost() {
    ULOG("connecting surface after move");
    // check if non-moved points are in hash
    // if so, set their coords to matching point
    // set pselect for newly moved points so remove_particles() will work
    int m,value,j,jwhich;
    double *p[3],*q;
    OnePoint3d key;

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

bigint FixFea::removeParticles() {
    ULOG("removing particles");

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

/* ---------------------------------------------------------------------- */

/**
 * makes surface using data from fea
*/
void FixFea::loadSurf() {
    char* fname = (char*)"";
    char* arr[1];
    int length = 0;
    std::string _temp;

    if (comm->me == 0) {
        ULOG("no surface detected, making surface from fea object");

        START_TRY
        _temp = this->fea->makeSpartaSurf();
        length = (int)(_temp.size());
        fname = (char*)_temp.c_str();
        END_TRY
    }
    MPI_Barrier(world);

    fprintf(screen, "-> %d process here\n", comm->me);
    MPI_Bcast(fname, length, MPI_CHAR, 0, world);
    arr[0] = fname;

    ULOG("running surface reader on resulting surface file: " + std::string(fname));

    if (true) {
        ReadSurf reader(sparta);
        fprintf(screen, "%d process here\n", comm->me);
        reader.command(1, arr);
    }
    // util::string_t cmd = "read_surf " + std::string(fname);
    // char* cmd_name = input->one(cmd.c_str());

    if (!(surf->exist))
        UERR("reading surf was unsuccessful");
    MPI_Barrier(world);
}

// void FixFea::readSurfFile(char* surf_file) {
//     surf->exist = 1;
//     if (surf->distributed) {
//         UERR("Fix fea does not work with distributed surfaces");
//     }

//     // if filename contains "*", search dir for latest surface file

//     if (comm->me == 0) {
//         if (screen) fprintf(screen,"Reading surface file ...\n");
//     }

//     MPI_Barrier(world);
//     double time1 = MPI_Wtime();

//     // -----------------------
//     // read surface data from file(s)
//     // -----------------------

//     int nsurf_total_old = surf->nsurf;
//     int nsurf_old = surf->nlocal;

//     // read file
//     read_file(file,LOCAL);

//     // surf counts, stats, error check

//       // set surf->nsurf based on single file header or base file
//   // set surf_new based on surf->nsurf

//     if (!multiproc) surf->nsurf = nsurf_total_old + nsurf_file;

//     nsurf_new = surf->nsurf;



//     if (surf->nlocal != surf->nsurf)
//       error->all(FLERR,"Surface element count does not match file");

//     delete [] file;

//     MPI_Barrier(world);
//     double time2 = MPI_Wtime();

//     // -----------------------
//     // transform and check surface elements
//     // -----------------------

//     // check for consecutive IDs

//     check_consecutive();

//     // process command-line args
//     // geometry transformations, group, type, etc

//     process_args(1,narg,arg);


//     // output extent of new surfs, tiny ones may have been created by clip

    
//     surf->output_extent(nsurf_old);

//     // compute normals of new surfs
//     surf->compute_tri_normal(nsurf_old);

//     // error check on new surfs
//     // all points must be inside or on surface of simulation box
//     surf->check_point_inside(nsurf_old);

//     // error checks that can be done before surfs are mapped to grid cells


//     surf->check_watertight_3d();
//     check_neighbor_norm_3d();

//     MPI_Barrier(world);
//     double time3 = MPI_Wtime();

//     // -----------------------
//     // map surfs to grid cells
//     // -----------------------

//     // sort particles

//     if (particle->exist) particle->sort();

//     // make list of surf elements I own
//     // clear grid of surf info including split cells

//     surf->setup_owned();
//     grid->unset_neighbors();
//     grid->remove_ghosts();

//     if (particle->exist && grid->nsplitlocal) {
//     Grid::ChildCell *cells = grid->cells;
//     int nglocal = grid->nlocal;
//     for (int icell = 0; icell < nglocal; icell++)
//         if (cells[icell].nsplit > 1)
//         grid->combine_split_cell_particles(icell,1);
//     }

//     grid->clear_surf();

//     MPI_Barrier(world);
//     double time4 = MPI_Wtime();

//     // assign surfs to grid cells

//     grid->surf2grid(1);

//     MPI_Barrier(world);
//     double time5 = MPI_Wtime();

//     // error check on any points too near other surfs
//     // done on per-grid-cell basis, expensive to do globally

//     if (dim == 2) surf->check_point_near_surf_2d();
//     else surf->check_point_near_surf_3d();

//     // re-setup grid ghosts and neighbors

//     grid->setup_owned();
//     grid->acquire_ghosts();
//     grid->reset_neighbors();
//     comm->reset_neighbors();

//     MPI_Barrier(world);
//     double time6 = MPI_Wtime();

//     // flag cells and corners as OUTSIDE or INSIDE

//     grid->set_inout();
//     grid->type_check();

//     // DEBUG
//     //grid->debug();

//     MPI_Barrier(world);
//     double time7 = MPI_Wtime();

//     // remove particles in any cell that is now INSIDE or has new surfs
//     // reassign particles in split cells to sub cell owner
//     // compress particles if any flagged for deletion

//     bigint ndeleted;
//     if (particle->exist) {
//     Grid::ChildCell *cells = grid->cells;
//     Grid::ChildInfo *cinfo = grid->cinfo;
//     int nglocal = grid->nlocal;
//     int delflag = 0;

//     for (int icell = 0; icell < nglocal; icell++) {
//         if (cinfo[icell].type == INSIDE) {
//         if (partflag == KEEP)
//             error->one(FLERR,"Particles are inside new surfaces");
//         if (cinfo[icell].count) delflag = 1;
//         particle->remove_all_from_cell(cinfo[icell].first);
//         cinfo[icell].count = 0;
//         cinfo[icell].first = -1;
//         continue;
//         }
//         if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
//         int nsurf = cells[icell].nsurf;
//         surfint *csurfs = cells[icell].csurfs;
//         int m;
//         if (dim == 2) {
//             Surf::Line *lines = surf->lines;
//             for (m = 0; m < nsurf; m++) {
//             if (lines[csurfs[m]].id > nsurf_total_old) break;
//             }
//         } else {
//             Surf::Tri *tris = surf->tris;
//             for (m = 0; m < nsurf; m++) {
//             if (tris[csurfs[m]].id >= nsurf_total_old) break;
//             }
//         }
//         if (m < nsurf && partflag == CHECK) {
//             if (cinfo[icell].count) delflag = 1;
//             particle->remove_all_from_cell(cinfo[icell].first);
//             cinfo[icell].count = 0;
//             cinfo[icell].first = -1;
//         }
//         }
//         if (cells[icell].nsplit > 1)
//         grid->assign_split_cell_particles(icell);
//     }
//     int nlocal_old = particle->nlocal;
//     if (delflag) particle->compress_rebalance();
//     bigint delta = nlocal_old - particle->nlocal;
//     MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
//     }

//     MPI_Barrier(world);
//     double time8 = MPI_Wtime();

//     // stats

//     double time_total = time6-time1;
//     double time_s2g = time5-time4;

//     if (comm->me == 0) {
//         if (screen) {
//             if (particle->exist)
//             fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
//             fprintf(screen,"  CPU time = %g secs\n",time_total);
//             fprintf(screen,"  read/check/sort/surf2grid/ghost/"
//                     "inout/particle percent = "
//                     "%g %g %g %g %g %g %g\n",
//                     100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
//                     100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
//                     100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
//                     100.0*(time8-time7)/time_total);
//             fprintf(screen,"  surf2grid time = %g secs\n",time_s2g);
//             fprintf(screen,"  map/comm1/comm2/comm3/comm4/split percent = "
//                     "%g %g %g %g %g %g\n",
//                     100.0*grid->tmap/time_s2g,100.0*grid->tcomm1/time_s2g,
//                     100.0*grid->tcomm2/time_s2g,100.0*grid->tcomm3/time_s2g,
//                     100.0*grid->tcomm4/time_s2g,100.0*grid->tsplit/time_s2g);
//         }

//         if (logfile) {
//             if (particle->exist)
//             fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
//             fprintf(logfile,"  CPU time = %g secs\n",time_total);
//             fprintf(logfile,"  read/check/sort/surf2grid/ghost/"
//                     "inout/particle percent = "
//                     "%g %g %g %g %g %g %g\n",
//                     100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
//                     100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
//                     100.0*(time6-time5)/time_total,100.0*(time7-time6)/time_total,
//                     100.0*(time8-time7)/time_total);
//             fprintf(logfile,"  surf2grid time = %g secs\n",time_s2g);
//             fprintf(logfile,"  map/comm1/comm2/comm3/comm4/split percent = "
//                     "%g %g %g %g %g %g\n",
//                     100.0*grid->tmap/time_s2g,100.0*grid->tcomm1/time_s2g,
//                     100.0*grid->tcomm2/time_s2g,100.0*grid->tcomm3/time_s2g,
//                     100.0*grid->tcomm4/time_s2g,100.0*grid->tsplit/time_s2g);
//         }
//     }
// }