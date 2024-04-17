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

#include "Elmer/elmer.h"
#include "fix_fea.h"

using namespace SPARTA_NS;

// DO NOT CHANGE THESE
enum{INT,DOUBLE};                     // several files
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP}; // several files

/* ---------------------------------------------------------------------- */

// DO NOT CHANGE THIS
#define INVOKED_PER_SURF 32

// THESE CAN CHANGE BUT BE VERY CAREFUL
#define SERR(_msg) error->all(FLERR, _msg)
#define ROOT 0
#define START_GUARD if (this->comm->me == ROOT) {
#define END_GUARD } MPI_Barrier(world);
#define START_TRY try {
#define END_TRY } catch (std::string _msg) { error->allNoFileAndLine(_msg.c_str()); } catch (std::exception& e) { error->all(FLERR, e.what()); } catch (...) { error->all(FLERR, "unidentified error occurred"); }

/* ---------------------------------------------------------------------- */

// for debugging, prints to terminal from all processes with the process identifier
void FixFea::processMsg(const char* msg) { fprintf(screen, "%d : %s\n", comm->me, msg); }

/*
Sends a char* to all other processes from the root process
*/
void FixFea::BcastString(char*& str) {
    int length;
    if (this->comm->me == ROOT) length = strlen(str)+1;
    MPI_Barrier(world);

    MPI_Bcast(&length,  1, MPI_INT,  ROOT, world);
    if (comm->me != ROOT) str = new char[length];
    MPI_Bcast(str, length, MPI_CHAR, ROOT, world);
    str[length] = '\0';
}

/*
Sends a char* to all other processes from the root process, then
splits the string at spaces into a char** (safer to do this then send
a char** directly)

returns: length of the char** 
*/
int FixFea::BcastStringIntoArrOfStrings(char* str, char**& arr) {
    BcastString(str);
    return util::splitStringInto(str, arr);
}

/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    /** temporary Variables used during construction **/
    char** arr;
    char* str = new char[1];
    int length;
    struct stat sb;
    char *compute_args, *surf_collide_args, *surf_modify_args, *customID;

    // setting pointers for nice printing
    util::_screen  = &*screen;
    util::_logfile = &*logfile;
    util::_me      = comm->me;

    // setting pointers so elmer can grab them
    this->dt       = &this->update->dt;
    this->timestep = &this->update->ntimestep;
    this->boxlo    = &this->domain->boxlo[0];
    this->boxhi    = &this->domain->boxhi[0];

    // got bored and made this
    START_GUARD
    // fprintf(screen, "-----------------------------\n                        ________\n ||         ||        // \n ||         ||       //\n ||         ||      ||\n ||         ||      ||\n ||         ||      ||\n ||         ||      ||\n  \\\\       //        \\\\\n   \\\\_____//          \\\\________\n\nFix fea by the ARL and the University of Cincinnati\n\n-----------------------------\n");
    END_GUARD

    // ULOG("Setting up fix fea");

    // checking if the number of args is correct
    if (narg > 3)
        SERR("Illegal fix fea command, too many inputs");
    else if (narg < 3)
        SERR("Illegal fix fea command, too few inputs");

    // making sure provided config path exists
    if (!(stat(arg[2],  &sb) == 0))
        SERR("Illegal fix fea command, python file does not exist");

    START_GUARD
    fprintf(screen, "  Loading config from: %s\n", arg[2]);
    fprintf(screen, "  Running: Python %s\n", Py_GetVersion());
    
    // ULOG("Loading config from: " + std::string(arg[2]));
    END_GUARD

    // wrapping in a try catch defined above
    START_GUARD
    START_TRY

    // reading python file
    this->python = new python::handler(arg[2]);

    // getting the sparta configuration from the python object
    PyObject* sparta_config = this->python->loadObjectWithSetupFromMain("sparta");
    python::loadAttrFromObjectAndConvert(sparta_config, "run_every",    this->run_every);
    python::loadAttrFromObjectAndConvert(sparta_config, "nevery",       this->nevery);
    python::loadAttrFromObjectAndConvert(sparta_config, "connect",      this->connectflag);
    python::loadAttrFromObjectAndConvert(sparta_config, "customID",     customID);
    python::loadAttrFromObjectAndConvert(sparta_config, "compute",      compute_args);
    python::loadAttrFromObjectAndConvert(sparta_config, "surf_collide", surf_collide_args);
    python::loadAttrFromObjectAndConvert(sparta_config, "surf_modify",  surf_modify_args);

    // making new fea class
    this->fea = new elmer::Elmer(this, comm->me, this->python);

    // dumping data on first timestep (like sparta does)
    this->fea->dump();

    this->should_update_surf = this->fea->shouldUpdateSurf();

    END_TRY
    END_GUARD

    MPI_Bcast(&this->run_every,   1, MPI_INT, ROOT, world);
    MPI_Bcast(&this->nevery,      1, MPI_INT, ROOT, world);
    MPI_Bcast(&this->connectflag, 1, MPI_INT, ROOT, world);
    MPI_Bcast(&this->should_update_surf, 1, MPI_INT, ROOT, world);

    // makes a surface file if none is provided
    if (!(surf->exist))
        this->loadSurf();
    else 
        SERR("can not use fix fea with existing surface, let sparta generate the surface, it ensures correct indicies");

    if (!grid->exist)
        error->all(FLERR,"Cannot use fix fea before grid is defined");

    // making sure the surface is the correct type
    if (surf->implicit)
        SERR("Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        SERR("Cannot use fix fea with distributed surfs");
    
    START_GUARD
    fprintf(screen, "Setting up FEA\n");
    END_GUARD
    
    // adding custom variable and getting the index of the surface variable added
    BcastString(customID);
    this->tindex = surf->add_custom(customID, DOUBLE, 0);
    START_GUARD
    fprintf(screen, "  Added surface variable:                %s\n", customID);
    END_GUARD

    // adding surface collision model
    length = BcastStringIntoArrOfStrings(surf_collide_args, arr);
    this->surf->add_collide(length, arr);
    START_GUARD
    fprintf(screen, "  Added surface collision model with id: %s\n", arr[0]);
    END_GUARD

    // modifying surface
    length = BcastStringIntoArrOfStrings(surf_modify_args, arr);
    this->surf->modify_params(length, arr);
    START_GUARD
    fprintf(screen, "  Modified surface with id:              %s\n", arr[0]);
    END_GUARD

    // adding compute
    length = BcastStringIntoArrOfStrings(compute_args, arr);
    modify->add_compute(length, arr);
    START_GUARD
    fprintf(screen, "  Added compute id:                      %s\n", arr[0]);
    END_GUARD

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;

    // error checks
    this->cqw = modify->compute[icompute];
    if (icompute < 0)
        SERR("Could not find fix fea compute ID");
    if (cqw->per_surf_flag == 0)
        SERR("Fix fea compute does not compute per-surf info");

    // getting the indicies for the pressures (stresses) from the compute
    START_GUARD
    energy_loc = util::findStringInArr((char*)"etot", arr, length);
    if (energy_loc == length)
        SERR("etot is not provided as compute arg");
    
    std::vector<long> _temp_indicies;
    if (this->should_update_surf) {
        force_locs[0] = util::findStringInArr((char*)"px", arr, length);
        if (force_locs[0] == length)
            SERR("px is not provided as compute arg");

        force_locs[1] = util::findStringInArr((char*)"py", arr, length);
        if (force_locs[1] == length)
            SERR("py is not provided as compute arg");

        force_locs[2] = util::findStringInArr((char*)"pz", arr, length);
        if (force_locs[2] == length)
            SERR("pz is not provided as compute arg");
        
        // getting the indicies for the shear stresses from the compute
        shear_locs[0] = util::findStringInArr((char*)"shx", arr, length);
        if (shear_locs[0] == length)
            SERR("shx is not provided as compute arg");

        shear_locs[1] = util::findStringInArr((char*)"shy", arr, length);
        if (shear_locs[1] == length)
            SERR("shy is not provided as compute arg");

        shear_locs[2] = util::findStringInArr((char*)"shz", arr, length);
        if (shear_locs[2] == length)
            SERR("shz is not provided as compute arg");

        // checking to make sure the indicies just detected do not go above the bounds of the compute array
        // they are shifted by 4 because that is the number of args before they are specified
        _temp_indicies = {
            energy_loc-4, force_locs[0]-4, force_locs[1]-4, force_locs[2]-4, shear_locs[0]-4, shear_locs[1]-4, shear_locs[2]-4
        };
    } else { 
        _temp_indicies = {energy_loc-4};
    }

    // if the indices exceed the bounds
    long max = util::max(_temp_indicies);
    if (max > 0 && max > cqw->size_per_surf_cols)
        SERR("Fix fea compute array is accessed out-of-range");
    
    END_GUARD

    // sending the data to all the processes
    MPI_Bcast(&energy_loc, 1, MPI_INT, ROOT, world);
    if (this->should_update_surf) {
        MPI_Bcast(force_locs,  3, MPI_INT, ROOT, world);
        MPI_Bcast(shear_locs,  3, MPI_INT, ROOT, world);
    }
    
    START_GUARD
    fprintf(screen, "  Added surf compute id:                 %s\n", modify->compute[icompute]->id);
    END_GUARD

    // getting number of processes
    MPI_Comm_size(world, &nprocs);

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
    delete str;

    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

/**
 * cleaning up the class and freeing up memory 
 */
FixFea::~FixFea() {
    // deleting the pointer
    START_GUARD
    delete this->fea;
    delete this->python;
    END_GUARD

    // freeing up memory
    memory->destroy(this->pselect);
    memory->destroy(this->qw_me);
    memory->destroy(this->qw);
    if (this->should_update_surf) {
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
    }

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
    this->nsurf = surf->nsurf;

    // creating the memory
    memory->create(this->qw,        this->nsurf, "fea:qw");
    memory->create(this->qw_me,     this->nsurf, "fea:qw_me");
    if (this->should_update_surf) {
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
    }
    memory->create(this->pselect, 3*surf->nsurf, "fea:pselect"); // needs "surf->nsurf"

    // setting all of the entries to zero
    memset(this->qw,      0,   this->nsurf*sizeof(double));
    memset(this->qw_me,   0,   this->nsurf*sizeof(double));
    if (this->should_update_surf) {
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
    }
    memset(this->pselect, 0, 3*this->nsurf*sizeof(int));

    // loading the boundary data
    START_TRY

    // sets the temperatures of the surface elements
    this->updateTemperatures();
    
    // no need to update surf here
    END_TRY

    START_GUARD
    START_TRY
    
    // dump the initial sif file
    this->fea->makeSif();

    END_TRY
    END_GUARD

    // waits for all processes to get here
    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step gets values from the compute and sends them
 * to elmer
 */
void FixFea::end_of_step() {
    long i, m;
    
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
    if (this->should_update_surf) {
        m = 0;
        for (i = comm->me; i < this->nsurf; i += this->nprocs) {
            this->qw_me[i] +=this->cqw->array_surf[m][this->energy_loc];
            this->px_me[i] +=this->cqw->array_surf[m][this->force_locs[0]];
            this->py_me[i] +=this->cqw->array_surf[m][this->force_locs[1]];
            this->pz_me[i] +=this->cqw->array_surf[m][this->force_locs[2]];
            this->shx_me[i]+=this->cqw->array_surf[m][this->shear_locs[0]];
            this->shy_me[i]+=this->cqw->array_surf[m][this->shear_locs[1]];
            this->shz_me[i]+=this->cqw->array_surf[m][this->shear_locs[2]];
            m++;
        }
    } else {
        m = 0;
        for (i = comm->me; i < this->nsurf; i += this->nprocs) {
            this->qw_me[i] +=this->cqw->array_surf[m][this->energy_loc];
            m++;
        }
    }
    MPI_Barrier(world);

    // if should run fea
    if (this->runCondition()) {
        // ULOG("got true run condition, running fea");
        double denominator = ((double)this->run_every)/((double)this->nevery);
        
        // averaging the values over time range
        if (this->should_update_surf) {
            for (i = comm->me; i < this->nsurf; i += this->nprocs) {
                this->qw_me[i]  /= denominator;
                this->px_me[i]  /= denominator;
                this->py_me[i]  /= denominator;
                this->pz_me[i]  /= denominator;
                this->shx_me[i] /= denominator;
                this->shy_me[i] /= denominator;
                this->shz_me[i] /= denominator;
            }

            // summing all of the per process arrays into shared arrays
            MPI_Allreduce(this->qw_me,  this->qw,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->px_me,  this->px,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->py_me,  this->py,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->pz_me,  this->pz,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->shx_me, this->shx, this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->shy_me, this->shy, this->nsurf, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(this->shz_me, this->shz, this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        } else {
            for (i = comm->me; i < this->nsurf; i += this->nprocs) {
                this->qw_me[i]  /= denominator;
            }
            // summing all of the per process arrays into shared arrays
            MPI_Allreduce(this->qw_me,  this->qw,  this->nsurf, MPI_DOUBLE, MPI_SUM, world);
        }

        // wrapping in try to catch if exception is,
        // fea throws exceptions if it encounters an error
        START_GUARD
        START_TRY
        // this->fea->dumpBefore();
        
        // runs the fea solver
        this->fea->run();

        // dumps the data obtained from the fea solver
        this->fea->dump();
        END_TRY
        END_GUARD

        // loading new temperatures and surface points from the fea result
        this->updateTemperatures();
        if (this->should_update_surf)
            this->updateSurf();

        // resetting arrays to all zeros
        memset(this->qw_me,  0, this->nsurf*sizeof(double));
        if (this->should_update_surf) {
            memset(this->px_me,  0, this->nsurf*sizeof(double));
            memset(this->py_me,  0, this->nsurf*sizeof(double));
            memset(this->pz_me,  0, this->nsurf*sizeof(double));
            memset(this->shx_me, 0, this->nsurf*sizeof(double));
            memset(this->shy_me, 0, this->nsurf*sizeof(double));
            memset(this->shz_me, 0, this->nsurf*sizeof(double));
        }
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
    START_GUARD
    this->fea->averageNodeTemperaturesInto(tvector, this->nsurf);

    // checking to make sure all values where set to something non zero
    for (int i = 0; i < this->nsurf; i++) {
        if (tvector[i] <= 0.0)
            SERR("wall temperature not set correctly");
    }
    END_GUARD

    MPI_Bcast(tvector, this->nsurf, MPI_DOUBLE, 0, world);
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

        START_GUARD
        this->fea->getNodePointAtIndex(i, 0, surf->tris[i].p1);
        END_GUARD

        this->pselect[3*i] = 1; // saying the first point of the triangle moved
        MPI_Bcast(surf->tris[i].p1, 3, MPI_DOUBLE, 0, world);

        START_GUARD
        this->fea->getNodePointAtIndex(i, 1, surf->tris[i].p2);
        END_GUARD

        this->pselect[3*i+1] = 1; // saying the second point of the triangle moved
        MPI_Bcast(surf->tris[i].p2, 3, MPI_DOUBLE, 0, world);

        START_GUARD
        this->fea->getNodePointAtIndex(i, 2, surf->tris[i].p3);
        END_GUARD

        this->pselect[3*i+2] = 1; // saying the third point of the triangle moved
        MPI_Bcast(surf->tris[i].p3, 3, MPI_DOUBLE, 0, world);
    }
    //this->fea->checkUpdatedAllBoundaryNodes();
    MPI_Bcast(this->pselect, 3*this->nsurf*sizeof(int), MPI_INT, 0, world);
    
    // connects the surfs, makes it water tight
    if (this->connectflag && this->groupbit != 1) this->connect3dPost();

    // ULOG("tri normal");
    surf->compute_tri_normal(0);

    // ULOG("point inside");
    // check that all points are still inside simulation box
    surf->check_point_inside(0);

    // assign split cell particles to parent split cell
    // assign surfs to grid cells
    // ULOG("unset neighbors");
    grid->unset_neighbors();
    // ULOG("remove ghosts");
    grid->remove_ghosts();

    if (grid->nsplitlocal) {
        Grid::ChildCell *cells = grid->cells;
        int nglocal = grid->nlocal;
        for (int icell = 0; icell < nglocal; icell++) {
            if (cells[icell].nsplit > 1)
                grid->combine_split_cell_particles(icell,1);
        }
    }
    // ULOG("clear surf");
    grid->clear_surf();
    // ULOG("surf to grid");
    grid->surf2grid(1,0);

    // checks
    // ULOG("near surf");
    surf->check_point_near_surf_3d();
    // ULOG("check water tight");
    surf->check_watertight_3d();

    // re-setup owned and ghost cell info
    // ULOG("setup owned");
    grid->setup_owned();
    // ULOG("acquire ghosts");
    grid->acquire_ghosts();
    // ULOG("grid reset neighbors");
    grid->reset_neighbors();
    // ULOG("comm reset neighbors");
    comm->reset_neighbors();

    // flag cells and corners as OUTSIDE or INSIDE
    // ULOG("set inout");
    grid->set_inout();
    // ULOG("type check");
    grid->type_check(0);

    // remove particles as needed due to surface move
    // set ndeleted for scalar output
    // ULOG("remove particles");
    if (particle->exist) {
        bigint ndeleted = this->removeParticles();
        ULOG("number of particles deleted: " + std::to_string(ndeleted));
    }

    // notify all classes that store per-grid data that grid may have changed
    // ULOG("notify change");
    grid->notify_changed();
}

/* ---------------------------------------------------------------------- */

void FixFea::connect3dPre() {
    // ULOG("connecting surface before move");
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
    // ULOG("connecting surface after move");
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
    // ULOG("removing particles");

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
    char **arr, *args;
    long length;
    std::string _temp;

    START_GUARD
    START_TRY
    // ULOG("no surface detected, making surface from fea object");
    PyObject* sparta_config = this->python->loadObjectWithSetupFromMain("sparta");
    python::loadAttrFromObjectAndConvert(sparta_config, "read_surf",  args);
    END_TRY    
    END_GUARD

    length = BcastStringIntoArrOfStrings(args, arr);

    START_GUARD
    START_TRY
    this->fea->makeSpartaSurf(arr[0]);
    END_TRY
    END_GUARD

    // ULOG("running surface reader on resulting surface file: " + std::string(arr[0]));

    ReadSurf reader(sparta);
    reader.command(length, arr);

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