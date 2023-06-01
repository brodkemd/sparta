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
   Contributing author: Arnaud Borner (NASA Ames)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_surf_temp_dynamic.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// Stefan-Boltzmann constants for different units and dimensions

#define SB_SI 5.670374419e-8
#define SB_CGS 5.670374419e-5

enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};

/* ---------------------------------------------------------------------- */
/**
 * 0 : ID, documented in fix command
 * 1 : surf/temp/dynamic, style name of this fix command
 * 2 : surf-ID, group ID for which surface elements to consider
 * 3 : Nevery, adjust surface temperature once every Nevery steps
 * 4 : source, computeID or fixID
        computeID = c_ID or c_ID[n] for a compute that calculates per surf values
        fixID = f_ID or f_ID[n] for a fix that calculates per surf values 
 * 5 : Tsurf_file, path to surface temperature file 
 * 6 : Tsurf_file_format,
        elmer = elmer data file format
        sparta = sparta surf data file format
 * 7 : emisurf = emissivity of the surface (unitless, 0 < emisurf <= 1)
 * 8 : customID = name of a custom per-surf variable to create
 */
FixSurfTempDynamic::FixSurfTempDynamic(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    if (narg != 9)
        error->all(FLERR,"Illegal fix surf/temp/dynamic command");
    if (surf->implicit)
        error->all(FLERR,"Cannot use fix surf/temp/dynamic with implicit surfs");
    if (surf->distributed)
        error->all(FLERR,"Cannot use fix surf/temp/dynamic with distributed surfs");

    int igroup = surf->find_group(arg[2]);
    if (igroup < 0)
        error->all(FLERR,"Fix surf/temp/dynamic group ID does not exist");
    
    groupbit = surf->bitmask[igroup];

    nevery = atoi(arg[3]);

    if (strncmp(arg[4],"c_",2) == 0) {
        source = COMPUTE;
        int n = strlen(arg[4]);
        id_qw = new char[n];
        strcpy(id_qw,&arg[4][2]);

        char *ptr = strchr(id_qw,'[');
        if (ptr) {
            if (id_qw[strlen(id_qw)-1] != ']')
                error->all(FLERR,"Invalid source in fix surf/temp/dynamic command");
            qwindex = atoi(ptr+1);
            *ptr = '\0';
        } else qwindex = 0;

        // error checks
        icompute = modify->find_compute(id_qw);
        cqw = modify->compute[icompute];
        if (icompute < 0)
            error->all(FLERR,"Could not find fix surf/temp/dynamic compute ID");
        if (cqw->per_surf_flag == 0)
            error->all(FLERR,"Fix surf/temp/dynamic compute does not compute per-surf info");
        if (qwindex == 0 && cqw->size_per_surf_cols > 0)
            error->all(FLERR,"Fix surf/temp/dynamic compute does not compute per-surf vector");
        if (qwindex > 0 && cqw->size_per_surf_cols == 0)
            error->all(FLERR,"Fix surf/temp/dynamic compute does not compute per-surf array");
        if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
            error->all(FLERR,"Fix surf/temp/dynamic compute array is accessed out-of-range");

    } else if (strncmp(arg[4],"f_",2) == 0) {
        source = FIX;
        int n = strlen(arg[4]);
        id_qw = new char[n];
        strcpy(id_qw,&arg[4][2]);

        char *ptr = strchr(id_qw,'[');
        if (ptr) {
            if (id_qw[strlen(id_qw)-1] != ']')
                error->all(FLERR,"Invalid source in fix surf/temp/dynamic command");
            qwindex = atoi(ptr+1);
            *ptr = '\0';
        } else qwindex = 0;

        // error checks

        ifix = modify->find_fix(id_qw);
        fqw = modify->fix[ifix];
        if (ifix < 0)
            error->all(FLERR,"Could not find fix surf/temp/dynamic fix ID");
        if (fqw->per_surf_flag == 0)
            error->all(FLERR,"Fix surf/temp/dynamic fix does not compute per-surf info");
        if (qwindex == 0 && fqw->size_per_surf_cols > 0)
            error->all(FLERR,"Fix surf/temp/dynamic fix does not compute per-surf vector");
        if (qwindex > 0 && fqw->size_per_surf_cols == 0)
            error->all(FLERR,"Fix surf/temp/dynamic fix does not compute per-surf array");
        if (qwindex > 0 && qwindex > fqw->size_per_surf_cols)
            error->all(FLERR,"Fix surf/temp/dynamic fix array is accessed out-of-range");
        if (nevery % fqw->per_surf_freq)
            error->all(FLERR,"Fix surf/temp/dynamic source not computed at compatible times");

    } else error->all(FLERR,"Invalid source in fix surf/temp/dynamic command");


    this->twall_file = arg[5];
    
    // Structure which would store the metadata
    struct stat sb;

    // if the twall file does not exist
    if (!(stat(twall_file,  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, twall_file does not exist");

    emi = input->numeric(FLERR,arg[7]);
    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR,"Fix surf/temp/dynamic emissivity must be > 0.0 and <= 1");

    // parsing the file format
    if (strcmp(arg[6],"elmer") == 0) {
        this->file_format_flag = 0;
    } else if (strcmp(arg[6],"sparta") == 0) {
        this->file_format_flag = 1;
    } else {
        error->all(FLERR, "Fix surf/temp/dynamic temperature file format must be elmer or sparta");
    }

    int n = strlen(arg[8]) + 1;
    char *id_custom = new char[n];
    strcpy(id_custom,arg[8]);

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
}

/* ---------------------------------------------------------------------- */

FixSurfTempDynamic::~FixSurfTempDynamic() {
    delete [] id_qw;
    memory->destroy(tvector_me);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

int FixSurfTempDynamic::setmask() { return 0 | END_OF_STEP | START_OF_STEP; }

/* ---------------------------------------------------------------------- */

void FixSurfTempDynamic::init() {
    // one-time initialization of temperature for all surfs in custom vector
    if (!firstflag) return;
    firstflag = 0;

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;

    for (int i = 0; i < nlocal; i++)
        tvector[i] = twall[i];

    // allocate per-surf vector for explicit all surfs
    memory->create(tvector_me,nlocal,"surf/temp/dynamic:tvector_me");
}


void FixSurfTempDynamic::load_surf_temps(int _nlocal) {
    // waits till all of the processes get here
    MPI_Barrier(world);

    // if this process is supposed to read the files
    if (this->file_handler) {

        // message
        std::cout << "--> " << comm->me << " reading surf file\n";
        
        // init vars for use later
        std::string line; std::vector<std::string> v;
        std::ifstream file;

        // opening the file
        file.open(this->twall_file);
        if (!file.is_open()) error->all(FLERR, "twall_file did not open");

        // boolean for reading file
        bool go = false;

        // reading the file
        while (file) {
            // getting the latest line
            std::getline(file, line);

            // trimming whitespaces off of the line
            boost::algorithm::trim(line);
            
            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));
            
            switch (this->file_format_flag) {
                case 0:
                    /* code */
                    break;
                
                case 1:
                    // if good to record data
                    if (go) {
                        if (v.size() == 2) {
                            // adding the data to the class array, position in array is first value in array
                            this->twall[std::stoi(v[0])-1] = std::stod(v[1]);
                        }
                    }
                    // sets a bool if based on if the 
                    if (!go) {
                        if (v.size() >= 2) go = (v[1] == (std::string)"SURFS");
                    }
            }
        }
        for (int i = 1; i < nprocs; i++) {
            MPI_Send(&twall, _nlocal, MPI_DOUBLE, i, 1, world);
        }
    } else {
        // master process is 0
        MPI_Recv(&twall, _nlocal, MPI_DOUBLE, 0, 1, world, MPI_STATUS_IGNORE);
    }

    // waits till all of the processes get here
    MPI_Barrier(world);
}


void FixSurfTempDynamic::start_of_step() {
    int nlocal = surf->nlocal;

    // reading the info from the file
    memset(twall, 0, nlocal*sizeof(double));

    // loads from the data file
    this->load_surf_temps(nlocal);
}


void FixSurfTempDynamic::write_elmer_file() {

}


/* ----------------------------------------------------------------------
   compute new surface element temperatures based on heat flux
   only invoked once every Nevery steps
------------------------------------------------------------------------- */
void FixSurfTempDynamic::end_of_step() {
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

    if (file_handler) {
        // writing elmer file
        write_elmer_file();
    }

    // Allreduce tvector_me with my owned surfs to tvector custom variable
    // so that all procs know new temperature of all surfs
    // NOTE: could possibly just Allreduce a vector size of surface group
    double *tvector = surf->edvec[surf->ewhich[tindex]];
    MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);
}
