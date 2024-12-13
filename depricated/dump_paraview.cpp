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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_paraview.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define START_ARGS 5

// customize by adding keyword

enum{ID,PROC,XLO,YLO,ZLO,XHI,YHI,ZHI,XC,YC,ZC,VOL,
     COMPUTE,FIX,VARIABLE};
enum{INT,DOUBLE,BIGINT,STRING};        // same as Dump

#define INVOKED_PER_GRID 16
#define CHUNK 8

/* ---------------------------------------------------------------------- */

DumpParaview::DumpParaview(SPARTA *sparta, int narg, char **arg) : Dump(sparta, narg, arg) {
    // checking input values
    if (narg == 4) error->all(FLERR,"No dump grid attributes specified");

    clearstep  = 1; // tells dump to run compute
    first_flag = 1; // dump on first timestep

    dimension = domain->dimension;

    // getting the grid cells we are looking
    int igroup = grid->find_group(arg[2]);
    if (igroup < 0) error->all(FLERR,"Dump grid group ID does not exist");
    groupbit = grid->bitmask[igroup];

    // how often to run
    nevery = atoi(arg[3]);

    char* prefix = arg[4];

    // computes, fixes, variables which the dump accesses
    ncompute = 0;
    id_compute = NULL;
    compute = NULL;

    nfix = 0;
    id_fix = NULL;
    fix = NULL;

    nvariable = 0;
    id_variable = NULL;
    variable = NULL;
    vbuf = NULL;

    parse_fields(narg, arg);
}

/* ---------------------------------------------------------------------- */

DumpParaview::~DumpParaview() {
    delete [] pack_func;
    delete [] vtype;
    delete [] field2index;
    delete [] argindex;

    for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
    memory->sfree(id_compute);
    delete [] compute;

    for (int i = 0; i < nfix; i++) delete [] id_fix[i];
    memory->sfree(id_fix);
    delete [] fix;

    for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
    memory->sfree(id_variable);
    delete [] variable;

    for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
    delete [] vbuf;

    memory->destroy(cpart);
}

/* ---------------------------------------------------------------------- */

void DumpParaview::parse_fields(int narg, char **arg) {
    // allocate field vectors
    pack_func   = new FnPtrPack[narg - START_ARGS];
    vtype       = new int[narg - START_ARGS];
    field2index = new int[narg - START_ARGS];
    argindex    = new int[narg - START_ARGS];

    // customize by adding to if statement
    for (int j = START_ARGS; j < narg; j++) {
        argindex[j-START_ARGS] = -1;
        if (strcmp(arg[j], "id") == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_id;
            if (sizeof(cellint) == sizeof(smallint))
                vtype[j-START_ARGS] = INT;
            else
                vtype[j-START_ARGS] = BIGINT;
        } else if (strcmp(arg[j], "idstr") == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_id;
            vtype[j-START_ARGS] = STRING;
        } else if (strcmp(arg[j], "split") == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_split;
            vtype[j-START_ARGS] = INT;
        } else if (strcmp(arg[j], "proc") == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_proc;
            vtype[j-START_ARGS] = INT;
        } else if (strcmp(arg[j], "vol") == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_vol;
            vtype[j-START_ARGS] = DOUBLE;

        // compute value = c_ID
        // if no trailing [], then index = 0, else index = int between []

        } else if (strncmp(arg[j], "c_", 2) == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_compute;
            vtype[j-START_ARGS] = DOUBLE;

            int n = strlen(arg[j]);
            char *suffix = new char[n];
            strcpy(suffix, &arg[j][2]);

            char *ptr = strchr(suffix,'[');
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']')
                error->all(FLERR,"Invalid attribute in dump grid command");
                argindex[j-START_ARGS] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j-START_ARGS] = 0;

            n = modify->find_compute(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump grid compute ID");
            if (modify->compute[n]->per_grid_flag == 0)
                error->all(FLERR,"Dump grid compute does not compute per-grid info");
            if (argindex[j-START_ARGS] == 0 && modify->compute[n]->size_per_grid_cols != 0)
                error->all(FLERR,"Dump grid compute does not calculate per-grid vector");
            if (argindex[j-START_ARGS] > 0 && modify->compute[n]->size_per_grid_cols == 0)
                error->all(FLERR,"Dump grid compute does not calculate per-grid array");
            if (argindex[j-START_ARGS] > 0 && argindex[j-START_ARGS] > modify->compute[n]->size_per_grid_cols)
                error->all(FLERR,"Dump grid compute array is accessed out-of-range");

            field2index[j-START_ARGS] = add_compute(suffix);
            delete [] suffix;

        // fix value = f_ID
        // if no trailing [], then index = 0, else index = int between []
        // if index = 0 and compute stores array, expand to one value per column
        } else if (strncmp(arg[j],"f_",2) == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_fix;
            vtype[j-START_ARGS] = DOUBLE;

            int n = strlen(arg[j]);
            char *suffix = new char[n];
            strcpy(suffix,&arg[j][2]);

            char *ptr = strchr(suffix,'[');
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']')
                error->all(FLERR,"Invalid attribute in dump grid command");
                argindex[j-START_ARGS] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[j-START_ARGS] = 0;

            n = modify->find_fix(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump grid fix ID");
            if (modify->fix[n]->per_grid_flag == 0)
                error->all(FLERR,"Dump grid fix does not compute per-grid info");
            if (argindex[j-START_ARGS] == 0 && modify->fix[n]->size_per_grid_cols != 0)
                error->all(FLERR,"Dump grid fix does not calculate a per-grid vector");
            if (argindex[j-START_ARGS] > 0 && modify->fix[n]->size_per_grid_cols == 0)
                error->all(FLERR,"Dump grid fix does not calculate per-grid array");
            if (argindex[j-START_ARGS] > 0 && argindex[j-START_ARGS] > modify->fix[n]->size_per_grid_cols)
                error->all(FLERR,"Dump grid fix array is accessed out-of-range");

            field2index[j-START_ARGS] = add_fix(suffix);
            delete [] suffix;

        // variable value = v_name
        } else if (strncmp(arg[j],"v_",2) == 0) {
            pack_func[j-START_ARGS] = &DumpParaview::pack_variable;
            vtype[j-START_ARGS] = DOUBLE;

            int n = strlen(arg[j]);
            char *suffix = new char[n];
            strcpy(suffix,&arg[j][2]);

            argindex[j-START_ARGS] = 0;

            n = input->variable->find(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump grid variable name");
            if (input->variable->grid_style(n) == 0)
                error->all(FLERR,"Dump grid variable is not grid-style variable");

            field2index[j-START_ARGS] = add_variable(suffix);
            delete [] suffix;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::init_style() {
    // find current ptr for each compute,fix,variable
    // check that fix frequency is acceptable
    int icompute;
    for (int i = 0; i < ncompute; i++) {
        icompute = modify->find_compute(id_compute[i]);
        if (icompute < 0)
            error->all(FLERR,"Could not find dump grid compute ID");
        compute[i] = modify->compute[icompute];
    }

    int ifix;
    for (int i = 0; i < nfix; i++) {
        ifix = modify->find_fix(id_fix[i]);
        if (ifix < 0) error->all(FLERR,"Could not find dump grid fix ID");
        fix[i] = modify->fix[ifix];
        if (nevery % modify->fix[ifix]->per_grid_freq)
            error->all(FLERR,"Dump grid and fix not computed at compatible times");
    }

    int ivariable;
    for (int i = 0; i < nvariable; i++) {
        ivariable = input->variable->find(id_variable[i]);
        if (ivariable < 0)
            error->all(FLERR,"Could not find dump grid variable name");
        variable[i] = ivariable;
    }

    // create cpart index of owned grid cells with particles in grid group
    reset_grid_count();
}

/* ---------------------------------------------------------------------- */

// need this for sparta but do nothing
void DumpParaview::write_header(bigint ndump) { return; }

/* ---------------------------------------------------------------------- */

int DumpParaview::count()
{
    // grow variable vbuf arrays if needed

    int nglocal = grid->nlocal;
    if (nglocal > maxgrid) {
        maxgrid = grid->maxlocal;
        for (int i = 0; i < nvariable; i++) {
            memory->destroy(vbuf[i]);
            memory->create(vbuf[i],maxgrid,"dump:vbuf");
        }
    }

    // invoke Computes for per-grid quantities
    if (ncompute) {
        for (int i = 0; i < ncompute; i++){
            if (!(compute[i]->invoked_flag & INVOKED_PER_GRID)) {
                compute[i]->compute_per_grid();
                compute[i]->invoked_flag |= INVOKED_PER_GRID;
            }
        }
    }

    // evaluate grid-style Variables for per-grid quantities
    if (nvariable) {
        for (int i = 0; i < nvariable; i++)
            input->variable->compute_grid(variable[i],vbuf[i],1,0);
    }

    // return # of grid cells with particles
    return ncpart;
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack() {
  for (int n = 0; n < size_one; n++) (this->*pack_func[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpParaview::write_data(int n, double *mybuf) {
  (this->*write_choice)(n,mybuf);
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParaview::add_compute(char *id) {
    int icompute;
    for (icompute = 0; icompute < ncompute; icompute++)
        if (strcmp(id,id_compute[icompute]) == 0) break;
    if (icompute < ncompute) return icompute;

    id_compute = (char **)
        memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
    delete [] compute;
    compute = new Compute*[ncompute+1];

    int n = strlen(id) + 1;
    id_compute[ncompute] = new char[n];
    strcpy(id_compute[ncompute],id);
    ncompute++;
    return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParaview::add_fix(char *id) {
    int ifix;
    for (ifix = 0; ifix < nfix; ifix++)
        if (strcmp(id,id_fix[ifix]) == 0) break;
    if (ifix < nfix) return ifix;

    id_fix = (char **)
        memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
    delete [] fix;
    fix = new Fix*[nfix+1];

    int n = strlen(id) + 1;
    id_fix[nfix] = new char[n];
    strcpy(id_fix[nfix],id);
    nfix++;
    return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParaview::add_variable(char *id) {
    int ivariable;
    for (ivariable = 0; ivariable < nvariable; ivariable++)
        if (strcmp(id,id_variable[ivariable]) == 0) break;
    if (ivariable < nvariable) return ivariable;

    id_variable = (char **)
        memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                        "dump:id_variable");
    delete [] variable;
    variable = new int[nvariable+1];
    delete [] vbuf;
    vbuf = new double*[nvariable+1];
    for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

    int n = strlen(id) + 1;
    id_variable[nvariable] = new char[n];
    strcpy(id_variable[nvariable],id);
    nvariable++;
    return nvariable-1;
}

/* ----------------------------------------------------------------------
   create cpart array to index owned grid cells with particles in grid group
   called by init or by any operation which changes the grid during a run
     e.g. fix balance, fix adapt, fix ablate
------------------------------------------------------------------------- */

void DumpParaview::reset_grid_count() {
    memory->destroy(cpart);
    int nglocal = grid->nlocal;
    memory->create(cpart,nglocal,"dump:cpart");

    Grid::ChildCell *cells = grid->cells;
    Grid::ChildInfo *cinfo = grid->cinfo;

    ncpart = 0;
    for (int i = 0; i < nglocal; i++) {
        if (!(cinfo[i].mask & groupbit)) continue;
        if (cells[i].nsplit > 1) continue;
        cpart[ncpart++] = i;
    }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint DumpParaview::memory_usage() {
    bigint bytes = Dump::memory_usage();
    bytes += memory->usage(cpart,grid->nlocal);
    return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpParaview::pack_compute(int n) {
    int index = argindex[n];
    Compute *c = compute[field2index[n]];

    // if one of post_process flags is set,
    //   invoke post_process_grid() or invoke post_process_tally()
    // else extract from compute's vector_grid and array_grid directly
    // dump buf only stores values for grid cells with particles
    //   use cpart indices to extract needed subset

    if (c->post_process_grid_flag)
        c->post_process_grid(index,1,NULL,NULL,NULL,1);
    else if (c->post_process_isurf_grid_flag)
        c->post_process_isurf_grid();

    if (index == 0 || c->post_process_grid_flag) {
        double *vector = c->vector_grid;
        for (int i = 0; i < ncpart; i++) {
            buf[n] = vector[cpart[i]];
            n += size_one;
        }
    } else {
        index--;
        double **array = c->array_grid;
        for (int i = 0; i < ncpart; i++) {
            buf[n] = array[cpart[i]][index];
            n += size_one;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack_fix(int n) {
    double *vector = fix[field2index[n]]->vector_grid;
    double **array = fix[field2index[n]]->array_grid;
    int index = argindex[n];

    if (index == 0) {
        for (int i = 0; i < ncpart; i++) {
            buf[n] = vector[cpart[i]];
            n += size_one;
        }
    } else {
        index--;
        for (int i = 0; i < ncpart; i++) {
            buf[n] = array[cpart[i]][index];
            n += size_one;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack_variable(int n) {
    double *vector = vbuf[field2index[n]];
    for (int i = 0; i < ncpart; i++) {
        buf[n] = vector[cpart[i]];
        n += size_one;
    }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump grid can output
   the grid property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpParaview::pack_id(int n) {
    Grid::ChildCell *cells = grid->cells;
    // NOTE: cellint (bigint) won't fit in double in some cases
    for (int i = 0; i < ncpart; i++) {
        buf[n] = cells[cpart[i]].id;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack_split(int n) {
    Grid::ChildCell *cells = grid->cells;
    for (int i = 0; i < ncpart; i++) {
        // convert to human readable format:
        //   split = 0: unsplit cell
        //   split = 1..N: split cell index + 1
        buf[n] = -cells[cpart[i]].nsplit + 1;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack_proc(int n) {
    Grid::ChildCell *cells = grid->cells;
    for (int i = 0; i < ncpart; i++) {
        buf[n] = cells[cpart[i]].proc;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParaview::pack_vol(int n) {
    Grid::ChildInfo *cinfo = grid->cinfo;
    for (int i = 0; i < ncpart; i++) {
        buf[n] = cinfo[cpart[i]].volume;
        n += size_one;
    }
}