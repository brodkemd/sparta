/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_CUT2D_H
#define SPARTA_CUT2D_H

#include "pointers.h"
#include "myvec.h"
#include "mylist.h"
#include "pool.h"

namespace SPARTA_NS {

class Cut2d : protected Pointers {
 public:
  // TEMP for VERBOSE output
  int icell;

  Cut2d(class SPARTA *);
  ~Cut2d();
  void surf2grid();
  void split();

 private:
  int nsurf;
  MyVec<int> used;
  MyVec<int> startpts;
  MyVec<int> endpts;
  
  struct PLone {
    int iline;
    PLone *prev,*next;
  };

  struct Cpt {
    int flag;
    double x[2];
    double dot;
    int ipl,oindex;
    Cpt *prev,*next;
  };

  struct Opt {
    int flag;
    double x[2];
    int cindex;
  };

  struct Entrypt {
    int iopt;
    int index;
    Entrypt *prev,*next;
  };
  
  struct Loop {
    int flag;
    double area;
    MyVec<int> lines;
  };

  struct PG {
    double area;
    MyVec<int> lines;
  };

  Pool<PLone> ppool;
  MyVec< MyList<PLone*> > pl;

  MyVec< MyVec<Opt> > opts;

  Pool<Cpt> cpool;
  MyList<Cpt*> cpts;
  Cpt *cindex[5];

  Pool<Entrypt> epool;
  MyList<Entrypt*> entrypts;

  MyVec<Loop> loops;
  MyVec<PG> pg;

  void line2pl(int, int *);
  void weiler_intersect(double *, double *, int *);
  void interleave(double *, double *, double *, double *, int, int, int);
  void weiler_walk(double *, double *, int *);
  void loop2pg(double *, double *, int *);
  void surf2pg(int, int *, int *);

  int cliptest(double *, double *, double *, double *);
  void clip(double *, double *, double *, double *, double *, double *);
  int sameborder(double *pt1, double *pt2, double *lo, double *hi);
  int corner(double *, double *, double *);
  int ptflag(double *pt, double *lo, double *hi);
};

}

#endif

/* ERROR/WARNING messages:

E: Bad grid of processors for create_grid

UNDOCUMENTED

E: Per-processor grid count is too big

UNDOCUMENTED

*/