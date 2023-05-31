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

#ifdef FIX_CLASS

FixStyle(surf/temp/dynamic,FixSurfTempDynamic)

#else

#ifndef SPARTA_FIX_SURF_TEMP_DYNAMIC_H
#define SPARTA_FIX_SURF_TEMP_DYNAMIC_H

#include "fix.h"
#include "surf.h"
#include <sys/stat.h>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iostream>

namespace SPARTA_NS {
    class FixSurfTempDynamic : public Fix {
        public:
            FixSurfTempDynamic(class SPARTA *, int, char **);
            virtual ~FixSurfTempDynamic();
            int setmask();
            virtual void init();
            virtual void end_of_step();
            virtual void start_of_step();

        protected:
            void load_surf_temps(int _nlocal);
            void write_elmer_file();

        private:
            int file_format_flag;

            bool file_handler;
            int source,icompute,ifix,firstflag;
            int groupbit;
            int nprocs;
            double emi;
            int tindex,qwindex;

            char *id_qw;
            char *twall_file;
            class Compute *cqw;
            class Fix *fqw;

            double prefactor,threshold;
            double *tvector_me;
            double *twall;
    };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use non-rcb fix balance with a grid cutoff

This is because the load-balancing will generate a partitioning
of cells to processors that is dispersed and which will not work
with a grid cutoff >= 0.0.

*/
