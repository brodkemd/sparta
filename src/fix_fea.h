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

FixStyle(fea,FixFea)

#else

#ifndef SPARTA_FIX_FEA_H
#define SPARTA_FIX_FEA_H

#include "fix.h"
#include "hash3.h"




namespace elmer {
    class Elmer;
}

namespace SPARTA_NS {
    /**
     * Class for this fix command
     */
    class FixFea : public Fix {
        public:
            FixFea(class SPARTA *, int, char **);
            FixFea(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
            ~FixFea();
            int setmask();
            virtual void init();
            virtual void end_of_step();

        private:
            void load_temperatures();
            void move_surf();
            bigint remove_particles();
            void connect_3d_pre();
            void connect_3d_post();
            void print(const char* str, int num_indent = 1, const char* end = "\n");
            bool run_condition();

#include "hashlittle.h"
#include "hash_options.h"
#ifdef SPARTA_MAP
  typedef std::map<OnePoint3d,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#else
  typedef std::tr1::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#endif

            class Compute *cqw;
            MyHash *hash;
            class elmer::Elmer* elmer;
            int groupbit, nprocs, icol, tindex, run_every, nsurf, dimension, connectflag, *pselect;
            double emi, threshold, *qw_avg_me, *qw_avg;
            bigint ndeleted;
                   
    };
}

#endif
#endif
