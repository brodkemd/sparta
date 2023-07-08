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
#include "hashlittle.h"

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
            virtual double compute_scalar();

        private:
            void update_temperatures();
            void update_surf();
            bigint remove_particles();
            void connect_3d_pre();
            void connect_3d_post();
            bool run_condition();
            bool checkVarSums(std::string& _name);

            #include "hash_options.h"

            #ifdef SPARTA_MAP
              typedef std::map<OnePoint3d,int> MyHash;
            #elif SPARTA_UNORDERED_MAP
              typedef std::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
            #else
              typedef std::tr1::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
            #endif

            MyHash *hash;
            class Compute *cqw;
            class elmer::Elmer* fea;

            int groupbit, nprocs, tindex, run_every, nsurf, dimension, connectflag, *pselect, shear_locs[3], force_locs[3], energy_loc;
            double energy_threshold, force_threshold, shear_threshold, *qw_avg_me, *qw_avg, *fx_avg_me, *fx_avg, *fy_avg_me, *fy_avg, *fz_avg_me, *fz_avg, *shx_avg_me, *shx_avg, *shy_avg_me, *shy_avg, *shz_avg_me, *shz_avg;
            bigint ndeleted;
                   
    };
}

#endif
#endif
