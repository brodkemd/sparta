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

#define DUMP_SURF_ARGS_SIZE 7
#define COMPUTE_SURF_ARGS_SIZE 5
#define SURF_ARGS_SIZE 1
#define DUMP_FEA_MODIFY_ARGS_SIZE 2

#include "fix.h"
#include "error.h"
#include "write_surf.h"
#include "surf.h"
#include "modify.h"
#include "output.h"
#include "dump.h"

#include <string>
// #include <iostream>

#include <sys/stat.h>

namespace SPARTA_NS {
    /**
     * Class for this fix command
     */
    class FixFea : public Fix {
        public:
            /***
             * inputs, format INDEX : DESCRIPTION or VALUE
             * 0  : id
             * 1  : fea
             * 2  : nevery, execute every this many time steps
             * 3  : elmer_exe, path to the executable for elmer
             * 4  : sif_file, path the sif file to run with elmer
             * 5  : surf_file, path or name of file to dump surface to
             * 6  : compute-id, id for compute 
             * 7  : group-id, group ID for which surface elements to perform calculation on
             * 8  : mix-ID, mixture ID for particles to perform calculation on
             * 9  : dump-id, 
             * 10 : dump-file, file to dump surf info to
             * 11 : select-ID, what surface elements to dump
             * 12 : dump compute format, this is c_[*]
             */
            FixFea(class SPARTA *, int, char **);
            FixFea(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
            virtual ~FixFea();
            int setmask();
            void end_of_step();

        protected:
            // std::string command = "";
            char* surf_args[SURF_ARGS_SIZE];
            class WriteSurf* writer;
    };
}

#endif
#endif
