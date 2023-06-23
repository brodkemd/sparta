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

#include <string>
#include <vector>
#include <array>

#include "UTIL/elmer_definitions.h"


namespace elmer {
    class Elmer;
}

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
            ~FixFea();
            int setmask();
            virtual void init();
            virtual void end_of_step();
            // virtual void start_of_step();

        private:
            void get_elmer(std::string& _buffer);
            void load_temperatures();
            void load_boundary();
            void setup_data_file();
            // void load_data();
            // void load_sif(std::string sif_file);
            void print(std::string str, int num_indent = 1, std::string end = "\n");
            bool run_condition();

            // under elmer
            class elmer::Elmer* elmer;

            // Structure which would store the metadata
            std::string meshDBstem, temperature_data_file;

            std::vector<std::array<int, elmer::boundary_size>> boundary_data;

            int groupbit, ngroup, nprocs, tindex, qwindex, run_every, last_nlocal, dimension, firstflag;
            
            double emi, prefactor, threshold, *tvector_me, *twall, *qw_avg;

            class Compute *cqw;
            class Fix *fqw;
    };
}

#endif
#endif
