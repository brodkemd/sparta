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

// #define DUMP_SURF_ARGS_SIZE 7
#define COMPUTE_SURF_ARGS_SIZE 5
// #define SURF_ARGS_SIZE 1
// #define DUMP_FEA_MODIFY_ARGS_SIZE 2
#define BOUNDARY_DATA_SIZE 4
#define BOUNDARY_LINE_SIZE 8

#include "fix.h"
#include "error.h"
#include "write_surf.h"
#include "surf.h"
#include "modify.h"
#include "output.h"
#include "dump.h"
#include "compute.h"
#include "comm.h"
#include "memory.h"
#include "domain.h"
#include "input.h"
#include "update.h"

#include <string>
#include <vector>
#include <array>

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
            virtual void start_of_step();

            // void debug_msg() {
            //     std::cout << "---> " << this->debug_num << "\n";
            //     this->debug_num++;
            // }

        private:
            void load_boundary();
            void load_data();
            void load_sif(std::string sif_file);
            void print(std::string str, bool indent = true, std::string end = "\n");

            // int debug_num = 1;
            // int compute_index;
            std::string command;

            std::string file_stem;
            std::string data_file;
            std::string sif_format;
            std::vector<double> data;
            std::vector<std::array<double, 4>> boundary_data;


            // surface temperature vars
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
            // std::string command = "";
            // char* surf_args[SURF_ARGS_SIZE];
            // class WriteSurf* writer;
    };
}

#endif
#endif
