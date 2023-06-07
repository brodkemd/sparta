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
#include <sstream>

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

        protected:
            class Elmer {
                public:
                    Elmer(std::string output_file, Error*& _error);
                    void write(std::string header = "");
                    void Boundary_Condition(int _n, std::vector<std::string> args);
                    void Initial_Condition(int _n, std::vector<std::string> args);


                private:
                    std::string sif_path = "";
                    std::string tab = "  ";
                    std::vector<std::vector<std::string>> commands;
                    std::vector<std::string> v;
                    Error* error;

                    void _add_section(std::string _name, std::string _n, std::vector<std::string> args);
            };

            class ConfigParser {
                public:
                    ConfigParser(std::string _file_name, Error*& _error, bool ignore_var_name_case = true);
                    std::size_t size();
                    std::pair<std::string, std::string> &operator[](int i);

                private:
                    Error* error;
                    std::string file_name;
                    std::vector<std::pair<std::string, std::string>> contents;
                    
                    void read_file(bool ignore_var_name_case = true);
            };

        private:
            void load_boundary();
            void load_data();
            void load_sif(std::string sif_file);
            void print(std::string str, int num_indent = 1, std::string end = "\n");
            bool run_condition();

            // int debug_num = 1;
            // int compute_index;

            std::string exe_path;
            std::string sif_path;

            std::string command;
            std::string tsurf_file;
            std::string meshDBstem;
            //std::string data_file;
            std::string customID;
            std::string groupID;
            std::string mixID;
            std::string sif_format;
            std::vector<double> data;
            std::vector<std::array<int, 4>> boundary_data;

            int run_every;
            int last_nlocal;

            // surface temperature vars
            bool file_handler;
            int source,icompute,ifix,firstflag;
            int groupbit;
            int nprocs;
            double emi;
            int tindex,qwindex;

            char *id_qw;
            class Compute *cqw;
            class Fix *fqw;

            double prefactor,threshold;
            double *tvector_me;
            double *twall;
            double *qw_avg;

            Elmer* elmer;

    };
}

#endif
#endif
