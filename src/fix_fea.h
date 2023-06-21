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

#include <string>
#include <vector>
#include <array>
#include <sys/stat.h>

namespace toml {
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
            virtual void start_of_step();

        protected:

            // template<typename dict>
            // void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

            // template<typename dict>
            // void run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

            void get_elmer(std::string& _buffer);

            // template<typename T>
            // void set(std::string _caller, std::string _name, toml::var_type_t _type, toml::node_t _val, T& _var);


            // toml::node_t resolve_name(toml::node_t _val, std::string _path);

            // toml::node_t get_var(toml::node_t _tbl, std::string _path, std::string _orig_path = "");
            // toml::node_t get_var(toml::table _tbl, std::string _path, std::string _orig_path = "");
                
        // public:
        //     void handle_emi(std::string _caller, toml::node_t val);
        //     void handle_tsurf_file(std::string _caller, toml::node_t val);
        //     void handle_both(std::string _caller, toml::node_t tbl);
        //     void handle_compute(std::string _caller, toml::node_t val);
        //     void handle_nevery(std::string _caller, toml::node_t val);
        //     void handle_groupID(std::string _caller, toml::node_t val);
        //     void handle_mixID(std::string _caller, toml::node_t val);
        //     void handle_customID(std::string _caller, toml::node_t val);
        //     void handle_sparta(std::string _caller, toml::node_t tbl);
        //     void handle_exe(std::string _caller, toml::node_t val);
        //     void handle_sif(std::string _caller, toml::node_t val);
        //     void handle_meshDBstem(std::string _caller, toml::node_t val);
        //     void handle_header(std::string _caller, toml::node_t tbl);
        //     void handle_simulation(std::string _caller, toml::node_t tbl);
        //     void handle_constants(std::string _caller, toml::node_t tbl);
        //     void handle_solver(std::string _caller, toml::node_t tbl);
        //     void handle_equation(std::string _caller, toml::node_t tbl);
        //     void handle_material(std::string _caller, toml::node_t tbl);
        //     void handle_body(std::string _caller, toml::node_t tbl);
        //     void handle_initial_condition(std::string _caller, toml::node_t tbl);
        //     void handle_boundary_condition(std::string _caller, toml::node_t tbl);
        //     void handle_elmer(std::string _caller, toml::node_t tbl);

        private:
            void load_boundary();
            void load_data();
            void load_sif(std::string sif_file);
            void print(std::string str, int num_indent = 1, std::string end = "\n");
            bool run_condition();


            // under "both"
            double emi;
            std::string tsurf_file;


            // under "sparta"
            int nevery;
            std::string groupID;
            std::string mixID;
            std::string customID;

            std::vector<std::string> compute_args;
            std::vector<std::string> surf_collide_args;
            std::vector<std::string> surf_modify_args;

            // under elmer
            class toml::Elmer* elmer;

            // Structure which would store the metadata
            std::string exe_path;
            std::string sif_path;
            std::string command;
            std::string meshDBstem;
            std::vector<double> data;
            std::vector<std::array<int, 4>> boundary_data;

            int run_every;
            int last_nlocal;
            int dimension;

            // surface temperature vars
            bool file_handler;
            int source,icompute,ifix,firstflag;
            int groupbit;
            int ngroup;
            int nprocs;

            int tindex,qwindex;

            char *id_qw;
            class Compute *cqw;
            class Fix *fqw;

            double prefactor,threshold;
            double *tvector_me;
            double *twall;
            double *qw_avg;

    };
}

// namespace toml {
//     typedef std::pair<std::string, void (SPARTA_NS::FixFea::*)(std::string, toml::value)> dict_item_t;
//     typedef std::vector<dict_item_t> dict_t;
// }


#endif
#endif
