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

FixStyle(paraview,FixParaview)

#else

#ifndef SPARTA_FIX_PARAVIEW_H
#define SPARTA_FIX_PARAVIEW_H

#include "fix.h"

namespace SPARTA_NS {
    /**
     * Class for this fix command
     */
    class FixParaview : public Fix {
        public:
            FixParaview(class SPARTA *, int, char **);
            FixParaview(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
            ~FixParaview();
            int setmask();
            void init();
            void end_of_step();
        
        private:
            Surf::Line *lines;
            Surf::Tri *tris;
            
            char **ids, **data_names, *prefix;

            typedef void (FixParaview::*FnPtrPack)(int, index_t*, void*);
            FnPtrPack *funcs;
            int *argindex, groupbit, *func_types, style;
            index_t num_expanded_args, *cell_index_map, num_cells, num_points, total_num_cells;
            double *vbuf;
            // tells the code to dump before the first timestep
            const bool dump_at_0 = true;
            bool there_is_variable, collect_into_one_file;

            void collectDataFromProcesses(index_t step, index_t num, double *&points, index_t &num_points, index_t *&connections, index_t &num_cells, index_t num_elements_per_cell_num, void *&data, index_t *&offsets, index_t &offsets_last_index);

            void handleGridArgs(int argc, char** argv);
            void handleSurfArgs(int argc, char** argv);

            void pointSwitch(index_t _index, index_t _element_index, double arr[3]);

            /**
             * The rest handle data
            */

            // per grid values
            void handleGridId(      int index_in_offsets, index_t* offsets, void* data);
            void handleGridSplit(   int index_in_offsets, index_t* offsets, void* data);
            void handleGridProc(    int index_in_offsets, index_t* offsets, void* data);
            void handleGridVol(     int index_in_offsets, index_t* offsets, void* data);
            void handleGridCompute( int index_in_offsets, index_t *offsets, void* data);
            void handleGridFix(     int index_in_offsets, index_t *offsets, void* data);
            void handleGridVariable(int index_in_offsets, index_t *offsets, void* data);

            // per surf values
            void handleSurfId(      int index_in_offsets, index_t* offsets, void* data);
            void handleSurfType(    int index_in_offsets, index_t* offsets, void* data);
            void handleSurfArea(    int index_in_offsets, index_t* offsets, void* data);
            void handleSurfCustom(  int index_in_offsets, index_t* offsets, void* data);
            void handleSurfCompute( int index_in_offsets, index_t* offsets, void* data);
            void handleSurfFix(     int index_in_offsets, index_t* offsets, void* data);
            void handleSurfVariable(int index_in_offsets, index_t* offsets, void* data);
    };
}

#endif
#endif
