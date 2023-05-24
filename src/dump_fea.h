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

#ifdef DUMP_CLASS

DumpStyle(fea,DumpFea)

#else

#ifndef SPARTA_DUMP_FEA_H
#define SPARTA_DUMP_FEA_H

#include "dump_surf.h"
#include "string.h"
#include "surf.h"
#include "compute.h"
#include "fix.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "input.h"

namespace SPARTA_NS {
    class DumpFea : public DumpSurf {
    public:
        DumpFea(class SPARTA *, int, char **);
        ~DumpFea();

    private:
        void write() override;
        // void openfile() override {};
    };
}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump surf attributes specified

Self-explanatory.

E: Invalid attribute in dump surf command

Self-explanatory.

E: Could not find dump surf compute ID

Self-explanatory.

E: Could not find dump surf fix ID

Self-explanatory.

E: Dump surf and fix not computed at compatible times

Fixes generate values on specific timesteps.  The dump surf output
does not match these timesteps.

E: Could not find dump surf variable name

Self-explanatory.

E: Invalid dump surf field for 2d simulation

Self-explanatory.

E: Dump surf compute does not compute per-surf info

Self-explanatory.

E: Dump surf compute does not calculate per-surf array

Self-explanatory.

E: Dump surf compute vector is accessed out-of-range

Self-explanatory.

E: Dump surf fix does not compute per-surf info

Self-explanatory.

E: Dump surf fix does not compute per-surf array

Self-explanatory.

E: Dump surf fix vector is accessed out-of-range

Self-explanatory.

E: Dump surf variable is not surf-style variable

Self-explanatory.

*/
