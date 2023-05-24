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
#include "fix_fea.h"

#include <iostream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <array>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h> 

struct CommandResult {
    std::string output;
    int exitstatus;
};

namespace SPARTA_NS {
    class DumpFea : public DumpSurf {
        public:
            DumpFea(class SPARTA *, int, char **);
            ~DumpFea();

            void modify_params(int narg, char ** arg) override;

        private:
            void write() override;
            std::string command = "";
    };
}

#endif
#endif