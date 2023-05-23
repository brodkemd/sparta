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
#include "write_surf.h"

#include <string>
#include <stdexcept>
#include <array>
#include <sys/stat.h>

struct CommandResult {
    std::string output;
    int exitstatus;
};

namespace SPARTA_NS {
    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return commandResult containing STDOUT (not stderr) output & exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    static CommandResult EXEC(const std::string &command) {
        int exitcode = 0;
        std::array<char, 1048576> buffer {};
        std::string result;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result += std::string(buffer.data(), bytesread);
            }
        } catch (...) {
            pclose(pipe);
            throw;
        }
        exitcode = WEXITSTATUS(pclose(pipe));
        return CommandResult{result, exitcode};
    }

    /**
     * Class for this fix command
    */
    class FixFea : public Fix {
        public:
            FixFea(class SPARTA *, int, char **);
            FixFea(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
            virtual ~FixFea();
            int setmask();
            void end_of_step();


        protected:
            std::string command = "";
            char* surf_path;
            WriteSurf* writer;
            char* surf_id;
    };

}

#endif
#endif
