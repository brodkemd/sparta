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

#include "fix_fea.h"
#include "error.h"
#include "surf.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg)
{
  /***
   * Format of this command
  */

  if (narg > 5)
    error->all(FLERR,"Illegal fix fea command, too many inputs, format \"fix id fea EXE FILE\"");
  else if (narg < 5)
    error->all(FLERR,"Illegal fix fea command, too few inputs, format \"fix id fea EXE FILE\"");

  if (!surf->exist)
    error->all(FLERR,"Illegal fix fea command, no surface to analyze");

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix fea command, nevery <= 0");

  // std::string msg  = "";
  // for (int i = 0; i < narg; i++) msg += (std::string(arg[i]) + " ");
  // error->message(FLERR, ((std::string)"FEA args: " + msg).c_str());
  // CommandResult command_result = EXEC("pwd" + this->command_end);

  const char* exe_path = arg[3];
  const char* file_path = arg[4];

  // Structure which would store the metadata
  struct stat sb;

  // Calls the function with path as argument
  // If the file/directory exists at the path returns 0
  // If block executes if path exists
  if (!(stat(exe_path, &sb) == 0)) error->all(FLERR,"Illegal fix fea command, exe path does not exist");
  if (!(stat(file_path, &sb) == 0)) error->all(FLERR,"Illegal fix fea command, file path does not exist");

  // must have " 2>&1" at end to pipe stderr to stdout
  this->command = std::string(exe_path) + " " + std::string(file_path) + " 2>&1";
}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea()
{
  return;
}

/* ---------------------------------------------------------------------- */

int FixFea::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


void FixFea::end_of_step()
{
  error->message(FLERR, "Running FEA");
  CommandResult command_result = EXEC(this->command);

  if (command_result.exitstatus) {
    fprintf(logfile, command_result.output.c_str());
    error->all(FLERR, "fix fea failed, see sparta log file");
  }
}
/* ---------------------------------------------------------------------- */

// void FixFea::init()
// {
//   if (maxion != particle->nspecies)
//     error->all(FLERR,"Number of particle species has changed since "
//                "fix fea was specified");
// }

// /* ----------------------------------------------------------------------
//    called when a particle with index is created
//     or when temperature dependent properties need to be updated
//    creation used temp_thermal and vstream to set particle velocity
//    if an ion, set ionambi and velambi for particle
// ------------------------------------------------------------------------- */

// void FixFea::update_custom(int index, double temp_thermal,
//                                 double, double,
//                                 double *vstream)
// {
//   int *ionambi = particle->eivec[particle->ewhich[ionindex]];
//   double **velambi = particle->edarray[particle->ewhich[velindex]];

//   // if species is not fea ion, set ionambi off and return

//   int ispecies = particle->particles[index].ispecies;

//   if (ions[ispecies] == 0) {
//     ionambi[index] = 0;
//     return;
//   }

//   // set velocity of electron
//   // based on electron mass, thermal temperature, and streaming velocity

//   ionambi[index] = 1;

//   double vscale = sqrt(2.0 * update->boltz * temp_thermal /
//                        particle->species[especies].mass);

//   double vn = vscale * sqrt(-log(random->uniform()));
//   double vr = vscale * sqrt(-log(random->uniform()));
//   double theta1 = MY_2PI * random->uniform();
//   double theta2 = MY_2PI * random->uniform();

//   velambi[index][0] = vstream[0] + vn*cos(theta1);
//   velambi[index][1] = vstream[1] + vr*cos(theta2);
//   velambi[index][2] = vstream[2] + vr*sin(theta2);
// }

// /* ----------------------------------------------------------------------
//    called when a surface reaction occurs
//    iorig = particle I before reaction
//    I,J = indices of two particles after reaction
//          either can be -1, meaning particle does not exist
// ------------------------------------------------------------------------- */

// void FixFea::surf_react(Particle::OnePart *iorig, int &i, int &j)
// {
//   int ispecies = iorig->ispecies;

//   // recombination reaction, just return

//   if (i < 0) return;

//   // exchange reaction
//   // if ion -> non-ion, unset ionambi flag

//   if (j < 0) {
//     if (ions[ispecies] == 0) return;
//     Particle::OnePart *particles = particle->particles;
//     if (ions[particles[i].ispecies] == 1) return;
//     int *ionambi = particle->eivec[particle->ewhich[ionindex]];
//     ionambi[i] = 0;
//   }

//   // dissociation reaction
//   // if non-ion -> ion + electron, create an fea ion
//   // use global temp_thermal and vstream for electron creation
//   // set j = -1 to delete electron that was just created by caller

//   else {
//     if (ions[ispecies] == 1) return;
//     Particle::OnePart *particles = particle->particles;
//     if (ions[particles[i].ispecies] == 0) return;
//     if (particles[j].ispecies != especies) return;
//     update_custom(i,update->temp_thermal,update->temp_thermal,
//                  update->temp_thermal,update->vstream);
//     j = -1;
//   }
// }
