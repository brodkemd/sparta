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

/* ----------------------------------------------------------------------
   This file was contributed by Marek Brodke, a student at the University
   of Cincinnati
------------------------------------------------------------------------- */

#include "fix_fea.h"
#include "stdlib.h"
#include "string.h"

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
#include "particle.h"
#include "mixture.h"
#include "update.h"

// #include "UTIL/toml.h"
#include "UTIL/elmer.h"

#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>

using namespace SPARTA_NS;

#define SB_SI  5.670374419e-8
#define SB_CGS 5.670374419e-5


enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};



/* ----------------------------------------------------------------------
    
    FixFea class methods

 ---------------------------------------------------------------------- */


/**
 * NOTE: Make sure that units are consistent, if you use si in sparta, make sure you set units
*/
FixFea::FixFea(SPARTA *sparta, int narg, char **arg) : Fix(sparta, narg, arg) {
    std::string groupID;
    std::string mixID;
    std::string customID;

    std::vector<std::string> compute_args;
    std::vector<std::string> surf_collide_args;
    std::vector<std::string> surf_modify_args;

    this->print("Setting up fix fea:", 0);

    // checking if the number of args is correct
    if (narg > 3)
        error->all(FLERR,"Illegal fix fea command, too many inputs");
    else if (narg < 3)
        error->all(FLERR,"Illegal fix fea command, too few inputs");

    // making sure there is a surface to analyze
    if (!surf->exist)
        error->all(FLERR,"Illegal fix fea command, no surface to analyze");

    // making sure the surface is the correct type
    if (surf->implicit)
        error->all(FLERR,"Cannot use fix fea with implicit surfs");
    if (surf->distributed)
        error->all(FLERR,"Cannot use fix fea with distributed surfs");

    // Structure which would store the metadata
    struct stat sb;

    if (!(stat(arg[2],  &sb) == 0))
        error->all(FLERR,"Illegal fix fea command, toml config file does not exist");
    this->print("Loading config from: " + std::string(arg[2]));

    this->elmer = new elmer::Elmer();

    try {
        toml::handler s(arg[2]);

        s.get_at_path(this->emi,               "both.emi");
        s.get_at_path(this->temperature_data_file,        "both.temperature_data_file");
        s.get_at_path(this->run_every,         "sparta.run_every");
        s.get_at_path(this->nevery,            "sparta.nevery");
        s.get_at_path(groupID,           "sparta.groupID");
        s.get_at_path(mixID,             "sparta.mixID");
        s.get_at_path(customID,          "sparta.customID");
        s.get_at_path(this->qwindex,           "sparta.qwindex");

        s.get_at_path(this->meshDBstem,        "both.meshDBstem");
        s.get_at_path(this->elmer->exe,        "elmer.exe");
        s.get_at_path(this->elmer->sif,        "elmer.sif");
        s.get_at_path(this->elmer->meshDBstem, "both.meshDBstem");
        s.get_at_path(this->elmer->base_temp,  "elmer.base_temp");

        s.get_at_path(compute_args,      "sparta.compute");
        s.get_at_path(surf_collide_args, "sparta.surf_collide");
        s.get_at_path(surf_modify_args,  "sparta.surf_modify");

        elmer->header.set(s);
        elmer->simulation.set(s);
        elmer->constants.set(s);
        elmer->equation.set(s);
        elmer->material.set(s);
        elmer->solver.set(s);
    
    } catch (std::string _msg) {
        error->all(FLERR, _msg.c_str()); 
    } catch (std::exception& e) {
        error->all(FLERR, e.what());
    } catch (...) {
        error->all(FLERR, "unidentified error occurred while setting variables");
    }

    if (emi <= 0.0 || emi > 1.0)
        error->all(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");

    this->print("emi = " + std::to_string(this->emi));


    if (!(stat(this->temperature_data_file.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, temperature_data_file path does not exist");

    this->print("temperature_data_file = " + this->temperature_data_file);


    if (this->run_every <= 0)
        error->all(FLERR,"Illegal fix fea command, run_every <= 0");

    this->print("run_every = " + std::to_string(this->run_every));

    if (this->nevery <= 0)
        error->all(FLERR,"Illegal fix fea command, nevery <= 0");
    
    if (this->run_every % this->nevery != 0)
        error->all(FLERR, "Illegal fix fea command, run_every must be a multiple of nevery");

    this->print("nevery = " + std::to_string(this->nevery));


    // get the surface group
    int igroup = surf->find_group(groupID.c_str());
    if (igroup < 0)
        error->all(FLERR,"Fix fea group ID does not exist");
    
    groupbit = surf->bitmask[igroup];
    this->print("groupID = " + groupID);


    int imix = particle->find_mixture((char*)mixID.c_str());
    if (imix < 0)
        error->all(FLERR,"Compute thermal/grid mixture ID does not exist");

    ngroup = particle->mixture[imix]->ngroup;
    this->print("mixID = " + mixID);
    

    this->tindex = surf->add_custom((char*)customID.c_str(), DOUBLE, 0);
    this->print("customID = " + customID);


    if (!(stat(this->elmer->exe.c_str(),  &sb) == 0))
        error->all(FLERR, "Illegal fix fea command, exe path does not exist");
    
    this->print("exe = " + this->elmer->exe);


    // if (!(stat(this->elmer->sif.c_str(),  &sb) == 0))
    //     error->all(FLERR, "Illegal fix fea command, sif path does not exist");
    
    // this->print("sif = " + this->elmer->sif);
    

    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer->meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error->all(FLERR, ("Illegal fix fea command, mesh database incomplete, " + (this->elmer->meshDBstem + "." + exts[i]) + " does not exist").c_str());
    }
    this->print("meshDBstem = " + this->elmer->meshDBstem);


    char** arr;
    int size;
    
    size = toml::vec_to_arr(compute_args, arr);

    // adding the needed compute
    modify->add_compute(size, arr);

    // the compute index, it was just made so it is the number of computes minus 1 because it is an index
    int icompute = modify->ncompute - 1;// modify->find_compute(id_qw);

    // error checks
    cqw = modify->compute[icompute];
    if (icompute < 0)
        error->all(FLERR,"Could not find fix fea compute ID");
    if (cqw->per_surf_flag == 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf info");
    if (qwindex == 0 && cqw->size_per_surf_cols > 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf vector");
    if (qwindex > 0 && cqw->size_per_surf_cols == 0)
        error->all(FLERR,"Fix fea compute does not compute per-surf array");
    if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
        error->all(FLERR,"Fix fea compute array is accessed out-of-range");

    this->print("added surf compute with id: " + std::string(modify->compute[icompute]->id));

    // prefactor and threshold in Stefan/Boltzmann equation
    // units of prefactor (SI) is K^4 / (watt - m^2)
    // same in 3d vs 2d, since SPARTA treats 2d cell volume as 1 m in z
    this->dimension = domain->dimension;

    // setting variables based on unit system
    if (strcmp(update->unit_style, "si") == 0) {
        prefactor = 1.0 / (emi * SB_SI);
        threshold = 1.0e-6;
    } else if (strcmp(update->unit_style, "cgs") == 0) {
        prefactor = 1.0 / (emi * SB_CGS);
        threshold = 1.0e-3;
    }

    // getting number of processes
    MPI_Comm_size(world, &nprocs);

    surf_collide_args[2] = "s_" + surf_collide_args[2];
    size = toml::vec_to_arr(surf_collide_args, arr);

    // adding surf collision model needed
    this->surf->add_collide(size, arr);

    size = toml::vec_to_arr(surf_modify_args, arr);

    // modifying params
    this->surf->modify_params(size, arr);
   
    // temporary
    delete [] arr;

    // try {
    //     std::string _temp;
    //     this->elmer->join(_temp);

    //     elmer::writeFile("test.out", _temp);
    // } catch (std::string _msg) {
    //     error->all(FLERR, _msg.c_str());
    // }

    this->elmer->simulation.Timestep_intervals = {this->run_every};
    this->elmer->simulation.Output_Intervals = {this->run_every};
    this->elmer->simulation.Output_File = this->temperature_data_file;

    firstflag = 1;

}

/* ---------------------------------------------------------------------- */

FixFea::~FixFea() {
    delete this->elmer;
    memory->destroy(tvector_me);
    memory->destroy(twall);
    memory->destroy(qw_avg);
    surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

/**
 * sets the mask to make the class run at the start and end of each timestep 
 */
int FixFea::setmask() {
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
 }// | START_OF_STEP; }

/* ---------------------------------------------------------------------- */

void FixFea::load_temperatures() {
    int nlocal = surf->nlocal;

    if (comm->me == 0) {

        std::vector<double> data;
        try {
            elmer::getLatestNodeData(this->temperature_data_file, data);
        } catch (std::string _msg) {
            error->all(FLERR, _msg.c_str()); 
        } catch (std::exception& e) {
            error->all(FLERR, e.what());
        }

        if (nlocal != this->boundary_data.size())
            error->all(FLERR, ("boundary data does not match required size, required size " + std::to_string(nlocal) + ", size " + std::to_string(this->boundary_data.size())).c_str());

        // used later
        double avg;
        
        // averages values for the nodes of a surface element and sets this average to the 
        // temperature of the surface element
        for (int i = 0; i < this->boundary_data.size(); i++) {
            // computes the average temperature of the nodes that make up the surface element
            // this value is used to set the surface element temperature
            avg = 0;
            for (int j = 1; j < elmer::boundary_size; j++) {
                // gets the data point corresponding to node id and adds it to the rolling sum
                avg += data[this->boundary_data[i][j]-1];
            }
            // computing the average by dividing the sum by the number of points and setting the
            // surface element
            this->twall[this->boundary_data[i][0]-1] = avg/(elmer::boundary_size - 1);
        }

        // checking to make sure all values where set
        for (int i = 0; i < nlocal; i++) {
            if (this->twall[i] == (double)0)
                error->all(FLERR, "wall temperature not set correctly");
        }

        this->elmer->node_temperature_data = data;
        // elmer::eraseFile(this->temperature_data_file);

        // // sending to other processes
        // for (int i = 1; i < nprocs; i++) {
        //     MPI_Send(&twall, nlocal, MPI_DOUBLE, i, 1, world);
        // }
        
    }
    // } else {
    //     // master process is 0, blocks until data is received
    //     MPI_Recv(&twall, nlocal, MPI_DOUBLE, 0, 1, world, MPI_STATUS_IGNORE);
    // }
    double *tvector = surf->edvec[surf->ewhich[tindex]];
    MPI_Barrier(world);
    MPI_Allreduce(this->twall, tvector, nlocal, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(this->twall, this->twall, nlocal, MPI_DOUBLE, MPI_SUM, world);
}


void FixFea::load_boundary() {
    if (comm->me == 0) {
        this->print("Loading boundary data");
        try {
            elmer::getBoundaryData(this->meshDBstem + ".boundary", this->boundary_data);
        } catch (std::string _msg) {
            error->all(FLERR, _msg.c_str()); 
        } catch (std::exception& e) {
            error->all(FLERR, e.what());
        }
    }
}



void FixFea::setup_data_file() {
    int count = elmer::count_lines_in_file(this->meshDBstem + ".nodes");

    std::ostringstream double_converter;
    double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << this->elmer->base_temp;

    std::string base_temp_str = double_converter.str();

    std::ofstream output(this->temperature_data_file);

    for (int i = 0; i < count; i++) 
        output << base_temp_str << "\n";

    output.close();
}

void FixFea::init() {
    if (!firstflag) return;
    firstflag = 0;

    double *tvector = surf->edvec[surf->ewhich[tindex]];
    int nlocal = surf->nlocal;
    this->last_nlocal = nlocal;

    memory->create(twall, nlocal,"fea:twall");
    memset(this->twall, 0, nlocal*sizeof(double));

    memory->create(this->qw_avg, nlocal, "fea:qw_avg");
    memset(this->qw_avg, 0, nlocal*sizeof(double));

    // allocate per-surf vector for explicit all surfs
    memory->create(tvector_me,nlocal,"fea:tvector_me");
    memset(tvector_me, 0, nlocal*sizeof(double));

    // loading the boundary data
    this->setup_data_file();
    this->load_boundary();
    this->elmer->load_elements();
    this->load_temperatures();

    // this->start_of_step(true);
    for (int i = 0; i < nlocal; i++) tvector[i] = twall[i];
    this->print("done initing");
}

/* ---------------------------------------------------------------------- */

bool FixFea::run_condition() {
    return update->ntimestep % this->run_every == 0;
}


/* ---------------------------------------------------------------------- */

/**
 * Runs at the end of each time step
 * writing the surface data to the file
 */
void FixFea::end_of_step() {
    int i,m,mask;
    bool go = false;
    double qw;

    // access source compute or fix
    // set new temperature via Stefan-Boltzmann eq for nown surfs I own
    // use Twall if surf is not in surf group or eng flux is too small
    // compute/fix output is just my nown surfs, indexed by M
    // store in tvector_me = all nlocal surfs, indexed by I
    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;

    int nlocal = surf->nlocal;
    if (this->last_nlocal != nlocal)
        error->all(FLERR, "detected surface change, this is not allowed with fix fea command");

    double **array;
    cqw->post_process_surf();
    array = cqw->array_surf;

    int icol = qwindex-1;

    m = 0;
    for (i = comm->me; i < nlocal; i += nprocs) {
        if (dimension == 3) mask = tris[i].mask;
        else mask = lines[i].mask;

        if (mask & groupbit) {
            qw = array[m][icol];
            if (qw > threshold) this->qw_avg[i] += qw/this->run_every;
        }
        m++;
    }

    if (this->run_condition()) {
        
        this->print("running");
        MPI_Allreduce(this->qw_avg, this->qw_avg, nlocal, MPI_DOUBLE, MPI_SUM, world);

        if (comm->me == 0) {
            
            int i;
            for (i = 0; i < nlocal; i++) {
                if (qw_avg[i] != (double)0.0) break;
            }
            go = (i == (nlocal - 1));
            if (go) {

                this->print("in if");
                this->elmer->simulation.Timestep_Sizes = {update->dt};

                this->elmer->bodys.clear();
                this->elmer->initial_conditions.clear();
                this->elmer->boundary_conditions.clear();

                this->print("after clear");

                for (int i = 0; i < nlocal; i++) {
                    elmer::Boundary_Condition* bc = new elmer::Boundary_Condition(i+1);
                    bc->Heat_Flux = this->qw_avg[i];
                    this->elmer->boundary_conditions.push_back(bc);
                }
                this->print("after for");
                

                try {
                    this->elmer->createInitCondsFromData();
                    this->elmer->makeSif();
                    this->print("done in if");
                    this->elmer->run();
                } catch (std::string& _msg ){
                    error->all(FLERR, _msg.c_str());
                } catch (std::exception& e) {
                    error->all(FLERR, e.what());
                }

            // error->all(FLERR, "done");
            } else {
                this->print("Skipping elmer, no flux detected");
            }

        }
        MPI_Barrier(world);
        // error->all(FLERR, "generated everything");



        // // m = 0;
        // // for (i = me; i < nlocal; i += nprocs) {
        // //     if (dimension == 3) mask = tris[i].mask;
        // //     else mask = lines[i].mask;
        // //     if (!(mask & groupbit)) tvector_me[i] = twall[i];
        // //     else {
        // //         qw = qw_avg[i];
        // //         if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
        // //         else tvector_me[i] = twall[i];
        // //     }
        // //     m++;
        // // }

        // // Allreduce tvector_me with my owned surfs to tvector custom variable
        // // so that all procs know new temperature of all surfs
        // // NOTE: could possibly just Allreduce a vector size of surface group
        // // NOTE: all data is put into tvector
        // double *tvector = surf->edvec[surf->ewhich[tindex]];
        // MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);

        // MPI_Barrier(world);
        // if (comm->me == 0) {
        //     this->print("Generating Boundary Conditions", 4);

        //     this->elmer->boundary_conditions.clear();
        //     this->elmer.initial_conditions.clear();

        //     for (int i = 0; i < nlocal; i++) {
        //         this->elmer->Boundary_Condition(i+1, {
        //             "Target Boundaries(1) = " + std::to_string(this->boundary_data[i][0]-1),
        //             "Temperature = " + std::to_string(tvector[this->boundary_data[i][0]-1])
        //         });
        //     }
        //     this->print("writing to elemer file", 4);
        //     this->elmer->write(this->sif_format);
        //     // this->print("done writing to elemer file", 4);


        //     // // running the command
        //     // error->all(FLERR, "would have run command");
        //     // CommandResult command_result = EXEC(this->command);
            
        //     // // if the command did not succeed
        //     // if (command_result.exitstatus) {
        //     //     // writing to the logfile/screen and erroring out
        //     //     if (logfile) {
        //     //         fprintf(logfile, command_result.output.c_str());
        //     //         error->all(FLERR, "fix fea failed, see sparta log file");
        //     //     } else if (screen) {
        //     //         fprintf(screen, command_result.output.c_str());
        //     //         error->all(FLERR, "fea exited with the above error");
        //     //     } else error->all(FLERR, "no place to output command error to");
        //     // }
            
        // }
        memset(this->qw_avg, 0, nlocal*sizeof(double));
        if (go) this->load_temperatures();
        MPI_Barrier(world);
    }
}

/* ---------------------------------------------------------------------- */


/**
 * custom printing for this class
 */
void FixFea::print(std::string str, int num_indent, std::string end) {
    std::string space = "";
    for (int i = 0; i < num_indent; i++)
        space += "  ";

    if (comm->me == 0) {
        if (screen)  fprintf(screen,  "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        if (logfile) fprintf(logfile, "%s%s%s", space.c_str(), str.c_str(), end.c_str());
    }
}


// void FixFea::get_elmer(std::string& _buffer) { this->elmer->join(_buffer); }

/*-----------------------------------------------------------


elmer namespace definitions


------------------------------------------------------------*/

// void elmer::Base::join(std::string& _buffer) {
//     // for good measure
//     if (id == INT_MIN)
//         _buffer = name + "\n";
//     else
//         _buffer = name + " "  + std::to_string(this->id) + "\n";

//     // for (auto it : *this) {
//     //     if (it.second.has_value())
//     //         _buffer+=(this->_tab + it.first + this->sep + it.second.value() + "\n");
//     // }

//     _buffer+=this->_end;
// }


// elmer::Header::Header()                         { this->name = "Header";             sep = " ";   }
// elmer::Constants::Constants()                   { this->name = "Constants";          sep = " = "; }
// elmer::Simulation::Simulation()                 { this->name = "Simulation";         sep = " = "; }
// elmer::Solver::Solver()                         { this->name = "Solver";             sep = " = "; }
// elmer::Equation::Equation()                     { this->name = "Equation";           sep = " = "; }
// elmer::Material::Material()                     { this->name = "Material";           sep = " = "; }
// elmer::Body::Body()                             { this->name = "Body";               sep = " = "; }
// elmer::Initial_Condition::Initial_Condition()   { this->name = "Initial Condition";  sep = " = "; }
// elmer::Boundary_Condition::Boundary_Condition() { this->name = "Boundary Condition"; sep = " = "; }


// elmer::Elmer::Elmer() {
//     this->header     = Header();
//     this->simulation = Simulation();
//     this->constants  = Constants();
// }


// void elmer::Elmer::join(std::string& _buffer) {
//     _buffer.clear();

//     std::string _temp;
    
//     this->header.join(_temp);
//     _buffer += (_temp + this->sep);
    
//     this->simulation.join(_temp);
//     _buffer += (_temp + this->sep);
    
//     this->constants.join(_temp);
//     _buffer += (_temp + this->sep);
    
//     for (int i = 0; i < this->solvers.size(); i++) {
//         solvers[i].join(_temp);
//         _buffer += (_temp + this->sep);
//     }

//     for (int i = 0; i < this->equations.size(); i++) {
//         equations[i].join(_temp);
//         _buffer += (_temp + this->sep);
//     }

//     for (int i = 0; i < this->materials.size(); i++) {
//         materials[i].join(_temp);
//         _buffer += (_temp + this->sep);
//     }

//     for (int i = 0; i < this->bodys.size(); i++) {
//         bodys[i].join(_temp);
//         _buffer += (_temp + this->sep);
//     }

//     for (int i = 0; i < this->initial_conditions.size(); i++) {
//         initial_conditions[i].join(_temp);
//         _buffer += (_temp + this->sep);
//     }

//     for (int i = 0; i < this->boundary_conditions.size(); i++) {
//         boundary_conditions[i].join(_temp);
//         _buffer += (_temp);
//     }
// }



// /**
//  * Loading wall temperature info on start of step
//  */
// void FixFea::start_of_step() {
//     // only running if a valid time step
//     if (!(this->run_condition())) return;

//     int nlocal = surf->nlocal;

//     // reading the info from the file
//     memset(twall, 0, nlocal*sizeof(double));

//     // waits till all of the processes get here
//     MPI_Barrier(world);

//     // if this process is supposed to read the files, reading the temperature file
//     if (comm->me == 0) {
//         // loading data from the temperature file, this is usually spit out from elmer
//         // this->load_data();

//         std::vector<double> data;
//         try {
//             elmer::getLatestNodeData(this->temperature_data_file, data);
//         } catch (std::string _msg) {
//             error->all(FLERR, _msg.c_str()); 
//         } catch (std::exception& e) {
//             error->all(FLERR, e.what());
//         }

//         if (nlocal != this->boundary_data.size())
//             error->all(FLERR, "boundary data does not match required size");

//         // used later
//         double avg;
        
//         // averages values for the nodes of a surface element and sets this average to the 
//         // temperature of the surface element
//         for (int i = 0; i < this->boundary_data.size(); i++) {
//             // computes the average temperature of the nodes that make up the surface element
//             // this value is used to set the surface element temperature
//             avg = 0;
//             for (int j = 1; j < elmer::boundary_size; j++) {
//                 // gets the data point corresponding to node id and adds it to the rolling sum
//                 avg += data[this->boundary_data[i][j]-1];
//             }
//             // computing the average by dividing the sum by the number of points and setting the
//             // surface element
//             this->twall[this->boundary_data[i][0]-1] = avg/(elmer::boundary_size - 1);
//         }

//         // checking to make sure all values where set
//         for (int i = 0; i < nlocal; i++) {
//             if (this->twall[i] == (double)0)
//                 error->all(FLERR, "wall temperature not set correctly");
//         }

//         // sending to other processes
//         for (int i = 1; i < nprocs; i++) {
//             MPI_Send(&twall, nlocal, MPI_DOUBLE, i, 1, world);
//         }
//     } else {
//         // master process is 0, blocks until data is received
//         MPI_Recv(&twall, nlocal, MPI_DOUBLE, 0, 1, world, MPI_STATUS_IGNORE);
//     }

//     // waits till all of the processes get here
//     MPI_Barrier(world);
// }




// /**
//  * loads data from the elmer output file
//  */
// void FixFea::load_sif(std::string sif_path) {
//     this->print("Loading sif format from: " + sif_path);

//     // Read from the text file
//     std::ifstream sif_file;

//     // opening the node file
//     sif_file.open(sif_path);

//     // temporary vector to store the split strings
//     std::vector<std::string> v;

//     // clearing data
//     //this->sif_format.clear();

//     if (sif_file.is_open()) {
//         //sif_file >> this->sif_format;
//     } else {
//         // catching if the file did not open
//         error->all(FLERR, ((std::string)"sif file did not open, " + sif_path).c_str());
//     }
//     // Close the file
//     sif_file.close();
// }


// template<typename First, typename... Args>
// void to_elmer_section(std::vector<std::string>& _buffer, bool _specify_size, First first, Args... args) {
//     _buffer.clear();

//     toml::value _tbl = toml::find(toml::data, first, args...);

//     if (!(_tbl.is_table()))
//         toml::error("can only create elmer section from table");


//     std::string param_name;
//     std::string name;

//     std::string opt_str;
//     int opt_int;
//     double opt_double;
//     bool opt_bool;
//     std::ostringstream double_converter;
//     std::string arr_str;
//     toml::array arr;

//     double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);

//     for (auto it : _tbl.as_table()) {
//         double_converter.str("");
//         name = it.first;
//         param_name = name;
//         std::replace(param_name.begin(), param_name.end(), '_', ' ');

//         if (it.second.is_array()) {
//             arr = it.second.as_array();
//             if (_specify_size)
//                 arr_str = param_name + "(" + std::to_string(arr.size()) + ") " + toml::elmer_separator;
//             else
//                 arr_str = param_name + " " + toml::elmer_separator;

//             for (auto elem : arr) {
//                 double_converter.str("");
//                 if (elem.is_string()) {
//                     opt_str = elem.as_string().str;
//                     arr_str += (" \"" + opt_str + "\"");
                    
//                 } else if (elem.is_integer()) {
//                     opt_int = elem.as_integer();
//                     arr_str += (" " + std::to_string(opt_int));

//                 } else if (elem.is_floating()) {
//                     opt_double = elem.as_floating();
//                     double_converter << opt_double;
//                     arr_str += (" " + double_converter.str());
                
//                 } else if (elem.is_boolean()) {
//                     opt_bool = elem.as_boolean();

//                     if (opt_bool) 
//                         arr_str += (" True");
//                     else 
//                         arr_str += (" False");

//                 } else {
//                     toml::error("invalid type in " + name);

//                 }
//             }

//             _buffer.push_back(arr_str);

//         } else if (it.second.is_string()) {
//             opt_str = it.second.as_string();
            
//             _buffer.push_back(
//                 param_name + " " + toml::elmer_separator + " \"" + opt_str + "\""
//             );
            
//         } else if (it.second.is_integer()) {
//             opt_int = it.second.as_integer();

//             _buffer.push_back(
//                 param_name + " " + toml::elmer_separator + " " + std::to_string(opt_int)
//             );

//         } else if (it.second.is_floating()) {
//             opt_double = it.second.as_floating();

//             double_converter << opt_double;
//             _buffer.push_back(
//                 param_name + " "  + toml::elmer_separator + " " + double_converter.str()
//             );


//         } else if (it.second.is_boolean()) {
//             opt_bool = it.second.as_boolean();

//             if (opt_bool) {
//                 _buffer.push_back(
//                     param_name + " "  + toml::elmer_separator + " " + "True"
//                 );
//             } else {
//                 _buffer.push_back(
//                     param_name + " "  + toml::elmer_separator + " " + "False"
//                 );
//             }

//         } else {
//             std::stringstream ss;
//             ss << "Invalid type for " << name << ", has type " << it.second.type();
//             toml::error(ss.str());
//         }
//     }
// }

// bool is_int(const std::string& s) {
//     std::string::const_iterator it = s.begin();
//     while (it != s.end() && std::isdigit(*it)) ++it;
//     return !s.empty() && it == s.end();
// }

// template<typename T, typename First, typename... Args>
// void to_elmer_section_with_id(std::vector<T>& _buffer, bool _specify_size, First first, Args... args) {
//     _buffer.clear();

//     toml::value _tbl = toml::find(toml::data, first, args...);

//     if (!(_tbl.is_table()))
//         toml::error("can only create elmer section with id from table");

//     for (auto it : _tbl.as_table()) {
//         if (is_int(it.first)) {
//             T _temp = T();
//             to_elmer_section(_temp.contents, _specify_size, first, args..., it.first);
//             _buffer.push_back(_temp);
//         } else {
//             toml::error("id for elmer section must be int");
//         }
//     }
// }

// int countFreq(std::string pat, std::string txt) {
//     int M = pat.length();
//     int N = txt.length();
//     int res = 0;
 
//     /* A loop to slide pat[] one by one */
//     for (int i = 0; i <= N - M; i++) {
//         /* For current index i, check for
//            pattern match */
//         int j;
//         for (j = 0; j < M; j++)
//             if (txt[i + j] != pat[j])
//                 break;
 
//         // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
//         if (j == M) {
//             res++;
//         }
//     }
//     return res;
// }

// toml::value toml::handle_arr(toml::value _data, std::string& _path, std::string _orig_path) {
//     toml::key _temp;

//     if (!(_data.is_table()))
//         error("can not resolve path, " + _orig_path + ", can only index an element of a table");


//     if (_path.find(toml::start_index) != std::string::npos) {
//         std::cout << "handling arr:" << _path << "\n";
        
//         _temp = _path.substr(0, _path.find(toml::start_index));
//         boost::algorithm::trim(_temp);
//         std::cout << "temp: " << _temp << "\n";
//         if (_path.find(toml::start_index) != 0) {
//             std::cout << "in if\n";
//             if (_data.contains(_temp)) {
//                 std::cout << "getting\n";
//                 // std::cout << _data << "\n";
//                 ;//.at(_temp);
//                 _data = toml::find(_data, _temp);
//                 std::cout << "_data: " << _data << "\n";
//             } else
//                 toml::error("can not resolve path, " + _orig_path);
//         }
//         std::cout << "HERERERE\n";
//         if (!(_data.is_array()))
//             toml::error("can not resolve path, " + _orig_path + ", tried indexing something that is not an array");


//         int index;
//         for (int i = 0; i < countFreq(_path, toml::start_index); i++) {
//             _temp = _path.substr(_path.find(toml::start_index) + toml::start_index.length(), _path.find(toml::end_index));
//             boost::algorithm::trim(_temp);
//             index = std::stoi(_temp);
//             if (index >= _data.as_array().size() || index < 0)
//                 toml::error("can not resolve path, " + _orig_path + ", index out of bounds");
//             std::cout << "index = " << index << "\n";
//             _data = _data.as_array()[index];
//             _path = _path.substr(_path.find(end_index) + end_index.length(), _path.length());
//         }

//     }
//     return _data;
// }


// toml::value toml::get_from(std::string _path, std::string _orig_path) {
//     if (_orig_path.length() == 0) {
//         _orig_path = _path;
//     }

//     toml::value _temp_val;
//     std::string _name;

    

//     // if (_path.find(toml::separator) != std::string::npos) {
//     //     _name = _path.substr(0, _path.find(toml::separator));
//     //     _data = toml::handle_arr(_data, _name, _orig_path);
//     //     if (_name.length() == 0) {
//     //         return get_from(_data, _path.substr(_path.find(toml::separator) + toml::separator.length(), _path.length()), _orig_path);
//     //     } else if (_data.contains(_name)) {
//     //         return get_from(_data[_name], _path.substr(_path.find(toml::separator) + toml::separator.length(), _path.length()), _orig_path);
//     //     } else {
//     //         error("can not resolve path, " + _orig_path);
//     //     }

//     // } else {
//     //     _data = handle_arr(_data, _path, _orig_path);

//     //     if (_path.length() == 0) {
//     //         if (data.is_string()) {
//     //             _data = resolve_name(_data.as_string());
//     //         }
//     //         return _data;
//     //     } else if (_data.contains(_path)) {
//     //         _temp_val = _data[_path];
//     //         if (_temp_val.is_string()) {
//     //             _data[_path] = resolve_name(_temp_val.as_string());
//     //         }
//     //         return _temp_val;
//     //     } else {
//     //         error("can not resolve path, " + _orig_path);
//     //     }
//     // }
//     return _temp_val;
// }


// toml::value toml::resolve_name(toml::string _name) {
//     if (_name.str.substr(0, indicator.length()) == indicator) {
//         if (num_recursions > max_recursions)
//             error("reached max allowed recursion resolving " + (_name.str) + " (there is probably a circular definition somewhere)");
//         std::cout << "resolving: " << _name << "\n";
//         num_recursions++;
//         return get_from(_name.str.substr(indicator.length(), _name.str.length()));
//     }
//     return _name;
// }


// template <typename First, typename... Args>
// toml::value toml::preprocess(First first, Args... args) {
//     toml::value _tbl = toml::find(toml::data, first, args...);
//     if (!(_tbl.is_table()))
//         toml::error("did not get a table in preprocess");

//     for (auto& it : _tbl) {
//         num_recursions = 0;
//         if (it.second.is_table()) {
//             toml::preprocess(first, args..., it.first);
//         } else {
//             if (it.second.is_array()) {
//                 toml::array _arr = it.second.as_array();
//                 for (auto& elem : _arr) {
//                     if (elem.is_string()){
//                         elem = resolve_name(elem.as_string());
//                     }
//                 }
//                 toml::data[it.first] = _arr;
//             } else if (it.second.is_string()) {
//                 toml::data[it.first] = resolve_name(it.second.as_string());
//             }
//         }
//     }
// }

// // template<typename T>
// // void inform_setting(std::string _name, T& _var) {
// //     std::cout << _name << " = " << _var << "\n";
// // }


// template<typename T, typename First, typename... Args>
// void set(T& _var, First first, Args... args) {
//     _var = toml::find<T>(toml::data, first, args...);
// }


// int vec_to_arr(std::vector<std::string>& _vec, char**& _arr) {
//     const int _size = _vec.size();

//     _arr = new char*[_size];

//     for (int i = 0; i < _size; i++)
//         _arr[i] = (char*)_vec[i].c_str();
    
//     return _size;
// }


// void toml::format_name(std::string_view _s, std::string& _out) {
//     _out.clear();
//     _out = _s;
//     std::replace(_out.begin(), _out.end(), '_', ' ');
// }


// void toml::table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<std::string>& _var, std::string _sep, bool _specify_size) {
    
//     if (__tbl.type() != toml::node_type::table)
//         error(_caller + " must be given a table");

//     // std::cout << "Running table for: " << _caller << "\n";

//     // converting to a table
//     toml::table* _tbl = __tbl.as<toml::table>();

//     std::string param_name;
//     std::string name;

//     std::optional<std::string> opt_str;
//     std::optional<int64_t> opt_int;
//     std::optional<double> opt_double;
//     std::optional<bool> opt_bool;
//     std::ostringstream double_converter;
//     std::string arr_str;
//     toml::array arr;

//     double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);

//     _var.clear();

//     for (auto it : (*_tbl)) {
//         double_converter.str("");
//         name = it.first.str();
//         name = _caller + "." + name;
//         format_name(it.first.str(), param_name);

//         //std::cout << "  running: " << name << "\n";

//         if (toml::array_t == it.second.type()) {
//             arr = *it.second.as_array();
//             if (_specify_size)
//                 arr_str = param_name + "(" + std::to_string(arr.size()) + ") " + _sep;
//             else
//                 arr_str = param_name + _sep;

//             for (auto&& elem : arr) {
//                 double_converter.str("");
//                 if (toml::string_t == elem.type()) {
//                     opt_str = elem.value<std::string>();
//                     if (opt_str.has_value()) {
//                         arr_str += (" \"" + opt_str.value() + "\"");
//                     } else
//                         error("could not get value for " + name);

                    
//                 } else if (toml::integer_t == elem.type()) {
//                     opt_int = elem.value<double>();
//                     if (opt_int.has_value()) {
//                         arr_str += (" " + std::to_string(opt_int.value()));

//                     } else
//                         error("could not get value for " + name);

//                 } else if (toml::double_t == elem.type()) {
//                     opt_double = elem.value<double>();
//                     if (opt_double.has_value()) {
//                         double_converter << opt_double.value();
//                         arr_str += (" " + double_converter.str());

//                     } else
//                         error("could not get value for " + name);
                
//                 } else if (toml::bool_t == elem.type()) {
//                     opt_bool = elem.value<bool>();
//                     if (opt_bool.has_value()) {
//                         if (opt_bool.value()) 
//                             arr_str += (" True");
//                         else 
//                             arr_str += (" False");
                        
//                     } else
//                         error("could not get value for " + name);
//                 } else {
//                     error("invalid type in " + name);

//                 }
//             }

//             _var.push_back(arr_str);

//         } else if (toml::string_t == it.second.type()) {
//             opt_str = it.second.value<std::string>();
//             if (opt_str.has_value()) {
//                 _var.push_back(
//                     param_name + " " + _sep + " \"" + opt_str.value() + "\""
//                 );
//             } else
//                 error("could not get value for " + name);

            
//         } else if (toml::integer_t == it.second.type()) {
//             opt_int = it.second.value<double>();
//             if (opt_int.has_value()) {
//                 _var.push_back(
//                     param_name + " " + _sep + " " + std::to_string(opt_int.value())
//                 );
//             } else
//                 error("could not get value for " + name);

//         } else if (toml::double_t == it.second.type()) {
//             opt_double = it.second.value<double>();
//             if (opt_double.has_value()) {
//                 double_converter << opt_double.value();
//                 _var.push_back(
//                     param_name + " "  + _sep + " " + double_converter.str()
//                 );
//             } else
//                 error("could not get value for " + name);

//         } else if (toml::bool_t == it.second.type()) {
//             opt_bool = it.second.value<bool>();
//             if (opt_bool.has_value()) {
//                 if (opt_bool.value()) {
//                     _var.push_back(
//                         param_name + " "  + _sep + " " + "True"
//                     );
//                 } else {
//                     _var.push_back(
//                         param_name + " "  + _sep + " " + "False"
//                     );
//                 }
//             } else
//                 error("could not get value for " + name);

//         } else {
//             std::stringstream ss;
//             ss << "Invalid type for " << name << ", has type " << it.second.type();
//             error(ss.str());
//         }

//         // std::cout << param_name << " : " << it.second.type() << "\n";
//     }

//     // for (auto it : _var) {
//     //     std::cout << "=> " << it << "\n";
//     // }
// }

// bool toml::is_int(const std::string& s) {
//     std::string::const_iterator it = s.begin();
//     while (it != s.end() && std::isdigit(*it)) ++it;
//     return !s.empty() && it == s.end();
// }

// template<typename T>
// void toml::id_table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<T>& _var, std::string _sep, bool _specify_size) {
//     if (__tbl.type() != toml::node_type::table)
//         error(_caller + " must be given a table");

//     // std::cout << "Running id table for: " << _caller << "\n";

//     // converting to a table
//     toml::table* _tbl = __tbl.as<toml::table>();

//     std::string _name;
//     std::string _key;

//     for (auto it : (*_tbl)) {
//         _key = it.first.str();
//         if (!(is_int(_key)))
//             error("elements for " + _caller + " must be a positive integer, not " + _key);

//         _name = _caller + "." + _key;

//         T inst = T();
//         inst.id = std::stoi(_key);

//         table_value_parser(_name, (*_tbl)[_key], inst.contents, _sep, _specify_size);
        

//         _var.push_back(inst);

//         // std::cout << " <> "<< it.first.str() << "\n";
//     }

// }


// void toml::join_array_sparta(std::string _caller, std::string _name, toml::node_t val, std::vector<std::string>& _buffer) {
//     _caller = _caller + "." + _name;
//     if (val.type() != toml::array_t)
//         error("can not join, " + _caller + ", it is not an array");

//     std::optional<std::string> opt_str;
//     std::optional<int64_t> opt_int;
//     std::optional<double> opt_double;
//     std::optional<bool> opt_bool;
//     std::ostringstream double_converter;
//     toml::array arr;

//     double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);

//     arr = *val.as_array();

//     _buffer.clear();

//     for (auto&& elem : arr) {
//         double_converter.str("");
//         if (toml::string_t == elem.type()) {
//             opt_str = elem.value<std::string>();
//             if (opt_str.has_value()) {
//                 _buffer.push_back(opt_str.value());
//             } else
//                 error("could not get value for " + _caller);

            
//         } else if (toml::integer_t == elem.type()) {
//             opt_int = elem.value<double>();
//             if (opt_int.has_value()) {
//                 _buffer.push_back(std::to_string(opt_int.value()));

//             } else
//                 error("could not get value for " + _caller);

//         } else if (toml::double_t == elem.type()) {
//             opt_double = elem.value<double>();
//             if (opt_double.has_value()) {
//                 double_converter << opt_double.value();
//                 _buffer.push_back(double_converter.str());

//             } else
//                 error("could not get value for " + _caller);
        
//         // } else if (toml::bool_t == elem.type()) {
//         //     opt_bool = elem.value<bool>();
//         //     if (opt_bool.has_value()) {
//         //         if (opt_bool.value()) 
//         //             arr_str += (" True");
//         //         else 
//         //             arr_str += (" False");
                
//         //     } else
//         //         error("could not get value for " + name);
//         } else {
//             error("invalid type in " + _caller);

//         }
//     }
// }


























// void FixFea::handle_both(std::string _caller, toml::node_t tbl) {
//     toml::dict_t options = {
//         std::make_pair("emi", &FixFea::handle_emi),
//         std::make_pair("tsurf_file", &FixFea::handle_tsurf_file)
//     };

//     this->run_table(_caller, "both", tbl, options);
// }


// void FixFea::handle_emi(std::string _caller, toml::node_t val) {
//     this->set(_caller, "emi", toml::double_t, val, this->emi);
//     if (emi <= 0.0 || emi > 1.0)
//         error->all(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");

//     this->print("emi = " + std::to_string(this->emi));
// }


// void FixFea::handle_tsurf_file(std::string _caller, toml::node_t val) {
//     this->set(_caller, "tsurf_file", toml::string_t, val, this->tsurf_file);

//     if (!(stat(this->tsurf_file.c_str(),  &sb) == 0))
//         error->all(FLERR, "Illegal fix fea command, tsurf_file path does not exist");

//     this->print("tsurf_file = " + this->tsurf_file);
// }



// void FixFea::handle_sparta(std::string _caller, toml::node_t tbl) {
//     toml::dict_t options = {
//         std::make_pair("nevery", &FixFea::handle_nevery),
//         std::make_pair("groupID", &FixFea::handle_groupID),
//         std::make_pair("mixID", &FixFea::handle_mixID),
//         std::make_pair("customID", &FixFea::handle_customID),
//         std::make_pair("compute", &FixFea::handle_compute)
//     };

//     this->run_table(_caller, "sparta", tbl, options);
// }


// void FixFea::handle_nevery(std::string _caller, toml::node_t val) {
//     this->set(_caller, "nevery", toml::integer_t, val, this->run_every);
//     if (this->run_every <= 0)
//         error->all(FLERR,"Illegal fix fea command, nevery <= 0");
//     this->print("nevery = " + std::to_string(this->run_every));
// }


// void FixFea::handle_groupID(std::string _caller, toml::node_t val) {
//     this->set(_caller, "groupID", toml::string_t, val, this->groupID);
//     // get the surface group
//     int igroup = surf->find_group(this->groupID.c_str());
//     if (igroup < 0)
//         error->all(FLERR,"Fix fea group ID does not exist");
    
//     groupbit = surf->bitmask[igroup];
//     this->print("groupID = " + this->groupID);
// }


// void FixFea::handle_mixID(std::string _caller, toml::node_t val) {
//     this->set(_caller, "mixID", toml::string_t, val, this->mixID);
//     int imix = particle->find_mixture((char*)this->mixID.c_str());
//     if (imix < 0)
//         error->all(FLERR,"Compute thermal/grid mixture ID does not exist");

//     ngroup = particle->mixture[imix]->ngroup;
//     this->print("mixID = " + this->mixID);
// }


// void FixFea::handle_customID(std::string _caller, toml::node_t val) {
//     this->set(_caller, "customID", toml::string_t, val, this->customID);
//     this->tindex = surf->add_custom((char*)this->customID.c_str(), DOUBLE, 0);
//     this->print("customID = " + this->customID);
// }


// void FixFea::handle_compute(std::string _caller, toml::node_t val) {
//     toml::join_array_sparta(_caller, "compute", val, this->compute_args);
// }


// void FixFea::handle_elmer(std::string _caller, toml::node_t tbl) {
//     toml::dict_t options = {
//         std::make_pair("exe", &FixFea::handle_exe),
//         std::make_pair("sif", &FixFea::handle_sif),
//         std::make_pair("meshDBstem", &FixFea::handle_meshDBstem),
//         std::make_pair("header", &FixFea::handle_header),
//         std::make_pair("simulation", &FixFea::handle_simulation),
//         std::make_pair("constants", &FixFea::handle_constants),
//         std::make_pair("solver", &FixFea::handle_solver),
//         std::make_pair("equation", &FixFea::handle_equation),
//         std::make_pair("material", &FixFea::handle_material),
//         std::make_pair("body", &FixFea::handle_body),
//         std::make_pair("initial_condition", &FixFea::handle_initial_condition),
//         std::make_pair("boundary_condition", &FixFea::handle_boundary_condition)
//     };

//     this->run_table(_caller, "elmer", tbl, options);
// }


// void FixFea::handle_exe(std::string _caller, toml::node_t val) {
//     this->set(_caller, "exe", toml::string_t, val, this->elmer.exe);

//     if (!(stat(this->elmer.exe.c_str(),  &sb) == 0))
//         error->all(FLERR, "Illegal fix fea command, exe path does not exist");
    
//     this->print("exe = " + this->elmer.exe);
// }

// void FixFea::handle_sif(std::string _caller, toml::node_t val) {
//     this->set(_caller, "sif", toml::string_t, val, this->elmer.sif);
    
//     if (!(stat(this->elmer.sif.c_str(),  &sb) == 0))
//         error->all(FLERR, "Illegal fix fea command, sif path does not exist");
    
//     this->print("sif = " + this->elmer.sif);
// }

// void FixFea::handle_meshDBstem(std::string _caller, toml::node_t val) {
//     this->set(_caller, "meshDBstem", toml::string_t, val, this->elmer.meshDBstem);

//     std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
//     for (int i = 0; i < 4; i++) {
//         if (!(stat((this->elmer.meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
//             error->all(FLERR, ("Illegal fix fea command, mesh database incomplete, " + (this->elmer.meshDBstem + "." + exts[i]) + " does not exist").c_str());
//     }
//     this->print("meshDBstem = " + this->elmer.meshDBstem);
// }


// void FixFea::handle_header(std::string _caller, toml::node_t __tbl) {
//     _caller = _caller + ".header";
//     toml::table_value_parser(_caller, __tbl, this->elmer.header.contents, "");
// }


// void FixFea::handle_simulation(std::string _caller, toml::node_t __tbl) {
//     _caller = _caller + ".simulation";
//     table_value_parser(_caller, __tbl, this->elmer.simulation.contents);
// }


// void FixFea::handle_constants(std::string _caller, toml::node_t __tbl) {
//     _caller = _caller + ".constants";
//     table_value_parser(_caller, __tbl, this->elmer.constants.contents);
// }


// void FixFea::handle_solver(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".solver";
//     id_table_value_parser(_caller, tbl, this->elmer.solvers);
// }


// void FixFea::handle_equation(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".equation";
//     id_table_value_parser(_caller, tbl, this->elmer.equations);
// }


// void FixFea::handle_material(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".material";
//     id_table_value_parser(_caller, tbl, this->elmer.materials);
// }


// void FixFea::handle_body(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".body";
//     id_table_value_parser(_caller, tbl, this->elmer.bodys);
// }


// void FixFea::handle_initial_condition(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".initial_condition";
//     id_table_value_parser(_caller, tbl, this->elmer.initial_conditions);
// }


// void FixFea::handle_boundary_condition(std::string _caller, toml::node_t tbl) {
//     _caller = _caller + ".boundary_condition";
//     id_table_value_parser(_caller, tbl, this->elmer.boundary_conditions);
// }




/*-------------------------------------------------


Utility functions


-------------------------------------------------*/

// /**
//  * executes a table
// */
// template<typename dict>
// void FixFea::run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options) {
//     if () 
//         error->all(FLERR, (_caller + " must be given a table").c_str());

//     if (_caller.length())
//         _caller = _caller + "." + _name;
//     else
//         _caller = _name;

//     for (toml::dict_item_t it : _options) {
//         toml::node_t val = _tbl[it.first];
//         (this->*it.second)(_caller, val);
//         _tbl.erase(it.first);
//     }

//     if (_tbl.size()) {
//         std::string msg;
//         if (_caller.length())
//             msg = "Invalid args in section \"" + _caller + "\":\n";
//         else
//             msg = "Invalid section:\n";
//         for (auto it : _tbl) {
//             msg+=("-> " + std::string(it.first.str()) + "\n");
//         }
//         error->all(FLERR, msg.c_str());
//     }
// }


// /**
//  * executes a table that wrapped as a node_t type
// */
// template<typename dict>
// void FixFea::run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options) {

//     if (__tbl.type() != toml::node_type::table) 
//         error->all(FLERR, (_caller + " must be given a table").c_str());

//     // converting to a table
//     toml::table* _tbl = __tbl.as<toml::table>();

//     // calling the other definition because there is nothing else unique needed
//     run_table(_caller, _name, (*_tbl), _options);
// }


















// template<typename T>
// void toml::set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val) {
//     // if the _val is none type, just set to default value immediately
//     if (_val.type() == none_t) {
//         _var = _default_val;
//         return;
//     }

//     if (_val.type() != _type) {
//         // adding the calling hierarchy to the name
//         _name = _caller + "." + _name;

//         // making message to user
//         std::stringstream ss;
//         ss << "\"" << _name << "\" has wrong type, " << _val.type();
//         ss << ", must be " << _type << " type";
//         error(ss.str());
//     }

//     std::optional<T> _var_opt = _val.value<T>();

//     if (_var_opt.has_value())
//         _var = _var_opt.value();
//     else
//         _var = _default_val;
// }

// toml::node_t FixFea::get_var(toml::node_t _tbl, std::string _path, std::string _orig_path) {
//     if (_tbl.type() == toml::table_t) {
//         get_var((*_tbl.as_table()), _path, _orig_path);
//     } else {

//     }
// }


// toml::node_t FixFea::get_var(toml::table _tbl, std::string _path, std::string _orig_path) {
//     if (_orig_path.length() == 0) {
//         _orig_path = _path;
//     }
//     std::string _name;
//     std::string _sep = ".";
//     toml::node_t _to_return;

//     if (_path.find(_sep) != std::string::npos) {
//         _name = _path.substr(0, _path.find(_sep));

//         if (_tbl.contains(std::string_view(_name))) {
//             _to_return = get_var(_tbl[_name], _path.substr(_path.find(_sep) + _sep.length(), _path.length()), _orig_path);
//         } else {
//             error->all(FLERR, ("Could not resolve path, " + _orig_path).c_str());
//         }   
//     } else {
//         if (_tbl.contains(_path)) {
//             _to_return = _tbl[_path];
//             if (_to_return.type() == toml::string_t) {
//                 _to_return = this->resolve_name(_to_return, _path);
//             }
//         } else {
//             error->all(FLERR, ("Could not resolve path, " + _orig_path).c_str());
//         }
//     }

//     return _to_return;

// }


// toml::node_t FixFea::resolve_name(toml::node_t _val, std::string _path) {
//     if (_val.type() == toml::string_t) {
//         std::string _temp;
//         std::optional<std::string> _var_opt_temp = _val.value<std::string>();
//         if (_var_opt_temp.has_value()) {
//             _temp = _var_opt_temp.value();
//             if (_temp.substr(0, toml::indicator.length()) == toml::indicator) {
//                 if (this->recursion_num > this->max_recursion) {
//                     error->all(FLERR, ("reached max allowed recursion resolving " + _path + " (there is probably a circular definition somewhere)").c_str());
//                 }

//                 this->recursion_num++;
//                 _temp = _temp.substr(toml::indicator.length(), _temp.length());

//                 return get_var(this->tbl, _temp);
//             }
//         } else error->all(FLERR, (char*)"no value gotten for optional var");
//     } else if (_val.type() == toml::none_t)
//         error->all(FLERR, ("could not resolve path, " + _path).c_str());

//     return _val;
// }


// // errors instead of setting to default
// template<typename T>
// void FixFea::set(std::string _caller, std::string _name, toml::var_type_t _type, toml::node_t _val, T& _var) {
//     _name = _caller + "." + _name;

//     // if the _val is none type, just set to default value immediately
//     if (_val.type() == toml::none_t) {
//         error->all(FLERR, ("no value provided for \"" + _name + "\"").c_str());
//     }

//     if (_val.type() == toml::string_t) {
//         _val = this->resolve_name(_val, _caller);
//     }

//     //std::cout << _val << "\n";
//     std::optional<T> _var_opt = _val.value<T>();

//     if (_val.type() != _type) {
//         std::stringstream ss;
//         ss << "\"" << _name << "\" has wrong type, " << _val.type() << ", must be " << _type << " type";
//         error->all(FLERR, ss.str().c_str());
//     }

//     _var_opt = _val.value<T>();

//     if (_var_opt.has_value()) {
//         _var = _var_opt.value();
//     } else {
//         error->all(FLERR, ("no value provided for \"" + _name + "\"").c_str());
//     }
// }