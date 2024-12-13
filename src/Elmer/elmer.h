#ifndef ELMER_H
#define ELMER_H

#include <array>
#include "python.h"
#include "util.h"
#include "../cVTU.h"

#define ELMER_LINE_START "  "
#define ELMER_ARRAY_SEP " "
#define ELMER_SECTION_END "End"

#define ELMER_DIMENSION 3

#define BOUNDARY_ELEMENT_START_READ 1 // index for a line in the boundary file to start reading data from (including it)
#define BOUNDARY_ELEMENT_LINE_SIZE 8
#define BOUNDARY_ELEMENT_ARRAY_SIZE (BOUNDARY_ELEMENT_LINE_SIZE-BOUNDARY_ELEMENT_START_READ)
#define BOUNDARY_ELEMENT_ARRAY_NODE_START 4 // is an index starting from second column
#define BOUNDARY_ELEMENT_NODE_COUNT (BOUNDARY_ELEMENT_LINE_SIZE-(BOUNDARY_ELEMENT_ARRAY_NODE_START+BOUNDARY_ELEMENT_START_READ))

// elements are not restricted to a size
#define ELEMENT_START_READ 1 // index for a line in the boundary file to start reading data from (including it)
#define ELEMENT_ARRAY_NODE_START 2 // is an index starting from second column

// difference of these two should equal ELMER_DIMENSION
#define NODE_START_READ 2 // index for a line in the boundary file to start reading data from (including it)
#define NODE_LINE_SIZE 5

// #define HEAT_FLUX_TOL 1e-20
// #define SHEAR_TOL 1e-20
// #define FORCE_TOL 1e-20

#define DOUBLE_TRUNCATE_PRECISION 20

// #define HEAT_FLUX_PRECISION 20
// #define SHEAR_PRECISION 20
// #define FORCE_PRECISION 20

/*
Notes:
    - On the elements file
        - the first column is removed when read
    - On the boundary file
        - the first column is removed when read
    - On the nodes file
        - the first and second column are remove when read

*/

/* ---------------------------------------------------------------------- */

// defining the os path separator based on the detected os
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    const char SEP = '\\';
#else
    const char SEP = '/';
#endif

/* ---------------------------------------------------------------------- */

// dummy classes
namespace SPARTA_NS { class FixFea;  }
namespace elmer     { class Server;  }

namespace elmer {
    
    /**
     * Base class that is used to construct the sections in a sif file
    */
    class Section {
        private:
            std::string _name, _sep;
            SPARTA_NS::index_t _id;
            std::vector<std::array<std::string, 2>> _content_pairs;

        public:
            Section() {}
            Section(std::string name, SPARTA_NS::index_t id = util::NO_INDEX_T, std::string sep = " = ");
            
            // writes the content of this class to the file buffer
            void joinInto(util::oFile& _buf);

            // adds a variable name and its value to the class contents
            void addEquality(std::string var, SPARTA_NS::index_t   val, bool include_count = false);
            void addEquality(std::string var, double val, bool include_count = false);
            void addEquality(std::string var, std::vector<SPARTA_NS::index_t>&   val, bool include_count = true);
            void addEquality(std::string var, std::vector<double>& val, bool include_count = true);
            void addEquality(std::string var, double*& arr, SPARTA_NS::index_t len, bool include_count = true);
    };


    class Elmer {
        private:
            SPARTA_NS::FixFea* sparta;
            python::handler*   python;
            Server*            server;          

            /* ---------------------------------------------------------------------- */

            SPARTA_NS::index_t   *cells, *offsets, *elements_additional_info, *boundary_additional_info, me, body_force_id, equation_id, material_id, num_nodes, num_elements, num_boundary;

            cVTU_cell_type_t *cell_types;

            bool   *allow_heat_flux, *allow_forces, *clamped,  gravity_on, has_been_initialized = false, elastic_solver_is_used;

            double *nodes, *node_velocities, *node_displacements,  *node_temperatures, base_temp;
            
            char   *sif, *meshDB, *node_data_file, *simulation_directory, *boundary_file, 
                   *element_file, *header_file, *node_file, *node_temperature_file_ext, 
                   *node_position_file_ext, *node_velocity_file_ext, *sif_file_config_str;

            /* ---------------------------------------------------------------------- */

            PyObject *header, *simulation, *equation, *body_force, *elmer, *material;
            std::vector<Section*> bodies, initial_conditions, boundary_conditions;

            /* ---------------------------------------------------------------------- */

        public:
            Elmer(SPARTA_NS::FixFea* _sparta, SPARTA_NS::index_t _me, python::handler* _h);
            ~Elmer();

            void makeSpartaSurf(char* fname);
            void averageNodeTemperaturesInto(double*& _temperatures, SPARTA_NS::index_t _length);
            void getNodePointAtIndex(SPARTA_NS::index_t _index, SPARTA_NS::index_t _boundary_index, double(&_point)[ELMER_DIMENSION]);
            void run();
            void dump();
            void makeSif();
            bool shouldUpdateSurf();

        protected:
            // void handleBoundaryFileSplitLine(std::vector<std::string>& split, SPARTA_NS::index_t line_number);
            // void handleElementFileSplitLine(std::vector<std::string>& split, SPARTA_NS::index_t line_number);
            void handleNodeFileSplitLine(std::vector<std::string>& split, SPARTA_NS::index_t line_number);
            void handleUserBoundaryConditions();
            void createInitialConditions();
            void createBoundaryConditions();
            void setupVectors();
            void updateNodes();
            void updateNodeFile();
            void loadNodeData();
            void loadNodes();
            void loadBoundaries();
            void loadElements();

            void dumpNodeTemperatures();
            // void dumpNodePositions();
            void dumpNodeVelocities();
    };
}

#endif

// bool shouldRun() { return true; }
// dummy class to allow the use of a pointer later
// namespace SPARTA_NS { class FixFea; }
// void join(util::oFile& _buf);
// void dumpNodeTemperaturesBefore();
// void dumpNodeVelocitiesBefore();
// void dumpNodePositionsBefore();
// void makeSif();
// void setupDumpDirectory();
// void dumpBefore();
// void addEquality(std::string var, unsigned long val);