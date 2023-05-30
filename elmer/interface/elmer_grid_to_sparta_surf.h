#ifndef ELMER_GRID_TO_SPARTA_SURF_H
#define ELMER_GRID_TO_SPARTA_SURF_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
// #include "tools/util.h"
#include <boost/algorithm/string.hpp>

class ElmerGridToSpartaSurf {
    public:
        ElmerGridToSpartaSurf(std::string _file_stem);

        void load_boundary();
        void load_nodes();
        void update_boundary_file();
        void make_sparta_surf();
    
    protected:
        std::string file_stem;
        std::vector<std::vector<std::string>> boundary_data;
        std::vector<std::pair<std::string, std::vector<std::string>>> node_data;
};

#endif