#ifndef ELMER_TO_SPARTA_TEMPERATURE_H
#define ELMER_TO_SPARTA_TEMPERATURE_H

#include <iostream>

#include "elmer_grid_to_sparta_surf.h"

constexpr auto max_precision {std::numeric_limits<long double>::digits10 + 1};

class ElmerToSpartaTemperature : public ElmerGridToSpartaSurf {
    public:
        ElmerToSpartaTemperature(std::string _data_file, std::string _file_stem);
        void make_sparta_surf_data();

    private:
        void load_data();
        std::string data_file;
        std::vector<std::string> data;

        
};

#endif