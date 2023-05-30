#ifndef ELMER_TO_SPARTA_TEMPERATURE_H
#define ELMER_TO_SPARTA_TEMPERATURE_H

#include <iostream>
#include "elmer_grid_to_sparta_surf.h"

class ElmerToSpartaTemperature : public ElmerGridToSpartaSurf {
    private:
        int temp;
    // ElmerToSpartaTemperature() {};
};

#endif