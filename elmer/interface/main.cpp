#include "main.h"

int main() {
    // std::cout << "\nStarting\n";
    std::string file_stem = "/home/marekbrodke/Documents/C++/sparta/elmer/test_data/mesh";

    ElmerGridToSpartaSurf inst(file_stem);
    inst.make_sparta_surf();
    
    std::string data_file = "/home/marekbrodke/Documents/C++/sparta/elmer/test_data/case.dat";

    ElmerToSpartaTemperature inst2(data_file, file_stem);
    inst2.make_sparta_surf_data();

    return 0;
}