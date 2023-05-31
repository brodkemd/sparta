#include "main.h"

int main() {
    std::cout << "\nStarting\n";
    std::string file_stem = "/home/marekbrodke/Documents/C++/sparta/elmer/test_mesh/mesh";

    ElmerGridToSpartaSurf inst(file_stem);
    inst.make_sparta_surf();

    return 0;
}