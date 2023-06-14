#include "fixfea.h"

// #include <iomanip>
// #include <cmath>
// #include <limits>
#include <fstream>

int main(int argc, char** argv) {

    // for (int i = 0; i < argc; i++) {
    //     std::cout << argv[i] << "\n";
    // }
    // return 0;

    // std::cout << std::to_string(INT64_MAX) << "\n";
    // std::cout << INT64_MAX << "\n";
    // std::cout << std::numeric_limits<double>::digits10 << "\n";
    // return 0;

    toml::table tbl;
    try {
        tbl = toml::parse_file(argv[1]);
    } catch (const toml::parse_error& err) {
        error(FLERR, "Parsing failed:\n" + std::string(err.description()));
    }

    toml::dict_t items = {
        std::make_pair("sparta", &FixFea::handle_sparta),
        std::make_pair("elmer", &FixFea::handle_elmer),
        std::make_pair("both", &FixFea::handle_both)
    };


    try {
        FixFea fix = FixFea();
        fix.run_table("", "", tbl, items);
        
        std::string buffer;
        fix.get_elmer(buffer);

        std::ofstream f("test.sif");
        if (!(f.is_open()))
            error(FLERR, "File did not open");
        
        f << buffer;
        f.close();


    } catch (...) {
        error(FLERR, "in running table");
    }

    return 0;
}