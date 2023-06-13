#include "toml.h"
// #include <iostream>
// using namespace toml;

void error(std::string file, int line, std::string msg) {
    std::cout << "ERROR: (" << file << ":" << line << ") " << msg << "\n";
    exit(1);
}


