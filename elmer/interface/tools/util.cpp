#include "util.h"

void error(const char *file, int line, const char *str) {
    std::cout << "ERROR: " << str << " " << "(" << file << ":" << line << ")\n";
    exit(1);
}
