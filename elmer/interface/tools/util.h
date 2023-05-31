#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <string>
 
template<typename T>
void print_vector(std::vector<T>& vec) {
    std::cout << "[";
    for (int i = 0; i < vec.size() - 1; i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << vec.back() << "]\n";
}

void error(const char *file, int line, const char *str);



#endif