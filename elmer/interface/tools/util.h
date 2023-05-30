#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>

void print() {
    std::cout << '\n';
}

template <typename T, typename... TAIL>
void print(const T &t, TAIL... tail) {
    std::cout << t << " ";
    print(tail...);
}


#endif