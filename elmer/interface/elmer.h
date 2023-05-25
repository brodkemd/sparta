#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/join.hpp>

class Elmer {
    public:
        Elmer(std::string output_file);
        void write();

        void add_line(std::string line);
        void add_section(std::string name, std::vector<std::string> content);

        /* elmer input sections */
        void Header(std::vector<std::string> content = {});
        void Simulation(std::vector<std::string> content);
        void Constants(std::vector<std::string> content = {});

        // n is the number of keyword commands
        void Body(int n, std::vector<std::string> content);
        void Material(int n, std::vector<std::string> content);
        void Body_Force(int n, std::vector<std::string> content);
        void Equation(int n, std::vector<std::string> content);
        void Solver(int n, std::vector<std::string> content);
        void Boundary_Conditions(int n, std::vector<std::string> content);
        void Component(int n, std::vector<std::string> content);

    protected:
        std::string sif_file = "";
        std::string tab = "  ";
        std::vector<std::string> data = {};

};

#endif