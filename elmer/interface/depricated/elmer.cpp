#include "elmer.h"


Elmer::Elmer(std::string output_file) { this->sif_file = output_file; }


void Elmer::add_section(std::string name, std::vector<std::string> content) {
    std::string line = name + "\n";

    for (std::string item : content) {
        line += this->tab + item + "\n";
    }
    line += "End";
    this->data.push_back(line);
}


void Elmer::write() {
    std::cout << "Writing to: " << this->sif_file << "\n";
    std::ofstream myfile(this->sif_file);

    for (std::string line : this->data) myfile << line << "\n\n";

    myfile.close();
}


/***
 * This begins the different types of input sections
 * 
*/

void Elmer::Header(std::vector<std::string> content) {
    if (content.size() == 0) {
        content = {
            "CHECK KEYWORDS Warn",
            "Mesh DB \".\" \".\"",
            "Include Path \"\"",
            "Results Directory \"\""
        };
    }
    this->add_section("Header", content);
}

void Elmer::Simulation(std::vector<std::string> content) {
    this->add_section("Simulation", content);
}


void Elmer::Constants(std::vector<std::string> content) {
    if (content.size() == 0) {
        content = {
            "Gravity(4) = 0 -1 0 9.82",
            "Stefan Boltzmann = 5.67e-08",
            "Permittivity of Vacuum = 8.8542e-12",
            "Permeability of Vacuum = 1.25663706e-6",
            "Boltzmann Constant = 1.3807e-23",
            "Unit Charge = 1.602e-19"
        };
    }
    this->add_section("Constants", content);
}

// n is the number of keyword commands
void Elmer::Body(int n, std::vector<std::string> content) {
    this->add_section("Body " + std::to_string(n), content);
}

void Elmer::Material(int n, std::vector<std::string> content) {
    this->add_section("Material " + std::to_string(n), content);
}

void Elmer::Body_Force(int n, std::vector<std::string> content) {
    this->add_section("Body Force " + std::to_string(n), content);
}

void Elmer::Equation(int n, std::vector<std::string> content) {
    this->add_section("Equation " + std::to_string(n), content);
}

void Elmer::Solver(int n, std::vector<std::string> content) {
    this->add_section("Solver " + std::to_string(n), content);
}

void Elmer::Boundary_Conditions(int n, std::vector<std::string> content) {
    this->add_section("Boundary Conditions " + std::to_string(n), content);
}

void Elmer::Component(int n, std::vector<std::string> content) {
    this->add_section("Component " + std::to_string(n), content);
}
