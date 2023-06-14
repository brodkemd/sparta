#include "classes.h"


void Base::join(std::string& _buffer) {
    // for good measure
    if (id == INT_MIN)
        _buffer = name + "\n";
    else
        _buffer = name + " "  + std::to_string(this->id) + "\n";

    for (std::string it : this->contents) {
        _buffer+=(this->tab + it + "\n");
    }
    _buffer+=this->end;
}


Header::Header() {
    this->name = "Header";
}


Constants::Constants() {
    this->name = "Constants";
}


Simulation::Simulation() {
    this->name = "Simulation";
}


Solver::Solver() {
    this->name = "Solver";
}


Equation::Equation() {
    this->name = "Equation";
}


Material::Material() {
    this->name = "Material";
}


Body::Body() {
    this->name = "Body";
}


Initial_Condition::Initial_Condition() {
    this->name = "Initial Condition";
}


Boundary_Condition::Boundary_Condition() {
    this->name = "Boundary Condition";
}


Elmer::Elmer() {
    this->header = Header();
    this->simulation = Simulation();
    this->constants = Constants();
}


void Elmer::join(std::string& _buffer) {
    _buffer.clear();

    std::string _temp;
    
    this->header.join(_temp);
    _buffer += (_temp + this->sep);
    
    this->simulation.join(_temp);
    _buffer += (_temp + this->sep);
    
    this->constants.join(_temp);
    _buffer += (_temp + this->sep);
    
    for (int i = 0; i < this->solvers.size(); i++) {
        solvers[i].join(_temp);
        _buffer += (_temp + this->sep);
    }

    for (int i = 0; i < this->equations.size(); i++) {
        equations[i].join(_temp);
        _buffer += (_temp + this->sep);
    }

    for (int i = 0; i < this->materials.size(); i++) {
        materials[i].join(_temp);
        _buffer += (_temp + this->sep);
    }

    for (int i = 0; i < this->bodys.size(); i++) {
        bodys[i].join(_temp);
        _buffer += (_temp + this->sep);
    }

    for (int i = 0; i < this->initial_conditions.size(); i++) {
        initial_conditions[i].join(_temp);
        _buffer += (_temp + this->sep);
    }

    for (int i = 0; i < this->boundary_conditions.size(); i++) {
        boundary_conditions[i].join(_temp);
        _buffer += (_temp);
    }
}