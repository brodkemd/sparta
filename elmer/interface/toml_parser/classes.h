#ifndef CLASSES_H
#define CLASSES_H

#include <vector>
#include <string>
#include <iostream>
#include <climits>

class Base {
    private:
        std::string tab = "  ";
        std::string end = "End";

    protected:
        std::string name;

    public:
        int id = INT_MIN;
        std::vector<std::string> contents;

        void join(std::string& _buffer);
};


class Header : public Base {
    public:
        Header();
};


class Constants : public Base {
    public:
        Constants();
};


class Simulation : public Base {
    public:
        Simulation();
};


class Solver : public Base {
    public:
        Solver();
};


class Equation : public Base {
    public:
        Equation();
};


class Material : public Base {
    public:
        Material();
};


class Body : public Base {
    public:
        Body();
};


class Initial_Condition : public Base {
    public:
        Initial_Condition();
};


class Boundary_Condition : public Base {
    public:
        Boundary_Condition();
};


class Elmer {
    private:
        std::string sep = "\n\n";

    public:
        Elmer();

        void join(std::string& _buffer);
        

        std::string name;

        Header     header;
        Simulation simulation;
        Constants  constants;
        std::vector<Solver> solvers;
        std::vector<Equation> equations;
        std::vector<Material> materials;
        std::vector<Body> bodys;
        std::vector<Initial_Condition> initial_conditions;
        std::vector<Boundary_Condition> boundary_conditions;


        std::string exe;
        std::string sif;
        std::string meshDBstem;

};

#endif