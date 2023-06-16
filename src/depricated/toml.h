#ifndef TOML_HPP
#define TOML_HPP

#include "toml_original.hpp"
#include <iomanip>
#include <algorithm>

#pragma once
namespace toml {
    void error(std::string _msg);

    typedef toml::v3::node_view<toml::v3::node> node_t;
    typedef toml::node_type var_type_t;

    const toml::v3::node_type string_t =  toml::node_type::string;
    const toml::v3::node_type integer_t =  toml::node_type::integer;
    const toml::v3::node_type table_t = toml::node_type::table;
    const toml::v3::node_type none_t = toml::node_type::none;
    const toml::v3::node_type array_t = toml::node_type::array;
    const toml::v3::node_type double_t =  toml::node_type::floating_point;
    const toml::v3::node_type bool_t = toml::node_type::boolean;

    const std::string indicator = "$";

    // template<typename T>
    // void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val);


    // errors instead of setting to default
    // template<typename T>
    // void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var);


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
            std::vector<Solver>             solvers;
            std::vector<Equation>           equations;
            std::vector<Material>           materials;
            std::vector<Body>               bodys;
            std::vector<Initial_Condition>  initial_conditions;
            std::vector<Boundary_Condition> boundary_conditions;


            std::string exe;
            std::string sif;
            std::string meshDBstem;

    };

    void format_name(std::string_view _s, std::string& _out);

    int vec_to_arr(std::vector<std::string>& _vec, char**& _arr);

    void join_array_sparta(std::string _caller, std::string _name, toml::node_t val, std::vector<std::string>& _buffer);



    void table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<std::string>& _var, std::string _sep = "=", bool _specify_size = false);

    bool is_int(const std::string& s);

    template<typename T>
    void id_table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<T>& _var, std::string _sep = "=", bool _specify_size = false);
}

#endif