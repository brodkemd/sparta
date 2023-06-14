#ifndef TOML_H
#define TOML_H

#include <iostream>
#include <string>
#include "toml.hpp"

#define FLERR __FILE__,__LINE__


void error(std::string file, int line, std::string msg);


namespace toml {
    typedef toml::v3::node_view<toml::v3::node> node_t;
    typedef toml::node_type var_type_t;

    const toml::v3::node_type string_t =  toml::node_type::string;
    const toml::v3::node_type integer_t =  toml::node_type::integer;
    const toml::v3::node_type table_t = toml::node_type::table;
    const toml::v3::node_type none_t = toml::node_type::none;
    const toml::v3::node_type array_t = toml::node_type::array;
    const toml::v3::node_type double_t =  toml::node_type::floating_point;
    const toml::v3::node_type bool_t = toml::node_type::boolean;

    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val) {
        // if the _val is none type, just set to default value immediately
        if (_val.type() == none_t) {
            _var = _default_val;
            return;
        }

        if (_val.type() != _type) {
            // adding the calling hierarchy to the name
            _name = _caller + "." + _name;

            // making message to user
            std::stringstream ss;
            ss << "\"" << _name << "\" has wrong type, " << _val.type();
            ss << ", must be " << _type << " type";
            error(FLERR, ss.str());
        }

        std::optional<T> _var_opt = _val.value<T>();

        if (_var_opt.has_value())
            _var = _var_opt.value();
        else
            _var = _default_val;
    }


    // errors instead of setting to default
    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var) {
        _name = _caller + "." + _name;

        // if the _val is none type, just set to default value immediately
        if (_val.type() == none_t) {
            error(FLERR, "no value provided for \"" + _name + "\"");
        }

        if (_val.type() != _type) {
            std::stringstream ss;
            ss << "\"" << _name << "\" has wrong type, " << _val.type() << ", must be " << _type << " type";
            error(FLERR, ss.str());
        }

        //std::cout << _val << "\n";
        std::optional<T> _var_opt = _val.value<T>();

        if (_var_opt.has_value()) {
            _var = _var_opt.value();
        } else {
            error(FLERR, "no value provided for \"" + _name + "\"");
        }
    }
}



#endif