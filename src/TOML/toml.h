#ifndef TOML_H
#define TOML_H

#ifndef FLERR
#define FLERR __FILE__,__LINE__
#endif

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "../error.h"
#include "../pointers.h"

#include "toml.hpp"

namespace config {
    typedef toml::v3::node_view<toml::v3::node> node_t;
    typedef std::pair<std::string, void (*)(std::string, node_t)> dict_item_t;
    typedef std::vector<dict_item_t> dict_t;
    typedef toml::node_type var_type_t; 

    const auto string_t =  toml::node_type::string;
    const auto integer_t =  toml::node_type::integer;
    const auto none_t = toml::node_type::none;

    Error* error;

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
            config::error->all(FLERR, ss.str().c_str());
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
            error->all(FLERR, ("no value provided for \"" + _name + "\"").c_str());
            return;
        }

        if (_val.type() != _type) {
            std::stringstream ss;
            ss << "\"" << _name << "\" has wrong type, " << _val.type() << ", must be " << _type << " type";
            config::error->all(FLERR, ss.str().c_str());
        }

        std::optional<T> _var_opt = _val.value<T>();

        if (_var_opt.has_value()) {
            _var = _var_opt.value();
        } else {
            config::error->all(FLERR, ("no value provided for \"" + _name + "\"").c_str());
        }
    }

    /**
     * executes a table
    */
    void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict_t& _options) {
        if (_tbl.type() != toml::node_type::table) 
            error->all(FLERR, (_caller + " must be given a table").c_str());

        if (_caller.length())
            _caller = _caller + "." + _name;
        else
            _caller = _name;

        for (dict_item_t it : _options) {
            node_t val = _tbl[it.first];
            it.second(_caller, val);
            _tbl.erase(it.first);
        }

        if (_tbl.size()) {
            std::string msg;
            if (_caller.length())
                msg = "Invalid args in section \"" + _caller + "\":\n";
            else
                msg = "Invalid section:\n";
            for (auto it : _tbl) {
                msg+=("-> " + std::string(it.first.str()) + "\n");
            }
            error->all(FLERR, msg.c_str());
        }
    }

    /**
     * executes a table that wrapped as a node_t type
    */
    void run_table(std::string _caller, std::string _name, node_t& __tbl, dict_t& _options) {

        if (__tbl.type() != toml::node_type::table) 
            error->all(FLERR, (_caller + " must be given a table").c_str());

        // converting to a table
        toml::table* _tbl = __tbl.as<toml::table>();

        // calling the other definition because there is nothing else unique needed
        run_table(_caller, _name, (*_tbl), _options);
    }
}
#endif