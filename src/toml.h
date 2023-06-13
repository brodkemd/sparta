#pragma once
#ifndef TOML_H
#define TOML_H

#include <vector>
#include <string>
#include <utility>

#include "./TOML/toml.hpp"

#pragma once
#ifndef TOML_NAMESPACE
#define TOML_NAMESPACE

namespace toml {
    typedef toml::v3::node_view<toml::v3::node> node_t;
    
    //typedef std::tuple<std::string,void (*)(std::string, node_t)> dict_item_t;
    //typedef std::vector<dict_item_t> dict_t;
    typedef toml::node_type var_type_t; 

    const toml::v3::node_type string_t =  toml::node_type::string;
    const toml::v3::node_type integer_t =  toml::node_type::integer;
    const toml::v3::node_type double_t =  toml::node_type::floating_point;
    const toml::v3::node_type none_t = toml::node_type::none;
    const std::string sep = ".";

    void error(std::string msg); // { throw msg; }


    // toml::table defaults;

    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val);

    // errors instead of setting to default
    template<typename T>
    void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var);

    /**
     * executes a table
    */
    // void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict_t& _options);

    /**
     * executes a table that wrapped as a node_t type
    */
    // void run_table(std::string _caller, std::string _name, node_t& __tbl, dict_t& _options);

}
#endif
#endif