#ifndef UTIL_HPP
#define UTIL_HPP

#include "toml.h"
#include "classes.h"

#include <algorithm>
#include <iomanip>


void format_name(std::string_view _s, std::string& _out) {
    _out.clear();
    _out = _s;
    std::replace(_out.begin(), _out.end(), '_', ' ');
}


void table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<std::string>& _var, std::string _sep = "=") {
    
    if (__tbl.type() != toml::node_type::table)
        error(FLERR, _caller + " must be given a table");

    std::cout << "Running table for: " << _caller << "\n";

    // converting to a table
    toml::table* _tbl = __tbl.as<toml::table>();

    std::string param_name;
    std::string name;

    std::optional<std::string> opt_str;
    std::optional<int64_t> opt_int;
    std::optional<double> opt_double;
    std::optional<bool> opt_bool;
    std::ostringstream double_converter;
    std::string arr_str;
    toml::array arr;

    double_converter << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);

    _var.clear();

    for (auto it : (*_tbl)) {
        double_converter.str("");
        name = it.first.str();
        name = _caller + "." + name;
        format_name(it.first.str(), param_name);

        //std::cout << "  running: " << name << "\n";

        if (toml::array_t == it.second.type()) {
            arr = *it.second.as_array();
            arr_str = param_name + "(" + std::to_string(arr.size()) + ") " + _sep;

            for (auto&& elem : arr) {
                double_converter.str("");
                if (toml::string_t == elem.type()) {
                    opt_str = elem.value<std::string>();
                    if (opt_str.has_value()) {
                        arr_str += (" \"" + opt_str.value() + "\"");
                    } else
                        error(FLERR, "could not get value for " + name);

                    
                } else if (toml::integer_t == elem.type()) {
                    opt_int = elem.value<double>();
                    if (opt_int.has_value()) {
                        arr_str += (" " + std::to_string(opt_int.value()));

                    } else
                        error(FLERR, "could not get value for " + name);

                } else if (toml::double_t == elem.type()) {
                    opt_double = elem.value<double>();
                    if (opt_double.has_value()) {
                        double_converter << opt_double.value();
                        arr_str += (" " + double_converter.str());

                    } else
                        error(FLERR, "could not get value for " + name);
                
                } else if (toml::bool_t == elem.type()) {
                    opt_bool = elem.value<bool>();
                    if (opt_bool.has_value()) {
                        if (opt_bool.value()) 
                            arr_str += (" True");
                        else 
                            arr_str += (" False");
                        
                    } else
                        error(FLERR, "could not get value for " + name);
                } else {
                    error(FLERR, "invalid type in " + name);

                }
            }

            _var.push_back(arr_str);

        } else if (toml::string_t == it.second.type()) {
            opt_str = it.second.value<std::string>();
            if (opt_str.has_value()) {
                _var.push_back(
                    param_name + " " + _sep + " \"" + opt_str.value() + "\""
                );
            } else
                error(FLERR, "could not get value for " + name);

            
        } else if (toml::integer_t == it.second.type()) {
            opt_int = it.second.value<double>();
            if (opt_int.has_value()) {
                _var.push_back(
                    param_name + " " + _sep + " " + std::to_string(opt_int.value())
                );
            } else
                error(FLERR, "could not get value for " + name);

        } else if (toml::double_t == it.second.type()) {
            opt_double = it.second.value<double>();
            if (opt_double.has_value()) {
                double_converter << opt_double.value();
                _var.push_back(
                    param_name + " "  + _sep + " " + double_converter.str()
                );
            } else
                error(FLERR, "could not get value for " + name);

        } else if (toml::bool_t == it.second.type()) {
            opt_bool = it.second.value<bool>();
            if (opt_bool.has_value()) {
                if (opt_bool.value()) {
                    _var.push_back(
                        param_name + " "  + _sep + " " + "True"
                    );
                } else {
                    _var.push_back(
                        param_name + " "  + _sep + " " + "False"
                    );
                }
            } else
                error(FLERR, "could not get value for " + name);

        } else {
            std::stringstream ss;
            ss << "Invalid type for " << name << ", has type " << it.second.type();
            error(FLERR, ss.str());
        }

        // std::cout << param_name << " : " << it.second.type() << "\n";
    }

    for (auto it : _var) {
        std::cout << "=> " << it << "\n";
    }
}

bool is_int(const std::string& s) {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

template<typename T>
void id_table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<T>& _var, std::string _sep = "=") {
    if (__tbl.type() != toml::node_type::table)
        error(FLERR, _caller + " must be given a table");

    std::cout << "Running id table for: " << _caller << "\n";

    // converting to a table
    toml::table* _tbl = __tbl.as<toml::table>();

    std::string _name;
    std::string _key;

    for (auto it : (*_tbl)) {
        _key = it.first.str();
        if (!(is_int(_key)))
            error(FLERR, "elements for " + _caller + " must be a positive integer, not " + _key);

        _name = _caller + "." + _key;

        T inst = T();
        inst.id = std::stoi(_key);

        table_value_parser(_name, (*_tbl)[_key], inst.contents);
        

        _var.push_back(inst);

        std::cout << " <> "<< it.first.str() << "\n";
    }

}

#endif