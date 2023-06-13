#include "fixfea.h"
#include <algorithm>
#include <sstream>
#include <iomanip>

/**
 * executes a table
*/
template<typename dict>
void FixFea::run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options) {
    if (_tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    if (_caller.length())
        _caller = _caller + "." + _name;
    else
        _caller = _name;

    for (toml::dict_item_t it : _options) {
        toml::node_t val = _tbl[it.first];
        (this->*it.second)(_caller, val);
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
        error(FLERR, msg);
    }
}


/**
 * executes a table that wrapped as a node_t type
*/
template<typename dict>
void FixFea::run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options) {

    if (__tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    // converting to a table
    toml::table* _tbl = __tbl.as<toml::table>();

    // calling the other definition because there is nothing else unique needed
    run_table(_caller, _name, (*_tbl), _options);
}


void FixFea::handle_both(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("emi", &FixFea::handle_emi),
        std::make_pair("tsurf_file", &FixFea::handle_tsurf_file)
    };

    this->run_table(_caller, "both", tbl, options);
}


void FixFea::handle_emi(std::string _caller, toml::node_t val) {
    toml::set(_caller, "emi", toml::double_t, val, this->emi);
    if (emi <= 0.0 || emi > 1.0)
        error(FLERR, "Fix fea emissivity must be > 0.0 and <= 1");
    this->print("emi = " + std::to_string(this->emi));
}


void FixFea::handle_tsurf_file(std::string _caller, toml::node_t val) {
    toml::set(_caller, "tsurf_file", toml::string_t, val, this->tsurf_file);

    if (!(stat(this->tsurf_file.c_str(),  &sb) == 0))
        error(FLERR, "Illegal fix fea command, tsurf_file path does not exist");

    this->print("tsurf_file = " + this->tsurf_file);
}


void FixFea::handle_sparta(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("nevery", &FixFea::handle_nevery),
        std::make_pair("groupID", &FixFea::handle_groupID),
        std::make_pair("mixID", &FixFea::handle_mixID),
        std::make_pair("customID", &FixFea::handle_customID)
    };

    this->run_table(_caller, "sparta", tbl, options);
}


void FixFea::handle_nevery(std::string _caller, toml::node_t val) {
    toml::set(_caller, "nevery", toml::integer_t, val, this->nevery);
    this->print("nevery = " + std::to_string(this->nevery));
}


void FixFea::handle_groupID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "groupID", toml::string_t, val, this->groupID);
    this->print("groupID = " + this->groupID);
}


void FixFea::handle_mixID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "mixID", toml::string_t, val, this->mixID);
    this->print("mixID = " + this->mixID);
}


void FixFea::handle_customID(std::string _caller, toml::node_t val) {
    toml::set(_caller, "customID", toml::string_t, val, this->customID);
    this->print("customID = " + this->customID);
}


void FixFea::handle_elmer(std::string _caller, toml::node_t tbl) {
    toml::dict_t options = {
        std::make_pair("exe", &FixFea::handle_exe),
        std::make_pair("sif", &FixFea::handle_sif),
        std::make_pair("meshDBstem", &FixFea::handle_meshDBstem),
        std::make_pair("header", &FixFea::handle_header),
        std::make_pair("simulation", &FixFea::handle_simulation),
        std::make_pair("constants", &FixFea::handle_constants),
        std::make_pair("solver", &FixFea::handle_solver),
        std::make_pair("equation", &FixFea::handle_equation),
        std::make_pair("material", &FixFea::handle_material),
        std::make_pair("body", &FixFea::handle_body),
        std::make_pair("initial_condition", &FixFea::handle_initial_condition),
        std::make_pair("boundary_condition", &FixFea::handle_boundary_condition)
    };

    this->run_table(_caller, "elmer", tbl, options);
}


void FixFea::handle_exe(std::string _caller, toml::node_t val) {
    toml::set(_caller, "exe", toml::string_t, val, this->elmer.exe);
    // Calls the function with path as argument
    // If the file/directory exists at the path returns 0
    // If block executes if path exists
    if (!(stat(this->elmer.exe.c_str(),  &sb) == 0))
        error(FLERR, "Illegal fix fea command, exe path does not exist");
    
    this->print("exe = " + this->elmer.exe);
}

void FixFea::handle_sif(std::string _caller, toml::node_t val) {
    toml::set(_caller, "sif", toml::string_t, val, this->elmer.sif);
    if (!(stat(this->elmer.sif.c_str(),  &sb) == 0))
        error(FLERR, "Illegal fix fea command, sif path does not exist");
    this->print("sif = " + this->elmer.sif);
}

void FixFea::handle_meshDBstem(std::string _caller, toml::node_t val) {
    toml::set(_caller, "meshDBstem", toml::string_t, val, this->elmer.meshDBstem);

    std::string exts[4] = {"boundary", "nodes", "header", "elements"}; // list of component file extensions
    for (int i = 0; i < 4; i++) {
        if (!(stat((this->elmer.meshDBstem + "." + exts[i]).c_str(),  &sb) == 0))
            error(FLERR, "Illegal fix fea command, mesh database incomplete, " + (this->elmer.meshDBstem + "." + exts[i]) + " does not exist");
    }
    this->print("meshDBstem = " + this->elmer.meshDBstem);
}

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

        std::cout << "  running: " << name << "\n";

        
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
            double_converter << opt_double.value();
            if (opt_double.has_value()) {
                _var.push_back(
                    param_name + " "  + _sep + " " + double_converter.str()
                );
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


void FixFea::handle_header(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".header";
    table_value_parser(_caller, __tbl, this->elmer.header.contents, "");
}


void FixFea::handle_simulation(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".simulation";
    table_value_parser(_caller, __tbl, this->elmer.simulation.contents);
}


void FixFea::handle_constants(std::string _caller, toml::node_t __tbl) {
    _caller = _caller + ".constants";
    table_value_parser(_caller, __tbl, this->elmer.constants.contents);
    return;
}

void FixFea::handle_solver(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".solver";

    //this->run_table(_caller, "solver", tbl, options);
    return;
}

void FixFea::handle_equation(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".equation";

    //this->run_table(_caller, "equation", tbl, options);
    return;
}

void FixFea::handle_material(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".material";

    //this->run_table(_caller, "material", tbl, options);
    return;
}


void FixFea::handle_body(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".body";

    //this->run_table(_caller, "body", tbl, options);
    return;
}

void FixFea::handle_initial_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".initial_condition";

    //this->run_table(_caller, "initial_condition", tbl, options);
    return;
}

void FixFea::handle_boundary_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".boundary_condition";

    //this->run_table(_caller, "boundary_condition", tbl, options);
    return;
}