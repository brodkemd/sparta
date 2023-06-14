#include "fixfea.h"
#include <sstream>

#include "util.hpp"


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
}


void FixFea::handle_solver(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".solver";
    id_table_value_parser(_caller, tbl, this->elmer.solvers);
}


void FixFea::handle_equation(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".equation";
    id_table_value_parser(_caller, tbl, this->elmer.equations);
}


void FixFea::handle_material(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".material";
    id_table_value_parser(_caller, tbl, this->elmer.materials);
}


void FixFea::handle_body(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".body";
    id_table_value_parser(_caller, tbl, this->elmer.bodys);
}


void FixFea::handle_initial_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".initial_condition";
    id_table_value_parser(_caller, tbl, this->elmer.initial_conditions);
}


void FixFea::handle_boundary_condition(std::string _caller, toml::node_t tbl) {
    _caller = _caller + ".boundary_condition";
    id_table_value_parser(_caller, tbl, this->elmer.boundary_conditions);
}




/*-------------------------------------------------


Utility functions


-------------------------------------------------*/

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


