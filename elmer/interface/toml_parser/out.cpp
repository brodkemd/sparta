#include "out.h"


template<typename T>
void toml::set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val) {
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
void toml::set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var) {
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

    std::cout << _val << "\n";
    std::optional<T> _var_opt = _val.value_exact<T>();

    if (_var_opt.has_value()) {
        _var = _var_opt.value();
    } else {
        error(FLERR, "no value provided for \"" + _name + "\"");
    }
}


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


/**
 * executes a table
*/
template<typename dict>
void FixFea::eval_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options) {
    if (_tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    if (_caller.length())
        _caller = _caller + "." + _name;
    else
        _caller = _name;

    for (toml::var_type_t it : _options) {
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
        error(FLERR, msg);
    }
}


/**
 * executes a table that wrapped as a node_t type
*/
template<typename dict>
void FixFea::eval_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options) {

    if (__tbl.type() != toml::node_type::table) 
        error(FLERR, _caller + " must be given a table");

    // converting to a table
    toml::table* _tbl = __tbl.as<toml::table>();

    // calling the other definition because there is nothing else unique needed
    run_table(_caller, _name, (*_tbl), _options);
}
