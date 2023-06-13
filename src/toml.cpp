// #include "toml.h"


// void toml::error(std::string msg) { throw msg; }

// template<typename T>
// void toml::set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val) {
//     // if the _val is none type, just set to default value immediately
//     if (_val.type() == none_t) {
//         _var = _default_val;
//         return;
//     }

//     if (_val.type() != _type) {
//         // adding the calling hierarchy to the name
//         _name = _caller + "." + _name;

//         // making message to user
//         std::stringstream ss;
//         ss << "\"" << _name << "\" has wrong type, " << _val.type();
//         ss << ", must be " << _type << " type";
//         error(ss.str());
//     }

//     std::optional<T> _var_opt = _val.value<T>();

//     if (_var_opt.has_value())
//         _var = _var_opt.value();
//     else 
//         _var = _default_val;
// }

// // errors instead of setting to default
// template<typename T>
// void toml::set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var) {
//     _name = _caller + "." + _name;

//     // if the _val is none type, just set to default value immediately
//     if (_val.type() == none_t) {
//         error("no value provided for \"" + _name + "\"");
//     }

//     if (_val.type() != _type) {
//         std::stringstream ss;
//         ss << "\"" << _name << "\" has wrong type, " << _val.type() << ", must be " << _type << " type";
//         error(ss.str());
//     }

//     // std::cout << _val << "\n";
//     std::optional<T> _var_opt = _val.value_exact<T>();

//     if (_var_opt.has_value()) {
//         _var = _var_opt.value();
//     } else {
//         error("no value provided for \"" + _name + "\"");
//     }
// }

// toml::table defaults;


// void _iter(toml::node_view<toml::node> _data, std::string _to_find, node_t& _var);

// void _iter(toml::table _data, std::string _to_find, node_t& _var);

// void _iter(toml::node_t _data, std::string _to_find, node_t& _var) {
//     if (_data.type() != toml::node_type::table) 
//         error(" must be given a table");

//     // converting to a table
//     toml::table* _tbl = (_data.as<toml::table>());
//     toml::_iter((*_tbl), _to_find, _var);
// }



// template<typename T>
// void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, bool _use_default_val) {
//     std::cout << "---> setting\n";
//     _name = _caller + "." + _name;

//     // if the _val is none type, just set to default value immediately
//     std::optional<T> temp;
//     node_t temp_var;

//     if (_use_default_val) {
//         std::cout << "---> using default\n";
//         if (_val.type() == none_t) {
//             std::cout << "---> caught none, using default\n";
//             _iter(toml::defaults, _name, temp_var);
            
//             if (temp_var.type() != _type)
//                 error(FLERR, "invalid type in defaults for "+_name);

//             std::cout << "---> getting default value\n";
//             temp = temp_var.value<T>();
//             if (temp.has_value()) {
//                 _var = temp.value();
//                 std::cout << _var << "\n";
//             } else
//                 error(FLERR, "no value retrieved from defaults for \"" + _name + "\"");

//             return;
//         }
//     }

//     if (_val.type() != _type) {
//         // adding the calling hierarchy to the name
        

//         // making message to user
//         std::stringstream ss;
//         ss << "\"" << _name << "\" has wrong type, " << _val.type();
//         ss << ", must be " << _type << " type";
//         error(FLERR, ss.str());
//     }

//     std::cout << "---> getting value\n";
//     temp = _val.value<T>();
//     if (temp.has_value()) {
//         std::cout << "---> setting value\n";
//         _var = temp.value();
//     } else { 
//         if (_use_default_val) {
//             std::cout << "---> getting default value\n";
//             _iter(toml::defaults, _name, temp_var);
//             std::cout << "after\n";
//             std::cout << temp_var << "\n";
//             //if (temp_var.type() != _type)
//             //   error(FLERR, "invalid type in defaults for "+_name);

//             std::cout << "---> setting with default value\n";
//             temp = temp_var.value<T>();
//             if (temp.has_value())
//                 _var = temp.value();
//             else
//                 error(FLERR, "no value retrieved from defaults for \"" + _name + "\"");

//             return;
//         } else {
//             error(FLERR, "no value provided for \"" + _name + "\"");
//         }
//     }
// }

// errors instead of setting to default
// template<typename T>
// void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var);

// /**
//  * executes a table
// */
// void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict_t& _options);

// /**
//  * executes a table that wrapped as a node_t type
// */
// void run_table(std::string _caller, std::string _name, node_t& __tbl, dict_t& _options);

// template<typename T>
// void set(std::string _caller, std::string _name, var_type_t _type, node_t _val, T& _var, T _default_val);
