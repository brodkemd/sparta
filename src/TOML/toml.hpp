#ifndef TOML_HPP
#define TOML_HPP

#include "toml_original.hpp"
#include <iomanip>

namespace toml {
    void error(std::string _msg) { throw _msg; }

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
            error(ss.str());
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
            error("no value provided for \"" + _name + "\"");
        }

        if (_val.type() != _type) {
            std::stringstream ss;
            ss << "\"" << _name << "\" has wrong type, " << _val.type() << ", must be " << _type << " type";
            error(ss.str());
        }

        //std::cout << _val << "\n";
        std::optional<T> _var_opt = _val.value<T>();

        if (_var_opt.has_value()) {
            _var = _var_opt.value();
        } else {
            error("no value provided for \"" + _name + "\"");
        }
    }




    class Base {
        private:
            std::string tab = "  ";
            std::string end = "End";

        protected:
            std::string name;

        public:
            int id = INT_MIN;
            std::vector<std::string> contents;

            void join(std::string& _buffer) {
                // for good measure
                if (id == INT_MIN)
                    _buffer = name + "\n";
                else
                    _buffer = name + " "  + std::to_string(this->id) + "\n";

                for (std::string it : this->contents) 
                    _buffer+=(this->tab + it + "\n");

                _buffer+=this->end;
            }
    };


    class Header : public Base {
        public:
            Header() {
                this->name = "Header";
            }
    };


    class Constants : public Base {
        public:
            Constants() {
                this->name = "Constants";
            }
    };


    class Simulation : public Base {
        public:
            Simulation() {
                this->name = "Simulation";
            }
    };


    class Solver : public Base {
        public:
            Solver() {
                this->name = "Solver";
            }
    };


    class Equation : public Base {
        public:
            Equation() {
                this->name = "Equation";
            }
    };


    class Material : public Base {
        public:
            Material() {
                this->name = "Material";
            }
    };


    class Body : public Base {
        public:
            Body() {
                this->name = "Body";
            }
    };


    class Initial_Condition : public Base {
        public:
            Initial_Condition() {
                this->name = "Initial Condition";
            }
    };


    class Boundary_Condition : public Base {
        public:
            Boundary_Condition() {
                this->name = "Boundary Condition";
            }
    };


    class Elmer {
        private:
            std::string sep = "\n\n";

        public:
            Elmer() {
                this->header = Header();
                this->simulation = Simulation();
                this->constants = Constants();
            }

            void join(std::string& _buffer) {
                _buffer.clear();

                std::string _temp;
                
                this->header.join(_temp);
                _buffer += (_temp + this->sep);
                
                this->simulation.join(_temp);
                _buffer += (_temp + this->sep);
                
                this->constants.join(_temp);
                _buffer += (_temp + this->sep);
                
                for (int i = 0; i < this->solvers.size(); i++) {
                    solvers[i].join(_temp);
                    _buffer += (_temp + this->sep);
                }

                for (int i = 0; i < this->equations.size(); i++) {
                    equations[i].join(_temp);
                    _buffer += (_temp + this->sep);
                }

                for (int i = 0; i < this->materials.size(); i++) {
                    materials[i].join(_temp);
                    _buffer += (_temp + this->sep);
                }

                for (int i = 0; i < this->bodys.size(); i++) {
                    bodys[i].join(_temp);
                    _buffer += (_temp + this->sep);
                }

                for (int i = 0; i < this->initial_conditions.size(); i++) {
                    initial_conditions[i].join(_temp);
                    _buffer += (_temp + this->sep);
                }

                for (int i = 0; i < this->boundary_conditions.size(); i++) {
                    boundary_conditions[i].join(_temp);
                    _buffer += (_temp);
                }
            }
            

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

    void format_name(std::string_view _s, std::string& _out) {
        _out.clear();
        _out = _s;
        std::replace(_out.begin(), _out.end(), '_', ' ');
    }


    void table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<std::string>& _var, std::string _sep = "=") {
        
        if (__tbl.type() != toml::node_type::table)
            error(_caller + " must be given a table");

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
                            error("could not get value for " + name);

                        
                    } else if (toml::integer_t == elem.type()) {
                        opt_int = elem.value<double>();
                        if (opt_int.has_value()) {
                            arr_str += (" " + std::to_string(opt_int.value()));

                        } else
                            error("could not get value for " + name);

                    } else if (toml::double_t == elem.type()) {
                        opt_double = elem.value<double>();
                        if (opt_double.has_value()) {
                            double_converter << opt_double.value();
                            arr_str += (" " + double_converter.str());

                        } else
                            error("could not get value for " + name);
                    
                    } else if (toml::bool_t == elem.type()) {
                        opt_bool = elem.value<bool>();
                        if (opt_bool.has_value()) {
                            if (opt_bool.value()) 
                                arr_str += (" True");
                            else 
                                arr_str += (" False");
                            
                        } else
                            error("could not get value for " + name);
                    } else {
                        error("invalid type in " + name);

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
                    error("could not get value for " + name);

                
            } else if (toml::integer_t == it.second.type()) {
                opt_int = it.second.value<double>();
                if (opt_int.has_value()) {
                    _var.push_back(
                        param_name + " " + _sep + " " + std::to_string(opt_int.value())
                    );
                } else
                    error("could not get value for " + name);

            } else if (toml::double_t == it.second.type()) {
                opt_double = it.second.value<double>();
                if (opt_double.has_value()) {
                    double_converter << opt_double.value();
                    _var.push_back(
                        param_name + " "  + _sep + " " + double_converter.str()
                    );
                } else
                    error("could not get value for " + name);

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
                    error("could not get value for " + name);

            } else {
                std::stringstream ss;
                ss << "Invalid type for " << name << ", has type " << it.second.type();
                error(ss.str());
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
            error(_caller + " must be given a table");

        std::cout << "Running id table for: " << _caller << "\n";

        // converting to a table
        toml::table* _tbl = __tbl.as<toml::table>();

        std::string _name;
        std::string _key;

        for (auto it : (*_tbl)) {
            _key = it.first.str();
            if (!(is_int(_key)))
                error("elements for " + _caller + " must be a positive integer, not " + _key);

            _name = _caller + "." + _key;

            T inst = T();
            inst.id = std::stoi(_key);

            table_value_parser(_name, (*_tbl)[_key], inst.contents);
            

            _var.push_back(inst);

            std::cout << " <> "<< it.first.str() << "\n";
        }

    }
}

#endif