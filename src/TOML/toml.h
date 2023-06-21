#ifndef TOML_H
#define TOML_H

#include <Python.h>
#include <string>
#include <vector>

#define PYTHON_STRING "from tomli import load\n\ndef get_from_file(_file) -> dict:\n    with open(_file, 'rb') as f:\n        return load(f)\n\ndata = {}\nrecurse_num = 0\nmax_recurse = 10\n\n\ndef error(_msg:str):\n    raise Exception(_msg)\n\n\ndef get_arr(_data, _name:str):\n    start_char = \"[\"\n    end_char   = \"]\"\n    \n    if start_char in _name:\n        _data = _data[_name[:_name.find(start_char)]]\n        for i in range(_name.count(start_char)):\n            index = int(_name[_name.find(start_char)+len(start_char): _name.find(end_char)].strip())\n            _data = _data[index]\n            _name = _name[_name.find(end_char)+len(end_char):]\n\n    return _data, _name\n\n\ndef get_from(_data:dict, _path:str, _orig_path:str=None):\n    if _orig_path is None:\n        _orig_path = _path\n\n    sep = \".\"\n    if sep in _path:\n        _name = _path[:_path.find(sep)]\n        _data, _name = get_arr(_data, _name)\n        # print(\"name:\", _name)\n        if len(_name) == 0:\n            return get_from(_data, _path[_path.find(sep)+len(sep):], _orig_path)\n        elif _name in _data:\n            return get_from(_data[_name], _path[_path.find(sep)+len(sep):], _orig_path)\n        else:\n            error(\"can not resolve path, \" + _orig_path)\n    else:\n        _data, _path = get_arr(_data, _path)\n        #print(\"path:\"+ _path +\"|\")\n        if len(_path) == 0: # this means is it an array value\n            return resolve_name(_data)\n        elif _path in _data:\n            if isinstance(_data[_path], str):\n                _data[_path] = resolve_name(_data[_path])\n\n            return _data[_path]\n        \n        else:\n            error(\"can not resolve path, \" + _orig_path)\n            \n    \n\ndef resolve_name(_name:str):\n    global data, recurse_num, max_recurse\n    \n    if recurse_num > max_recurse:\n        error(f\"reached max allowed recursion resolving {_name} (there is probably a circular definition somewhere)\")\n    \n    indicator = \"$\"\n    if _name.startswith(indicator):\n        recurse_num+=1;\n        _name = get_from(data, _name[len(indicator):])\n    return _name\n        \n    \n\ndef iter(_data):\n    global recurse_num\n    for item in _data:\n        recurse_num = 0\n        if isinstance(_data[item], dict):\n            iter(_data[item])\n        else:\n            if isinstance(_data[item], list):\n                str_indicies = []\n                for i, val in enumerate(_data[item]):\n                    if isinstance(val, str):\n                        str_indicies.append(i)\n\n                for i in range(len(_data[item])):\n                    val = _data[item][i]\n                    if isinstance(val, str):\n                        _data[item][i] = resolve_name(val)\n                        if _data[item][i] != val:\n                            if i in str_indicies:\n                                str_indicies.remove(i)\n                \n                if len(str_indicies):\n                    for i, val in enumerate(_data[item]):\n                        _data[item][i] = str(val)\n            \n            elif isinstance(_data[item], str):\n                _data[item] = resolve_name(_data[item])\n\ndef get_at_path(_path):\n    global data\n    return get_from(data, _path)\n\n\ndef section_to_string(_path, _include_count, _sep = \"= \"):\n    global data\n    _section = get_from(data, _path)\n\n    s = \"\"\n    for item in _section:   \n        if isinstance(_section[item], dict): error(\"Can not handle nested section\")\n        elif isinstance(_section[item], list):\n            if _include_count:\n                s += (str(item) + f\"({len(_section[item])}) {_sep}\" + \" \".join([str(val) for val in _section[item]]) + \"\\n\")\n            else:\n                s += (str(item) + f\" {_sep}\" + \" \".join([str(val) for val in _section[item]]) + \"\\n\")\n        else: s += (str(item) + f\" {_sep}\" + str(_section[item]) + \"\\n\")\n    \n    return s\n\n\ndef Main(_file):\n    global data\n    data = get_from_file(_file)\n    iter(data)\n    \n    #print(section_to_string(\"elmer.constants\", False, \"\"))\n    # print(get_at_path(\"sparta.compute[0]\"))\n\n# Main(\"/home/marekbrodke/Documents/C++/play/fea.toml\")\n\n# from json import dumps\n# with open(\"out.json\", 'w') as f: f.write(dumps(data, indent=4))\n"


namespace toml {
    void error(std::string _msg);
    // toml::value handle_arr(toml::value _data, std::string& _path, std::string _orig_path);
    // template <typename First, typename... Args>
    // toml::value preprocess(First first, Args... args);
    // toml::value resolve_name(toml::string _name);
    // toml::value get_from(std::string _orig_path = "");

    class handler {
        public:
            handler(std::string _file) {
                Py_Initialize();

                PyRun_SimpleString(PYTHON_STRING);
                PyObject* u_name = PyUnicode_FromString("__main__");
                PyObject* m = PyImport_GetModule(u_name);
                PyObject* Main_func = PyObject_GetAttrString(m, "Main");

                this->pArgs = PyTuple_New(1);

                if (Main_func && PyCallable_Check(Main_func)) {
                    pValue = PyUnicode_FromString(_file.c_str());

                    PyTuple_SetItem(pArgs, 0, pValue);
                    PyObject_CallObject(Main_func, pArgs);
                    this->_err();

                    this->get_at = PyObject_GetAttrString(m, "get_at_path");
                    this->_err();

                    if (!(this->get_at && PyCallable_Check(this->get_at))) {  
                        error("get_at_path function does not exists or is not callable");
                    }

                    this->to_string = PyObject_GetAttrString(m, "section_to_string");
                    this->_err();

                    if (!(this->to_string && PyCallable_Check(this->to_string))) {  
                        error("section_to_string function does not exists or is not callable");
                    }
                }
            }

            ~handler() {
                Py_Finalize();
                Py_DECREF(pValue);
            }

            template<typename T>
            void get_at_path(T& _var, std::string _path) {
                this->pArgs = PyTuple_New(1);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();
                this->_handle_var(_var, out, _path);
                Py_DECREF(out);
            }

            std::string section_to_string(std::string _path, bool _include_count, std::string _sep = "= ") {
                this->pArgs = PyTuple_New(3);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyTuple_SetItem(pArgs, 1, PyBool_FromLong(_include_count));
                PyTuple_SetItem(pArgs, 2, PyUnicode_FromString(_sep.c_str()));

                PyObject* out = PyObject_CallObject(this->to_string, pArgs);
                this->_err();
                
                std::string _to_return;
                _handle_var(_to_return, out, _path);
                return _to_return;
            }

        private:
            PyObject *get_at, *to_string, *pValue, *pArgs, *type, *value, *traceback;

            void _err() {
                if (PyErr_Occurred()) {
                    PyErr_Fetch(&this->type, &this->value, &this->traceback);
                    std::string _value = _PyUnicode_AsString(PyObject_Repr(this->value));
                    // std::cout << _value  << "\n";
                    
                    error(_value);
                }
            }



            void _handle_var(double& _var, PyObject *_src, std::string _path) {
                if (!(PyFloat_Check(_src)))
                    error(_path + " is not the correct type, must be float");
                _var = PyFloat_AsDouble(_src);
            }

            void _handle_var(std::string& _var, PyObject *_src, std::string _path) {
                if (!(PyUnicode_Check(_src)))
                    error(_path + " is not the correct type, must be string");
                _var = _PyUnicode_AsString(_src);
            }

            void _handle_var(int& _var, PyObject *_src, std::string _path) {
                if (!(PyLong_Check(_src)))
                    error(_path + " is not the correct type, must be int");
                _var = PyLong_AsLong(_src);
            }


            template<typename T>
            void _list(std::vector<T>& _var, PyObject *_src, std::string _path) {
                _var.clear();
                T _temp;
                for (int i = 0; i < PyList_Size(_src); i++) {
                    _handle_var(_temp, PyList_GetItem(_src, i), _path + "[" + std::to_string(i) + "]");
                    _var.push_back(_temp);
                }
            }

            void _handle_var(std::vector<std::string>& _var, PyObject *_src, std::string _path) {
                if (!(PyList_Check(_src)))
                    error(_path + " is not the correct type, must be array of strings");

                _list(_var, _src, _path);
            }

            void _handle_var(std::vector<double>& _var, PyObject *_src, std::string _path) {
                if (!(PyList_Check(_src)))
                    error(_path + " is not the correct type, must be array of floats");
                
                _list(_var, _src, _path);
            }

            void _handle_var(std::vector<int>& _var, PyObject *_src, std::string _path) {
                if (!(PyList_Check(_src)))
                    error(_path + " is not the correct type, must be array of ints");
                _list(_var, _src, _path);
            }
    };

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
}


#endif