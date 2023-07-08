#ifndef TOML_H
#define TOML_H

#include <string>
#include <vector>

#include "python_config.h"


namespace toml {
    void error(std::string _msg) { throw _msg; }

    const double noDouble   = std::numeric_limits<double>::max();
    const int noInt         = std::numeric_limits<int>::max();
    const char noString     = std::char_traits<char>::eof();
    const bool noBool       = 2;

    int vec_to_arr(std::vector<std::string>& _vec, char**& _arr) {
        const int _size = _vec.size();
        _arr = new char*[_size];
        for (int i = 0; i < _size; i++) _arr[i] = (char*)_vec[i].c_str();
        return _size;
    }

    class handler {
        public:
            handler(std::string _file) {
                Py_Initialize();

                PyRun_SimpleString(PYTHON_STRING); this->_err();

                PyObject* u_name    = PyUnicode_FromString("__main__");
                PyObject* m         = PyImport_GetModule(u_name);
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
                }
            }

            ~handler() {
                Py_Finalize();
                Py_DECREF(pValue);
            }

            template<typename T>
            void get_at_path(T& _var, std::string _path, bool _strict) {
                this->pArgs = PyTuple_New(2);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyTuple_SetItem(pArgs, 1, PyBool_FromLong((long)_strict));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();
                this->_handle_var(_var, out, _path);
                Py_DECREF(out);
            }


        private:
            PyObject *get_at, *pValue, *pArgs, *type, *value, *traceback;

            void _err() {
                if (PyErr_Occurred()) {
                    PyErr_Fetch(&this->type, &this->value, &this->traceback);
                    std::string _value = _PyUnicode_AsString(PyObject_Repr(this->value));
                    error(_value);
                }
            }

            void _handle_var(double& _var, PyObject *_src, std::string _path) {
                if (Py_IsNone(_src))
                    _var = noDouble;

                if (!(PyFloat_Check(_src)))
                    error(_path + " is not the correct type, must be float");
                _var = PyFloat_AsDouble(_src);
            }

            void _handle_var(std::string& _var, PyObject *_src, std::string _path) {
                if (Py_IsNone(_src))
                    _var = noString;

                if (!(PyUnicode_Check(_src)))
                    error(_path + " is not the correct type, must be string");
                _var = _PyUnicode_AsString(_src);
            }

            void _handle_var(int& _var, PyObject *_src, std::string _path) {
                if (Py_IsNone(_src))
                    _var = noInt;

                if (!(PyLong_Check(_src)))
                    error(_path + " is not the correct type, must be int");
                _var = PyLong_AsLong(_src);
            }

            void _handle_var(bool& _var, PyObject *_src, std::string _path) {
                if (Py_IsNone(_src))
                    _var = 2;

                if (!(PyBool_Check(_src)))
                    error(_path + " is not the correct type, must be int");
                
                if (PyObject_IsTrue(_src))
                    _var = true;
                else
                    _var = false;
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

            void _handle_var(std::vector<bool>& _var, PyObject *_src, std::string _path) {
                if (!(PyList_Check(_src)))
                    error(_path + " is not the correct type, must be array of ints");
                _list(_var, _src, _path);
            }
            
    };
}

#endif