#ifndef TOML_H
#define TOML_H

#include "python_config.h"
#include "util.hpp"


namespace toml {
    const util::double_t noDouble       = std::numeric_limits<util::double_t>::min();
    const util::int_t noInt             = std::numeric_limits<int>::min();
    const util::string_t noString  = std::to_string((char)std::char_traits<char>::eof());
    const util::bool_t noBool           = 2;

    /* ---------------------------------------------------------------------- */

    class handler {
        public:
            handler(util::string_t _file) {
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
                        UERR("get_at_path function does not exists or is not callable");
                    }
                }
            }

            /* ---------------------------------------------------------------------- */

            ~handler() {
                Py_Finalize();
                Py_DECREF(pValue);
            }

            /* ---------------------------------------------------------------------- */

            template<typename T>
            void getAtPath(T& _var, util::string_t _path, util::bool_t _strict) {
                this->pArgs = PyTuple_New(2);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyTuple_SetItem(pArgs, 1, PyBool_FromLong((long)_strict));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();
                this->_handleVar(_var, out, _path);
                Py_DECREF(out);
            }

            // template<typename T>
            // void getAtPathOrRef(T& _var, T& _ref, util::string_t _path) {
            //     this->pArgs = PyTuple_New(2);
            //     PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
            //     PyTuple_SetItem(pArgs, 1, PyBool_FromLong((long)_strict));
            //     PyObject* out = PyObject_CallObject(this->get_at, pArgs);
            //     this->_err();
            //     this->_handleVar(_var, out, _path);
            //     Py_DECREF(out);
            // }


        private:
            PyObject *get_at, *pValue, *pArgs, *type, *value, *traceback;

            /* ---------------------------------------------------------------------- */

            void _err() {
                if (PyErr_Occurred()) {
                    PyErr_Fetch(&this->type, &this->value, &this->traceback);
                    util::string_t _value = _PyUnicode_AsString(PyObject_Repr(this->value));
                    UERR(_value);
                }
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(util::double_t& _var, PyObject *_src, util::string_t _path) {
                if (Py_IsNone(_src)) {
                    _var = noDouble;
                    return;
                }

                if (!(PyFloat_Check(_src)))
                    UERR(_path + " is not the correct type, must be float");
                _var = PyFloat_AsDouble(_src);
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(util::string_t& _var, PyObject *_src, util::string_t _path) {
                if (Py_IsNone(_src)) {
                    _var = noString;
                    return;
                }

                if (!(PyUnicode_Check(_src)))
                    UERR(_path + " is not the correct type, must be string");
                _var = _PyUnicode_AsString(_src);
            }

            // /* ---------------------------------------------------------------------- */

            void _handleVar(int& _var, PyObject *_src, util::string_t _path) {
                if (Py_IsNone(_src)) {
                    _var = noInt;
                    return;
                }

                if (!(PyLong_Check(_src)))
                    UERR(_path + " is not the correct type, must be int");
                _var = PyLong_AsLong(_src);
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(util::int_t& _var, PyObject *_src, util::string_t _path) {
                if (Py_IsNone(_src)) {
                    _var = noInt;
                    return;
                }

                if (!(PyLong_Check(_src)))
                    UERR(_path + " is not the correct type, must be int");
                _var = PyLong_AsLong(_src);
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(util::bool_t& _var, PyObject *_src, util::string_t _path) {
                if (Py_IsNone(_src)) {
                    _var = noBool;
                    return;
                }

                if (!(PyBool_Check(_src)))
                    UERR(_path + " is not the correct type, must be bool");
                
                if (PyObject_IsTrue(_src))
                    _var = true;
                else
                    _var = false;
            }

            // /* ---------------------------------------------------------------------- */

            // void _handleVarOr(util::double_t& _var, util::double_t& _or, PyObject *_src, util::string_t _path) {
            //     if (Py_IsNone(_src)) {
            //         _var = &_or;
            //         return;
            //     }

            //     if (!(PyFloat_Check(_src)))
            //         UERR(_path + " is not the correct type, must be float");
            //     _var = PyFloat_AsDouble(_src);
            // }

            // /* ---------------------------------------------------------------------- */

            // void _handleVarOr(util::string_t& _var, PyObject *_src, util::string_t _path) {
            //     if (Py_IsNone(_src)) {
            //         _var = noString;
            //         return;
            //     }

            //     if (!(PyUnicode_Check(_src)))
            //         UERR(_path + " is not the correct type, must be string");
            //     _var = _PyUnicode_AsString(_src);
            // }

            // // /* ---------------------------------------------------------------------- */

            // void _handleVarOr(int& _var, PyObject *_src, util::string_t _path) {
            //     if (Py_IsNone(_src)) {
            //         _var = noInt;
            //         return;
            //     }

            //     if (!(PyLong_Check(_src)))
            //         UERR(_path + " is not the correct type, must be int");
            //     _var = PyLong_AsLong(_src);
            // }

            // /* ---------------------------------------------------------------------- */

            // void _handleVarOr(util::int_t& _var, PyObject *_src, util::string_t _path) {
            //     if (Py_IsNone(_src)) {
            //         _var = noInt;
            //         return;
            //     }

            //     if (!(PyLong_Check(_src)))
            //         UERR(_path + " is not the correct type, must be int");
            //     _var = PyLong_AsLong(_src);
            // }

            // /* ---------------------------------------------------------------------- */

            // void _handleVarOr(util::bool_t& _var, PyObject *_src, util::string_t _path) {
            //     if (Py_IsNone(_src)) {
            //         _var = noBool;
            //         return;
            //     }

            //     if (!(PyBool_Check(_src)))
            //         UERR(_path + " is not the correct type, must be bool");
                
            //     if (PyObject_IsTrue(_src))
            //         _var = true;
            //     else
            //         _var = false;
            // }

            /* ---------------------------------------------------------------------- */

            template<typename T>
            void _list(std::vector<T>& _var, PyObject *_src, util::string_t _path) {
                _var.clear();
                T _temp;
                for (util::int_t i = 0; i < PyList_Size(_src); i++) {
                    _handleVar(_temp, PyList_GetItem(_src, i), _path + "[" + std::to_string(i) + "]");
                    _var.push_back(_temp);
                }
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(std::vector<util::string_t>& _var, PyObject *_src, util::string_t _path) {
                if (!(PyList_Check(_src)))
                    UERR(_path + " is not the correct type, must be array of strings");

                _list(_var, _src, _path);
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(std::vector<util::double_t>& _var, PyObject *_src, util::string_t _path) {
                if (!(PyList_Check(_src)))
                    UERR(_path + " is not the correct type, must be array of floats");
                
                _list(_var, _src, _path);
            }

            /* ---------------------------------------------------------------------- */

            void _handleVar(std::vector<util::int_t>& _var, PyObject *_src, util::string_t _path) {
                if (!(PyList_Check(_src)))
                    UERR(_path + " is not the correct type, must be array of ints");
                
                _list(_var, _src, _path);
            }

            // /* ---------------------------------------------------------------------- */

            // void _handleVar(std::vector<long>& _var, PyObject *_src, util::string_t _path) {
            //     if (!(PyList_Check(_src)))
            //         UERR(_path + " is not the correct type, must be array of ints");
            //     _list(_var, _src, _path);
            // }

            /* ---------------------------------------------------------------------- */

            void _handleVar(std::vector<util::bool_t>& _var, PyObject *_src, util::string_t _path) {
                if (!(PyList_Check(_src)))
                    UERR(_path + " is not the correct type, must be array of ints");
                _list(_var, _src, _path);
            }            
    };
}

#endif