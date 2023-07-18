#ifndef TOML_H
#define TOML_H

#include "python_config.h"
#include "util.hpp"


namespace toml {
    const util::double_t noDouble       = std::numeric_limits<util::double_t>::min();
    const util::int_t noInt             = std::numeric_limits<int>::min();
    const util::string_t noString       = std::to_string((char)std::char_traits<char>::eof());
    const util::bool_t noBool           = 2;

    /* ---------------------------------------------------------------------- */
    enum{INT,DOUBLE,STRING,BOOL,LIST_INT,LIST};

    class Item_t {
        private:
            std::vector<Item_t> _list;
            std::string _src;
            int type = -1;
        public:
            Item_t() {}
            Item_t(util::double_t _key) { *this = _key; }
            // Item_t(int _key)         { *this = _key; }
            Item_t(util::int_t _key)    { *this = _key; }
            Item_t(util::string_t _key) { *this = _key; }
            Item_t(util::bool_t _key)   { *this = _key; }
            Item_t(std::vector<util::double_t> _key) { *this = _key; }
            Item_t(std::vector<util::int_t> _key)    { *this = _key; }
            Item_t(std::vector<util::string_t> _key) { *this = _key; }
            Item_t(std::vector<util::bool_t> _key)   { *this = _key; }
            Item_t(char* _key)       { *this = _key; }
            
            // Item_t& operator+(const Item_t& other) {
            //     if (this->type != other.type)
            //         UERR("Can not add types that are not the same");
                
            //     switch (this->type) {
            //         case INT:
            //             this->_src = std::to_string(std::stoi(this->_src) + std::stoi(other._src));
            //             break;
                    
            //         case DOUBLE:
            //             this->_src = util::dtos(std::stod(this->_src) + std::stod(other._src));
            //             break;
                    
            //         case STRING:
            //             this->_src = this->_src + other._src;
            //             break;
                    
            //         case BOOL:
            //             UERR("Can not add boolean values");
            //             break;
                    
            //         case LIST:
            //             for (auto it : other._list) this->_list.push_back(it);
            //             break;
            //     }
            // } 

            void operator = (const util::double_t _val) { type = toml::DOUBLE; this->_list.clear(); _src = util::dtos(_val); }
            void operator = (const util::int_t _val)    { type = toml::INT;    this->_list.clear(); _src = std::to_string(_val); }
            // void operator = (const long _val)        { type = toml::LONG; _src = std::to_string(_val); }
            void operator = (const util::string_t _val) { type = toml::STRING; this->_list.clear(); _src = _val; }
            void operator = (const util::bool_t _val)   { type = toml::BOOL;   this->_list.clear(); if (_val) _src = "True"; else _src = "False"; }
            void operator = (const char* _val)          { type = toml::STRING; this->_list.clear(); _src = std::string(_val); }
            template<typename T>
            void operator = (const std::vector<T> _val) {
                type = toml::LIST;
                this->_src.clear();
                this->_list.clear();
                Item_t _temp = Item_t();
                for (const auto& it : _val) {
                    _temp = it;
                    _list.push_back(_temp);
                }
            }

            bool operator==(const Item_t& other) { 
                if (this->type != other.type) return false;
                if (this->type == toml::LIST) {
                    if (this->_list.size() != other._list.size())
                        return false;
                    for (std::size_t i = 0; i < this->_list.size(); i++) {
                        if (!(this->_list[i] == other._list[i])) return false;
                    }
                }

                return false; 
            }

            Item_t& operator[](util::int_t _index) {
                if (this->type != LIST) UERR("Can not index non list object");
                if (_index >= (util::int_t)this->_list.size()) UERR("Index out of bounds for indexing object");
                return this->_list[_index];
            }

            util::string_t toString(std::string _sep = " ") { 
                if (this->type == LIST) {
                    util::string_t to_return;
                    for (auto it : this->_list) to_return+=(it.toString() + _sep);
                    return to_return;
                }
                return this->_src;
            }

            util::double_t toDouble() {
                if (this->type != DOUBLE) UERR("Can not convert non double to double");
                return std::stod(this->_src);
            }

            util::int_t toInt() {
                if (this->type != INT) UERR("Can not convert non int to int");
                return std::stoi(this->_src);
            }

            util::bool_t toBool() {
                if (this->type != BOOL) UERR("Can not convert non bool to bool");
                return (util::bool_t)std::stoi(this->_src);
            }

            void clear() { this->_list.clear(); this->_src.clear(); }

            friend class OrderedDict_t;
    };

    class OrderedDict_t {
        private:
            std::vector<Item_t> keys, vals;
        public:
            OrderedDict_t() { keys.clear(); vals.clear(); }

            Item_t& operator[](Item_t _key) {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == _key) return vals[i];
                }
                keys.push_back(_key);
                Item_t temp = Item_t();
                vals.push_back(temp);
                return vals[vals.size()-1];
            }

            void toString(std::string& _buf, std::string _start_entry="", std::string _sep_entry=" : ", std::string _end_entry="\n") {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    _buf += (_start_entry + keys[i]._src + _sep_entry + vals[i]._src + _end_entry);
                }
            }

            void toFile(util::oFile& _buf, std::string _start_entry="", std::string _sep_entry=" : ", std::string _end_entry="\n") {
                for (std::size_t i = 0; i < this->keys.size(); i++) {
                    _buf << (_start_entry + keys[i]._src + _sep_entry + vals[i]._src + _end_entry);
                }
            }

            void clear() { this->keys.clear(); this->vals.clear(); }
    };

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


            void getDictAtPath(OrderedDict_t& _d, util::string_t _path) {
                this->pArgs = PyTuple_New(2);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyTuple_SetItem(pArgs, 1, PyBool_FromLong((long)(true)));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();

                if (!(PyDict_Check(out))) {
                    UERR("Did not get dictionary at path: " + _path);
                }
                
                PyObject *key, *value;
                Item_t _key, _value;
                Py_ssize_t pos = 0;
                while (PyDict_Next(out, &pos, &key, &value)) {
                    if (PyUnicode_Check(key)) _key = _PyUnicode_AsString(key);
                    else if (PyLong_Check(key)) _key = PyLong_AsLong(key);
                    else if (PyFloat_Check(key)) _key = PyFloat_AsDouble(key);
                    else UERR("got invalid key type in dictionary at path: " + _path);

                    if (PyUnicode_Check(value)) _value = _PyUnicode_AsString(value);
                    else if (PyLong_Check(value)) _value = PyLong_AsLong(value);
                    else if (PyFloat_Check(value)) _value = PyFloat_AsDouble(value);
                    else if (PyBool_Check(value)) {
                        if (PyObject_IsTrue(value)) _value = true;
                        else _value = false;
                    }
                    else if (PyList_Check(value)) {
                        if (PyObject_Length(value)){
                            PyList_As
                        } else {
                            _value = std::vector<util::int_t>({});
                        }
                    }
                    else UERR("got invalid type for the value of: " + _key.toString());
                    _d[_key] = _value;
                }
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