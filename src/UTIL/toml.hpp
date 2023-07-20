#ifndef TOML_H
#define TOML_H

#include "python_config.h"
#include "util.hpp"

#define TERR(_msg) UERR(_msg + ", while handling key: " + toml::cur_path)
#define TLOG(_msg) ULOG(_msg + ", handling key: " + toml::cur_path)

namespace toml {
    const util::double_t noDouble       = std::numeric_limits<util::double_t>::min();
    const util::int_t noInt             = std::numeric_limits<util::int_t>::min();
    const util::string_t noString       = std::to_string((char)std::char_traits<char>::eof());
    const util::bool_t noBool           = 2;
    static std::string cur_path         = "";

    /* ---------------------------------------------------------------------- */
    enum{NONE,INT,DOUBLE,STRING,BOOL,LIST};

    class Item_t {
        private:
            int type = NONE;

            std::vector<Item_t> _list;
            util::double_t _d_src;
            util::int_t _i_src;
            util::string_t _s_src;
            util::bool_t _b_src = 0;
        public:
            Item_t() { this->type = NONE; }
            Item_t(util::double_t _key) { *this = _key; }
            // Item_t(int _key)         { *this = _key; }
            Item_t(util::int_t _key)    { *this = _key; }
            // Item_t(Item_t _key)         { *this = _key; }
            Item_t(util::string_t _key) { *this = _key; }
            Item_t(util::bool_t _key)   { *this = _key; }
            Item_t(char* _key)          { *this = _key; }
            Item_t(const char* _key)    { *this = _key; }
            Item_t(std::vector<util::double_t> _key) { *this = _key; }
            Item_t(std::vector<util::int_t> _key)    { *this = _key; }
            Item_t(std::vector<util::string_t> _key) { *this = _key; }
            Item_t(std::vector<util::bool_t> _key)   { *this = _key; }
            Item_t(std::vector<Item_t> _key)         { *this = _key; }

            void operator = (util::double_t _val) {
                //ULOG("In double equal");
                type = toml::DOUBLE;
                _d_src = _val;
            }

            void operator = (util::int_t _val) {
                //ULOG("In int equal");
                type = toml::INT;
                _i_src = _val;
            }

            // void operator = (const long _val)        { type = toml::LONG; _src = std::to_string(_val); }
            void operator = (util::string_t _val) {
                //ULOG("In string equal");
                type = toml::STRING;
                _s_src = _val;
            }

            void operator = (util::bool_t _val) {
                //ULOG("In bool equal");
                type = toml::BOOL;
                _b_src = _val;
            }

            void operator = (char* _val) {
                //ULOG("In char* equal");
                type = toml::STRING;
                _s_src = std::string(_val);
            }

            void operator = (const char* _val) {
                //ULOG("In const char* equal");
                type = toml::STRING;
                _s_src = std::string(_val);
            }

            template<typename T>
            void operator = (std::vector<T> _val) {
                //ULOG("In list equal");
                type = toml::LIST;
                this->_list.clear();
                Item_t _temp = Item_t();
                for (T it : _val) {
                    _temp = it;
                    _list.push_back(_temp);
                }
            }

            bool operator==(const Item_t& other) { 
                if (this->type != other.type) return false;
                switch (this->type) {
                    case NONE:
                        return true;
                    
                    case INT:
                        return this->_i_src == other._i_src;
                    
                    case DOUBLE:
                        return this->_d_src == other._d_src;
                    
                    case STRING:
                        return this->_s_src == other._s_src;
                    
                    case BOOL:
                        return this->_b_src == other._b_src;

                    case LIST:
                        if (this->_list.size() != other._list.size())
                            return false;
                        for (std::size_t i = 0; i < this->_list.size(); i++) {
                            if (!(this->_list[i] == other._list[i])) return false;
                        }
                        break;
                }
                return false; 
            }

            Item_t operator[](util::int_t _index) {
                Item_t to_return = Item_t(0.0);
                switch (this->type) {
                    case NONE:
                        UERR("can not index none object");
                    
                    case INT:
                        UERR("can not index int object");
                    
                    case DOUBLE:
                        UERR("can not index double object");
                    
                    case STRING:
                        if (_index >= (util::int_t)this->_s_src.length() || _index < 0)
                            TERR("Index out of bounds for indexing object");
                        return this->_list[_index];
                    
                    case BOOL:
                        UERR("can not index bool object");

                    case LIST:
                        if (_index >= (util::int_t)this->_list.size() || _index < 0)
                            TERR("Index out of bounds for indexing object");
                        return this->_list[_index];
                }

                return to_return;
                
            }

            util::string_t toString(std::string _sep = " ", bool _include_string_parenthesis = false, util::int_t wrap_every = toml::noInt) {
                switch (this->type) {
                    case NONE:
                        return std::string("");
                    
                    case INT:
                        return std::to_string(this->_i_src);
                    
                    case DOUBLE:
                        return util::dtos(this->_d_src);
                    
                    case STRING:
                        if (_include_string_parenthesis)
                            return "\"" + _s_src + "\"";
                        else
                            return _s_src;
                    
                    case BOOL:
                        if (_b_src) return "True";
                        else return "False";

                    case LIST:
                        util::string_t to_return = "";
                        if (wrap_every == toml::noInt) {
                            for (auto it : this->_list)
                                to_return+=(it.toString(_sep, _include_string_parenthesis) + _sep);
                        } else {
                            util::int_t count = 1;
                            for (auto it : this->_list) {
                                if (count >= wrap_every) {
                                    to_return += "\n";
                                    count = 0;
                                }
                                to_return+=(it.toString(_sep, _include_string_parenthesis) + _sep);
                                count++;
                            }
                        }
                        return to_return ;
                }
                return std::string("");
            }


            util::double_t toDouble() {
                switch (this->type) {
                    case NONE:
                        TERR("Can not convert None to double");
                    
                    case INT:
                        return (util::double_t)this->_i_src;
                    
                    case DOUBLE:
                        return this->_d_src;
                    
                    case STRING:
                        TERR("Can not convert string to double");
                    
                    case BOOL:
                        TERR("Can not convert Bool to double");

                    case LIST:
                        TERR("Can not convert List to double");
                }
                return 0.0;
            }


            util::int_t toInt() {
                switch (this->type) {
                    case NONE:
                        TERR("Can not convert None to int");
                    
                    case INT:
                        return this->_i_src;
                    
                    case DOUBLE:
                        TERR("Can not convert double to int");
                    
                    case STRING:
                        TERR("Can not convert string to int");
                    
                    case BOOL:
                        TERR("Can not convert bool to int");

                    case LIST:
                        TERR("Can not convert list to int");
                }
                return 0;
            }


            util::bool_t toBool() {
                switch (this->type) {
                    case NONE:
                        TERR("Can not convert None to bool");
                    
                    case INT:
                        TERR("Can not convert int to bool");
                    
                    case DOUBLE:
                        TERR("Can not convert double to bool");
                    
                    case STRING:
                        TERR("Can not convert string to bool");
                    
                    case BOOL:
                        return this->_b_src;

                    case LIST:
                        TERR("Can not convert list to bool");
                }
                return false;
            }


            void append(Item_t _to_add) {
                switch (this->type) {
                    case NONE:
                        TERR("can not append anything to None type");
                    
                    case INT:
                        TERR("can not append anything to int type");
                    
                    case DOUBLE:
                        TERR("can not append anything to double type");
                    
                    case STRING:
                        TERR("can not append anything to string type");
                    
                    case BOOL:
                        TERR("can not append anything to bool type");

                    case LIST:
                        this->_list.push_back(_to_add);
                }
            }

            void clear() {
                _list.clear();
                _d_src = 0.0;
                _i_src = 0;
                _s_src.clear();
                _b_src = false;
            }

            void setType(int _type) {
                if (_type != this->type)
                    this->clear();
                this->type = _type;
            }

            util::int_t length() {
                switch (this->type) {
                    case NONE:
                        TERR("None does not have length");
                    
                    case INT:
                        TERR("int does not have length");
                    
                    case DOUBLE:
                        TERR("double does not have length");
                    
                    case STRING:
                        return (util::int_t)this->_s_src.length();
                    
                    case BOOL:
                        TERR("bool does not have length");

                    case LIST:
                        return (util::int_t)this->_list.size();
                }
                return -1;
            }
            int getType() { return this->type; }

            friend class OrderedDict_t;
    };

    class OrderedDict_t {
        private:
            std::vector<Item_t> keys, vals;
        public:
            OrderedDict_t() { keys.clear(); vals.clear(); }

            // Item_t& operator[](Item_t _key) {
            //     for (std::size_t  i = 0; i < this->keys.size(); i++) {
            //         if (keys[i] == _key) return vals.at(i);
            //     }
            //     keys.push_back(_key);
            //     Item_t temp = Item_t();
            //     vals.push_back(temp);
            //     return vals.at(vals.size()-1);
            // }

            Item_t getItem(Item_t key) {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == key) return vals.at(i);
                }
                UERR(key.toString() + " is not a valid key");
                return vals.at(0);
            }

            void setItem(Item_t key, Item_t val) {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == key) {
                        vals[i] == val;
                        return;
                    }
                }
                keys.push_back(key);
                vals.push_back(val);
            }

            // void toString(std::string& _buf, std::string _start_entry="", std::string _sep_entry=" : ", std::string _end_entry="\n") {
            //     for (std::size_t  i = 0; i < this->keys.size(); i++) {
            //         _buf += (_start_entry + keys[i].toString() + _sep_entry + vals[i].toString() + _end_entry);
            //     }
            // }

            // void toFile(util::oFile& _buf, std::string _start_entry="", std::string _sep_entry=" : ", std::string _end_entry="\n") {
            //     for (std::size_t i = 0; i < this->keys.size(); i++) {
            //         _buf << (_start_entry + keys[i].toString() + _sep_entry + vals[i].toString() + _end_entry);
            //     }
            // }

            Item_t getKeys() {
                Item_t to_return = this->keys;
                return to_return;
            }

            Item_t getValues() {
                Item_t to_return = this->vals;
                return to_return;
            }

            void clear() { this->keys.clear(); this->vals.clear(); }
            util::int_t length() { return (util::int_t)this->vals.size(); }
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

                cur_path = _path;
                if (!(PyDict_Check(out)))
                    TERR("Did not get dictionary");

                PyObject *key, *value;// , *_temp;
                Item_t _key, _value;
                Py_ssize_t pos = 0;
                while (PyDict_Next(out, &pos, &key, &value)) {
                    // ULOG(_PyUnicode_AsString(key));
                    // _temp = PyObject_Type(key);
                    // _temp = PyObject_Str(key);
                    // char* _out = (char*)_PyUnicode_AsString(_temp);
                    // ULOG(_out);
                    // _temp = PyObject_Type(value);
                    // _temp = PyObject_Str(_temp);
                    // _out = (char*)_PyUnicode_AsString(_temp);
                    // ULOG(_out);
                    //ULOG("Key");
                    this->_setVar(key, _key);
                    //ULOG("Get type: " + std::to_string(_key.getType()));
                    cur_path = (_path + "." + _key.toString());
                    
                    //ULOG("Value");
                    this->_setVar(value, _value);
                    //ULOG("Get type: " + std::to_string(_value.getType()));
                    //ULOG("setting");
                    _d.setItem(_key, _value);
                }
                // Py_DecRef(key); Py_DecRef(value);
            }

            void getDictKeysAtPath(std::vector<Item_t>& _d, util::string_t _path) {
                _d.clear();
                this->pArgs = PyTuple_New(2);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyTuple_SetItem(pArgs, 1, PyBool_FromLong((long)(true)));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();

                cur_path = _path;
                if (!(PyDict_Check(out)))
                    TERR("Did not get dictionary keys");

                PyObject *key, *value;// , *_temp;
                Item_t _key;
                Py_ssize_t pos = 0;
                while (PyDict_Next(out, &pos, &key, &value)) {
                    // ULOG(_PyUnicode_AsString(key));
                    // _temp = PyObject_Type(key);
                    // _temp = PyObject_Str(key);
                    // char* _out = (char*)_PyUnicode_AsString(_temp);
                    // ULOG(_out);
                    // _temp = PyObject_Type(value);
                    // _temp = PyObject_Str(_temp);
                    // _out = (char*)_PyUnicode_AsString(_temp);
                    // ULOG(_out);
                    //ULOG("Key");
                    this->_setVar(key, _key);
                    _d.push_back(_key);
                    //ULOG("Get type: " + std::to_string(_key.getType()));
                    // cur_path = (_path + "." + _key.toString());
                }
                // Py_DecRef(key); Py_DecRef(value);
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

            void _setVar(PyObject*& _in, Item_t& _out) {
                if (PyUnicode_Check(_in)) {
                    // ULOG("In string");
                    _out = _PyUnicode_AsString(_in);
                } else if (PyLong_Check(_in)) {
                    // ULOG("In int");
                    _out = PyLong_AsLong(_in);
                } else if (PyFloat_Check(_in)) {
                    // ULOG("In float");
                    _out = PyFloat_AsDouble(_in);
                } else if (PyBool_Check(_in)) {
                    // ULOG("In bool");
                    if (PyObject_IsTrue(_in)) _out = true;
                    else _out = false;
                } else if (PyList_Check(_in)) {
                    // ULOG("In list");
                    _out = std::vector<Item_t>({});
                    // _out.clear();
                    if (PyObject_Length(_in)) {
                        Item_t _temp;
                        PyObject* _py_temp;
                        for (Py_ssize_t i = 0; i < PyList_Size(_in); i++) {
                            _py_temp = PyList_GetItem(_in, i);
                            this->_setVar(_py_temp, _temp);
                            _out.append(_temp);
                        }
                    }
                } else { TERR("got invalid type"); }
                //ULOG("Return type:" + std::to_string(_out.getType()));
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