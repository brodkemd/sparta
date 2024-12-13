#ifndef TOML_H
#define TOML_H

#include "python_config.h"
#include "util.hpp"

#define TERR(_msg) UERR(_msg + ", while handling key: " + toml::cur_path)
#define TLOG(_msg) ULOG(_msg + ", handling key: " + toml::cur_path)

namespace elmer{
    class Section;
}

namespace toml {
    const double noDouble       = std::numeric_limits<double>::min();
    const long noInt             = std::numeric_limits<long>::min();
    const std::string noString       = std::to_string((char)std::char_traits<char>::eof());
    const bool noBool           = 2;
    static std::string cur_path         = "";

    /* ---------------------------------------------------------------------- */
    enum{NONE,INT,DOUBLE,STRING,BOOL,LIST};

    std::string getTypeName(long _type) {
        switch (_type) {
            case NONE:
                return "None";
            
            case INT:
                return "Int";
            
            case DOUBLE:
                return "Double";
            
            case STRING:
                return "String";
            
            case BOOL:
                return "Bool";

            case LIST:
                return "List";
        }
        return "Invalid";
    }

    class Item_t {
        private:
            int type = NONE;

            std::vector<Item_t> _list = {};
            double _d_src = 0.0;
            long _i_src = 0;
            std::string _s_src = "";
            bool _b_src = 0;

        public:
            Item_t()                    { this->type = NONE; }
            Item_t(double _key) { *this = _key; }
            Item_t(long _key)    { *this = _key; }
            Item_t(std::string _key) { *this = _key; }
            Item_t(bool _key)   { *this = _key; }
            Item_t(char* _key)          { *this = _key; }
            Item_t(const char* _key)    { *this = _key; }
            Item_t(std::vector<double> _key) { *this = _key; }
            Item_t(std::vector<long> _key)    { *this = _key; }
            Item_t(std::vector<std::string> _key) { *this = _key; }
            Item_t(std::vector<bool> _key)   { *this = _key; }
            Item_t(std::vector<Item_t> _key)         { *this = _key; }

            void operator = (double _val) { type = toml::DOUBLE; _d_src = _val; }
            void operator = (long _val)    { type = toml::INT; _i_src = _val; }
            void operator = (std::string _val) { type = toml::STRING; _s_src = _val; }
            void operator = (bool _val)   { type = toml::BOOL; _b_src = _val; }
            void operator = (char* _val)          { type = toml::STRING; _s_src = std::string(_val); }
            void operator = (const char* _val)    { type = toml::STRING; _s_src = std::string(_val); }

            template<typename T>
            void operator = (std::vector<T> _val) {
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
                        return true;
                }
                return false; 
            }


            Item_t operator[](long _index) {
                Item_t to_return = Item_t(0.0);
                switch (this->type) {
                    case NONE:
                        UERR("can not index none object");
                    
                    case INT:
                        UERR("can not index int object");
                    
                    case DOUBLE:
                        UERR("can not index double object");
                    
                    case STRING:
                        if (_index >= (long)this->_s_src.length() || _index < 0)
                            TERR("Index out of bounds for indexing object");
                        return this->_list[_index];
                    
                    case BOOL:
                        UERR("can not index bool object");

                    case LIST:
                        if (_index >= (long)this->_list.size() || _index < 0)
                            TERR("Index out of bounds for indexing object");
                        return this->_list[_index];
                }
                return to_return;
            }


            std::string toString(std::string _sep = " ", bool _include_string_parenthesis = false, long wrap_every = toml::noInt) {
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
                        std::string to_return = "";
                        if (wrap_every == toml::noInt) {
                            for (auto it : this->_list)
                                to_return+=(it.toString(_sep, _include_string_parenthesis) + _sep);
                        } else {
                            long count = 1;
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


            double toDouble() {
                switch (this->type) {
                    case NONE:
                        TERR("Can not convert None to double");
                    
                    case INT:
                        return (double)this->_i_src;
                    
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


            long toInt() {
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


            bool toBool() {
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


            std::vector<Item_t> toVector() {
                switch (this->type) {
                    case NONE:
                        TERR("Can not convert None to vector");
                    
                    case INT:
                        TERR("Can not convert None to vector");
                    
                    case DOUBLE:
                        TERR("Can not convert None to vector");
                    
                    case STRING:
                        TERR("Can not convert None to vector");
                    
                    case BOOL:
                        TERR("Can not convert None to vector");

                    case LIST:
                        return this->_list;
                }
                return this->_list;
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


            void clear() { _list.clear(); _d_src = 0.0;  _i_src = 0; _s_src.clear(); _b_src = false; }

            void setType(int _type) {
                if (_type != this->type)
                    this->clear();
                this->type = _type;
            }

            long length() {
                switch (this->type) {
                    case NONE:
                        TERR("None does not have length");
                    
                    case INT:
                        TERR("int does not have length");
                    
                    case DOUBLE:
                        TERR("double does not have length");
                    
                    case STRING:
                        return (long)this->_s_src.length();
                    
                    case BOOL:
                        TERR("bool does not have length");

                    case LIST:
                        return (long)this->_list.size();
                }
                return -1;
            }
            int getType() { return this->type; }

            friend class OrderedDict_t;
    };


    long listToCharArray(Item_t _list, char**& _arr) {
        std::vector<std::string> _temp;
        _temp.clear();
        for (auto it : _list.toVector())
            _temp.push_back(it.toString());
        return util::vecToArr(_temp, _arr);
    }

    template<typename T>
    class Dict_t {
        protected:
            std::vector<T> keys, vals;

            virtual std::string keyToString(T key) {
                Item_t _temp = key;
                return _temp.toString();
            }
        
        public:
            Dict_t() { this->keys.clear(); this->vals.clear(); }

            virtual T getKey(T key) {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == key) return vals.at(i);
                }
                UERR(this->keyToString(key) + " is not a valid key");
                return vals.at(0);
            }

            virtual void setKey(T key, T val) {
                for (std::size_t  i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == key) {
                        vals[i] = val;
                        return;
                    }
                }
                keys.push_back(key);
                vals.push_back(val);
            }

            virtual bool hasKey(T key) { return util::find(this->keys, key) != util::npos; }

            virtual void removeKey(T key) {
                std::size_t i;
                for (i = 0; i < this->keys.size(); i++) {
                    if (keys[i] == key) {
                        break;
                    }
                }
                this->keys.erase(this->keys.begin() + i);
                this->vals.erase(this->vals.begin() + i);
            }

            void clear() { this->keys.clear(); this->vals.clear(); }
            long length() { return (long)this->vals.size(); }
            
    };

    class OrderedDict_t : public Dict_t<Item_t> {
        public:
            OrderedDict_t() : Dict_t() {}

            Item_t getKeys() {
                Item_t to_return = this->keys;
                return to_return;
            }

            Item_t getValues() {
                Item_t to_return = this->vals;
                return to_return;
            }            
    };

    class handler {
        public:
            handler() {}
            handler(std::string _file) {
                Py_Initialize();

                // FILE* file_buf = fopen(_file.c_str(), "r");

                // if (file_buf != NULL) {
                //     PyRun_SimpleFile(file_buf, _file.c_str());
                //     fclose(file_buf);
                // } else {
                //     UERR("could not open fea config file");
                // }


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

            ~handler() { Py_Finalize(); Py_DECREF(pValue); }

            /* ---------------------------------------------------------------------- */

            void getAtPath(Item_t& _var, std::string _path, long _type) {
                cur_path = _path;
                this->pArgs = PyTuple_New(1);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();
            
                this->_setVar(out, _var);
                if (_var.getType() != _type)
                    TERR("Invalid type \"" + getTypeName(_var.getType()) + "\"");

                Py_DECREF(out);
            }

            template<typename T> // for elmer::section class
            void getDictAtPath(T& _s, std::string _path) { getDictAtPath(_s.contents, _path); }


            void getDictAtPath(OrderedDict_t& _d, std::string _path) {
                cur_path = _path;
                this->pArgs = PyTuple_New(1);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();

                if (!(PyDict_Check(out)))
                    TERR("Did not get dictionary");

                PyObject *key, *value;
                Item_t _key, _value;
                Py_ssize_t pos = 0;
                while (PyDict_Next(out, &pos, &key, &value)) {
                    this->_setVar(key, _key);
                    cur_path = (_path + "." + _key.toString());
                    
                    this->_setVar(value, _value);
                    _d.setKey(_key, _value);
                }
            }

            void getDictKeysAtPath(std::vector<Item_t>& _d, std::string _path) {
                cur_path = _path;
                _d.clear();
                this->pArgs = PyTuple_New(1);
                PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(_path.c_str()));
                PyObject* out = PyObject_CallObject(this->get_at, pArgs);
                this->_err();

                cur_path = _path;
                if (!(PyDict_Check(out)))
                    TERR("Did not get dictionary keys");

                PyObject *key, *value;
                Item_t _key;
                Py_ssize_t pos = 0;
                while (PyDict_Next(out, &pos, &key, &value)) {
                    this->_setVar(key, _key);
                    _d.push_back(_key);
                }
                // Py_DecRef(key); Py_DecRef(value);
            }


        private:
            PyObject *get_at, *pValue, *pArgs, *type, *value, *traceback;

            /* ---------------------------------------------------------------------- */

            void _err() {
                if (PyErr_Occurred()) {
                    PyErr_Fetch(&this->type, &this->value, &this->traceback);
                    std::string _value = _PyUnicode_AsString(PyObject_Repr(this->value));
                    TERR(_value);
                }
            }

            void _setVar(PyObject*& _in, Item_t& _out) {
                if (PyUnicode_Check(_in)) {
                    _out = _PyUnicode_AsString(_in);
                } else if (PyBool_Check(_in)) { // bool must be before int
                    if (PyObject_IsTrue(_in)) _out = true;
                    else _out = false;
                } else if (PyLong_Check(_in)) {
                    _out = PyLong_AsLong(_in);
                } else if (PyFloat_Check(_in)) {
                    _out = PyFloat_AsDouble(_in);
                } else if (PyList_Check(_in)) {
                    _out = std::vector<Item_t>({});

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
            }     
    };
}

#endif