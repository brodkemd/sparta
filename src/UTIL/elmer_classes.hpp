#ifndef ELMER_CLASSES_H
#define ELMER_CLASSES_H

#include <algorithm>

#include "toml.hpp"

namespace elmer {
    util::string_t makeVarName(util::string_t _var_name) {
        std::replace(_var_name.begin(), _var_name.end(), '_', ' ');
        return _var_name;
    }

    /**
     * Base class that all elmer sections inherit from
    */
    class Section {
        private:
            
            util::string_t _name, _sep;

            // tab in the file
            util::string_t _tab = "  ";

            // ends the section
            util::string_t _end = "End";
            util::int_t _id;
            util::bool_t _include_count;

        public:
            toml::OrderedDict_t contents;
            Section() {}
            Section(util::string_t name, util::int_t id = toml::noInt, util::string_t sep = " = ", util::bool_t include_count = true) {
                this->_name = name;
                this->_sep = sep;
                this->_id = id;
                this->contents = toml::OrderedDict_t();
                this->_include_count = include_count;
            }

            void join(util::oFile& _buf) {
                if (this->_id == toml::noInt)
                    _buf << this->_name << "\n";
                else
                    _buf << this->_name << " " << _id << "\n";
                
                toml::Item_t keys = this->contents.getKeys();
                toml::Item_t vals = this->contents.getValues();
                for (util::int_t i = 0; i < keys.length(); i++) {
                    _buf << (_tab + makeVarName(keys[i].toString()));
                    if (vals[i].getType() == toml::LIST) {
                        if (this->_include_count)
                            _buf << "(" << std::to_string(vals[i].length()) << ")";
                        _buf << _sep << vals[i].toString(" ", true, 1000);
                    } else {
                        _buf << _sep << vals[i].toString(" ", true);
                    }
                    _buf << "\n";
                }
                // this->contents.toFile(_buf, _tab, this->_sep);
                _buf << this->_end << "\n";
            }

            /**
             * Return id for this class instance.
             *
             * @return integer id for this class instance.
             */
            util::int_t getId() { return this->_id; }

            /**
             * Sum numbers in a vector.
             *
             * @param values Container whose values are summed.
             * @return sum of `values`, or 0.0 if `values` is empty.
             */
            toml::Item_t getItem(toml::Item_t key) { return this->contents.getItem(key); }

            bool hasKey(toml::Item_t key) {
                return this->contents.hasKey(key);
            }

            void setItem(toml::Item_t key, toml::Item_t val) {
                this->contents.setItem(key, val);
            }
            void clear() { this->contents.clear(); }
            util::int_t length() { return this->contents.length(); }
    };
}

#endif