#include "elmer.h"

namespace elmer{
    Section::Section(std::string name, long id, std::string sep, bool include_count) {
        this->_name = name; this->_sep = sep; this->_id = id;
    }

    /* ---------------------------------------------------------------------- */

    void Section::joinInto(util::oFile& _buf) {
        if (this->_id == util::NO_INT)
            _buf << this->_name << "\n";
        else
            _buf << this->_name << " " << _id << "\n";

        for (auto it : this->_content_pairs)
            _buf << ELMER_LINE_START << it[0] << this->_sep << it[1] << "\n";
        _buf << ELMER_SECTION_END << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Section::addEquality(std::string var, long val, bool include_count) {
        if (include_count)
            this->_content_pairs.push_back({var+"(1)", std::to_string(val)});
        else
            this->_content_pairs.push_back({var,       std::to_string(val)});
    }

    /* ---------------------------------------------------------------------- */

    void Section::addEquality(std::string var, double val, bool include_count) {
        if (include_count)
            this->_content_pairs.push_back({var+"(1)", util::dtos(val)});
        else
            this->_content_pairs.push_back({var,       util::dtos(val)});
    }

    /* ---------------------------------------------------------------------- */

    void Section::addEquality(std::string var, std::vector<long>& val, bool include_count) {
        std::array<std::string, 2> buf;
        if (include_count)
            buf[0] = var + "(" + std::to_string(val.size()) + ")";
        
        buf[1] = "";
        if (val.size()) {
            buf[1] = std::to_string(val[0]);
            for (std::size_t i = 1; i < val.size(); i++)
                buf[1] += " " + std::to_string(val[i]);
        } else
            UERR("can not interpret vector with no length");
        
        this->_content_pairs.push_back(buf);
    }

    /* ---------------------------------------------------------------------- */

    void Section::addEquality(std::string var, std::vector<double>& val, bool include_count) {
        std::array<std::string, 2> buf;
        if (include_count)
            buf[0] = var + "(" + std::to_string(val.size()) + ")";
        else 
            buf[0] = var;

        buf[1] = "";
        if (val.size()) {
            buf[1] = util::dtos(val[0]);
            for (std::size_t i = 1; i < val.size(); i++)
                buf[1] += " " + util::dtos(val[i]);
        } else
            UERR("can not interpret vector with no length");
        this->_content_pairs.push_back(buf);
    }

    /* ---------------------------------------------------------------------- */

    void Section::addEquality(std::string var, double*& arr, long len, bool include_count) {
        std::array<std::string, 2> buf;
        if (include_count)
            buf[0] = var + "(" + std::to_string(len) + ")";
        else 
            buf[0] = var;

        buf[1] = "";
        if (len) {
            buf[1] = util::dtos(arr[0]);
            for (long i = 1; i < len; i++)
                buf[1] += " " + util::dtos(arr[i]);
        } else
            UERR("can not interpret array with no length");
        this->_content_pairs.push_back(buf);
    }
}