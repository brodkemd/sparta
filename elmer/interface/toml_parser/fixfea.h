#ifndef FIXFEA_H
#define FIXFEA_H

#include "classes.h"
#include "toml.h"

#include <sys/stat.h>


class FixFea {
    protected:
        // under "both"
        double emi;
        std::string tsurf_file;

        // under "sparta"
        int nevery;
        std::string groupID;
        std::string mixID;
        std::string customID;

        // under elmer
        Elmer elmer = Elmer();

        // Structure which would store the metadata
        struct stat sb;     

    public:
        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        template<typename dict>
        void run_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

        void print(std::string _msg) { std::cout << _msg << "\n"; }

        void get_elmer(std::string& _buffer) {
            this->elmer.join(_buffer);
        }

        // /**
        //  * executes a table
        // */
        // template<typename dict>
        // void FixFea::eval_table(std::string _caller, std::string _name, toml::table& _tbl, dict& _options);

        // /**
        //  * executes a table that wrapped as a node_t type
        // */
        // template<typename dict>
        // void FixFea::eval_table(std::string _caller, std::string _name, toml::node_t& __tbl, dict& _options);

    public:
        void handle_emi(std::string _caller, toml::node_t val);
        void handle_tsurf_file(std::string _caller, toml::node_t val);
        void handle_both(std::string _caller, toml::node_t tbl);
        void handle_nevery(std::string _caller, toml::node_t val);
        void handle_groupID(std::string _caller, toml::node_t val);
        void handle_mixID(std::string _caller, toml::node_t val);
        void handle_customID(std::string _caller, toml::node_t val);
        void handle_sparta(std::string _caller, toml::node_t tbl);
        void handle_exe(std::string _caller, toml::node_t val);
        void handle_sif(std::string _caller, toml::node_t val);
        void handle_meshDBstem(std::string _caller, toml::node_t val);
        void handle_header(std::string _caller, toml::node_t tbl);
        void handle_simulation(std::string _caller, toml::node_t tbl);
        void handle_constants(std::string _caller, toml::node_t tbl);
        void handle_solver(std::string _caller, toml::node_t tbl);
        void handle_equation(std::string _caller, toml::node_t tbl);
        void handle_material(std::string _caller, toml::node_t tbl);
        void handle_body(std::string _caller, toml::node_t tbl);
        void handle_initial_condition(std::string _caller, toml::node_t tbl);
        void handle_boundary_condition(std::string _caller, toml::node_t tbl);
        void handle_elmer(std::string _caller, toml::node_t tbl);
};

namespace toml {
    typedef std::pair<std::string, void (FixFea::*)(std::string, toml::node_t)> dict_item_t;
    typedef std::vector<dict_item_t> dict_t;
}

#endif