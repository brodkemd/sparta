#ifndef ELMER_H
#define ELMER_H

#include <vector>
#include <optional>
#include <limits>
#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

#include "elmer_definitions.h"


namespace elmer {
    void error(std::string _msg) { throw _msg; }

    // used in the EXEC function below, represents data returned from a system command
    struct CommandResult {
        std::string output;
        int exitstatus;
    };


    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return commandResult containing STDOUT (not stderr) output & exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    static CommandResult EXEC(const std::string &command) {
        int exitcode = 0;
        std::array<char, 1048576> buffer {};
        std::string result;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result += std::string(buffer.data(), bytesread);
            }
        } catch (...) {
            pclose(pipe);
            throw;
        }
        exitcode = WEXITSTATUS(pclose(pipe));
        return CommandResult{result, exitcode};
    }


    void readFile(std::string fileName, std::vector<std::string>& lines) {
        std::ifstream myfile(fileName);
        std::string line;
        while (std::getline(myfile, line)) {
            if (line.length() >= indicator.length()) {
                if (line.substr(0, indicator.length()) == indicator){
                    lines.clear();
                    continue;
                }
            }
            lines.push_back(line + "\n");
        }
        myfile.close();
    }


    void readFileNoCheck(std::string fileName, std::vector<std::string>& lines) {
        std::ifstream myfile(fileName);
        std::string line;
        while (std::getline(myfile, line)) {
            lines.push_back(line + "\n");
        }
        myfile.close();
    }


    void getLatestNodeData(std::string filename, std::vector<double>& data) {
        std::vector<std::string> lines;

        data.clear();

        readFile(filename, lines);

        for (std::string it : lines) {
            boost::algorithm::trim(it);
            if (it.find(" ") == std::string::npos) {
                data.push_back(std::stod(it));
            }
        }
    }


    void getBoundaryData(std::string filename, std::vector<std::array<int, boundary_size>>& data) {
        std::vector<std::string> lines, split;
        std::array<int, boundary_size> arr;
        data.clear();

        readFileNoCheck(filename, lines);

        int count = 1;
        for (std::string it : lines) {
            boost::algorithm::trim(it);
            if (it.length() == 0) continue;

            // splitting the line at spaces
            boost::split(split, it, boost::is_any_of(" "));

            if (split[4] != (std::string)"303") 
                error("element is not a triangle in boundary file at line: " + std::to_string(count));

            // getting rid of stuff I do not need
            split.erase(split.begin()+1, split.begin() + 5);

            // catching errors
            if (split.size() != boundary_size)
                error("too many boundary elements to assign at line: " + std::to_string(count));

            // adding the data
            for (int i = 0; i < boundary_size; i++)
                arr[i] = std::stoi(split[i]);

            // adding the data to the class vector
            data.push_back(arr);
            count++;
        }
    }


    template<typename T, typename S>
    class dict_t {
        public:
            dict_t() {};

            typedef typename std::pair<T, std::optional<S>> item_type;
            typedef typename std::vector<std::pair<T, std::optional<S>>> vec_type;
            typedef typename vec_type::iterator iterator;
            typedef typename vec_type::const_iterator const_iterator;

            std::optional<S>& operator[](T _key) {
                size_t _ind = this->_find(_key);
                if (_ind == std::string::npos) {
                    this->_src.push_back(std::make_pair(_key, std::optional<S>{}));
                    return this->_src.back().second;
                } 
                return this->_src[_ind].second;
            }

            std::size_t length() { return this->_src.size(); }


            inline iterator begin() noexcept { return _src.begin(); }
            inline const_iterator cbegin() const noexcept { return _src.cbegin(); }
            inline iterator end() noexcept { return _src.end(); }
            inline const_iterator cend() const noexcept { return _src.cend(); }

        private:
            size_t _find(T _key) {
                for (size_t i = 0; i < this->_src.size(); i++) {
                    if (_src[i].first == _key) {
                        return i;
                    }
                }
                return std::string::npos;
            }

            vec_type _src;



    };


    class Base : public dict_t<std::string, std::string> {
        private:
            std::string _tab = "  ";
            std::string _end = "End";

        protected:
            std::string name, sep;

        public:
            int id = INT_MIN;
            // dict_t<std::string, std::string> contents;

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
            const std::string sep = "\n\n";

        public:
            Elmer();

            void join(std::string& _buffer);
            void run();

            std::string name, exe, sif, meshDBstem;

            Header     header;
            Simulation simulation;
            Constants  constants;
            std::vector<Solver>             solvers;
            std::vector<Equation>           equations;
            std::vector<Material>           materials;
            std::vector<Body>               bodys;
            std::vector<Initial_Condition>  initial_conditions;
            std::vector<Boundary_Condition> boundary_conditions;
    };
}


#endif