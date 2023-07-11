#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <array>
#include <vector>
#include <limits>
#include <iomanip>
#include <string>
#include "float.h"

namespace util {
    std::ostringstream makeDoubleConverter() {
        std::ostringstream __double_converter;
        __double_converter << std::scientific<<std::setprecision(std::numeric_limits<double>::digits10+2);
        return __double_converter;
    }

    FILE* _screen;
    FILE* _logfile;
    int _me;

    static std::ostringstream _double_converter = makeDoubleConverter();

    // error command
    // NOTE: all functions that using this command must be wrapped in a try-catch statement that catches strings
    void error(std::string _msg) { throw _msg; }

    /**
     * custom printing for this class
    */
    void print(std::string str, int num_indent=1, std::string end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (int i = 0; i < num_indent; i++) space += "  ";
            if (_screen)  fprintf(_screen,  "%s%s%s", space.c_str(), str.c_str(), end.c_str());
            if (_logfile) fprintf(_logfile, "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    void printToFile(std::string str, int num_indent=0, std::string end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (int i = 0; i < num_indent; i++) space += "  ";
            if (_logfile) fprintf(_logfile, "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    void printToScreen(std::string str, int num_indent=1, std::string end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (int i = 0; i < num_indent; i++) space += "  ";
            if (_screen) fprintf(_screen,  "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    // // used in the EXEC function below, represents data returned from a system command
    // struct CommandResult { std::string output; int exitstatus; };

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return commandResult containing STDOUT (not stderr) output & exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // static CommandResult EXEC(const std::string &command) {
    //     int exitcode = 0;
    //     std::array<char, 1048576> buffer {};
    //     std::string result;

    //     FILE *pipe = popen(command.c_str(), "r");
    //     if (pipe == nullptr) {
    //         throw std::runtime_error("popen() failed!");
    //     }
    //     try {
    //         std::size_t bytesread;
    //         while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
    //             result += std::string(buffer.data(), bytesread);
    //         }
    //     } catch (...) {
    //         pclose(pipe);
    //         error("unhandled exception occured");
    //     }
    //     exitcode = pclose(pipe);
    //     return CommandResult{result, exitcode};
    // }


    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return  exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    int EXEC(const std::string &command) {
        std::array<char, 128> buffer {};
        std::string result;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result = std::string(buffer.data(), bytesread);
                print(result, 0, "");
            }
        } catch (...) {
            pclose(pipe);
            error("unhandled exception occured");
        }
        return pclose(pipe);
    }

    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return  exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    int EXEC(const std::string &command, std::string _start, std::string _end) {
        std::array<char, 128> buffer {};
        std::string result = "";
        std::size_t _start_ind, _end_ind;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result += std::string(buffer.data(), bytesread);
                _start_ind = result.find(_start);
                _end_ind = result.find(_end, _start_ind+_start.length()+1);
                if (_start_ind != std::string::npos && _end_ind != std::string::npos) {
                    util::print(result.substr(_start_ind+_start.length(), _end_ind - _start_ind - _start.length()), 0);
                    result.clear();
                }
            }
        } catch (...) {
            pclose(pipe);
            error("unhandled exception occured");
        }
        return pclose(pipe);
    }

    std::string dtos(double _val) {
        _double_converter.str("");
        _double_converter << _val;
        return _double_converter.str();
    }

    void copyFile(std::string from, std::string to) {
        std::ifstream ini_file{from};
        std::ofstream out_file{to};
        if (ini_file && out_file)
            out_file << ini_file.rdbuf();
        else
            error("can not copy "  + from + " to " + to);

        ini_file.close(); out_file.close();
    }

    // left trims a string
    void ltrim(std::string& _s) {
        std::size_t j;
        for (j = 0; j < _s.length(); j++) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(j, _s.length() - j);
    }

    // right trims a string
    void rtrim(std::string& _s) {
        std::size_t j;
        for (j = _s.length()-1; j >= 0; j--) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(0,  j+1);
    }

    // trims a string from right and left
    void trim(std::string& _s) { rtrim(_s); ltrim(_s); }


    void split(std::string& _s, std::vector<std::string>& _v, char sep) {
        _v.clear();
        std::size_t last, j;
        last = 0;
        for (j = 0; j <= _s.length(); j++) {
            if (j == _s.length()) {
                _v.push_back(_s.substr(last, j-last));
            } else if (_s[j] == sep) {
                if (j != last)
                    _v.push_back(_s.substr(last, j-last));
                last = j+1;
            }
        }
    }

    template<typename T>
    int find(std::vector<T> _v, T _find) {
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] == _find) return (int)i;
        }
        return -1;
    }

    std::size_t find(std::string _buf, std::string _to_find) {
        for (std::size_t i = 0; i < (_buf.length() - _to_find.length()); i++) {
            if (_buf.substr(i, _to_find.length()) == _to_find) return i;
        }
        return std::string::npos;
    }


    int max(std::vector<int> _v) {
        int _max = INT_MIN;
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }


    double max(std::vector<double> _v) {
        double _max = -DBL_MAX;
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }


    double max(double* _v, int _size) {
        double _max = -DBL_MAX;
        for (int i = 0; i < _size; i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    int max(int* _v, int _size) {
        int _max = INT_MIN;
        for (int i = 0; i < _size; i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    /**
     * erases the file at the path inputted
    */ 
    void eraseFile(std::string filename) {
        std::ofstream ofs;
        ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
        if (!(ofs.is_open())) error(filename + " did not open");
        ofs.close();
    }

    /**
     * writes the inputted string to file at the path inputted
    */ 
    void writeFile(std::string filename, std::string& lines) {
        std::ofstream out(filename);
        if (!(out.is_open())) error(filename + " did not open");
        out << lines;
        out.close();
    }

    /**
     * counts the number of lines in the file at the path inputted
    */ 
    int count_lines_in_file(std::string _file) {
        int number_of_lines = 0;
        std::string line;
        std::ifstream myfile(_file);

        while (std::getline(myfile, line)) ++number_of_lines;
        return number_of_lines;
    }

    /**
     * reads the contents of the file at the path inputted into the inputted
     * vector
    */ 
    void readFile(std::string fileName, std::vector<std::string>& lines) {
        std::ifstream myfile(fileName);
        std::string line;
        while (std::getline(myfile, line)) {
            lines.push_back(line + "\n");
        }
        myfile.close();
    }
}

#endif