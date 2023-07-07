#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <array>
#include <vector>

namespace util {
    // error command
    // NOTE: all functions that using this command must be wrapped in a try-catch statement that catches strings
    void error(std::string _msg) { throw _msg; }

    // used in the EXEC function below, represents data returned from a system command
    struct CommandResult { std::string output; int exitstatus; };

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
            error("unhandled exception occured");
        }
        exitcode = pclose(pipe);
        return CommandResult{result, exitcode};
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
        long unsigned int j;
        for (j = 0; j < _s.length(); j++) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(j, _s.length() - j);
    }

    // right trims a string
    void rtrim(std::string& _s) {
        long unsigned int j;
        for (j = _s.length()-1; j >= 0; j--) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(0,  j+1);
    }

    // trims a string from right and left
    void trim(std::string& _s) { rtrim(_s); ltrim(_s); }


    void split(std::string& _s, std::vector<std::string>& _v, char sep) {
        _v.clear();
        long unsigned int last, j;
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
        for (long unsigned int i = 0; i < _v.size(); i++) {
            if (_v[i] == _find) return (int)i;
        }
        return -1;
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