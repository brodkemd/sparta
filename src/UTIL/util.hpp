#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <array>
#include <vector>
#include <iomanip>
#include <string>
#include "float.h"

/* ---------------------------------------------------------------------- */

#define STARTCOLOR "\033[0;32m"
#define ENDCOLOR "\033[0m"

/* ---------------------------------------------------------------------- */

#define UERR(_msg) util::error("\n  FILE: " + std::string(__FILE__) + "\n  FUNCTION: " + std::string(__PRETTY_FUNCTION__) + "\n  LINE: " + std::to_string(__LINE__) + "\n  MESSAGE: " + _msg)
#define ULOG(_msg) util::printColor("[" + util::getTime() + "] [" + util::formatFunc(__PRETTY_FUNCTION__) + "]", std::string(" ") + _msg)

/* ---------------------------------------------------------------------- */

namespace util {
    typedef long int_t;
    typedef double double_t;
    typedef std::string string_t;
    typedef bool bool_t;
    
    /* ---------------------------------------------------------------------- */

    std::ostringstream makeDoubleConverter() {
        std::ostringstream __double_converter;
        __double_converter << std::scientific << std::setprecision(std::numeric_limits<double_t>::digits10+2);
        return __double_converter;
    }

    /* ---------------------------------------------------------------------- */

    FILE* _screen;
    FILE* _logfile;
    int_t _me;
    int_t npos = -1;
    std::ostringstream _double_converter = makeDoubleConverter();

    /* ---------------------------------------------------------------------- */

    // error command
    // NOTE: all functions that using this command must be wrapped in a try-catch statement that catches strings
    void error(string_t _msg) { throw _msg; }

    /* ---------------------------------------------------------------------- */

    /**
     * gets current time as a string
    */
    string_t getTime() {
        timeval curTime;
        gettimeofday(&curTime, NULL);
        int milli = curTime.tv_usec / 1000;

        char buffer [80];
        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&curTime.tv_sec));

        char currentTime[84] = "";
        sprintf(currentTime, "%s:%03d", buffer, milli);
        return string_t(currentTime);
    }

    /* ---------------------------------------------------------------------- */

    /**
     * custom printing for this class
    */
    void print(string_t str, int_t num_indent=1, string_t end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            string_t space = "";
            for (int_t i = 0; i < num_indent; i++) space += "  ";
            if (_screen)  fprintf(_screen,  "%s", (space+str+end).c_str());
            if (_logfile) fprintf(_logfile, "%s", (space+str+end).c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    /**
     * custom printing for this class
    */
    void printColor(string_t color_string, string_t str, string_t end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            if (_screen)  fprintf(_screen,  "%s%s%s%s", STARTCOLOR, color_string.c_str(), ENDCOLOR, (str+end).c_str());
            if (_logfile) fprintf(_logfile, "%s", (color_string+str+end).c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    string_t formatFunc(string_t _func) {
        if (_func.find("::") == string_t::npos) return string_t("");
        std::size_t _start = _func.find(' ');
        std::size_t _end   = _func.find('(');
        if (_end <= _start || _start == string_t::npos)
            _func = _func.substr(0, _end);
        else {
            if (_func.substr(0, _start) == "virtual")
                _func = _func.substr(_start+1, _func.length()-_start-1);
            _start = _func.find(' ');
            _end   = _func.find('(');
            _func = _func.substr(_start+1, _end-_start-1);
        }
        return _func;
    }

    /* ---------------------------------------------------------------------- */

    void printToFile(string_t str, int_t num_indent=0, string_t end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            string_t space = "";
            for (int_t i = 0; i < num_indent; i++) space += "  ";
            if (_logfile) fprintf(_logfile, "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    void printToScreen(string_t str, int_t num_indent=1, string_t end = "\n") {
        // only prints on the main process
        if (_me == 0) {
            string_t space = "";
            for (int_t i = 0; i < num_indent; i++) space += "  ";
            if (_screen) fprintf(_screen,  "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    // // used in the EXEC function below, represents data returned from a system command
    // struct CommandResult { string_t output; int_t exitstatus; };

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return commandResult containing STDOUT (not stderr) output & exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // static CommandResult EXEC(const string_t &command) {
    //     int_t exitcode = 0;
    //     std::array<char, 1048576> buffer {};
    //     string_t result;

    //     FILE *pipe = popen(command.c_str(), "r");
    //     if (pipe == nullptr) {
    //         throw std::runtime_error("popen() failed!");
    //     }
    //     try {
    //         std::size_t bytesread;
    //         while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
    //             result += string_t(buffer.data(), bytesread);
    //         }
    //     } catch (...) {
    //         pclose(pipe);
    //         error("unhandled exception occured");
    //     }
    //     exitcode = pclose(pipe);
    //     return CommandResult{result, exitcode};
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return  exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    int_t verboseExec(const string_t &command) {
        std::array<char, 128> buffer {};

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                print(string_t(buffer.data(), bytesread), 0, "");
            }
        } catch (...) {
            pclose(pipe);
            error("unhandled exception occurred");
        }
        return pclose(pipe);
    }

    /* ---------------------------------------------------------------------- */

    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return  exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    int_t limitedExec(const string_t &command, string_t _start, string_t _end) {
        std::array<char, 128> buffer {};
        string_t result, temp;
        result = temp = "";
        std::size_t _start_ind, _end_ind;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            std::size_t bytesread;
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                temp    = string_t(buffer.data(), bytesread);
                //data   += temp;
                result += temp;
                _start_ind = result.find(_start);
                _end_ind = result.find(_end, _start_ind+_start.length()+1);
                if (_start_ind != string_t::npos && _end_ind != string_t::npos) {
                    // printing only the part wanted
                    ULOG(result.substr(_start_ind+_start.length(), _end_ind - _start_ind - _start.length()));
                    result.clear();
                }
                printToFile(temp, 0, "");
            }
        } catch (...) { pclose(pipe); error("unhandled exception occured"); }
        return pclose(pipe);
    }

    /* ---------------------------------------------------------------------- */

    /**
     * Execute system command and get STDOUT result.
     * Regular system() only gives back exit status, this gives back output as well.
     * @param command system command to execute
     * @return  exitstatus
     * of command. Empty if command failed (or has no output). If you want stderr,
     * use shell redirection (2&>1).
     */
    int_t quietExec(const string_t &command) {
        std::array<char, 128> buffer {};
        string_t result = "";
        std::size_t bytesread;

        FILE *pipe = popen(command.c_str(), "r");
        if (pipe == nullptr) {
            throw std::runtime_error("popen() failed!");
        }
        try {
            while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                result += string_t(buffer.data(), bytesread);
            }
        } catch (...) {
            pclose(pipe);
            error("unhandled exception occured");
        }
        printToFile(result, 0, "");
        return pclose(pipe);
    }

    /* ---------------------------------------------------------------------- */

    string_t dtos(double_t _val) {
        _double_converter.str("");
        _double_converter << _val;
        return _double_converter.str();
    }

    /* ---------------------------------------------------------------------- */

    int_t vecToArr(std::vector<string_t>& _vec, char**& _arr) {
        const int_t _size = _vec.size();
        _arr = new char*[_size];
        for (int_t i = 0; i < _size; i++) _arr[i] = (char*)_vec[i].c_str();
        return _size;
    }

    /* ---------------------------------------------------------------------- */

    class oFile {
        private:
            std::ofstream _out;
            string_t _f_name;
        public:
            oFile(string_t _file_name) {
                _f_name = _file_name;
                ULOG("opening: " + _f_name);
                _out.open(_file_name);
                if (!(_out.is_open())) {
                    UERR("could not open: " + _f_name);
                }
            }

            ~oFile() {
                if (_out.is_open()) this->close();
            }

            void close() {
                ULOG("closing: " + _f_name);
                if (_out.is_open()) {
                    _out.close();
                    if (_out.is_open())
                        UERR("could not close: " + _f_name);
                }
                
            }

            friend oFile& operator<<(oFile& _f, const double_t& _val) { 
                _f._out << util::dtos(_val);
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const int_t& _val) { 
                _f._out << _val;
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const int& _val) { 
                _f._out << _val;
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const std::size_t& _val) { 
                _f._out << _val;
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const string_t& _val) { 
                _f._out << _val;
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const char*& _val) { 
                _f._out << _val;
                return _f;
            }

            friend oFile& operator<<(oFile& _f, const char& _val) { 
                _f._out << _val;
                return _f;
            }
    };

    /* ---------------------------------------------------------------------- */

    class iFile {
        private:
            std::ifstream _out;
            string_t _f_name;
        public:
            iFile(string_t _file_name) {
                _f_name = _file_name;
                ULOG("opening: " + _f_name);
                _out.open(_file_name);
                if (!(_out.is_open())) {
                    UERR("could not open: " + _f_name);
                }
            }

            ~iFile() {
                if (_out.is_open()) this->close();
            }

            void close() {
                ULOG("closing: " + _f_name);
                if (_out.is_open()) {
                    _out.close();
                    if (_out.is_open())
                        UERR("could not close: " + _f_name);
                }
            }

            bool_t getLine(string_t& _line) {
                if (std::getline(_out, _line))
                    return true;
                return false;
            }

            bool_t getLines(std::vector<string_t>& _lines, bool_t clear = true, string_t end = "\n") {
                string_t _line;
                if (clear) _lines.clear();
                while (std::getline(_out, _line)) {
                    _lines.push_back(_line + end);
                }
                return false;
            }
    };

    /* ---------------------------------------------------------------------- */

    void copyFile(string_t from, string_t to) {
        std::ifstream ini_file{from};
        std::ofstream out_file{to};
        if (ini_file && out_file) {
            ULOG("Copying: " + from + " --> " + to);
            out_file << ini_file.rdbuf();
        }
        else error("can not copy:"  + from + " --> " + to);
        ini_file.close(); out_file.close();
    }

    /* ---------------------------------------------------------------------- */

    // left trims a string
    void ltrim(string_t& _s) {
        std::size_t j;
        for (j = 0; j < _s.length(); j++) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(j, _s.length() - j);
    }

    /* ---------------------------------------------------------------------- */

    // right trims a string
    void rtrim(string_t& _s) {
        std::size_t j;
        for (j = _s.length()-1; j >= 0; j--) { if (_s[j] != ' ' && _s[j] != '\n') break; }
        _s = _s.substr(0,  j+1);
    }

    /* ---------------------------------------------------------------------- */

    // trims a string from right and left
    void trim(string_t& _s) { rtrim(_s); ltrim(_s); }

    /* ---------------------------------------------------------------------- */

    void split(string_t& _s, std::vector<string_t>& _v, char sep) {
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

    /* ---------------------------------------------------------------------- */

    template<typename T>
    int_t find(std::vector<T> _v, T _find) {
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] == _find) return i;
        }
        return npos;
    }

    /* ---------------------------------------------------------------------- */

    int_t max(std::vector<int_t> _v) {
        int_t _max = INT_MIN;
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    /* ---------------------------------------------------------------------- */

    double_t max(std::vector<double_t> _v) {
        double_t _max = DBL_MIN;
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    /* ---------------------------------------------------------------------- */

    double_t max(double_t* _v, int_t _size) {
        double_t _max = DBL_MIN;
        for (int_t i = 0; i < _size; i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    /* ---------------------------------------------------------------------- */

    int_t max(int_t* _v, int_t _size) {
        int_t _max = INT_MIN;
        for (int_t i = 0; i < _size; i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }

    

    /* ---------------------------------------------------------------------- */

    // /**
    //  * erases the file at the path inputted
    // */ 
    // void eraseFile(string_t filename) {
    //     std::ofstream ofs;
    //     ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    //     if (!(ofs.is_open())) error(filename + " did not open");
    //     ofs.close();
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * writes the inputted string to file at the path inputted
    */ 
    void writeFile(string_t filename, string_t& lines) {
        util::oFile out{filename};
        out << lines;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * counts the number of lines in the file at the path inputted
    */ 
    int_t countLinesInFile(string_t _file) {
        int_t number_of_lines = 0;
        string_t line;
        util::iFile myfile(_file);

        while (myfile.getLine(line)) number_of_lines++;
        return number_of_lines;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * reads the contents of the file at the path inputted into the inputted
     * vector
    */ 
    void readFile(string_t fileName, std::vector<string_t>& lines) {
        util::iFile myfile(fileName);
        myfile.getLines(lines, true);
        myfile.close();
    }
}

#endif