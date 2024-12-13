#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include <limits>
#include "limits.h"

#include "../spatype.h"


/* ---------------------------------------------------------------------- */

#define STARTCOLOR "\033[0;32m"
#define ENDCOLOR "\033[0m"

/* ---------------------------------------------------------------------- */

#define UERR(_msg) util::error("\n  FILE: " + std::string(__FILE__) + "\n  FUNCTION: " + std::string(__PRETTY_FUNCTION__) + "\n  LINE: " + std::to_string(__LINE__) + "\n  MESSAGE: " + _msg)
#define ULOG(_msg) util::printColor("[" + util::formatFunc(__PRETTY_FUNCTION__) + "]", std::string(" ") + _msg)

/* ---------------------------------------------------------------------- */

#define ITERATE_OVER_NONEMPTY_LINES_FROM_FILE(file_name, func) std::string __line; SPARTA_NS::index_t __line_number = 1; std::ifstream __buf(file_name); if (__buf.is_open()) { while (std::getline(__buf,__line)) { util::trim(__line); if (__line.length() == 0) continue; func(__line, __line_number); __line_number++; } } else UERR("Failed to open: " + std::string(file_name)); __buf.close();

#define ITERATE_OVER_NONEMPTY_LINES_FROM_FILE_AND_SPLIT(file_name, func) std::vector<std::string> __split; std::string __line; SPARTA_NS::index_t __line_number = 1; std::ifstream __buf(file_name); if (__buf.is_open()) { while (std::getline(__buf,__line)) { util::trim(__line); if (__line.length() == 0) continue; util::splitStringAtWhiteSpace(__line, __split); func(__split, __line_number); __line_number++; } } else UERR("Failed to open: " + std::string(file_name)); __buf.close();

/* ---------------------------------------------------------------------- */

namespace util {
    // parameters set from sparta
    inline FILE* _screen;
    inline FILE* _logfile;
    inline SPARTA_NS::index_t _me;
    const  int npos = -1;

    // variables used a lot, so defining them here
    const int NO_INT = INT_MIN;
    const SPARTA_NS::index_t NO_INDEX_T = INDEX_T_MAX;

    /*
     * error command
     * NOTE: all functions that using this command must be wrapped in a 
     * try-catch statement that catches strings
    */
    void error(std::string _msg);

    /*
     * misc
    */
    // long max(std::vector<long> _v);

    int max(std::vector<int> _v);

    /*
     * Double to string stuff
    */
    std::string dtos(double _val);

    /*
     * Converting double and double arrays to bytes and string and bytes
     * for hashing purposes
    */
    // void doubleToByteArray(double val, char*& bytes);
    // std::string byteArrayToString(char*& bytes, long size);
    // std::string hashDouble(double val);
    std::string hashDoubleArray(double* arr, SPARTA_NS::index_t size);

    /**
     * gets current time as a string
    */
    // std::string getTime();
    
    /*
     * string and char* handling functions
    */
    void splitStringAtWhiteSpace(const std::string& input, std::vector<std::string>& out);
    std::string joinBySpaces(const std::vector<std::string>& words);
    void trim(std::string& str);
    int splitStringInto(char*& s, char**& arr);
    int findStringInArr(char* s, char** arr, int length);
    int countNumChar(char* s, char c);
    int find(char* s, char c, int start = 0);

    /**
     * custom printing
    */
    void print(std::string str, int num_indent=1, std::string end = "\n");
    void printColor(std::string color_string, std::string str, std::string end = "\n");
    std::string formatFunc(std::string _func);
    void printToFile(std::string str, int num_indent=0, std::string end = "\n");
    void printToScreen(std::string str, int num_indent=1, std::string end = "\n");
    

    /*
     * file stuff
    */
    bool fileExists(std::string filename);
    void copyFile(std::string from, std::string to);    

    /* 
     * class for an out file stream
    */
    class oFile {
        private:
            std::ofstream _out;
            std::string _f_name;
        public:
            oFile(std::string _file_name);
            ~oFile() { if (_out.is_open()) this->close(); }
            void close();
            friend oFile& operator<<(oFile& _f, const double& _val);
            friend oFile& operator<<(oFile& _f, const long& _val);
            friend oFile& operator<<(oFile& _f, const SPARTA_NS::index_t& _val);
            friend oFile& operator<<(oFile& _f, const int& _val);
            friend oFile& operator<<(oFile& _f, const std::size_t& _val);
            friend oFile& operator<<(oFile& _f, const std::string& _val);
            friend oFile& operator<<(oFile& _f, const char*& _val);
            friend oFile& operator<<(oFile& _f, const char& _val);
    };

    /*
     * class for an in file stream
    */
    class iFile {
        private:
            std::ifstream _out;
            std::string _f_name;
        public:
            iFile(std::string _file_name);
            ~iFile() { if (_out.is_open()) this->close(); }
            void close();
            bool getLine(std::string& _line);
            bool getLines(std::vector<std::string>& _lines, bool clear = true, std::string end = "\n");
    };
}

#endif


// // left trims a string
// void ltrim(std::string& _s);

// // right trims a string
// void rtrim(std::string& _s);

// // trims a string from right and left
// void trim(std::string& _s);


// double max(std::vector<double> _v);
// double max(double* _v, long _size);
// long   max(long* _v, long _size);

/**
 * writes the inputted string to file at the path inputted
*/ 
// void writeFile(std::string filename, std::string& lines);

/**
 * counts the number of lines in the file at the path inputted
*/ 
// long countLinesInFile(std::string _file);


/**
 * reads the contents of the file at the path inputted into the inputted
 * vector
*/ 
// void readFile(std::string fileName, std::vector<std::string>& lines);

// void split(std::string& _s, std::vector<std::string>& _v, char sep);


// long vecToArr(std::vector<std::string>& _vec, char**& _arr);

// void trim(std::string& _s);

// template<typename T>
// long find(std::vector<T> _v, T _find);