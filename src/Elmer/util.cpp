#include "util.h"

// #include <iomanip>
#include <cstring>
#include <sys/time.h>

const char* char_fmt =  "%.16e"; // ("%." + std::to_string(std::numeric_limits<double>::digits10) + "e").c_str();
const std::regex whitespace_pattern("\\s+");
const std::sregex_token_iterator regex_end;

/* ---------------------------------------------------------------------- */

namespace util {
    const char b64chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    void error(std::string _msg) { throw _msg; }

    /* ---------------------------------------------------------------------- */

    void splitStringAtWhiteSpace(const std::string& input, std::vector<std::string>& out) {
        std::sregex_token_iterator iter(input.begin(), input.end(), whitespace_pattern, -1);
        out.clear();
        for (; iter != regex_end; ++iter) out.push_back(*iter);
    }

    /* ---------------------------------------------------------------------- */

    // void doubleToByteArray(double val, char*& bytes) {
    //     bytes = new char[sizeof(double)];
    //     memcpy(bytes, &val, sizeof(double));
    // }

    // /* ---------------------------------------------------------------------- */

    // std::string byteArrayToString(char*& bytes, long size) {
    //     return std::string(bytes, bytes+size);
    // }

    // /* ---------------------------------------------------------------------- */

    // std::string hashDouble(double val) {
    //     char* buf;
    //     doubleToByteArray(val, buf);
    //     std::string to_return = byteArrayToString(buf, sizeof(double));
    //     delete [] buf;
    //     return to_return;
    // }

    /* ---------------------------------------------------------------------- */

    // std::string hashDoubleArray(double* arr, long size) {
    //     std::string to_return = "";
    //     if (size) {
    //         char* buf;
    //         for (long i = 0; i < size; i++) {
    //             doubleToByteArray(arr[i], buf);
    //             to_return += byteArrayToString(buf, sizeof(double));
    //         }
    //         delete [] buf;
    //     }
    //     return to_return;
    // }

    SPARTA_NS::index_t b64EncodedSize(SPARTA_NS::index_t inlen) {
        SPARTA_NS::index_t ret = inlen;
        if (inlen % 3 != 0) ret += 3 - (inlen % 3);
        ret /= 3;
        ret *= 4;
        return ret;
    }

    void b64Encode(const unsigned char *in, SPARTA_NS::index_t len, std::string& out) {
        SPARTA_NS::index_t elen, i, j, v;

        if (in == NULL || len == 0) { out = ""; return; }
        
        elen = b64EncodedSize(len);
        out.clear();
        out.resize(elen);

        for (i=0, j=0; i<len; i+=3, j+=4) {
            v = in[i];
            v = i+1 < len ? v << 8 | in[i+1] : v << 8;
            v = i+2 < len ? v << 8 | in[i+2] : v << 8;
            out[j]   = b64chars[(v >> 18) & 0x3F];
            out[j+1] = b64chars[(v >> 12) & 0x3F];

            if (i+1 < len) out[j+2] = b64chars[(v >> 6) & 0x3F];
            else out[j+2] = '=';
            
            if (i+2 < len) out[j+3] = b64chars[v & 0x3F];
            else out[j+3] = '=';
        }
    }


    std::string hashDoubleArray(double* val, const SPARTA_NS::index_t len) {
        unsigned char* bytes = new unsigned char[len*sizeof(double)];
        memcpy(bytes, val, len*sizeof(double));
        std::string to_return;
        b64Encode(bytes, sizeof(double)*len, to_return);
        delete [] bytes;
        return to_return;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * gets current time as a string
    */
    // std::string getTime() {
    //     timeval curTime;
    //     gettimeofday(&curTime, NULL);
    //     long milli = curTime.tv_usec / 1000;

    //     char buffer [80];
    //     strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&curTime.tv_sec));

    //     char currentTime[84] = "";
    //     sprintf(currentTime, "%s:%03ld", buffer, milli);
    //     return std::string(currentTime);
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * custom printing for this class
    */
    void print(std::string str, long num_indent, std::string end) {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (long i = 0; i < num_indent; i++) space += "  ";
            if (_screen)  fprintf(_screen,  "%s", (space+str+end).c_str());
            if (_logfile) fprintf(_logfile, "%s", (space+str+end).c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    int countNumChar(char* s, char c) {
        int length = (int)strlen(s);
        int count = 0;
        for (int i = 0; i < length; i++) {
            if (s[i] == c) count++;
        }
        return count;
    }

    /* ---------------------------------------------------------------------- */

    int find(char* s, char c, int start) {
        int length = (int)strlen(s);
        for (int i = start; i < length; i++) {
            if (s[i] == c) return i;
        }
        return length;
    }

    /* ---------------------------------------------------------------------- */

    int splitStringInto(char*& s, char**& arr) {
        int count = countNumChar(s, ' ') + 1;
        int next, last = 0;
        arr = new char*[count];
        for (int i = 0; i < count; i++) {
            next = find(s, ' ', last);
            //print(std::to_string(next-last));
            arr[i] = new char[next - last];
            for (int j = 0; j < next - last; j++)
                arr[i][j] = s[j+last];
            
            arr[i][next-last] = '\0';
            last = next+1;
        }
        return count;
    }

    /* ---------------------------------------------------------------------- */

    int findStringInArr(char* s, char** arr, int length) {
        for (int i = 0; i < length; i++) {
            if (strcmp(s, arr[i]) == 0) return i;
        }
        return length;
    }

    /* ---------------------------------------------------------------------- */

    /**
     * custom printing for this class
    */
    void printColor(std::string color_string, std::string str, std::string end) {
        // only prints on the main process
        if (_me == 0) {
            if (_screen)  fprintf(_screen,  "%s%s%s%s", STARTCOLOR, color_string.c_str(), ENDCOLOR, (str+end).c_str());
            if (_logfile) fprintf(_logfile, "%s", (color_string+str+end).c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    std::string formatFunc(std::string _func) {
        if (_func.find("::") == std::string::npos) return std::string("");
        std::size_t _start = _func.find(' ');
        std::size_t _end   = _func.find('(');
        if (_end <= _start || _start == std::string::npos) {
            _func = _func.substr(0, _end);
        } else {
            if (_func.substr(0, _start) == "virtual")
                _func = _func.substr(_start+1, _func.length()-_start-1);
            _start = _func.find(' ');
            _end   = _func.find('(');
            _func = _func.substr(_start+1, _end-_start-1);
        }
        return _func;
    }

    /* ---------------------------------------------------------------------- */

    void printToFile(std::string str, long num_indent, std::string end) {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (long i = 0; i < num_indent; i++) space += "  ";
            if (_logfile) fprintf(_logfile, "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    void printToScreen(std::string str, long num_indent, std::string end) {
        // only prints on the main process
        if (_me == 0) {
            std::string space = "";
            for (long i = 0; i < num_indent; i++) space += "  ";
            if (_screen) fprintf(_screen,  "%s%s%s", space.c_str(), str.c_str(), end.c_str());
        }
    }

    /* ---------------------------------------------------------------------- */

    std::string dtos(double val) {
        char buf[100];
        sprintf(buf, char_fmt, val); 
        return std::string(buf);
    }

    /* ---------------------------------------------------------------------- */

    bool fileExists(std::string filename) {
        std::ifstream infile(filename);
        return infile.good();
    }

    /* ---------------------------------------------------------------------- */

    oFile::oFile(std::string _file_name) {
        _f_name = _file_name;
        // ULOG("opening: " + _f_name);
        _out.open(_file_name);
        if (!(_out.is_open()))
            UERR("could not open: " + _f_name);
    }

    void oFile::close() {
        // ULOG("closing: " + _f_name);
        if (_out.is_open()) {
            _out.close();
            if (_out.is_open())
                UERR("could not close: " + _f_name);
        }
    }

    oFile& operator<<(oFile& _f, const double& _val) { 
        _f._out << util::dtos(_val);
        return _f;
    }

    oFile& operator<<(oFile& _f, const long& _val) { 
        _f._out << _val;
        return _f;
    }

    oFile& operator<<(oFile& _f, const int& _val) { 
        _f._out << _val;
        return _f;
    }

    oFile& operator<<(oFile& _f, const std::size_t& _val) { 
        _f._out << _val;
        return _f;
    }

    oFile& operator<<(oFile& _f, const std::string& _val) { 
        _f._out << _val;
        return _f;
    }

    oFile& operator<<(oFile& _f, const char*& _val) { 
        _f._out << _val;
        return _f;
    }

    oFile& operator<<(oFile& _f, const char& _val) { 
        _f._out << _val;
        return _f;
    }

    /* ---------------------------------------------------------------------- */
    
    iFile::iFile(std::string _file_name) {
        _f_name = _file_name;
        // ULOG("opening: " + _f_name);
        _out.open(_file_name);
        if (!(_out.is_open()))
            UERR("could not open: " + _f_name);
    }


    void iFile::close() {
        // ULOG("closing: " + _f_name);
        if (_out.is_open()) {
            _out.close();
            if (_out.is_open())
                UERR("could not close: " + _f_name);
        }
    }

    bool iFile::getLine(std::string& _line) {
        if (std::getline(_out, _line)) return true;
        return false;
    }

    bool iFile::getLines(std::vector<std::string>& _lines, bool clear, std::string end) {
        std::string _line;
        if (clear) _lines.clear();
        while (std::getline(_out, _line))
            _lines.push_back(_line + end);
        return false;
    }

    /* ---------------------------------------------------------------------- */

    void copyFile(std::string from, std::string to) {
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

    void trim(std::string& str) {
        str.erase(str.begin(), std::find_if(str.begin(), str.end(), 
                    [](unsigned char c) { return !std::isspace(c); }));
        str.erase(std::find_if(str.rbegin(), str.rend(), 
                    [](unsigned char c) { return !std::isspace(c); }).base(), str.end());
    }

    /* ---------------------------------------------------------------------- */

    std::string joinBySpaces(const std::vector<std::string>& words) {    
        std::string result;
        for (size_t i = 0; i < words.size(); ++i) {
            result += words[i];
            if (i != words.size() - 1)  // Don't add a space after the last word
                result += " ";
        }
        return result;
    }

    /* ---------------------------------------------------------------------- */

    int max(std::vector<int> _v) {
        int _max = NO_INT;
        for (std::size_t i = 0; i < _v.size(); i++) {
            if (_v[i] > _max) _max = _v[i];
        }
        return _max;
    }
}


/* --------------- hashing functions ---------------------- */

    // std::string sha256_raw(const unsigned char* data, size_t len) {
    //     unsigned char hash[SHA256_DIGEST_LENGTH];
    //     SHA256_CTX sha256;
    //     SHA256_Init(&sha256);
    //     SHA256_Update(&sha256, data, len);
    //     SHA256_Final(hash, &sha256);

    //     // Return the raw bytes as a string 
    //     return std::string(reinterpret_cast<char*>(hash), SHA256_DIGEST_LENGTH); 
    // }

    // std::string longToHash(long& val) {
    //     unsigned char *bytes = reinterpret_cast<unsigned char*>(&val);
    //     return sha256_raw(bytes, sizeof(val));
    // }

    // std::string stringToHash(std::string& val) {
    //     unsigned char *bytes = (unsigned char*)val.c_str();
    //     return sha256_raw(bytes, val.length());
    // }

    // std::string doubleToHash(double& val) {
    //     unsigned char *bytes = reinterpret_cast<unsigned char*>(&val);
    //     return sha256_raw(bytes, sizeof(double));
    // }

    // std::string doubleArrayToHash(double* arr, long length) {
    //     std::string out = "";
    //     for (int i = 0; i < length; i++) out+=doubleToHash(arr[i]);
    //     return out;
    // }

        // /* ---------------------------------------------------------------------- */

    // double max(std::vector<double> _v) {
    //     double _max = DBL_MIN;
    //     for (std::size_t i = 0; i < _v.size(); i++) {
    //         if (_v[i] > _max) _max = _v[i];
    //     }
    //     return _max;
    // }

    // /* ---------------------------------------------------------------------- */

    // double max(double* _v, long _size) {
    //     double _max = DBL_MIN;
    //     for (long i = 0; i < _size; i++) {
    //         if (_v[i] > _max) _max = _v[i];
    //     }
    //     return _max;
    // }

    // /* ---------------------------------------------------------------------- */

    // long max(long* _v, long _size) {
    //     long _max = INT_MIN;
    //     for (long i = 0; i < _size; i++) {
    //         if (_v[i] > _max) _max = _v[i];
    //     }
    //     return _max;
    // }

    

    /* ---------------------------------------------------------------------- */

    // /**
    //  * erases the file at the path inputted
    // */ 
    // void eraseFile(std::string filename) {
    //     std::ofstream ofs;
    //     ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    //     if (!(ofs.is_open())) error(filename + " did not open");
    //     ofs.close();
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * writes the inputted string to file at the path inputted
    */ 
    // void writeFile(std::string filename, std::string& lines) {
    //     util::oFile out{filename};
    //     out << lines;
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * counts the number of lines in the file at the path inputted
    */ 
    // long countLinesInFile(std::string _file) {
    //     long number_of_lines = 0;
    //     std::string line;
    //     util::iFile myfile(_file);

    //     while (myfile.getLine(line)) number_of_lines++;
    //     return number_of_lines;
    // }

    /* ---------------------------------------------------------------------- */

    /**
     * reads the contents of the file at the path inputted into the inputted
     * vector
    */ 
    // void readFile(std::string fileName, std::vector<std::string>& lines) {
    //     util::iFile myfile(fileName);
    //     myfile.getLines(lines, true);
    //     myfile.close();
    // }

    // void split(std::string& _s, std::vector<std::string>& _v, char sep) {
    //     _v.clear();
    //     std::size_t last, j;
    //     last = 0;
    //     for (j = 0; j <= _s.length(); j++) {
    //         if (j == _s.length()) {
    //             _v.push_back(_s.substr(last, j-last));
    //         } else if (_s[j] == sep) {
    //             if (j != last)
    //                 _v.push_back(_s.substr(last, j-last));
    //             last = j+1;
    //         }
    //     }
    // }

    // /* ---------------------------------------------------------------------- */

    // template<typename T>
    // long find(std::vector<T> _v, T _find) {
    //     for (std::size_t i = 0; i < _v.size(); i++) {
    //         if (_v[i] == _find) return i;
    //     }
    //     return npos;
    // }

    // /* ---------------------------------------------------------------------- */

    // // left trims a string
    // void ltrim(std::string& _s) {
    //     std::size_t j;
    //     for (j = 0; j < _s.length(); j++) { if (_s[j] != ' ' && _s[j] != '\n') break; }
    //     _s = _s.substr(j, _s.length() - j);
    // }

    // /* ---------------------------------------------------------------------- */

    // // right trims a string
    // void rtrim(std::string& _s) {
    //     std::size_t j;
    //     for (j = _s.length()-1; j >= 0; j--) { if (_s[j] != ' ' && _s[j] != '\n') break; }
    //     _s = _s.substr(0,  j+1);
    // }

    // void trim(std::string& _s) { rtrim(_s); ltrim(_s); }
    
    /* ---------------------------------------------------------------------- */

    // long vecToArr(std::vector<std::string>& _vec, char**& _arr) {
    //     const long _size = _vec.size();
    //     char* pc;
    //     _arr = new char*[_size];
    //     for (long i = 0; i < _size; i++) {
    //         pc = new char[_vec[i].size() + 1];
    //         std::strcpy(pc, _vec[i].c_str());
    //         _arr[i] = pc;
    //     }
    //     return _size;
    // }

    /* ---------------------------------------------------------------------- */

    // // used in the EXEC function below, represents data returned from a system command
    // struct CommandResult { std::string output; long exitstatus; };

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return commandResult containing STDOUT (not stderr) output & exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // static CommandResult EXEC(const std::string &command) {
    //     long exitcode = 0;
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

    /* ---------------------------------------------------------------------- */

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return  exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // long verboseExec(const std::string &command) {
    //     std::array<char, 128> buffer {};

    //     FILE *pipe = popen(command.c_str(), "r");
    //     if (pipe == nullptr) {
    //         throw std::runtime_error("popen() failed!");
    //     }
    //     try {
    //         std::size_t bytesread;
    //         while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
    //             print(std::string(buffer.data(), bytesread), 0, "");
    //         }
    //     } catch (...) {
    //         pclose(pipe);
    //         error("unhandled exception occurred");
    //     }
    //     return pclose(pipe);
    // }

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return  exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // long limitedExec(const std::string &command, std::string _start, std::string _end) {
    //     std::array<char, 128> buffer {};
    //     std::string result, temp;
    //     result = temp = "";
    //     std::size_t _start_ind, _end_ind;
    //     // ULOG("right before popen");

    //     FILE *pipe = popen(command.c_str(), "r");
    //     // ULOG("right after popen");
    //     if (pipe == nullptr) {
    //         throw std::runtime_error("popen() failed!");
    //     }
    //     try {
    //         std::size_t bytesread;
    //         // ULOG("before read");
    //         while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
    //             // ULOG("in loop");
    //             temp    = std::string(buffer.data(), bytesread);
    //             //data   += temp;
    //             result += temp;
    //             _start_ind = result.find(_start);
    //             _end_ind = result.find(_end, _start_ind+_start.length()+1);
    //             if (_start_ind != std::string::npos && _end_ind != std::string::npos) {
    //                 // printing only the part wanted
    //                 ULOG(result.substr(_start_ind+_start.length(), _end_ind - _start_ind - _start.length()));
    //                 result.clear();
    //             }
    //             printToFile(temp, 0, "");
    //         }
    //     } catch (...) { pclose(pipe); error("unhandled exception occured"); }
    //     return pclose(pipe);
    // }

    // /* ---------------------------------------------------------------------- */

    // /**
    //  * Execute system command and get STDOUT result.
    //  * Regular system() only gives back exit status, this gives back output as well.
    //  * @param command system command to execute
    //  * @return  exitstatus
    //  * of command. Empty if command failed (or has no output). If you want stderr,
    //  * use shell redirection (2&>1).
    //  */
    // long quietExec(const std::string &command) {
    //     std::array<char, 128> buffer {};
    //     std::string result = "";
    //     std::size_t bytesread;

    //     FILE *pipe = popen(command.c_str(), "r");
    //     if (pipe == nullptr) {
    //         throw std::runtime_error("popen() failed!");
    //     }
    //     try {
    //         while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
    //             result += std::string(buffer.data(), bytesread);
    //         }
    //     } catch (...) {
    //         pclose(pipe);
    //         error("unhandled exception occured");
    //     }
    //     printToFile(result, 0, "");
    //     return pclose(pipe);
    // }