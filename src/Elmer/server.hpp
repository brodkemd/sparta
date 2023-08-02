#ifndef SERVER_HPP
#define SERVER_HPP

#include "elmer_classes.hpp"
#include "shell_server_config.h"

namespace elmer {
    class Server {
        private:
            toml::Item_t simulation_directory, start_file, end_file, done_file, run_file, server_code_output_file;
            util::string_t contents;
        
        public:
            Server() {}
            ~Server() { this->deleteServerFile(); }
            void set(toml::handler& _h);
            void makeServerFile();
            void start();
            void end();
            void runCommand(util::string_t cmd);
            void waitForDoneFile();
            void deleteServerFile();
    };

    void Server::set(toml::handler& _h) {
        ULOG("Setting Elmer Server");
        _h.getAtPath(this->start_file,              "server.start_file",    toml::STRING);
        _h.getAtPath(this->end_file,                "server.end_file",      toml::STRING);
        _h.getAtPath(this->run_file,                "server.run_file",      toml::STRING);
        _h.getAtPath(this->done_file,               "server.done_file",     toml::STRING);
        _h.getAtPath(this->server_code_output_file, "server.server_file",   toml::STRING);
    }

    void Server::makeServerFile() {
        util::oFile out(this->server_code_output_file.toString());
        out << "START_FILE=" << this->start_file.toString() << "\n";
        out << "END_FILE="   << this->end_file.toString()   << "\n";
        out << "RUN_FILE="   << this->run_file.toString()   << "\n";
        out << "DONE_FILE="  << this->done_file.toString()  << "\n";
        out << SHELL_SERVER_STRING;
    }

    void Server::start() {
        if (util::fileExists(this->done_file.toString()))
            remove(this->done_file.toString().c_str());
        if (util::fileExists(this->run_file.toString()))
            remove(this->run_file.toString().c_str());
        if (util::fileExists(this->end_file.toString()))
            remove(this->end_file.toString().c_str());
        if (util::fileExists(this->start_file.toString()))
            remove(this->start_file.toString().c_str());

        ULOG("Starting Server");
        util::oFile out(this->start_file.toString());
        out.close();
        ULOG("Waiting for server to confirm");
        ULOG("Run the following command in a seperate terminal: bash " + this->server_code_output_file.toString());
        while (util::fileExists(this->start_file.toString())) { sleep(1); }
        ULOG("Server Confirmed, continuing");
    }

    void Server::end() {
        ULOG("Shutting Down Server");
        util::oFile out(this->end_file.toString());
        out.close();

        ULOG("Waiting for server to confirm");
        while (util::fileExists(this->end_file.toString())) { sleep(1); }
        ULOG("Server Confirmed, continuing");
    }

    void Server::runCommand(util::string_t cmd) {
        ULOG("Sending command to command server");
        util::oFile out(this->run_file.toString());
        out << cmd << "\n";
    }

    void Server::deleteServerFile() {
        remove(this->server_code_output_file.toString().c_str());
    }

    void Server::waitForDoneFile() {
        ULOG("Waiting for Server Command Completion");
        util::string_t line;
        while (!(util::fileExists(done_file.toString()))) { sleep(1); }
        util::iFile read_this(this->done_file.toString());
        read_this.getLine(line);
        read_this.close();
        util::trim(line);
        remove(done_file.toString().c_str());
        if (std::stoi(line)) {
            this->deleteServerFile();
            UERR("Server returned non-zero exit code: " + line);
        }
        
        
        ULOG("Got server command completion, continuing");
    }

} // namespace elmer



#endif