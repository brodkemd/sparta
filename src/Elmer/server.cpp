#include "server.h"
#include "util.h"
#include "shell_server_config/shell_server_config.h"

namespace elmer {
    Server::Server(python::handler* _h) {
        // ULOG("Setting up Elmer Server");
        PyObject* server_config = _h->loadObjectWithSetupFromMain("server");
        python::loadAttrFromObjectAndConvert(server_config, "start_file",  this->start_file);
        python::loadAttrFromObjectAndConvert(server_config, "end_file",    this->end_file);
        python::loadAttrFromObjectAndConvert(server_config, "run_file",    this->run_file);
        python::loadAttrFromObjectAndConvert(server_config, "done_file",   this->done_file);
        python::loadAttrFromObjectAndConvert(server_config, "server_file", this->server_code_output_file);
        python::loadAttrFromObjectAndConvert(server_config, "command",     this->command);
        python::loadAttrFromObjectAndConvert(server_config, "head",        this->head);
        python::loadAttrFromObjectAndConvert(server_config, "tail",        this->tail);
    }

    /* ---------------------------------------------------------------------- */

    Server::~Server() {
        // ULOG("Calling destructor");
        this->end();
        this->deleteServerFile(); 
    }

    /* ---------------------------------------------------------------------- */

    void Server::makeServerFile() {
        util::oFile out(this->server_code_output_file);
        out << this->head << "\n\n";
        out << "START_FILE=" << this->start_file << "\n";
        out << "END_FILE="   << this->end_file   << "\n";
        out << "RUN_FILE="   << this->run_file   << "\n";
        out << "DONE_FILE="  << this->done_file  << "\n";
        out << SHELL_SERVER_STRING;
        out << "\n\n" << this->tail;
    }

    /* ---------------------------------------------------------------------- */

    void Server::start() {
        if (util::fileExists(this->done_file))
            remove(this->done_file);
        if (util::fileExists(this->run_file))
            remove(this->run_file);
        if (util::fileExists(this->end_file))
            remove(this->end_file);
        if (util::fileExists(this->start_file))
            remove(this->start_file);

        ULOG("Starting Server");
        util::oFile out(this->start_file);
        out.close();
        // ULOG("Waiting for server to confirm");
        ULOG("Run the following command in a seperate terminal: bash " + this->server_code_output_file);
        while (util::fileExists(this->start_file)) { sleep(1); }
        // ULOG("Server Confirmed, continuing");
    }

    /* ---------------------------------------------------------------------- */

    void Server::end() {
        ULOG("Shutting Down Server");
        util::oFile out(this->end_file);
        out.close();

        // ULOG("Waiting for server to confirm");
        while (util::fileExists(this->end_file)) { sleep(1); }
        // ULOG("Server Confirmed, continuing");
    }

    /* ---------------------------------------------------------------------- */

    void Server::runCommand() {
        ULOG("Sending command to command server");
        util::oFile out(this->run_file);
        out << command << "\n";
    }

    /* ---------------------------------------------------------------------- */

    void Server::deleteServerFile() { remove(this->server_code_output_file); }

    /* ---------------------------------------------------------------------- */

    void Server::waitForDoneFile() {
        ULOG("Waiting for Server Command Completion");
        std::string line;
        while (!(util::fileExists(done_file))) { sleep(1); }
        util::iFile read_this(this->done_file);
        read_this.getLine(line);
        read_this.close();
        util::trim(line);
        remove(done_file);
        if (std::stoi(line)) {
            this->deleteServerFile();
            UERR("Server returned non-zero exit code: " + line);
        }   
        ULOG("Got server command completion, continuing");
    }
} // namespace elmer