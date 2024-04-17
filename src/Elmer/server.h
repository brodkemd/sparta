#ifndef SERVER_H
#define SERVER_H

#include "python.h"

namespace elmer {
    class Server {
        private:
            char *simulation_directory, *start_file, *end_file, *done_file, *run_file, *server_code_output_file, *command, *head, *tail;
        
        public:
            Server() {}
            Server(python::handler* _h);
            ~Server();
            // void set(python::handler& _h);
            void makeServerFile();
            void start();
            void end();
            void runCommand();
            void waitForDoneFile();
            void deleteServerFile();
    };
} // namespace elmer


#endif