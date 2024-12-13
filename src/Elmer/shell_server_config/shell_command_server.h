#ifndef SHELL_COMMAND_SERVER_H
#define SHELL_COMMAND_SERVER_H

#define SHELL_SERVER_STRING "\nstty -echoctl # hide ^C\nto_exit() {\n    echo \"1\" > $DONE_FILE\n    touch $END_FILE\n}\n\ntrap 'to_exit' SIGINT\n\nwhile [ 1 ]; do\n    if [ -f $START_FILE ]; then\n        rm $START_FILE\n    fi\n    if [ -f $RUN_FILE ]; then\n        while IFS= read -r line; do\n            eval \"$line\"\n            echo $? > $DONE_FILE\n            break\n        done < $RUN_FILE\n        rm $RUN_FILE\n    fi\n    if [ -f $END_FILE ]; then\n        echo \"Caught exit file, exiting ...\"\n        rm $END_FILE\n        break\n    fi\n    sleep 1\ndone"

#endif