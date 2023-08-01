#ifndef SHELL_SERVER_CONFIG_H
#define SHELL_SERVER_CONFIG_H

#define SHELL_SERVER_STRING "RUN_FILE=$1\nEXIT_FILE=$2\n\nwhile [ 1 ]; do\n    if [ -f $RUN_FILE ]; then\n        echo $RUN_FILE\" exists. Waiting\"\n        sleep 2\n        while IFS= read -r line; do\n            eval \"$line\"\n            break\n        done < $RUN_FILE\n        rm $RUN_FILE\n    fi\n    if [ -f $EXIT_FILE ]; then\n        echo \"Caught exit file, exiting ...\"\n        rm $EXIT_FILE\n        break\n    fi\n    sleep 1\ndone"

#endif