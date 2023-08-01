RUN_FILE=$1
EXIT_FILE=$2

while [ 1 ]; do
    if [ -f $RUN_FILE ]; then
        echo $RUN_FILE" exists. Waiting"
        sleep 2
        while IFS= read -r line; do
            eval "$line"
            break
        done < $RUN_FILE
        rm $RUN_FILE
    fi
    if [ -f $EXIT_FILE ]; then
        echo "Caught exit file, exiting ..."
        rm $EXIT_FILE
        break
    fi
    sleep 1
done