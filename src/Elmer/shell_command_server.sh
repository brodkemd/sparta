
stty -echoctl # hide ^C
to_exit() {
    echo "1" > $DONE_FILE
    touch $END_FILE
}

trap 'to_exit' SIGINT

while [ 1 ]; do
    if [ -f $START_FILE ]; then
        rm $START_FILE
    fi
    if [ -f $RUN_FILE ]; then
        echo $RUN_FILE" exists. Waiting"
        sleep 2
        while IFS= read -r line; do
            eval "$line" && echo $? > $DONE_FILE
            break
        done < $RUN_FILE
        rm $RUN_FILE
    fi
    if [ -f $END_FILE ]; then
        echo "Caught exit file, exiting ..."
        rm $END_FILE
        break
    fi
    sleep 1
done