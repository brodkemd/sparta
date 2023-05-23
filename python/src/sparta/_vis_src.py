import glob, multiprocessing, time, traceback
from types import FunctionType
from queue import Empty

from ._src import error

# handles the kwargs for the functions below
def _handle_kwargs(input_kwargs:dict, allowed_kwargs:dict, strict:bool=False) -> dict:
    for item in allowed_kwargs:
        # if there is a key error then reverts to default
        try:
            # setting the class object from the input dict
            allowed_kwargs[item] = input_kwargs[item]
            # removing key from the inputted dict
            del input_kwargs[item]
        except KeyError:
            pass
    
    # if there is strict enforcement, then caught if any setting are left over
    if strict:
        if len(list(input_kwargs.keys())):
            for item in list(input_kwargs.keys()):
                raise KeyError(f"\"{item}\" is not a valid argument")
    
    return allowed_kwargs


# handles file list and glob paths needed for the video generating function below
def handle_file_directive(**kwargs):
    max_len = 0
    directed = {}
    should_sort = {}
    
    # going through the inputted parameters
    for arg in kwargs:
        # the files (list or str or glob pattern)
        files = kwargs[arg]
        should_sort[arg] = False # if these file(s) need sorted\
        if files is None: continue
        
        # parsing the glob pattern and keeping it if it is a list
        if isinstance(files, str): files = glob.glob(files)
        elif not isinstance(files, list[str]): raise TypeError(f"invalid type for \"{arg}\"")
        
        # getting the largest length
        if len(files) > max_len:  max_len = len(files)
        
        # adding to the list of args that have a file list that is not None
        directed[arg] = files
    
    # checking to make sure each arg has the correct number of files with it
    # seeing if it needs to be sorted later
    for arg in directed:
        files = directed[arg]

        # checking to make sure the files list is a correct size
        if len(files) != max_len and len(files) != 1:
            error(f"not enough files supplied for \"{arg}\"")
        else:
            # making the list the correct length if it is only one element
            if len(files) == 1:
                directed[arg] = files*max_len
            else:
                # if the list is more than one element, then it needs to be sorted
                should_sort[arg] = True
    
    return directed, max_len, should_sort


'''
Function to run on threads

'''
def thread_func(
        func:FunctionType,
        queue_in:multiprocessing.Queue, # queue to get data from
        queue_out:multiprocessing.Queue, # queue to get data from
        stop,      # shared mem integer to indicate for this function to stop
    ) -> None:
    print(multiprocessing.current_process().name, "starting up")
    
    # max time allowed for threads to wait for data
    timeout_time = 2.
    
    try:
        last_update = time.time()
        
        # runs until the queue is empty
        while not queue_in.empty():
            # trying to get data from the queue, if none is available then it checks to make sure
            # the max time to wait for data has not been reached.
            # If the max wait time is reached then this process exits
            try:
                data = queue_in.get(block=False) # getting latest object from the queue
            except Empty:
                if time.time() - last_update > timeout_time:
                    print(multiprocessing.current_process().name, "caught timeout shutting down")
                    break
                continue
            
            
            print(multiprocessing.current_process().name, "making id", f"{data[0]}")
            
            # making an image based on the file names from the queue
            out_file_name = func(**data[1], **data[-1])
            
            # setting the last time an action was perform to the current time (just made an image)
            last_update = time.time()
            
            # outputing the data
            queue_out.put([data[0], out_file_name])

            # if this function should stop (if the int is set to 1 from the main thread)
            if stop.value:
                print(multiprocessing.current_process().name, " caught end signal shutting down")
                break

        print(multiprocessing.current_process().name, "done")
    
    # catches exception (primarily for keyboard interrupt)
    except Exception as e:
        print(multiprocessing.current_process().name, " caught error:", traceback.format_exc())

