import inspect, os

def _checkType(_msg:str, _var:any, *args) -> None:
    # just going through the inputted types and checking them
    for _type in args:
        if isinstance(_var, _type):
            return
    
    # if here then none were a match
    raise TypeError(f"From {inspect.stack()[1][3]}: {_msg}")

def error(
        msg:str, # message to send to user
        exeception=True # whether to raise exception, msg will be exception message
    ) -> None:
    if exeception: raise Exception(f"From \'{inspect.stack()[1][3]}\': {msg}")
    else:
        print(f"Error from {inspect.stack()[1][3]}: {msg}")
        exit(1)
        
def get_cur_file_dir() -> str:
    """
    - **Description**: returns the name of the directory where the file calling this function is located
    - **inputs**: (nothing)
    - **returns**:
        - `file_dir:str`: the name of the directory where the file calling this function is located
    - **Example Usage** (some ways to use it, not all)
        ```python
        import my
        my.get_cur_file_dir()
        ```
    """
    return os.path.sep.join(inspect.stack()[1].filename.split(os.path.sep)[:-1])


def get_cur_filename() -> str:
    """
    - **Description**: returns the name of the file that this function is called from
    - **inputs**: (nothing)
    - **returns**:
        - `filename:str`: he name of the file that this function is called from
    - **Example Usage** (some ways to use it, not all)
        ```python
        import my
        my.get_cur_filename()
        ```
    """
    return inspect.stack()[1].filename

def local_path(*args) -> str:
    """
    - **Description**: prepends the inputted file path with the directory of the file this function is called from
    - **inputs**:
        - `args`: items to join together as a local path
    - **returns**:
        - `new_file_path:str`: the inputted file path with the directory of the file this function is called from added to the beginning
    - **Example Usage** (some ways to use it, not all)
        ```python
        import my
        filename = "local_image.png"
        my.local_path("out", filename)
        ```
    """
    _path = os.path.sep.join(args)
    return os.path.join(os.path.sep.join(inspect.stack()[1].filename.split(os.path.sep)[:-1]), _path)