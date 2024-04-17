from tomli import load
import os

def get_from_file(_file) -> dict:
    with open(_file, 'rb') as f: return load(f)

data = {}
recurse_num = 0
max_recurse = 10
strict = True
keys = {
    "HOME" : os.path.expanduser('~'),
    "PWD"  : os.getcwd()
}

def error(_msg:str):
    if strict: raise Exception(_msg)

def get_arr(_data, _name:str):
    start_char = "["
    end_char   = "]"
    if start_char in _name:
        _data = _data[_name[:_name.find(start_char)]]
        for i in range(_name.count(start_char)):
            index = int(_name[_name.find(start_char)+len(start_char): _name.find(end_char)].strip())
            if not isinstance(_data, list):
                error("can not index list: " + _name)
            try: _data = _data[index]
            except IndexError as e: error(str(e))
            _name = _name[_name.find(end_char)+len(end_char):]
    return _data, _name

def get_from(_data:dict, _path:str, _orig_path:str=None):
    
    if _orig_path is None: _orig_path = _path
    sep = "."
    
    if _path in keys:
        return keys[_path]
    
    if sep in _path:
        _name = _path[:_path.find(sep)]
        _data, _name = get_arr(_data, _name)
        # print("name:", _name)
        if len(_name) == 0: return get_from(_data, _path[_path.find(sep)+len(sep):], _orig_path)
        elif _name in _data: return get_from(_data[_name], _path[_path.find(sep)+len(sep):], _orig_path)
        else: error("can not resolve path, " + _orig_path)
    else:
        _data, _path = get_arr(_data, _path)
        if len(_path) == 0: # this means is it an array value
            return resolve_name(_data)
        elif _path in _data:
            if isinstance(_data[_path], str): _data[_path] = resolve_name(_data[_path])
            return _data[_path]
        else: error("can not resolve path, " + _orig_path)

def resolve_name(_name:str):
    global data, recurse_num, max_recurse
    if recurse_num > max_recurse:
        error(f"reached max allowed recursion resolving {_name} (there is probably a circular definition somewhere)")
    start = "$("
    end   = ")"
    while _name.count(start):
        recurse_num+=1
        key = _name[_name.find(start)+len(start):_name.find(end, _name.find(start)+len(start))]
        _name = _name[:_name.find(start)] + str(get_from(data, key)) + _name[_name.find(end, _name.find(start)+len(start)) + len(end):]
    return _name

def iter(_data):
    global recurse_num
    for item in _data:
        recurse_num = 0
        if isinstance(_data[item], dict): iter(_data[item])
        else:
            if isinstance(_data[item], list):
                str_indicies = []
                for i, val in enumerate(_data[item]):
                    if isinstance(val, str): str_indicies.append(i)
                for i in range(len(_data[item])):
                    val = _data[item][i]
                    if isinstance(val, str):
                        _data[item][i] = resolve_name(val)
                        if _data[item][i] != val:
                            if i in str_indicies: str_indicies.remove(i)
                if len(str_indicies):
                    for i, val in enumerate(_data[item]): _data[item][i] = str(val)
            elif isinstance(_data[item], str): _data[item] = resolve_name(_data[item])

def get_at_path(_path):
    global data
    return get_from(data, _path)

def Main(_file):
    global data
    data = get_from_file(_file)
    iter(data)