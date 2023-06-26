import sys, os
from ._src import get_cur_file_dir

def dump2cfg(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "dump2cfg.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def dump2xyz(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "dump2xyz.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def dumpsort(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "dumpsort.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def grid_refine(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "grid_refine.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def implicit_grid(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "implicit_grid.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def jagged2d(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "jagged2d.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def jagged3d(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "jagged3d.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def log2txt(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "log2txt.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def logplot(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "logplot.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def stl2surf(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "stl2surf.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def surf_create(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "surf_create.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def surf_transform(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "surf_transform.py"), *[str(arg) for arg in args]])
    os.system(_cmd)
