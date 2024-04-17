import sys, os
from ._src import get_cur_file_dir

def surf2paraview(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "paraview", "surf2paraview.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def grid2paraview(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "paraview", "grid2paraview.py"), *[str(arg) for arg in args]])
    os.system(_cmd)

def grid2paraview_cells(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", "paraview", "grid2paraview_cells.py"), *[str(arg) for arg in args]])
    os.system(_cmd)
