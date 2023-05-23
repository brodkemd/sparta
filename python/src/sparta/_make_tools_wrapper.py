import os, sys
try:
    from ._src import get_cur_file_dir
except:
    from _src import get_cur_file_dir

look_at = os.path.join(get_cur_file_dir(), "tools")

header = """import sys, os
from ._src import get_cur_file_dir
"""

function_format = """
def {}(*args):
    _cmd = " ".join([sys.executable, os.path.join(get_cur_file_dir(), "tools", \"{}\"), *[str(arg) for arg in args]])
    os.system(_cmd)
"""

with open(os.path.join(get_cur_file_dir(), "_tools_wrapper.py"), 'w') as f:
    f.write(header)
    for item in os.listdir(look_at):
        if item.endswith(".py"):
            func = item.strip(".py")
            f.write(function_format.format(func, item))

del look_at, function_format, f, func, item, header
        