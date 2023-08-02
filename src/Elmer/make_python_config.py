
header_name = "python_config"

header_format = f"""
#ifndef {header_name.upper()}_H
#define {header_name.upper()}_H

#include <Python.h>

#define PYTHON_STRING "{{contents}}"

#endif
"""

with open("toml.py", 'r') as f:
    header_format = header_format.format(
        contents=f.read().replace("\n", "\\n").replace("'", "\'").replace("\"", "\\\"")
    )

with open(header_name + ".h", 'w') as f:
    f.write(header_format.strip())