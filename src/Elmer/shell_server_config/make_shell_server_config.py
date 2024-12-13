import os

print("Making Shell Server Config")

base_file_name = "shell_command_server"

file = os.path.abspath(os.path.join(os.path.dirname(__file__), base_file_name))

header_format = f"""
#ifndef {base_file_name.upper()}_H
#define {base_file_name.upper()}_H

#define SHELL_SERVER_STRING "{{contents}}"

#endif
"""

print("  Writing:", file+".sh")
with open(file+".sh", 'r') as f:
    header_format = header_format.format(
        contents=f.read().replace("\n", "\\n").replace("'", "\'").replace("\"", "\\\"")
    )

print("  Writing:", file+".h")
with open(file+".h", 'w') as f:
    f.write(header_format.strip())