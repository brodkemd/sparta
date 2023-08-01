
header_name = "shell_server_config"

header_format = f"""
#ifndef {header_name.upper()}_H
#define {header_name.upper()}_H

#define SHELL_SERVER_STRING "{{contents}}"

#endif
"""

with open("shell_command_server.sh", 'r') as f:
    header_format = header_format.format(
        contents=f.read().replace("\n", "\\n").replace("'", "\'").replace("\"", "\\\"")
    )

with open(header_name + ".h", 'w') as f:
    f.write(header_format.strip())