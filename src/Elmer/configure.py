import os, glob

print("Updating style headers")

cwd = os.path.dirname(os.path.abspath(__file__))
src = os.path.dirname(cwd)


names = {
    "collide"      : ["collide_*.h"],
    "command"      : [
                        "adapt_*.h",
                        "balance_*.h",
                        "create_*.h",
                        "move_*.h",
                        "read_*.h",
                        "remove_*.h",
                        "run.h",
                        "scale_*.h",
                        "write_*.h"
                     ],
    "compute"      : ["compute_*.h"],
    "dump"         : ["dump_*.h"],
    "fix"          : ["fix_*.h"],
    "react"        : ["react_*.h"],
    "region"       : ["region_*.h"],
    "surf_collide" : ["surf_collide_*.h"],
    "surf_react"   : ["surf_react_*.h"]
}

ignore = {
    "collide"      : [],
    "command"      : [],
    "compute"      : [],
    "dump"         : [],
    "fix"          : ["fix_emit.h"],
    "react"        : ["react_bird.h"],
    "region"       : [],
    "surf_collide" : [],
    "surf_react"   : []
}

for name in names:
    files = []
    for arg in names[name]:
        files += glob.glob(os.path.join(src, arg))

    print("  Writing:", os.path.join(os.path.split(src)[-1], f"style_{name}.h"))
    with open(os.path.join(src, f"style_{name}.h"), 'w') as f:
        for file in files:
            file = os.path.split(file)[-1]
            if file not in ignore[name]:
                f.write(f"#include \"{file}\"\n")
            else:
                print("    Ignoring header:", file)



