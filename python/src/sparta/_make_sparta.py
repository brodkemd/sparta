import os
cwd = os.path.sep.join(__file__.split(os.path.sep)[:-1])
in_f = os.path.join(cwd, "_sparta_template.py")
out_f = os.path.join(cwd, "sparta.py")
doc_src = "doc"

commands = {
    "adapt_grid.html": "adapt_grid",
    "balance_grid.html": "balance_grid",
    "boundary.html": "boundary",
    "bound_modify.html": "bound_modify",
    "clear.html": "clear",
    "collide.html": "collide",
    "collide_modify.html": "collide_modify",
    "compute.html": "compute",
    "create_box.html": "create_box",
    "create_grid.html": "create_grid",
    "create_particles.html": "create_particles",
    "dimension.html": "dimension",
    "dump.html": "dump",
    "dump_modify.html": "dump_modify",
    "echo.html": "echo",
    "fix.html": "fix",
    "global.html": "global",
    "group.html": "group",
    "if.html": "if",
    "include.html": "include",
    "jump.html": "jump",
    "label.html": "label",
    "log.html": "log",
    "mixture.html": "mixture",
    "move_surf.html": "move_surf",
    "next.html": "next",
    "package.html": "package",
    "partition.html": "partition",
    "print.html": "print",
    "quit.html": "quit",
    "react.html": "react",
    "react_modify.html": "react_modify",
    "read_grid.html": "read_grid",
    "read_isurf.html": "read_isurf",
    "read_particles.html": "read_particles",
    "read_restart.html": "read_restart",
    "read_surf.html": "read_surf",
    "region.html": "region",
    "remove_surf.html": "remove_surf",
    "reset_timestep.html": "reset_timestep",
    "restart.html": "restart",
    "run.html": "run",
    "scale_particles.html": "scale_particles",
    "seed.html": "seed",
    "shell.html": "shell",
    "species.html": "species",
    "species_modify.html": "species_modify",
    "stats.html": "stats",
    "stats_modify.html": "stats_modify",
    "stats_style.html": "stats_style",
    "suffix.html": "suffix",
    "surf_collide.html": "surf_collide",
    "surf_react.html": "surf_react",
    "surf_modify.html": "surf_modify",
    "timestep.html": "timestep",
    "uncompute.html": "uncompute",
    "undump.html": "undump",
    "unfix.html": "unfix",
    "units.html": "units",
    "variable.html": "variable",
    "write_grid.html": "write_grid",
    "write_isurf.html": "write_isurf",
    "write_restart.html": "write_restart",
    "write_surf.html": "write_surf"
}

indent = 4
tab = " " * indent

to_add = ""
command_docs_dict = []
for command in commands:
    docs = command
    command = commands[command]
    command = command[0].upper() + command[1:]
    command_docs_dict.append(f"\"{command}\" : f\"{{cwd}}{os.path.sep}{doc_src}{os.path.sep}{docs}\"")

    # adding function format
    to_add+=f"""
    \"\"\"
    run  
        \'class_instance.man(class_instance.{command})\'
    to see more info on this command.
    (class_instance is the class that you are calling this function from)
    \"\"\"
    def {command}(self, args:list|str) -> None:
        if isinstance(args, list):
            self.commands.append([\"{command.lower()}\"]+args)
        elif isinstance(args, str):
            self.commands.append([\"{command.lower()} \"+args])
        else:
            raise TypeError("invalid input type for args")
"""

print("reading template from:", in_f)
with open(in_f, "r") as f: lines = f.read()

lines = lines.strip()
lines = lines[:lines.rfind("\n")]

# print("inserting docs")
# lines = lines.replace(" _docs:dict[str, str] = {}", " _docs:dict[str, str] = {" + (",\n"+tab+tab+tab).join(command_docs_dict) + "}")

print("adding commands")
lines = lines + to_add

print("unloading file to:", out_f)
with open(out_f, 'w') as f:
    f.write(lines + to_add)


del cwd, in_f, out_f, doc_src, commands, indent, tab, to_add, command_docs_dict, command, docs, lines, f