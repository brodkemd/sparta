from collections import OrderedDict

# surf_file: surface file used in the simulation
# output_file: file to output to
def make_surf_temperature_file(surf_file:str, output_file:str, temperature:float, change_temperatures:dict[int, float]={}):
    
    # reading the surf file to get the surface element ids
    with open(surf_file, 'r') as f:
        lines = [line.strip().split() for line in f.readlines() if len(line.strip())]

    # finding where the surface elements start
    start = None
    for i in range(len(lines)):
        if lines[i][0].lower() == "triangles":
            start = i
            break

    # just getting the portion that contains the surface elements
    lines = lines[start+1:]
    
    # setting the temperatures of the surface elements
    surf_elements = OrderedDict()
    for item in lines:
        id = int(item[0])
        
        # if the surface element does not need a custom temperature
        if id not in change_temperatures:
            surf_elements[id] = float(temperature)
        else:
            surf_elements[id] = float(change_temperatures[id])
    
    # the header for the file
    header = "\n".join([
        "ITEM: TIMESTEP", 
        "NAN",
        "ITEM: NUMBER OF SURFS",
        "NAN",
        "ITEM: BOX BOUNDS NAN NAN NAN",
        "-NAN NAN",
        "-NAN NAN",
        "-NAN NAN",
        "ITEM: SURFS id c"
    ])
    
    # writing the result to the output file
    with open(output_file, 'w') as f:
        f.write(header + "\n")
        
        # writing the temperature file
        for item in surf_elements: f.write(f"{item} {surf_elements[item]}\n")