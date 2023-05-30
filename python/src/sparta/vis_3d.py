import os, ffmpeg, multiprocessing, time
from collections import OrderedDict
from types import FunctionType

from mayavi import mlab
import matplotlib as mpl
import numpy as np

from ._vis_src import thread_func, handle_file_directive, _handle_kwargs
from ._src import error

_TOTAL_KWARGS_3D = {
    # arg                      default value
    # "surf_plot_style"        : "b-",
    # "surf_plot_args"         : dict(),
    # "particle_edge_color"    : "black",
    # "particle_line_width"    : 0.1,
    # "legend_args"            : dict(loc='upper right', ncol=1, fontsize=16),
    # "legend_particle_size"   : 80,
    "max_cell_size"          : None,
    "include_particle_every" : 1,
    "surf_color"             : (0.0,0.5,0.5),
    "add_bounding_box"       : False,
    "use_existing_figure"    : False,
    "particle_size"          : 0.01,
    "line_width"             : 1.0,
    "figsize"                : (600, 600),
    "image_size"             : (600, 600),
    "view"                   : None
}

"""
Kwargs options:
    - surf_plot_style               => style of the surface plot
    - surf_plot_args:dict[str, any] => arguements to pass to the plot of the surface function,
                                        this is any valid argument for plt.plot
    - particle_edge_color:str       => color for the edge of the particle
    - include_particle_every:int    => only include a particle after skipping over this many (used to reduce render time)
    - particle_line_width:float     => width of the border of the particles
    - legend_args:dict[str, any]    => args to pass to the plt.legend command in matplotlib
    - legend_particle_size:float    => size of the particle in the legend
    - particle_size:float           => size of the particle in the plot
    - box_line_size:float           => thickness of bounding box line
    - scale_fig_by:float            => fraction to scale the file by
    - figsize:tuple[int]            => size to make the image
    
"""

def _add_boundary_box(x_lims, y_lims, z_lims, line_width):   
    x = []
    y = []
    z = []
    connections = []
    N = 2
    k = 1
    for i in range(2):
        for j in range(2):
            x.append(np.linspace(x_lims[j], x_lims[j], N))
            y.append(np.linspace(y_lims[0], y_lims[1], N))
            z.append(np.linspace(z_lims[i], z_lims[i], N))
            connections.append(np.vstack((np.arange((k-1)*N, k*N-1, 1), np.arange((k-1)*N + 1, k*N, 1))).T)
            k+=1
    
    for i in range(2):
        for j in range(2):
            x.append(np.linspace(x_lims[0], x_lims[1], N))
            y.append(np.linspace(y_lims[j], y_lims[j], N))
            z.append(np.linspace(z_lims[i], z_lims[i], N))
            connections.append(np.vstack((np.arange((k-1)*N, k*N-1, 1), np.arange((k-1)*N + 1, k*N, 1))).T)
            k+=1
    
    for i in range(2):
        for j in range(2):
            x.append(np.linspace(x_lims[j], x_lims[j], N))
            y.append(np.linspace(y_lims[i], y_lims[i], N))
            z.append(np.linspace(z_lims[0], z_lims[1], N))
            connections.append(np.vstack((np.arange((k-1)*N, k*N-1, 1), np.arange((k-1)*N + 1, k*N, 1))).T)
            k+=1
    
    connections = np.vstack(connections)
    x = np.hstack(x)
    y = np.hstack(y)
    z = np.hstack(z)
    src = mlab.pipeline.scalar_scatter(x, y, z, np.ones(len(x)))   
    src.mlab_source.dataset.lines = connections
    src.update()

    # The stripper filter cleans up connected lines
    lines = mlab.pipeline.stripper(src)
    
    mlab.pipeline.surface(lines, colormap='gray', line_width=line_width)



def grid_3d(grid_file:str, output_file:str=None, **kwargs) -> None:
    # kwargs for this function
    Kwargs = _TOTAL_KWARGS_3D
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)
    
    # reading the file
    with open(grid_file, 'r') as f: lines = [line.strip().split() for line in f.readlines()]

    #ids = []
    centers = []
    #bounds = []
    params = []
    # loading the data from the lines
    for i, line in enumerate(lines):
        # if the line is not empty
        if len(line) > 1:
            # getting the box dimensions
            if line[1] == "BOX":
                x_lims = [float(item) for item in lines[i+1]]
                y_lims = [float(item) for item in lines[i+2]]
                z_lims = [float(item) for item in lines[i+3]]
            
            # loading the grid cells
            elif line[1] == "CELLS":
                # getting the indices of box parameters
                _x_lo_ind = line.index("xlo") - 2
                _y_lo_ind = line.index("ylo") - 2
                _z_lo_ind = line.index("zlo") - 2
                _x_hi_ind = line.index("xhi") - 2
                _y_hi_ind = line.index("yhi") - 2
                _z_hi_ind = line.index("zhi") - 2
                
                _x_c_ind = line.index("xc") - 2
                _y_c_ind = line.index("yc") - 2
                _z_c_ind = line.index("zc") - 2
                
                # loading the data from the rest of the lines
                for j in range(i+1, len(lines)):
                    # if (j - i - 1) == 1:
                    #     print(lines[j])
                    #ids.append(lines[j][0])
                    centers.append([
                        float(lines[j][_x_c_ind]),
                        float(lines[j][_y_c_ind]),
                        float(lines[j][_z_c_ind])
                    ])
                    # getting the bounds
                    params.append([
                        abs(float(lines[j][_x_hi_ind]) - float(lines[j][_x_lo_ind])),
                        abs(float(lines[j][_y_hi_ind]) - float(lines[j][_y_lo_ind])),
                        abs(float(lines[j][_z_hi_ind]) - float(lines[j][_z_lo_ind]))
                    ])
                    # bounds.append([
                    #     [float(lines[j][_x_lo_ind]), float(lines[j][_x_hi_ind])],
                    #     [float(lines[j][_y_lo_ind]), float(lines[j][_y_hi_ind])],
                    #     [float(lines[j][_z_lo_ind]), float(lines[j][_z_hi_ind])]
                    # ])
                break
    
    if not Kwargs["use_existing_figure"]:
        # creating figure using the size detected in the file
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=Kwargs["figsize"])

    x = []
    y = []
    z = []
    connections = []
    N = 2
    k = 1
    if Kwargs["max_cell_size"] is not None:
        first_len = len(centers)
        centers, params = zip(*((center, param) for center, param in zip(centers, params) if max(param[:2]) < Kwargs["max_cell_size"]))
        centers = list(centers)
        params = list(params)
        
        if first_len > len(centers):
            _add_boundary_box(x_lims, y_lims, z_lims, 1.)
    
    # for box in bounds:
    for i in range(len(centers)):
        for j in range(3):
            x.append(centers[i][0] - params[i][0]/2)
            y.append(centers[i][1] - params[i][1]/2)
            z.append(centers[i][2] - params[i][2]/2)

            x.append(centers[i][0] + params[i][0]/2 if j == 0 else centers[i][0] - params[i][0]/2)
            y.append(centers[i][1] + params[i][1]/2 if j == 1 else centers[i][1] - params[i][1]/2)
            z.append(centers[i][2] + params[i][2]/2 if j == 2 else centers[i][2] - params[i][2]/2)
            connections.append(np.vstack((np.arange((k-1)*N, k*N-1, 1), np.arange((k-1)*N + 1, k*N, 1))).T)
            k+=1

    connections = np.vstack(connections)
    
    src = mlab.pipeline.scalar_scatter(x, y, z)
    src.mlab_source.dataset.lines = connections
    src.update()
    
    # The stripper filter cleans up connected lines
    lines = mlab.pipeline.stripper(src)
    
    mlab.pipeline.surface(lines, colormap='gray', line_width=1.)


    if not Kwargs["use_existing_figure"]:
        if Kwargs["view"] is not None: mlab.view(*Kwargs["view"])
        if output_file is None:
            mlab.show()
        else:
            mlab.savefig(output_file, Kwargs["image_size"])
            mlab.close(mlab.gcf())    

def surf_3d(
        surf_file:str, # file to read the surface from
        output_file:str         = None,
        surf_file_is_dump:bool  = False, # if the surface is from a dump
        **kwargs
    ):
    # kwargs for this function
    Kwargs = _TOTAL_KWARGS_3D.copy()
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)

    surf_coords:list[list[float]] = [] # x,y,z coords for surface
    tris:list[list[float]] = [] # triangles for surface

    # loads lines from surface file
    lines = []
    with open(surf_file, 'r') as f: lines = f.readlines()

    # if the surface is from a dump from sparta
    if surf_file_is_dump:
        surfs_start:int = 0 # where the surfs start
        x_ind:int = 0  # x-coord index in the file
        y_ind:int = 0  # y-coord index in the file
        z_ind:int = 0  # z-coord index in the file
        for i, line in enumerate(lines):
            new_line = line.strip().split() # splits line into list
            if len(new_line)>=2:
                if new_line[1] == "SURFS":
                    # where the surface data starts
                    surfs_start = i+1
                    for j, item in enumerate(new_line):
                        if item == "v1x":   x_ind = j - 2 # the x coord of the surface
                        elif item == "v1y": y_ind = j - 2 # the y coord of the surface
                        elif item == "v1z": z_ind = j - 2 # the z coord of the surface
                    break # no longer need anything so exit the loop
    # if the surface is NOT from a dump from sparta
    # i.e. it is an input file into sparta
    else:
        surfs_start:int = 0 # where the surfs start
        tris_start:int = 0 # where the triangle data starts
        x_ind:int = 1 # x-coord index in the file
        y_ind:int = 2 # y-coord index in the file
        z_ind:int = 3  # z-coord index in the file
        for i, line in enumerate(lines):
            new_line = line.strip().split() # splitting the line into a list
            if len(new_line): # if there is data in the line
                if new_line[0] == "Points":
                    # where the surface data starts
                    surfs_start = i+2
                    
                if new_line[0] == "Triangles":
                    tris_start = i+2
                    break # no longer need any info so exits loop

    # iterates through the loaded lines
    for i in range(surfs_start, len(lines)):
        data = lines[i].strip().split() # splits lines into list
        # if the length of the line is 0 then there is no more data, so exit loop
        if not len(data): break
        # adds list to list (will contain surface order pairs)
        surf_coords.append([float(data[x_ind].strip()), float(data[y_ind].strip()), float(data[z_ind].strip())]) # grabs x,y,z coord of surface point

    # iterates through the loaded lines
    for i in range(tris_start, len(lines)):
        data = lines[i].strip().split() # splits lines into list
        # if the length of the line is 0 then there is no more data, so exit loop
        if not len(data): break
        # adds list to list (will contain surface order pairs)
        tris.append([int(data[1].strip())-1, int(data[2].strip())-1, int(data[3].strip())-1]) # grabs verticies of triangle

    del lines # deletes, no longer need var

    # getting the surface data from the lists
    x = [i[0] for i in surf_coords]
    y = [i[1] for i in surf_coords]
    z = [i[2] for i in surf_coords]

    # adding the first coords to close the surface
    x.append(surf_coords[0][0])
    y.append(surf_coords[0][1])
    z.append(surf_coords[0][2])

    if not Kwargs["use_existing_figure"]:
        # creating figure using the size detected in the file
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=Kwargs["figsize"])

    # plotting the surface, using the inputted style and the inputted arguements for it
    mlab.triangular_mesh(x, y, z, tris, color=Kwargs["surf_color"])
    
    if not Kwargs["use_existing_figure"]:
        if Kwargs["view"] is not None: mlab.view(*Kwargs["view"])
        if output_file is None:
            mlab.show()
        else:
            mlab.savefig(output_file, Kwargs["image_size"])
            mlab.close(mlab.gcf())


# class for 3d particle image
def particle_3d(
        particle_file:str, # file to read data from
        output_file:str              = None, # file to output video to
        colors:OrderedDict[str, str] = {}, # structure = {species name : "color"}
        **kwargs
    ):
    # kwargs for this function
    Kwargs = _TOTAL_KWARGS_3D
    Kwargs = _handle_kwargs(kwargs, Kwargs)

    # reading the lines from the input file
    lines = []
    with open(particle_file, 'r') as f: lines = f.readlines()
    
    atoms_start:int = 0 # where the atoms start
    x_ind:int = 0  # x-coord index in the file
    y_ind:int = 0  # y-coord index in the file
    z_ind:int = 0  # y-coord index in the file
    id_ind:int = 0 # species id of the particle
    x_lims:list[float] = [] # upper and lower bound for x
    y_lims:list[float] = [] # upper and lower bound for y
    z_lims:list[float] = [] # upper and lower bound for y

    # iterating through the lines loaded from the input file
    for i, line in enumerate(lines):
        new_line = line.strip().split() # turning line into list by splitting at whitespace
        if len(new_line)>=2:
            # if lines contain info on atoms
            if new_line[1] == "ATOMS":
                atoms_start = i+1 # the data on the atoms start on the next line
                # setting the indicies from the position of column indicators, 
                # they are shifted by -2 to account start strings in the line
                for j, item in enumerate(new_line):
                    if item == "x":   x_ind = j - 2
                    elif item == "y": y_ind = j - 2
                    elif item == "z": z_ind = j - 2
                    elif item == "type": id_ind = j - 2
                break # last thing that we are looking for so breaks loop
            elif new_line[1] == "BOX":
                # gets the bounds of the simulation box
                x_lims = [float(item) for item in lines[i+1].strip().split()]
                y_lims = [float(item) for item in lines[i+2].strip().split()]
                z_lims = [float(item) for item in lines[i+3].strip().split()]

    # creating an ordered dict to contain info on each molecule id
    particle_coords:OrderedDict[str, list[float]] = OrderedDict() # x,y,z coords for each molecule
    for i, item in enumerate(colors):
        # molecule id is i+1
        particle_coords[i+1] = []
    
    # loading particle data from loaded lines, only loading a particle every "include_particle_every" increment
    for i in range(atoms_start, len(lines), Kwargs["include_particle_every"]):
        data = lines[i].strip().split() # turning line into list
        particle_coords[int(data[id_ind])].append([]) # adding list to the list of this particle id
        particle_coords[int(data[id_ind])][-1].append(float(data[x_ind].strip())) # adding x coord
        particle_coords[int(data[id_ind])][-1].append(float(data[y_ind].strip())) # adding y coord
        particle_coords[int(data[id_ind])][-1].append(float(data[z_ind].strip())) # adding z coord

    del lines # no longer needed so it is deleted to save memory
    
    if not Kwargs["use_existing_figure"]:
        # creating figure using the size detected in the file
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=Kwargs["figsize"])
        
    if Kwargs["add_bounding_box"]:
        #print("adding bounding box")
        _add_boundary_box(x_lims, y_lims, z_lims, line_width=Kwargs["line_width"])

    # # iterating through the particle data
    for item in particle_coords:
        try:
                # getting the color for the particle from the dict
            color = colors[list(colors.keys())[item-1]]
        except KeyError:
            # defaults to red (if the dict does not contain particle info)
            color = "red"
        
        # scatter plotting the particles with id "item" with the suface
        mlab.points3d(
            # plotting x, y, and z
            [i[0] for i in particle_coords[item]], [i[1] for i in particle_coords[item]], [i[2] for i in particle_coords[item]],
            scale_factor = Kwargs["particle_size"],
            color=mpl.colors.to_rgb(color)
        )
    
    if not Kwargs["use_existing_figure"]:
        if Kwargs["view"] is not None: mlab.view(*Kwargs["view"])
        
        if output_file is None:
            mlab.show()
        else:
            mlab.savefig(output_file, Kwargs["image_size"])
            mlab.close()


def surf_temperature_3d(
        surf_temperature_file:str, # file to read data from
        surf_file:str,
        output_file:str              = None, # file to output to
        temperature_data_index:int   = 1,
        **kwargs
    ) -> str|None:

    Kwargs = _TOTAL_KWARGS_3D.copy()
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)
    
    if not Kwargs["use_existing_figure"]:
        # creating figure using the size detected in the file
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=Kwargs["figsize"])

    surf_coords:list[list[float]] = [] # x,y,z coords for surface
    tris:list[list[float]] = [] # triangles for surface

    data:dict = {}

    # loads lines from surface file
    lines = []
    with open(surf_file, 'r') as f: lines = [line.strip().split() for line in f.readlines()]

    surfs_start:int = 0 # where the surfs start
    tris_start:int = 0 # where the triangle data starts
    x_ind:int = 1 # x-coord index in the file
    y_ind:int = 2 # y-coord index in the file
    z_ind:int = 3  # z-coord index in the file
    for i, line in enumerate(lines):
        new_line = line # splitting the line into a list
        if len(new_line): # if there is data in the line
            if new_line[0] == "Points":
                # where the surface data starts
                surfs_start = i+2
                
            if new_line[0] == "Triangles":
                tris_start = i+2
                break # no longer need any info so exits loop

    # iterates through the loaded lines
    for i in range(tris_start, len(lines)):
        data = lines[i] # splits lines into list
        # if the length of the line is 0 then there is no more data, so exit loop
        if not len(data): break
        # adds list to list (will contain surface order pairs)
        tris.append([int(data[1].strip())-1, int(data[2].strip())-1, int(data[3].strip())-1]) # grabs verticies of triangle


    # iterates through the loaded lines
    for i in range(surfs_start, len(lines)):
        data = lines[i] # splits lines into list
        # if the length of the line is 0 then there is no more data, so exit loop
        if not len(data): break
        # adds list to list (will contain surface order pairs)
        surf_coords.append([float(data[x_ind].strip()), float(data[y_ind].strip()), float(data[z_ind].strip())]) # grabs x,y,z coord of surface point

    del lines # deletes, no longer need var

    # getting the surface data from the lists
    x = [i[0] for i in surf_coords]
    y = [i[1] for i in surf_coords]
    z = [i[2] for i in surf_coords]

    # adding the first coords to close the surface
    x.append(surf_coords[0][0])
    y.append(surf_coords[0][1])
    z.append(surf_coords[0][2])
    
    with open(surf_temperature_file, 'r') as f: lines = [line.strip().split() for line in f.readlines()]

    surfs_start = 0
    for i, line in enumerate(lines):
        if len(line)>=2:
            if line[1] == "SURFS":
                # where the surface data starts
                surfs_start = i+1
                break # no longer need anything so exit the loop

    t = []
    for i, line in enumerate(lines[surfs_start:]):
        t.append(float(line[temperature_data_index]))
    
    mesh = mlab.triangular_mesh(x, y, z, tris)
    mesh.mlab_source.dataset.cell_data.scalars = t
    mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data' 
    mesh.mlab_source.update() 
    mesh.parent.update()

    mesh2 = mlab.pipeline.set_active_attribute(mesh, cell_scalars='Cell data')
    s2 = mlab.pipeline.surface(mesh2) 

    mlab.scalarbar(s2)
    
    if not Kwargs["use_existing_figure"]:
        if Kwargs["view"] is not None: mlab.view(*Kwargs["view"])
        
        if output_file is None:
            mlab.show()
        else:
            mlab.savefig(output_file, Kwargs["image_size"])
            mlab.close()


def sim_3d(
        particle_file:str            = None, # file to read data from
        surf_file:str                = None,
        grid_file:str                = None,
        output_file:str              = None, # file to output to
        colors:OrderedDict[str, str] = {}, # structure = {species name : "color"}
        **kwargs
    ) -> str|None:
    # kwargs for this function
    
    if "add_bounding_box" not in kwargs:
        kwargs["add_bounding_box"] = grid_file is None
    
    Kwargs = _TOTAL_KWARGS_3D.copy()
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)

    Kwargs["use_existing_figure"] = True
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=Kwargs["figsize"])
    
    #print(Kwargs)
    # plotting each component if they are provided
    if grid_file is not None: grid_3d(grid_file, **Kwargs)
    if surf_file is not None: surf_3d(surf_file, **Kwargs)
    if particle_file is not None: particle_3d(particle_file, colors=colors, **Kwargs)

    if Kwargs["view"] is not None:
        mlab.view(*Kwargs["view"])

    if output_file is None:
        mlab.show()
    else:
        #print("saving to:", output_file)
        mlab.savefig(output_file, Kwargs["image_size"])
        mlab.close(mlab.gcf())
        return output_file


# makes a video from the output of a 2d simulation
def sim_3d_video(
        particle_files:list[str]|str    = None,
        surf_files:list[str]|str        = None,
        grid_files:list[str]|str        = None,
        frame_files:list[str]|str       = None,
        sort_key_particles:FunctionType = lambda x:int(x[x.rfind("/")+1:x.rfind(".")]),
        sort_key_surfs:FunctionType     = lambda x:int(x[x.rfind("/")+1:x.rfind(".")]),
        sort_key_grids:FunctionType     = lambda x:int(x[x.rfind("/")+1:x.rfind(".")]),
        sort_key_frames:FunctionType    = lambda x:int(x[x.rfind("/")+1:x.rfind(".")]),
        colors:dict[str, str]           = {}, # format = {molecule id : color in the plot}
        output_file:str                 = "out.mp4", # path of output file
        frame_rate:int                  = 10, # frames per second
        num_cores:int                   = multiprocessing.cpu_count()//4,
        **kwargs
    ):
    mlab.options.offscreen = True
    Kwargs = _TOTAL_KWARGS_3D.copy()
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)
    
    imgs = []
    
    # if no frame files were provided
    if frame_files is None:
        # using the above function to make appropriate file lists
        files, num_files, should_sort = handle_file_directive(
            particle_files=particle_files,
            surf_files=surf_files,
            grid_files=grid_files
        )

        # there should be something to do
        if not len(files): error("Nothing to make")
        
        # if particle files were provided, sorting them if they need to be sorted
        # if they were not provided, then just making a list of None the same size as the other file lists
        output_path = ""
        if particle_files is not None:
            if should_sort["particle_files"]: particle_files = sorted(files["particle_files"], key=sort_key_particles)
            else:                             particle_files = files["particle_files"]
            if not len(output_path):
                output_path = particle_files[0][:particle_files[0].rfind(os.path.sep)]
        else:                                 particle_files = [None]*num_files

        # if surf files were provided, sorting them if they need to be sorted
        # if they were not provided, then just making a list of None the same size as the other file lists
        if surf_files is not None:
            if should_sort["surf_files"]: surf_files = sorted(files["surf_files"], key=sort_key_surfs)
            else:                         surf_files = files["surf_files"]
            if not len(output_path):
                output_path = surf_files[0][:surf_files[0].rfind(os.path.sep)]
        else:                             surf_files = [None]*num_files

        # if grid files were provided, sorting them if they need to be sorted
        # if they were not provided, then just making a list of None the same size as the other file lists
        if grid_files is not None:
            if should_sort["grid_files"]: grid_files = sorted(files["grid_files"], key=sort_key_grids)
            else:                         grid_files = files["grid_files"]
            
            if not len(output_path):
                output_path = grid_files[0][:grid_files[0].rfind(os.path.sep)]
        else:                             grid_files = [None]*num_files
        
        '''
        this is the multithreading code
        
        '''

        # creating queues for threads
        q_in = multiprocessing.Queue() # queue to the threads
        q_out = multiprocessing.Queue() # queue from the threads
        
        # adding the data to the queue in the order below (must be in this order currently)
        for i in range(num_files): 
            q_in.put([
                i, {
                    "particle_file" : particle_files[i],
                    "surf_file"     : surf_files[i],
                    "grid_file"     : grid_files[i],
                    "output_file"   : os.path.join(output_path, f"{i}.png"),
                    "colors"        : colors
                },
                Kwargs
            ])

        # creates a process manages
        m = multiprocessing.Manager()
        procs = m.dict() # dict of process
        stop = multiprocessing.Value('i', 0) # creating shared mem integer
        try:
            processes = []
            # creates number of processes corresponding inputted number of cores to run on
            for i in range(num_cores-1): # -1 is because the main thread is also used
                # creating a process
                processes.append(multiprocessing.Process(target=thread_func, args=(sim_3d, q_in, q_out, stop)))
                # starting the above process
                processes[-1].start()
            
            # runs the thread function on the main thread
            thread_func(sim_3d, q_in, q_out, stop)
            
            for i in range(len(processes)): processes[i].join() # waits for the process to end
        
        # catches keyboard interrupt 
        except Exception as e:
            print("caught error", e, "shutting down cores")
            # waits for process to stop
            while True:
                stop.value = 1
                pids = [pid for pid, running in procs.items() if running]
                if not len(pids): break
                print('running jobs:', pids)
                time.sleep(1)
            # exits code
            exit()
        
        # getting the files that the threads created and their ids from the output queue
        _to_sort = []
        while not q_out.empty(): _to_sort.append(q_out.get())
        
        # sorting these files based on their id
        _to_sort = sorted(_to_sort, key=lambda item: item[0])
        
        # getting just the image path from the sorted list
        imgs = [item[1] for item in _to_sort]

    else:
        # using the above function to make appropriate file lists
        files, num_files, should_sort = handle_file_directive(frame_files=frame_files)

        # sorting the frames list if the above function says it needs it
        if should_sort["frame_files"]: frame_files = sorted(files["frame_files"], key=sort_key_frames)
        else:                          frame_files = files["frame_files"]

        # creates list of images if the inputted list if they are images
        for file in frame_files: imgs.append(file)
    
    # creates file for ffmpeg to read info from
    ffmpeg_file = os.path.join(os.getcwd(), ".sparta_temp")
    with open(ffmpeg_file, 'w') as f:
        for file in imgs: f.write(f"file \'{file}\'\n")
    del imgs # deletes no longer need object (frees mem)
    
    # running ffmpeg
    ffmpeg.input(
        ffmpeg_file, # file to load all of the files from
        r=str(frame_rate),       # framerate
        f='concat',  # do not know
        safe='0'     # do not know
    ).output(output_file, vcodec='libx264').run()
    
    # deletes ffmpeg file that is no longer needed
    os.remove(ffmpeg_file)
    
    mlab.options.offscreen = False