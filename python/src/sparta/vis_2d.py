import os, ffmpeg, glob, multiprocessing, time, traceback
from collections import OrderedDict
import matplotlib.pyplot as plt
from types import FunctionType
from queue import Empty
import numpy as np

from ._src import error
from ._vis_src import _handle_kwargs, handle_file_directive, thread_func

'''
Must modify:

File "/home/marekbrodke/anaconda3/envs/qiskit/lib/python3.10/site-packages/matplotlib/collections.py", in function named "set_sizes"

with the following (replace the code with this) ->

    if sizes is None:
            self._sizes = np.array([])
            self._transforms = np.empty((0, 3, 3))
        else:
            \'\'\'
            My modifications
            \'\'\'

            if isinstance(sizes, int) or isinstance(sizes, float):
                self._sizes = np.array([sizes])
            else:
                self._sizes = np.asarray(sizes)

            self._transforms = np.zeros((len(self._sizes), 3, 3))
            scale = np.sqrt(self._sizes) * dpi / 72.0 * self._factor
            self._transforms[:, 0, 0] = scale
            self._transforms[:, 1, 1] = scale
            self._transforms[:, 2, 2] = 1.0
        self.stale = True

'''

# kwargs used as reference for the functions in this file
_TOTAL_KWARGS_2D = {
    # arg                      default value
    "surf_plot_style"        : "b-",
    "surf_plot_args"         : dict(),
    "particle_edge_color"    : "black",
    "include_particle_every" : 1,
    "particle_line_width"    : 0.1,
    "legend_args"            : dict(loc='upper right', ncol=1, fontsize=16),
    "legend_particle_size"   : 80,
    "particle_size"          : 4,
    "figsize"                : (5, 5)
}


def grid_2d(grid_file:str, output_file:str=None, ax:plt.Axes=None, **kwargs) -> None:
    # kwargs for this function
    Kwargs = {"figsize" : (5, 5)}
    Kwargs = _handle_kwargs(kwargs, Kwargs)

    # reading the file
    with open(grid_file, 'r') as f: lines = [line.strip().split() for line in f.readlines()]

    #ids = []
    # centers = []
    bounds = []
    
    # loading the data from the lines
    for i, line in enumerate(lines):
        # if the line is not empty
        if len(line) > 1:
            # getting the box dimensions
            if line[1] == "BOX":
                x_lims = [float(item) for item in lines[i+1].strip().split()]
                y_lims = [float(item) for item in lines[i+2].strip().split()]
            
            # loading the grid cells
            elif line[1] == "CELLS":
                # getting the indices of box parameters
                _x_lo_ind = line.index("xlo") - 2
                _y_lo_ind = line.index("xlo") - 2
                _x_hi_ind = line.index("xhi") - 2
                _y_hi_ind = line.index("yhi") - 2
                
                # loading the data from the rest of the lines
                for j in range(i+1, len(lines)):
                    #ids.append(lines[j][0])
                    # centers.append([float(lines[j][-2]), float(lines[j][-1])])
                    # getting the bounds
                    bounds.append([
                        [float(lines[j][_x_lo_ind]), float(lines[j][_x_hi_ind])],
                        [float(lines[j][_y_lo_ind]), float(lines[j][_y_hi_ind])]
                    ])
                break
    
    
    #centers = np.array(centers)
    
    # if an axis is provided then one is not made
    provided_ax = True
    if ax is None:
        provided_ax = False
    
        fig = plt.figure(figsize=Kwargs["figsize"])
        ax = fig.add_subplot(111)
    
    # setting the axis limits (the grid determines the geometry so this is ok)
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    
    # used to make the grid lines
    num_points = 2
    template = np.ones(num_points)

    for box in bounds:
        # plotting only the inner grid lines in the x direction
        if box[0][0] != x_lims[0]:
            ax.plot(box[0][0]*template, np.linspace(box[1][0], box[1][1], num_points), 'k-', linewidth = ax.spines["left"].get_linewidth())
        
        # plotting only the inner grid lines in the y direction
        if box[1][0] != y_lims[0]:
            ax.plot(np.linspace(box[0][0], box[0][1], num_points), box[1][0]*template, 'k-', linewidth = ax.spines["left"].get_linewidth())

    #ax.plot(centers[:,0], centers[:, 1], 'r.', markersize=1)

    # finishing touches on the figure if no ax was provided
    if not provided_ax:
        ax.set_aspect(1)

        # removing ticks from axes
        ax.set_xticks([])
        ax.set_yticks([]) 
        
        fig.tight_layout()

        if output_file is not None:
            # creating an image from figure
            fig.savefig(output_file)
            # deleting the figure to free up memory
            plt.close(fig)
        else:
            plt.show()
            
# class for 2d particle image
def surf_2d(
        surf_file:str, # file to read the surface from
        surf_file_is_dump:bool  = False, # if the surface is from a dump
        output_file:str         = None,
        ax:plt.Axes             = None,
        **kwargs
    ) -> None:
    # kwargs for this function
    Kwargs = {
        # arg                      default value
        "surf_plot_style"        : "b-",
        "surf_plot_args"         : dict(),
        "figsize"                : (5, 5)
    }

    """
    Kwargs options:
        - surf_plot_style               => style of the surface plot
        - surf_plot_args:dict[str, any] => arguements to pass to the plot of the surface function,
                                           this is any valid argument for plt.plot     
    """
    Kwargs = _handle_kwargs(kwargs, Kwargs)
        
    surf_coords:list[list[float]] = [] # x,y,z coords for surface

    # loads lines from surface file
    lines = []
    with open(surf_file, 'r') as f: lines = f.readlines()

    # if the surface is from a dump from sparta
    if surf_file_is_dump:
        surfs_start:int = 0 # where the surfs start
        x_ind:int = 0  # x-coord index in the file
        y_ind:int = 0  # y-coord index in the file
        for i, line in enumerate(lines):
            new_line = line.strip().split() # splits line into list
            if len(new_line)>=2:
                if new_line[1] == "SURFS":
                    # where the surface data starts
                    surfs_start = i+1
                    for j, item in enumerate(new_line):
                        if item == "v1x":   x_ind = j - 2 # the x coord of the surface
                        elif item == "v1y": y_ind = j - 2 # the y coord of the surface
                    break # no longer need anything so exit the loop
    # if the surface is NOT from a dump from sparta
    # i.e. it is an input file into sparta
    else:
        surfs_start:int = 0 # where the surfs start
        x_ind:int = 1 # x-coord index in the file
        y_ind:int = 2 # y-coord index in the file
        for i, line in enumerate(lines):
            new_line = line.strip().split() # splitting the line into a list
            if len(new_line): # if there is data in the line
                if new_line[0] == "Points":
                    # where the surface data starts
                    surfs_start = i+2
                    break # no longer need any info so exits loop
    
    # iterates through the loaded lines
    for i in range(surfs_start, len(lines)):
        data = lines[i].strip().split() # splits lines into list
        # if the length of the line is 0 then there is no more data, so exit loop
        if not len(data): break
        # adds list to list (will contain surface order pairs)
        surf_coords.append([])
        surf_coords[-1].append(float(data[x_ind].strip())) # grabs x coord of surface point
        surf_coords[-1].append(float(data[y_ind].strip())) # grabs y coord of surface point

    del lines # deletes no longer need var

    # if an axis is provided then one is not made
    provided_ax = True
    if ax is None:
        provided_ax = False
        # creating figure using the size detected in the file
        fig = plt.figure(figsize=Kwargs["figsize"])#plt.figure(figsize=(max(x_lims)-min(x_lims), max(y_lims)-min(y_lims)))
        ax = fig.add_subplot(111) # creating axis

        # removing ticks from axes
        ax.set_xticks([])
        ax.set_yticks([])

    # getting the surface data from the lists
    x = [i[0] for i in surf_coords]
    y = [i[1] for i in surf_coords]

    # adding the first coords to close the surface
    x.append(surf_coords[0][0])
    y.append(surf_coords[0][1])
    
    # plotting the surface, using the inputted style and the inputted arguements for it
    ax.plot(x, y, Kwargs["surf_plot_style"], **Kwargs["surf_plot_args"])

    if not provided_ax:
        ax.set_aspect(1)

        # setting the figure to have a tight layout
        fig.tight_layout()
        
        if output_file is not None:
            # creating an image from figure
            fig.savefig(output_file)
            # deleting the figure to free up memory
            plt.close(fig)
        else:
            plt.show()


def particle_2d(
        particle_file:str, # file to read data from
        output_file:str              = None, # file to output video to
        colors:OrderedDict[str, str] = {}, # structure = {species name : "color"}
        ax:plt.Axes                  = None,
        **kwargs
    ) -> None:
    # kwargs for this function
    Kwargs = {
        # arg                      default value
        "particle_edge_color"    : "black",
        "include_particle_every" : 1,
        "particle_line_width"    : 0.1,
        "legend_args"            : dict(loc='upper right', ncol=1, fontsize=16),
        "legend_particle_size"   : 80,
        "particle_size"          : 4,
        "figsize"                : (5, 5)
    }

    """
    Kwargs options:
        - particle_edge_color:str       => color for the edge of the particle
        - include_particle_every:int    => only include a particle after skipping over this many (used to reduce render time)
        - particle_line_width:float     => width of the border of the particles
        - legend_args:dict[str, any]    => args to pass to the plt.legend command in matplotlib
        - legend_particle_size:float    => size of the particle in the legend
        - particle_size:float           => size of the particle in the plot     
    """
    Kwargs = _handle_kwargs(kwargs, Kwargs)

    # reading the lines from the input file
    lines = []
    with open(particle_file, 'r') as f: lines = f.readlines()
    
    atoms_start:int = 0 # where the atoms start
    x_ind:int = 0  # x-coord index in the file
    y_ind:int = 0  # y-coord index in the file
    id_ind:int = 0 # species id of the particle
    x_lims:list[float] = [] # upper and lower bound for x
    y_lims:list[float] = [] # upper and lower bound for y

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
                    elif item == "type": id_ind = j - 2
                break # last thing that we are looking for so breaks loop
            elif new_line[1] == "BOX":
                # gets the bounds of the simulation box
                x_lims = [float(item) for item in lines[i+1].strip().split()]
                y_lims = [float(item) for item in lines[i+2].strip().split()]
    
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

    del lines # no longer needed so it is deleted to save memory


    # if an axis is provided then one is not made
    provided_ax = True
    if ax is None:
        provided_ax = False
        # creating figure using the size detected in the file
        fig = plt.figure(figsize=Kwargs["figsize"])#plt.figure(figsize=(max(x_lims)-min(x_lims), max(y_lims)-min(y_lims)))
        ax = fig.add_subplot(111) # creating axis
        
        # removing ticks from axes
        ax.set_xticks([])
        ax.set_yticks([])

    # setting axes limits from detected bounds
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)

    s=[]
    labels = []

    # iterating through the particle data
    for item in particle_coords:
        try:
            # getting the color for the particle from the dict
            color = colors[list(colors.keys())[item-1]]
        except KeyError:
            # defaults to red (if the dict does not contain particle info)
            color = "red"
        
        # scatter plotting the particles with id "item" with the suface
        scatter = ax.scatter(
            # plotting x and y
            [i[0] for i in particle_coords[item]], 
            [i[1] for i in particle_coords[item]],
            s = Kwargs["particle_size"], # the marker size
            c=color, # the color of the marker
            edgecolors=Kwargs["particle_edge_color"], # edge color of the marker
            linewidths=Kwargs["particle_line_width"]  # width of the border line of the marker
        )
        # adding plot to list, to be used later by legend
        s.append(scatter)
        # adding the label to the list
        labels.append(list(colors.keys())[item-1])

    # adding a legend to the plot
    legend = ax.legend(
        s, # the scatter plot objects
        labels, # the labels for the respective scatter plots
        **Kwargs["legend_args"] # adding args to legend from input
    )

    # sets the size for the markers in the legend
    for i in range(len(legend.legendHandles)):
        legend.legendHandles[i]._sizes = Kwargs["legend_particle_size"]

    if not provided_ax:
        ax.set_aspect(1)

        # setting the figure to have a tight layout
        fig.tight_layout()
        
        if output_file is not None:
            # creating an image from figure
            fig.savefig(output_file)
            # deleting the figure to free up memory
            plt.close(fig)
        else:
            plt.show()


# class for 2d particle image
def sim_2d(
        particle_file:str            = None, # file to read data from
        surf_file:str                = None,
        grid_file:str                = None,
        output_file:str              = None, # file to output video to
        colors:OrderedDict[str, str] = {}, # structure = {species name : "color"}
        **kwargs
    ):
    # kwargs for this function
    Kwargs = _TOTAL_KWARGS_2D.copy()

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
    """
    Kwargs = _handle_kwargs(kwargs, Kwargs, strict=True)

    fig = plt.figure(figsize=Kwargs["figsize"])# figsize=(20, 20) plt.figure(figsize=(max(x_lims)-min(x_lims), max(y_lims)-min(y_lims)))
    ax = fig.add_subplot(111) # creating axis
    
    # removing ticks from axes
    ax.set_xticks([])
    ax.set_yticks([])
    
    # plotting each component if they are provided
    if grid_file is not None: grid_2d(grid_file, ax=ax, **Kwargs)
    if surf_file is not None: surf_2d(surf_file, ax=ax, **Kwargs)
    if particle_file is not None: particle_2d(particle_file, colors=colors, ax=ax, **Kwargs)
    
    ax.set_aspect(1)

    # setting the figure to have a tight layout
    fig.tight_layout()
    
    if output_file is not None:
        # creating an image from figure
        fig.savefig(output_file)
        # deleting the figure to free up memory
        plt.close(fig)
    else:
        plt.show()

'''
Function to run on threads

'''
# def thread_func_2d(
#         queue_in:multiprocessing.Queue, # queue to get data from
#         queue_out:multiprocessing.Queue, # queue to get data from
#         stop,      # shared mem integer to indicate for this function to stop
#     ) -> None:
#     print(multiprocessing.current_process().name, "starting up")
    
#     # max time allowed for threads to wait for data
#     timeout_time = 2.
    
#     try:
#         last_update = time.time()
        
#         # runs until the queue is empty
#         while not queue_in.empty():
#             # trying to get data from the queue, if none is available then it checks to make sure
#             # the max time to wait for data has not been reached.
#             # If the max wait time is reached then this process exits
#             try:
#                 data = queue_in.get(block=False) # getting latest object from the queue
#             except Empty:
#                 if time.time() - last_update > timeout_time:
#                     print(multiprocessing.current_process().name, "caught timeout shutting down")
#                     break
#                 continue
            
#             # creating the output file name from the id of the data and determines a path to place it
#             # in from the other provided files in the data.
#             # if no file is found then the cwd is used (default for mpl)
#             out_file_name = f"{data[0]}.png"
#             for i in range(1, 4):
#                 if data[i] is not None:
#                     out_file_name = os.path.join(data[i][:data[i].rfind(os.path.sep)], out_file_name)
#                     break
#             else:
#                 print(f"Can not find path for data with id: {data[0]} (using local instead)")

            
#             print(multiprocessing.current_process().name, "making", f"{out_file_name}")
            
#             # making an image based on the file names from the queue
#             sim_2d(
#                 particle_file = data[1],
#                 surf_file     = data[2],
#                 grid_file     = data[3],
#                 output_file   = out_file_name,
#                 colors        = data[4],
#                 **data[5]
#             )
            
#             # setting the last time an action was perform to the current time (just made an image)
#             last_update = time.time()
            
#             # outputing the data
#             queue_out.put([data[0], out_file_name])

#             # if this function should stop (if the int is set to 1 from the main thread)
#             if stop.value:
#                 print(multiprocessing.current_process().name, " caught end signal shutting down")
#                 break

#         print(multiprocessing.current_process().name, "done")
    
#     # catches exception (primarily for keyboard interrupt)
#     except Exception as e:
#         print(multiprocessing.current_process().name, " caught error:", traceback.format_exc())


# makes a video from the output of a 2d simulation
def sim_2d_video(
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
        num_cores:int                   = multiprocessing.cpu_count()//2,
        **kwargs
    ):
    Kwargs = _TOTAL_KWARGS_2D.copy()
    Kwargs = _handle_kwargs(kwargs, Kwargs)
    
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
                processes.append(multiprocessing.Process(target=thread_func, args=(sim_2d, q_in, q_out, stop)))
                # starting the above process
                processes[-1].start()
            
            # runs the thread function on the main thread
            thread_func(sim_2d, q_in, q_out, stop)
            
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