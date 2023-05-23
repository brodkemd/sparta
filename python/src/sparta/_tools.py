import numpy as np

def surf_create_2d(file:str, functions:list):
    # Examples input for functions:
    # funcs = [
    #     (lambda var: [None, 1/4*var], [0., 1.], 10),
    #     (lambda var: [1, None], [None, None], 10),
    #     (lambda var: [None, -1/4*var], [1., 0.], 10)
    # ]
    #
    # None in the lambda function results get replaced by the array generated from the range (the second element in the tuple) and the number of points (the last element in the tuple).
    # None in the range just allows this function to insert something that works
    
    # list to hold the points
    points = []
    count = 1
    
    # evaluating the functions
    for i, func in enumerate(functions):
        if func[1][0] is None:
            if functions[i-1][1][1] is not None:
                result = functions[i-1][0](functions[i-1][1][1])
                for val in result:
                    if val is not None:
                        functions[i][1][0] = val
                        break
                
        
        if func[1][1] is None:
            #functions[i][1][1] = functions[i+1][0](functions[i-1][1][0])
            if functions[i+1][1][0] is not None:
                result = functions[i+1][0](functions[i+1][1][0])
                for val in result:
                    if val is not None:
                        functions[i][1][1] = val
                        break
        
        func = functions[i]
        input_var = np.linspace(func[1][0], func[1][1], func[2])[:-1]
        output_var = func[0](input_var)
        
        if output_var[0] is not None:
            if not isinstance(output_var[0], np.ndarray):
                output_var[0] = output_var[0] * np.ones(len(input_var))
            output_var[1] = input_var
        else:
            if not isinstance(output_var[1], np.ndarray):
                output_var[1] = output_var[1] * np.ones(len(input_var))
            output_var[0] = input_var
        
        for i in range(len(output_var[0])):
            points.append([count, output_var[0][i], output_var[1][i]])
            count+=1
    
    # writing the info to the file
    with open(file, "w") as f:
        header = f"surf file from surf file from Pizza.py\n\n{len(points)} points\n{len(points)} lines\n\nPoints\n\n"
        f.write(header)
        for point in points:
            f.write(" ".join([str(val) for val in point]) + "\n")
        
        f.write("\nLines\n\n")
        for i in range(1, len(points)):
            f.write(f"{i} {i} {i+1}\n")
        f.write(f"{i+1} {i+1} {1}\n")