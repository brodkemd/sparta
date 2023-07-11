import os, subprocess

cwd = os.path.dirname(__file__)

## converts the input_file (markdown) to html
def markdown2html(input_file:str, style:str="light"):
    global cwd
    # the file to read the format from
    format_file = os.path.join(cwd, "markdown_format.html")
    
    # creating the output file path from the input file path
    output_file = input_path[:input_path.rfind(".")+1] + "html"

    # not case sensitive
    style = style.lower()

    # valid styles and their corresponding css files
    style_files = {
        "normal" : "github-markdown.css",
        "light"  : "github-markdown-light.css",
        "dark"   : "github-markdown-dark.css"
    }

    # error checking
    if style not in style_files:
        print(f"Error: {style} is not a valid style, valid styles are:", ", ".join(list(style_files.keys())))
        exit(1)

    # running the command to compile markdown into html
    result = subprocess.run([
        "marked", "-i", input_file, "-o", output_file
    ], stdout=subprocess.PIPE)
    
    # if the command did not succeed
    if result.returncode:
        print("Error:", result.stderr)
        print("DONE")
        exit(result.returncode)
    
    # Open the file for reading and read the input to a temp variable
    # the previous command already wrote to the output file, so reading
    # from it
    with open(output_file, 'r') as f: md = f.read()

    # getting the format
    with open(format_file, 'r') as f: html_format = f.read()

    # the key is replaced by the value of the key
    replacers = {
        "####" : md,
        "&&&&" : style_files[style]
    }

    for _replace in replacers:
        html_format = html_format.replace(_replace, replacers[_replace])



    # If necessary, could print or edit the results at this point.
    # Open the HTML file and write the output.
    with open(output_file, 'w') as f: f.write(html_format)


for item in os.listdir(cwd):
    if item.endswith(".md"):
        input_path = os.path.join(cwd, item)
        markdown2html(input_path, style="Dark")