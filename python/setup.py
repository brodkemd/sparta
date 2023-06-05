from setuptools import setup
from src.sparta import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

### to install mayavi, run pip install git+https://github.com/enthought/mayavi.git

import src.sparta._make_sparta
import src.sparta._make_tools_wrapper
import src.sparta._make_paraview_wrapper

s = setup(
    install_requires=[
        'matplotlib',
        "numpy",
        "ffmpeg",
        "PyQt5"
    ],
    name='sparta',
    version=__version__,
    packages=["sparta"],
    author='Marek Brodke, with support from the University of Cincinnati',
    description='Provides functionaliy for sparta DSMC',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='development',
    python_requires='>=3.8',
    license="MIT",
    author_email="brodkemd@mail.uc.edu",
    url="https://github.com/brodkemd/sparta",
    package_dir={"sparta" : "src/sparta"},
    package_data={"sparta" : ["doc/**", "pizza/**", "paraview/**", "db.txt", "tools/**", "elmer_doc/**"]},
)

#installation_path = s.command_obj['install'].install_lib