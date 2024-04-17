# Compiling
## Sparta
### Prerequisites
Must create an environment for anaconda called ```sparta```. The command below can make the environment.
```shell
conda create -n sparta python=3
```
### Directions
1. Allocate a job, to compile the code:
    ```shell
    srun --nodes 1 --exclusive --time 010:00:00 --pty bash
    ```

2. Download the code:
    ```shell
    git clone https://github.com/brodkemd/sparta.git; cd sparta/src;
    ```

3. Load the required modules:
    ```shell
    module purge; module load gnu/8.3.0 openmpi/4.1.4 anaconda/3.0; module list;
    ```

4. Activate the python environment to get access to python:
    ```shell
    conda activate sparta
    ```

5. Clean and build:
    ```shell
    make clean-all; make mpi -j 30;
    ```

6. Load the library for the executable, this must be done every time it is run:
    ```shell
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/anaconda3/envs/sparta/lib
    ```

7. Install the python library:
    ```shell
    cd ../python; pip install .
    ```

8. Relinquish resources and go back to the top directory:
    ```shell
    conda deactivate; exit;
    ```

## Elmer
1. Allocating resources to compile the code:
    ```shell
    srun --nodes 1 --exclusive --time 010:00:00 --pty bash
    ```

2. Load modules:
    ```shell
    module purge; module load gnu/8.3.0 openmpi/4.1.4 cmake; module list
    ```

2. Cloning the files:
    ```shell
    mkdir elmer; cd elmer; git clone https://github.com/ElmerCSC/elmerfem.git
    ```

3. Building:
    - ```shell
      mkdir build install; cd build;
      ```
    
    - ```shell
      cmake -DWITH_ELMERGUI:BOOL=FALSE -DWITH_MPI:BOOL=FALSE -DCMAKE_C_COMPILER="gcc" -DCMAKE_CXX_COMPILER="g++" -DCMAKE_Fortran_COMPILER="gfortran" -DCMAKE_INSTALL_PREFIX=../install ../elmerfem
      ``` 
    - ```shell
      make install -j 30
      ```

4. Relinquish resources and go back to the top directory:
    ```shell
    exit
    ```

# Running
Here is a slurm batch script. It runs on 5 nodes with 20 tasks each nodes each task has 1 core, totalling in 100 tasks, for 36 hours. Replace `SPARTA_COMMAND` with the command to run sparta.
```shell
#!/bin/bash
#SBATCH --nodes=5
#SBATCH --ntasks=100
#SBATCH --tasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=36:00:00

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load gnu/8.3.0 openmpi/4.1.4 anaconda/3.0

# adds path to environment vairable to load the python shared library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/anaconda3/envs/sparta/lib/
export PMIX_MCA_gds=hash

# runs with mpi, '--mca mpi_warn_on_fork 0' allows system calls
mpirun --mca mpi_warn_on_fork 0 SPARTA_COMMAND
```