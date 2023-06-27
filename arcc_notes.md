# Compiling
## Sparta
1. Allocate a job:
    ```shell
    srun --nodes 1 --exclusive --time 010:00:00 --pty bash
    ```

2. Load the required modules:
    ```shell
    module unload autotools prun/1.3 ohpc singularity/3.7.1
    module load gnu/8.3.0 openmpi/4.1.4 anaconda/3.0
    ```

3. Activate the python environment to get access to python:
    ```shell
    conda activate sparta
    ```

4. Clean and build:
    ```shell
    make clean-all
    make mpi -j 30
    ```

5. Load the library for the executable:
    ```shell
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.conda/envs/sparta/lib/
    ```

## Elmer
1. Cloning the files:
    ```shell
    mkdir elmer
    cd elmer
    git clone https://github.com/ElmerCSC/elmerfem.git
    ```

2. Allocating resources:
    ```shell
    srun --nodes 1 --exclusive --time 010:00:00 --pty bash
    module load gnu/12.1.0 openmpi/4.1.4
    ```


3. Building:
    ```shell
    mkdir build install
    cmake -DWITH_ELMERGUI:BOOL=FALSE -DWITH_MPI:BOOL=FALSE -DCMAKE_C_COMPILER="gcc" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../install ../elmerfem
    ```