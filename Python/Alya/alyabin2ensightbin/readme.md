# Running on marenostrum
```
module load python/3.6.1
mpirun -n 30 python alya2ensight-mpi.py case_name case_path output_path 
```
# You might need to run it on a highmem dedicated node if the mesh is large
# Only runs with MPI. Minimum: mpirun -np 2
```
usage: alya2ensight-mpi.py [-h] [--format FORMAT]
                           task_name input_folder output_folder

positional arguments:
  task_name        Name of the alya task
  input_folder     Folder with input alyabins
  output_folder    Folder for the output ensight case

optional arguments:
  -h, --help       show this help message and exit
  --format FORMAT  Format of the data to expect: mpio, alyabin, auto(default)
```