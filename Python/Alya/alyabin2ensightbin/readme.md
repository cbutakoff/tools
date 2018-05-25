# Running on marenostrum

module load python/3.6.1
mpirun -n 30 python alya2ensight-mpi.py case_name case_path output_path

# You might need to run it on a highmem dedicated node if the mesh is large
