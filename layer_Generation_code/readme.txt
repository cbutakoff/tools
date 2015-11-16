To run the code, you need to have admesh installed and in the path
also for windows, find 
	sprintf(fixcmd,"admesh 
and replace with 
	sprintf(fixcmd,"admesh.exe 

Alternatively you can set the variable targetarea=-1 to disable the usage of admesh, but it will remove uniform remeshing.

You need 
ITK not really used internally, i think
VTK important
VMTK for uniform remeshing
admesh (hopefully) https://sites.google.com/a/varlog.com/www/admesh-htm



Usage:
GenerateAtrialLayers -endo endo.vtk -epi epi.vtk -layers N_layers -o output_prefix -nfaces -1 -voxel_mc 1 -voxel_lap 0.2 -iter 1500 -smooth_rad 3 

Comments:
-endo --- mask image (short): 0 - background, nonzero - voxels insode endocardium
-epi --- mask image (short, same size as endo): 0 - background, nonzero - voxels inside epicardium
-layers --- number of layers to generate, equally spaced
-o output_prefix --- the prefix for the files to generate layers. Can include path. The number of the layer will be appended in the end
-nfaces --- leave it at -1
-voxel_mc --- voxel size for marching cubes. The bigger it is the smoother the mesh, but also less detailed. 1 is good
-voxel_lap --- voxel size for solving laplace equation. The smaller the higher the accuracy of the algorithm, but the number of iterations has to be increased, the algorithm will take longer and the memory requirements scale exponentially. 0.2 is good for atria. 
-iter --- number of iterations for solving the laplace equation, 1500 for voxel size 0.2 is ok
-smooth_rad --- standard deviation for the gaussian smoothing, 3 is ok. 


For every generated mesh layer smooth it using the following:
SmoothMeshThroughImage -i layer.vtk -o layer.vtk -v 0.5 -f -1 -s 3

-i --- input mesh
-o --- output mesh
-v --- voxel size for smoothing
-f --- leave at -1
-s --- standard deviation for the gaussian smoothing
