# Introduction and licensing
Everything in this repository is licensed under the Attribution-NonCommercial-ShareAlike 4.0 International License.

The collection is in a bit of a disarray and documentation need to be improved. 


#Contents of subfolders


##VolumetricLVMesh 
Generates a volumetric mesh of cardiac Left Ventricle (half-ellipsoid like). 
The elements are wedges aligned radially. Might have some additional requirements.
To be verified. This was ported from VTK 5.x and needs testing. An example mesh with associate local coordinates can be seen in the following figure.

Please cite B. Paun, B. Bijnens, T. Iles, P. A. Iaizzo, C. Butakoff, "Patient Independent Representation of the Detailed
Cardiac Ventricular Anatomy", Medical Image Analysis, Accepted 2016.

![Volumetric Meshing Image](https://raw.githubusercontent.com/cbutakoff/tools/master/Pics/volmeshbump.gif)

Todo:
- Test

##layer_Generation_code
Generates a volumetric mesh of any structure with two boundaries. 
The elements are wedges aligned radially. 
This was ported from VTK 5.x and needs testing. 

Please cite B. Paun, B. Bijnens, T. Iles, P. A. Iaizzo, C. Butakoff, "Patient Independent Representation of the Detailed
Cardiac Ventricular Anatomy", Medical Image Analysis, Accepted 2016.

![Volumetric Meshing Image](https://raw.githubusercontent.com/cbutakoff/tools/master/Pics/layergeneration.gif)

Todo:
- Test


