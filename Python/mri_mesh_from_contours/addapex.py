#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:58:49 2013
Needs: 

Takes the endo and epi meshes created by the other scripts and adds the apex
Assumes that the apex is in the end of the point list of the mesh
Ignores any scalars stored with the mesh

Again assumes that the slices are along X axis. But this code can be made universal by fitting a plane to the last 
slice and introducing a new coordinate system. Something to do one day.

@author: Constantine Butakoff
"""
import numpy as np
import mytools
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
#from mayavi import mlab
#import matplotlib.pyplot as plt
#from mayavi.mlab import triangular_mesh
import vtk




#
#  Parameters, adapt to your taste
#

#the cap generation is not particularly beautiful, so to avoid creating vertices very close to the
#old vertices, set the minimum distance for the new ones. If your new vertices are too far, reduce this number
mindist = 0.5;

#apply loop subdivision to make the mesh nicer
#the way I generate the cap, you loose point correspondence for the cap anyway
#but with this at least the whole mesh will look nicer. Be careful since this also smoothes the mesh a 
#little. I put only 1 subdivision.
make_beautiful = True; 

#input/output meshes
input_filename_epi = 'epi_mesh.vtk'
output_filename_epi = 'epi_mesh_capped.vtk'
input_filename_endo = 'endo_mesh.vtk'
output_filename_endo = 'endo_mesh_capped.vtk'



#load the number of samples per slice
npzfile = np.load('points_final.npz')
sampling_endo = npzfile['sampling_endo']
sampling_epi = npzfile['sampling_epi']

#------------------------------------------------------------
#
#
#
#
#     Functions
#
#
#
#------------------------------------------------------------
    


def GenerateApex(input_filename, output_filename, sampling):
    """
    Generate apex for the mesh. Assumes that the apex is in the end of the point list of the mesh
    
    input_filename string: vtk mesh to add apex to
    output_filename string: vtk file to store the result
    sampling int: number of vertices per layer

    Return nothing
    """        
    #input_filename = 'endo_mesh.vtk'
    #output_filename = 'endo_mesh_capped.vtk'
    #sampling = sampling_epi
    
    #read the mesh
    rdr = vtk.vtkPolyDataReader()
    rdr.SetFileName(input_filename)
    rdr.Update()
    mesh = rdr.GetOutput()
    meshpts = mytools.ExtractVTKPoints(mesh)
    meshfaces = mytools.ExtractVTKTriFaces(mesh)
    m, n = meshpts.shape
    
    #get the last 4 layers (after resampling they do not coincide with the MRI)
    pts = meshpts[-4*sampling:,:]
    #mlab.points3d(pts[:,0],pts[:,1],pts[:,2], scale_mode='none', scale_factor=0.2 )
    
    
    #assuming slices are along X axis
    #project points into YZ plane and use griddata to interpolate the third coordinate
    miny = pts[:,1].min()
    maxy = pts[:,1].max()
    minz = pts[:,2].min()
    maxz = pts[:,2].max() 
    
    gx, gy = np.mgrid[miny:maxy:10j, minz:maxz:10j]
    grid = griddata(pts[:,1:3],pts[:,0],(gx,gy),method = 'cubic')
    #finished interpolating
    
    #get just the last slice for triangulation. We will add these to the new points and run triangulation
    pts_last = meshpts[-sampling:,:]

    #we also need to know if the coordinate is increasing or decreasing to keep only 
    #those new points, that lie beyond the original mesh. Don't want any overlaps
    slice_coord_first = meshpts[0,0]
    slice_coord_last = pts_last[0,0]
    
    slices_increase = False
    if slice_coord_last > slice_coord_first:
        slices_increase = True
        
    #create a Nx3 array of grid points where we interpolated the x coordinate    
    gridpts = np.column_stack( (grid.flatten(),gx.flatten(),gy.flatten()) )
                              
    #leave only points beyond the current mesh
    if slices_increase:
        cover_pts = gridpts[gridpts[:,0]>slice_coord_last+mindist,:]
    else :
        cover_pts = gridpts[gridpts[:,0]<slice_coord_last-mindist,:]
        
    #combine generated points with the last slice and tirangulate
    tri = Delaunay( np.concatenate( (pts_last[:,1:3], cover_pts[:,1:3]) ) )
    mergedpts = np.concatenate( (meshpts, cover_pts) )
    
    capfaces = tri.simplices.copy()
    
    #eliminate faces that have only boundary points 0:sampling
    #sometimes the edges in the slice are curvy and there are triangles connecting 
    #just the boundary points. We don't want to have them in the mesh, they will ruin everything.
    n1 = capfaces<sampling
    n2 = n1.all(axis=1)
    capfaces1 = capfaces[ n2==False, : ]
    
    #the new points will be added to the end of the mesh points
    #so we need to update the faces to account for the new point ids.
    #Also we need to take into account that the first "sampling" points of the 
    #generated apex are the last points of the mesh
    capfaces1 = capfaces1 + m - sampling
    allfaces = np.concatenate( (meshfaces, capfaces1)) #concatenate faces
    
    
    #create polydata and save it
    pd = mytools.CreatePolyData(mergedpts,allfaces)
    
    
    if make_beautiful:
        ls = vtk.vtkLoopSubdivisionFilter()
        ls.SetInput(pd)
        ls.SetNumberOfSubdivisions(1)
        ls.Update()
        mytools.SaveVTKPolyData(ls.GetOutput(),output_filename)
    else:
        mytools.SaveVTKPolyData(pd,output_filename)
        
        
        
#------------------------------------------------------------
#
#
#
#
#     The code
#
#
#
#------------------------------------------------------------
            
             
        
GenerateApex(input_filename_endo, output_filename_endo, sampling_endo)
GenerateApex(input_filename_epi, output_filename_epi, sampling_epi)
