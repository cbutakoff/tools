#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:58:49 2013
Needs: VTK, Numpy, mayavi optionally(disable it in the code)

Creates the mesh from ladmanrked slices. If necessary it uses a reference point, which is a point that determines the starting point 
of the new landmarks. Basically if the reference point is chosen consistently across a population of
segmentations, then the reference point would would help maintain point correspondence.

@author: Constantine Butakoff
"""

import numpy as np
from scipy import interpolate
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from mayavi.mlab import triangular_mesh, points3d, figure, show
import mytools

#
#  Parameters, adapt to your taste
#


input_filename = 'points_corrected.npz'  #numpy file to read the corrections
output_endo_filename = 'endo_mesh.vtk'
output_epi_filename  = 'epi_mesh.vtk'
output_points_filename = 'points_final.npz'

remove_0_slice = True; #this just removes the 0 slice if you don't like it

sampling_endo = 40;  #number of points to put on endo
sampling_epi = 60; #number of points to put on epi
slice_sampling = 25; #number of vertical slices to generate (there will be 1 extra)
ax_slice = 0; #axis along which slices are placed
ax_1 = 1; #first inslice axis
ax_2 = 2; #second inslice axis. Make sure the triple ax_slice, ax_1, ax_2 form right-handed coordinate system. The order should be (0,1,2) or (2,0,1) or (1,2,0) 

inplane_smoothing = 7; #set by trial and error, depends on number of sampling points
betweenplane_smoothing = 5; #set by trial and error, depends on number of sampling points


#reference point. If you have it sotred somewhere, copy it here
pt_ref = np.array([0,0,0]);



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
    




def ResampleContoursInplane( points, sliceids, nsamples ):
    """
    Resamples every contour
    
    points Nx3 np.array: the points
    sliceids 1xN numpy array: stores the slice number for every point
    nsamples int: number of samples to take for each contour

    Return Mx3 numpy array of new points
    """    
    
    nslices = int(sliceids.max())+1
    points_new = np.zeros( (nsamples*nslices,3) )
    
    for i in range(nslices):
        pts = points[i==sliceids, :]
        cntr = pts.mean(axis=0)
        coord1 = pts[:,ax_1] #get inplane coordinates
        coord2 = pts[:,ax_2]
        
        ref_coord1 = pt_ref[ax_1] - cntr[ax_1]
        ref_coord2 = pt_ref[ax_2] - cntr[ax_2] 
        angle = np.arctan2( ref_coord2, ref_coord1 )
        
        coord1_cent = coord1 - cntr[ax_1]
        coord2_cent = coord2 - cntr[ax_2]
        angles = np.arctan2( coord2_cent, coord1_cent );
        
        #find closest point to the reference point according to the angle
        d = abs(angles - angle);
        ind = np.argmin(d)
        
        #shoft the points so they start at the reference point
        pts = np.roll(pts, shift=-ind, axis=0);
        m, n = pts.shape
        
        #fit spline and get new points
        tck, u = interpolate.splprep( [pts[:,ax_1], pts[:,ax_2]], s=inplane_smoothing)
        c = interpolate.splev(np.linspace(0,1,nsamples+1), tck)
        c[0] = c[0][0:-1]
        c[1] = c[1][0:-1]
        
        first_el = i*nsamples
        last_el = (i+1)*nsamples
        points_new[first_el:last_el,ax_slice] = pts[0,ax_slice]
        points_new[first_el:last_el,ax_1] = c[0]
        points_new[first_el:last_el,ax_2] = c[1]
        
    return points_new
    



def ResampleContoursLongitudinal( points, slice_sampling, contour_sampling ): 
    """
    Resamples the contours longitudinally using spline
    !Call after ResampleContoursInplane, bacause this assumes there is equal number of 
    points in every contour
   
    points Nx3 np.array: the points
    contour_sampling int: numer of samples each contour has
    slice_sampling int: number of samples to take 

    Return Mx3 numpy array of new points
    """    
    
    pts1_new=np.zeros( (contour_sampling*slice_sampling,3) );
    
    for i in range(contour_sampling):
        pts = points[ i::contour_sampling, : ];
        m, n = pts.shape
    
        tck, u = interpolate.splprep( [pts[:,0], pts[:,1], pts[:,2]], s=betweenplane_smoothing)
        c = interpolate.splev(np.linspace(0,1,slice_sampling), tck)
        
    
        pts1_new[i::contour_sampling, 0] = c[0]
        pts1_new[i::contour_sampling, 1] = c[1]
        pts1_new[i::contour_sampling, 2] = c[2]
    
    return pts1_new



def TriangulateGrid(samples_horizontal, samples_vertical):
    """
    Creates triangulation for a cylinder of samples_horizontal x samples_vertical points

    Return Mx3 numpy array of triangles
    """    

    x = np.arange(0, samples_horizontal+1)
    y = np.arange(0, samples_vertical)
    pts_eq_x, pts_eq_y = np.meshgrid(x, y);
    last_x = pts_eq_x[:,-1];
    last_y = pts_eq_y[:,-1];
    pts_eq_x1 = pts_eq_x[:,0:-1];
    pts_eq_y1 = pts_eq_y[:,0:-1];
    xx = np.concatenate( (pts_eq_x1.flatten(), last_x.flatten()) )
    yy = np.concatenate( (pts_eq_y1.flatten(),last_y.flatten()) )
    tri = Delaunay( np.column_stack( (xx, yy) ) )
    
    faces = tri.simplices
    for i in range(pts_eq_x1.size, pts_eq_x.size):
        faces[faces==i] = (samples_horizontal)*(i-pts_eq_x1.size) 
        
    return faces


#------------------------------------------------------------
#
#
#
#
#     MAIN CODE
#
#
#
#------------------------------------------------------------
    

#load the points
npzfile = np.load(input_filename)
pts_endo = npzfile['pts_endo']
pts_epi = npzfile['pts_epi']
sliceids_endo = npzfile['sliceids_endo']
sliceids_epi = npzfile['sliceids_epi']

#remove the 0th slice
if remove_0_slice:
    idx_endo = sliceids_endo>0
    idx_epi = sliceids_epi>0
    pts_endo = pts_endo[idx_endo, :]
    pts_epi = pts_epi[idx_epi, :]
    sliceids_endo = sliceids_endo[idx_endo]
    sliceids_epi = sliceids_epi[idx_epi]
    sliceids_endo = sliceids_endo-1
    sliceids_epi = sliceids_epi-1

#continue normally
nslices = int(sliceids_endo.max())
    
    
print 'Resampling inplane endo'
pts_endo_new = ResampleContoursInplane( pts_endo, sliceids_endo, sampling_endo )
print 'Resampling inplane epi'
pts_epi_new = ResampleContoursInplane( pts_epi, sliceids_epi, sampling_epi )

    
pts_endo_final = ResampleContoursLongitudinal( pts_endo_new, slice_sampling, sampling_endo)
pts_epi_final = ResampleContoursLongitudinal( pts_epi_new, slice_sampling, sampling_epi)


faces_endo = TriangulateGrid(sampling_endo, slice_sampling)
faces_epi = TriangulateGrid(sampling_epi, slice_sampling)

#plot for verification and save
figure()
triangular_mesh(pts_endo_final[:,0], pts_endo_final[:,1],pts_endo_final[:,2], faces_endo)
points3d(pts_endo[:,0], pts_endo[:,1], pts_endo[:,2], scale_mode='none', scale_factor=0.2)
pd_endo = mytools.CreatePolyData(pts_endo_final, faces_endo)
mytools.SaveVTKPolyData( pd_endo, output_endo_filename )


figure()
triangular_mesh(pts_epi_final[:,0], pts_epi_final[:,1],pts_epi_final[:,2], faces_epi)
points3d(pts_epi[:,0], pts_epi[:,1], pts_epi[:,2], scale_mode='none', scale_factor=0.2)
pd_epi = mytools.CreatePolyData(pts_epi_final, faces_epi)
mytools.SaveVTKPolyData( pd_epi, output_epi_filename )


np.savez(output_points_filename , pts_epi_final=pts_epi_final, pts_endo_final=pts_endo_final, sampling_endo = sampling_endo, sampling_epi = sampling_epi, slice_sampling = slice_sampling)


show()
