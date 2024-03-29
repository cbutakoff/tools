#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 09:58:49 2013
Needs OpenCV, Numpy, ITK

Extracts the contours from masks of epi and endo, stored in nrrd files (generated by seg3d)
The MRI slices are assumed to be along X axis (if not - modify the code). 

@author:  Constantine Butakoff
"""
import itk
import numpy as np
import cv2


#
#  Parameters, adapt to your taste
#
endo_filename = 'endo.nrrd'
epi_filename = 'epi.nrrd'
output_filename = 'points.npz'  #numpy file with epi and endo

#tested only when slices are along X axis. In other situations double check that the points still fit the image. There could be some orientation problem
#rearrange according to your images
ax_slice = 0  #slice axis -- slices along X axis
ax_1 = 1  #inplane axis 1
ax_2 = 2  #inplane axis 2


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
    


def ExtractMRIContoursAlongXAxis(image, origin, spacing):
    """
    Extract the contours from the mask
    
    image itkImage: the image
    origin 1x3 numpy array: image origin
    spacing 1x3 numpy array: image spacing

    Return (pts, sliceids) -- Mx3 numpy array of new points and 1xM array of values of slice number for every point
    """        
    
    #first create a list of slice contours, then concatenate them into a matrix
    slices = list() 
    npts = 0 #total number of points
    for slice_no in range(xsize):
        slice_data = np.zeros((ysize, zsize),np.uint8)

        for i in range(ysize):
            for j in range(zsize):
                slice_data[i,j] = image.GetPixel([slice_no,i,j])
                
        #trace the contour
        ret,thresh = cv2.threshold(slice_data*10,5,255,cv2.THRESH_BINARY)
        contours, hierarchy = cv2.findContours(thresh,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
        
        #save the contour vertices
        pts2d = np.zeros( (len(contours[0]), 3) )
        pts2d[:,ax_slice] = slice_no
        for i in range(len(contours[0])):
            pts2d[i,ax_1] = contours[0][i][0][1]-1 #swap dimensions! and correct for (0,0)
            pts2d[i,ax_2] = contours[0][i][0][0]-1
            
            
        #transform the vertices into physical coordinates
        pts3d = origin + pts2d * spacing
        slices.append(pts3d)
        
        #calculate the number of points
        m,n = pts3d.shape
        npts = npts + m
    #end for slice_no    
        
    #store all the points in a single matrix
    pts = np.zeros((npts,3))
    sliceids = np.zeros(npts)
    
    ind=0
    for i in range(len(slices)):
        m,n = slices[i].shape
        pts[ind:ind+m,:] = slices[i]
        sliceids[ind:ind+m] = i #indicator array, for every point stores its slice number
        ind = ind+m
    #end for i
    
    return (pts, sliceids)




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
    

#
#
#  The remianing processing
#


pixelType = itk.UC
imageType = itk.Image[pixelType, 3]
readerType = itk.ImageFileReader[imageType]
reader = readerType.New()
reader.SetFileName(endo_filename)
reader.Update()
image = reader.GetOutput()

#extract a slice
size = image.GetLargestPossibleRegion().GetSize()
xsize = size.GetElement(0)
ysize = size.GetElement(1)
zsize = size.GetElement(2)

itkorigin = image.GetOrigin()
itkspacing = image.GetSpacing()
origin = np.zeros(3);
spacing = np.zeros(3);
origin[0] = itkorigin.GetElement(0)
origin[1] = itkorigin.GetElement(1)
origin[2] = itkorigin.GetElement(2)
spacing[0] = itkspacing.GetElement(0)
spacing[1] = itkspacing.GetElement(1)
spacing[2] = itkspacing.GetElement(2)

print 'Extracting endo...'
pts_endo, sliceids_endo = ExtractMRIContoursAlongXAxis(reader.GetOutput(), origin, spacing)
reader = readerType.New()
reader.SetFileName(epi_filename)
reader.Update()

print 'Extracting epi...'
pts_epi, sliceids_epi = ExtractMRIContoursAlongXAxis(reader.GetOutput(), origin, spacing)

np.savez(output_filename , pts_epi=pts_epi, pts_endo=pts_endo, sliceids_endo=sliceids_endo, sliceids_epi=sliceids_epi)
