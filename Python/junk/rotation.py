#!/usr/bin/python
#parameters: imagefile outputimagefile normalsfile
#file with normals has 2 rows, each with x,y,z coordinates of the 2 planes aligned with the long axis

import sys

# coding: utf-8

# In[28]:

import vtk
import numpy as np


# # Parameters

# In[29]:




imagefilename = str(sys.argv[1])
outputfilename = str(sys.argv[2])
normalsfile = str(sys.argv[3])

normals = np.loadtxt(normalsfile)

#plane perpendicular to both
#its normal aligned with the long axis
n = np.cross(normals[0,:], normals[1,:])

# In[30]:

#two points of the line to align with axis Z
#p1 = np.array([351.664963853608,176.121879872494,245.396158561977]); #basal
#p2 = np.array([514.891869997372,226.322051596017,621.038567084717]); #apical

output_padding = 100; #how much to expand the output image to fit the rotated version


# # Internal stuff. Do not modify

# In[31]:

#  interanl stuff
#
#
#
#l = p2-p1;
l = n
l = l/np.linalg.norm(l);


axis = np.cross([0.0,0.0,0.1],l) #axis of rotation
axis = axis/np.linalg.norm(axis);

t = np.arccos( l.dot([0.0,0.0,1.0]) ); #angle fo rotation


# In[121]:




# In[122]:

def CreateRotation(axis, t):
    u = axis[0]
    v = axis[1]
    w = axis[2]
    
    Txz = np.matrix([[u/np.sqrt(u*u+v*v), v/np.sqrt(u*u+v*v), .0, .0],                 [-v/np.sqrt(u*u+v*v), u/np.sqrt(u*u+v*v), .0, .0],                 [.0,.0,1.0,.0],[.0,.0,.0,1.0]])
    Tz = np.matrix([[w/np.sqrt(u*u+v*v+w*w), 0.0, -np.sqrt(u*u+v*v)/np.sqrt(u*u+v*v+w*w),  .0],               [.0,1.0,0.0,.0],                [np.sqrt(u*u+v*v)/np.sqrt(u*u+v*v+w*w), 0.0, w/np.sqrt(u*u+v*v+w*w), .0],                [.0,.0,.0,1.0]])
    Rz = np.matrix([[np.cos(t), -np.sin(t), .0, .0],                [np.sin(t), np.cos(t), .0, .0],                [.0,.0,1.0,.0],[.0,.0,.0,1.0]])
    T = np.linalg.inv(Txz).dot(np.linalg.inv(Tz).dot(Rz.dot(Tz.dot(Txz))))
    
    return T

def TransformPoints(pts, T):
    pts_a = np.append(pts, np.ones((1, pts.shape[1])), axis=0)
    pts1_a = T.dot(pts_a)
    pts1 = np.array(pts1_a[0:3,:])
    return pts1


# In[123]:

#read he image
reader = vtk.vtkDataSetReader()
reader.SetFileName(imagefilename);
reader.Update()
image = reader.GetOutput()


# In[124]:

T = CreateRotation(axis,t)

#add translation
#Compute the center of the image
bounds = image.GetBounds()
center = np.zeros(3);
center[0] = (bounds[1] + bounds[0]) / 2.0;
center[1] = (bounds[3] + bounds[2]) / 2.0;
center[2] = (bounds[5] + bounds[4]) / 2.0;

Trans = np.matrix([ [1.0, 0.0, 0.0, -center[0]],                     [0.0, 1.0, 0.0, -center[1]],                     [0.0, 0.0, 1.0, -center[2]],                    [0.0, 0.0, 0.0, 1.0]])
#T = np.linalg.inv(Trans).dot(np.linalg.inv(T).dot(Trans))


# In[125]:

changeinfo = vtk.vtkImageChangeInformation()
changeinfo.SetInput(reader.GetOutput())
changeinfo.CenterImageOn()
changeinfo.Update()
image = changeinfo.GetOutput()


# In[126]:

#pts1 = TransformPoints(pts,T)


# In[127]:

matr = vtk.vtkMatrix4x4()
for i in range(4):
    for j in range(4):
        matr.SetElement(i,j, T[i,j]);


# In[128]:

# Rotate about the center of the image
transform = vtk.vtkTransform();
transform.SetMatrix(matr);


# In[129]:

#Reslice does all of the work
reslice = vtk.vtkImageReslice();
reslice.SetInput(image);
reslice.SetResliceTransform(transform);
reslice.SetInterpolationModeToCubic();
reslice.SetOutputSpacing( image.GetSpacing() );
reslice.SetOutputOrigin( image.GetOrigin() );
reslice.SetOutputExtent( image.GetExtent()[0]-output_padding, image.GetExtent()[1]+output_padding,                         image.GetExtent()[2]-output_padding, image.GetExtent()[3]+output_padding,                         image.GetExtent()[4]-output_padding, image.GetExtent()[5]+output_padding); 
reslice.Update()


# In[130]:

#write
wr = vtk.vtkDataSetWriter()
wr.SetFileName(outputfilename)
wr.SetInput(reslice.GetOutput())
wr.SetFileTypeToBinary()
wr.Write()


# In[130]:


