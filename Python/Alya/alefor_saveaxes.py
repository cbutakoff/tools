#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import progressbar
import io
import os



def write_vector_mpio(filename, vector, time):
    ndim = vector.shape[1]
    npts = vector.shape[0]

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(b'XFIEL00\0')
        f.write(b'VECTO00\0')
        f.write(b'NPOIN00\0')
        f.write(b'REAL000\0')
        f.write(b'8BYTE00\0')
        f.write(b'SEQUE00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array([0], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array([1], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array([0], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag2
        np.array([time], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
        f.write(b'0000000\0')
        f.write(b'NONE000\0') #1
        f.write(b'NONE000\0') #2
        f.write(b'NONE000\0') #3
        f.write(b'NONE000\0') #4
        f.write(b'NONE000\0') #5
        f.write(b'NONE000\0') #6
        f.write(b'NONE000\0') #7
        f.write(b'NONE000\0') #8
        f.write(b'NONE000\0') #9
        f.write(b'NONE000\0') #10
        vector.astype('float64').tofile(f)
# In[3]:



surface_filename = "moving_boundary.vtk" #only the surfaec that is to be moving. Contains pointdata: V1(normal), V2(tangent), V3=V1xV2. 
problem_name = "piece"
volume_filename = "piece/piece_0_0.vtu"
v1_name = "V1"
v2_name = "V2"
v3_name = "V3"


print("Reading volume")
vr = vtk.vtkXMLUnstructuredGridReader()
vr.SetFileName(volume_filename)
vr.Update()
vol = vr.GetOutput()

print("Reading surface")
rdr = vtk.vtkPolyDataReader()
rdr.SetFileName(surface_filename)
rdr.Update()
mesh = rdr.GetOutput()

v1a = mesh.GetPointData().GetArray( v1_name )
v2a = mesh.GetPointData().GetArray( v2_name )
v3a = mesh.GetPointData().GetArray( v3_name )


matrix = np.zeros( (vol.GetNumberOfPoints(), 9) )
loc = vtk.vtkPointLocator()
loc.SetDataSet(vol)
loc.BuildLocator()

for i in progressbar.progressbar(range(mesh.GetNumberOfPoints())):
    pt = mesh.GetPoint(i)
    ptid = loc.FindClosestPoint(pt) 
    matrix[ptid, 0:3] = v1a.GetTuple(i)
    matrix[ptid, 3:6] = v2a.GetTuple(i)
    matrix[ptid, 6:9] = v3a.GetTuple(i)
    
print("1st row: ",matrix[0,:])

filename = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format(problem_name, 1, 1)
print('Saving ', filename)
write_vector_mpio(  filename, matrix, 0 )




