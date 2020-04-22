#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import progressbar
import io
import os
import sys


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



surface_filename = sys.argv[1] #"moving_boundary.vtk" #only the surface that is to be moving. 
volume_filename = sys.argv[2] #"piece/piece_0_0.vtu"
problem_name = sys.argv[3] #"piece"

surface_label = 5
savevtk =  False
volmesh_scale = 0.01 #if volmesh has a scale different to the surface mesh, put here the scaling factor for the transform
boundary_array = 'BoundaryId'


print("Reading volume")
vr = vtk.vtkXMLUnstructuredGridReader()
vr.SetFileName(volume_filename)
vr.Update()


tff = vtk.vtkTransform()
tff.Scale(volmesh_scale,volmesh_scale,volmesh_scale)

print("Rescaling vol mesh")
tf = vtk.vtkTransformFilter()
tf.SetInputData(vr.GetOutput())
tf.SetTransform(tff)
tf.Update()

vol = tf.GetOutput()


print("Reading surface")
rdr = vtk.vtkPolyDataReader()
rdr.SetFileName(surface_filename)
rdr.Update()


print(f"Extracting suurface {surface_label}")
th = vtk.vtkThreshold()
th.SetInputData(rdr.GetOutput())
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, boundary_array)
th.ThresholdBetween(surface_label-0.1, surface_label+0.1)
th.Update()

sf = vtk.vtkDataSetSurfaceFilter()
sf.SetInputData(th.GetOutput())
sf.Update()



print("Calculating normals")
normalgen = vtk.vtkPolyDataNormals()
normalgen.SetInputData(sf.GetOutput())
normalgen.SplittingOff()
normalgen.Update()
mesh = normalgen.GetOutput()

normals = mesh.GetPointData().GetNormals()


matrix = np.zeros( (vol.GetNumberOfPoints(), 9) )

normal_coeff = 1
#check if normals need to be flipped
#the normals will point inside, for now
d = 0.001

loc = vtk.vtkCellLocator()
loc.SetDataSet(vol)
loc.BuildLocator()

for i in progressbar.progressbar(range(mesh.GetNumberOfPoints())):
    c = np.array(mesh.GetPoint(i))
    n = np.array(normals.GetTuple(i))
    pt1 = c+n*d
    pt2 = c-n*d

    GenCell = vtk.vtkGenericCell()
    pcoords = [0,0,0]
    weights = [0,0,0,0,0,0,0,0,0,0,0,0]
    cellid1 = loc.FindCell ( pt1, 1e-10, GenCell, pcoords, weights) 
    cellid2 = loc.FindCell ( pt2, 1e-10, GenCell, pcoords, weights) 
    if ((cellid1<0) and (cellid2<0)) or ((cellid1>=0) and (cellid2>=0)):
        continue
    else:
        if cellid1>0:
           normal_coeff = 1
           break
        else:
           normal_coeff = -1
           break


loc = vtk.vtkPointLocator()
loc.SetDataSet(vol)
loc.BuildLocator()

for i in progressbar.progressbar(range(mesh.GetNumberOfPoints())):
    pt = mesh.GetPoint(i)
    ptid = loc.FindClosestPoint(pt) 

    n = normals.GetTuple(i)*normal_coeff
    t1 = np.roll(n,1) #create random vector != n
    t2 = np.cross(n, t1) #find orthogonal vector
    t1 = t2/np.linalg.norm(t2) #normalise
    t2 = np.cross(n,t1) #find the 3rd vector

    matrix[ptid, 0:3] = n 
    matrix[ptid, 3:6] = t1
    matrix[ptid, 6:9] = t2
    
print("1st row: ",matrix[0,:])

filename = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format(problem_name, 1, 1)
print('Saving ', filename)
write_vector_mpio(  filename, matrix, 0 )

if savevtk:
    normals = vtk.vtkFloatArray()
    normals.SetName('MyNormal')
    normals.SetNumberOfComponents(3)
    normals.SetNumberOfTuples(vol.GetNumberOfPoints())
    for i in progressbar.progressbar(range(vol.GetNumberOfPoints())):
        normals.SetTuple(i, matrix[i,  0:3])

    vol.GetPointData().AddArray(normals)

    wwr = vtk.vtkDataSetWriter()
    wwr.SetFileTypeToBinary()
    wwr.SetFileName('volnormals.vtk')
    wwr.SetInputData(vol)
    wwr.Write()





