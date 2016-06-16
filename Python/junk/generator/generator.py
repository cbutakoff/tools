#!/usr/bin/env python
    
import vtk
import random
from array import array


reader = vtk.vtkPolyDataReader()
reader.SetFileName('./sphere.vtk')
reader.Update()

mesh = reader.GetOutput()

npts = reader.GetOutput().GetNumberOfPoints()

ntimes = 10  #number of timesteps

times = vtk.vtkFloatArray()
times.SetNumberOfValues( npts )

array = vtk.vtkFloatArray()
array.SetNumberOfValues( npts )

array.SetName('ecg')
times.SetName('times')


mesh.GetPointData().AddArray( array )
mesh.GetPointData().AddArray( times )


for i in range(ntimes):
    for j in range(npts):
        array.SetValue(j, random.random())
        times.SetValue(j, float(i)/10.0 )
        

    name = 'mesh{:02d}.vtk'.format(i)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName( name )
    writer.SetInput( mesh )
    writer.Update()
    
    #mesh.GetPointData().RemoveArray( 'ecg' )
    #mesh.GetPointData().RemoveArray( 'times' )
    
