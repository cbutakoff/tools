#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 21:01:41 2017

@author: costa
"""

mesh_temp_transmural = 'temperature_transmural.vtk'
simulationmesh_filename = 'final_mesh.vtk'
output_filename = 'transmural.txt'
output_vtk = 'transmural.vtk'


import vtk

print("Reading mesh")
rd2 = vtk.vtkDataSetReader()
rd2.SetFileName(mesh_temp_transmural)
rd2.Update()

print("Reading mesh 2" )
rd = vtk.vtkDataSetReader()
rd.SetFileName(simulationmesh_filename)
rd.Update()


mesh_simulation = rd.GetOutput()

t_mesh = rd2.GetOutput()
t = rd2.GetOutput().GetPointData().GetArray('temperature')


loc = vtk.vtkPointLocator()
loc.SetDataSet(mesh_simulation)
loc.BuildLocator()

from TextProgressBar import TextProgressBar


import numpy as np

tt = np.zeros(mesh_simulation.GetNumberOfPoints())

pb = TextProgressBar(t_mesh.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 40)
for i in range(t.GetNumberOfTuples()):
    pb.UpdateProgress(i + 1)
    ptid = loc.FindClosestPoint( t_mesh.GetPoint(i) )
    tt[ptid] = t.GetTuple1(i)


ttt = vtk.vtkFloatArray()
ttt.SetName('temperature')
ttt.SetNumberOfComponents(1)
ttt.SetNumberOfTuples(tt.shape[0])

with open(output_filename, 'w') as f:
    for i in range(tt.shape[0]):
        f.write(f'{i+1}  {tt[i]}\n')
        ttt.SetTuple1(i, tt[i])
        
mesh_simulation.GetPointData().AddArray(ttt)

print("Saving")
wr = vtk.vtkDataSetWriter()
wr.SetFileName(output_vtk)
wr.SetFileTypeToBinary()
wr.SetInputData(mesh_simulation)
wr.Write()

    
