#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 19:59:54 2017

@author: costa
"""

mesh_streeter_filename = "for_fibers_fib.vtk"
mesh_temp_longitudinal = "temperature_longitudinal.vtk"
output_filename = "final_mesh.vtk"

import vtk
import numpy
from TextProgressBar import TextProgressBar

print('Reading one mesh')
rd1 = vtk.vtkDataSetReader()
rd1.SetFileName(mesh_streeter_filename)
rd1.Update()

mesh_streeter = rd1.GetOutput()

print("Reading another mesh")
rd2 = vtk.vtkDataSetReader()
rd2.SetFileName(mesh_temp_longitudinal)
rd2.Update()

print("Gradient filter")
grad_filter = vtk.vtkGradientFilter()
grad_filter.SetInputData(rd2.GetOutput())
grad_filter.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,'temperature')
grad_filter.Update()

mesh_full = rd2.GetOutput() #this is full mesh
grad = grad_filter.GetOutput().GetPointData().GetArray('Gradients')

fibers_full = vtk.vtkFloatArray()
fibers_full.SetName('Fibers')
fibers_full.SetNumberOfComponents(3)
fibers_full.SetNumberOfTuples(mesh_full.GetNumberOfPoints())

fibers_mask = numpy.zeros(mesh_full.GetNumberOfPoints(), dtype=numpy.int8)
fibers_mask.fill(0)

loc = vtk.vtkPointLocator()
loc.SetDataSet(mesh_full)
loc.BuildLocator()

streeter_fibers = mesh_streeter.GetPointData().GetArray('Fibers')

print("Copy the streeter fibers from the submesh")
pb = TextProgressBar(mesh_streeter.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 40)
for i in range(mesh_streeter.GetNumberOfPoints()):
    ptid = loc.FindClosestPoint( mesh_streeter.GetPoint(i) )
    fibers_full.SetTuple(ptid, streeter_fibers.GetTuple(i))
    fibers_mask[ptid] = 1
    pb.UpdateProgress(i + 1)
    
print("Append the gradients as fibers")
pb = TextProgressBar(mesh_full.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 40)
for i in range(mesh_full.GetNumberOfPoints()):
    pb.UpdateProgress(i + 1)
    if fibers_mask[i]==0:
        v = numpy.array(grad.GetTuple3(i))
        v = v/numpy.linalg.norm(v)
        fibers_full.SetTuple3(i, v[0], v[1], v[2])
        
mesh_full.GetPointData().AddArray(fibers_full)

print("Saving")
wr = vtk.vtkDataSetWriter()
wr.SetFileName(output_filename)
wr.SetFileTypeToBinary()
wr.SetInputData(mesh_full)
wr.Write()



