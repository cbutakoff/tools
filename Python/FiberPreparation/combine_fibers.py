#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 19:59:54 2017

@author: costa
"""

mesh_streeter_filename = "for_fibers_fib.vtk"
mesh_temp_longitudinal = "temperature_longitudinal.vtk"
mesh_complete = "myo.vtk"
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

print("Reading yet another mesh")
rd3 = vtk.vtkDataSetReader()
rd3.SetFileName(mesh_complete)
rd3.Update()

mesh_full = rd3.GetOutput()


print("Gradient filter")
grad_filter = vtk.vtkGradientFilter()
grad_filter.SetInputData(rd2.GetOutput())
grad_filter.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,'temperature')
grad_filter.Update()

grad_mesh = grad_filter.GetOutput()
grad = grad_mesh.GetPointData().GetArray('Gradients')

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
pb = TextProgressBar(grad_mesh.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 40)
for i in range(grad_mesh.GetNumberOfPoints()):
    pb.UpdateProgress(i + 1)

    ptid = loc.FindClosestPoint( grad_mesh.GetPoint(i) )
    if fibers_mask[ptid]==0:
        v = numpy.array(grad.GetTuple3(i))
        v = v/numpy.linalg.norm(v)
        fibers_full.SetTuple3(ptid, v[0], v[1], v[2])
        
mesh_full.GetPointData().AddArray(fibers_full)

print("Saving")
wr = vtk.vtkDataSetWriter()
wr.SetFileName(output_filename)
wr.SetFileTypeToBinary()
wr.SetInputData(mesh_full)
wr.Write()



