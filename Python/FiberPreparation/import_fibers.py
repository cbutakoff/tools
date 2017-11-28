#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 19:21:45 2017

@author: costa
"""

mesh_filename = "for_fibers.vtk"
output_filename = "for_fibers_fib.vtk"
fibers_filename = "fibers_field.csv"

import vtk
import numpy

fibers = numpy.loadtxt(fibers_filename, usecols=[1,2,3])

rd = vtk.vtkDataSetReader()
rd.SetFileName(mesh_filename)
rd.Update()

mesh = rd.GetOutput()

vtk_fibers = vtk.vtkFloatArray()
vtk_fibers.SetName('Fibers')
vtk_fibers.SetNumberOfComponents(3)
vtk_fibers.SetNumberOfTuples(mesh.GetNumberOfPoints())

for i in range(mesh.GetNumberOfPoints()):
    vtk_fibers.SetTuple3(i, fibers[i,0], fibers[i,1], fibers[i,2])
    
mesh.GetPointData().AddArray(vtk_fibers);

wr = vtk.vtkDataSetWriter()
wr.SetFileTypeToBinary()
wr.SetFileName(output_filename)
wr.SetInputData(mesh)
wr.Write()
    