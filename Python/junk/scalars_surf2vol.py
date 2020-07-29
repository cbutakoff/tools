import vtk
from progressbar import progressbar
import numpy as np

volmesh_filename = 'aha.vtu'
surfmesh_filename = 'surf_0_0.vtp'
outsurf_filename = 'surf_aha.vtp'
cellarray = 'AHA'
volmesh_rescale = 0.75


print("Reading surf")
rd_surf = vtk.vtkXMLPolyDataReader()
rd_surf.SetFileName(surfmesh_filename)
rd_surf.Update()
surfmesh = rd_surf.GetOutput()

print("Getting centroids")
cc = vtk.vtkCellCenters()
cc.SetInputConnection( rd_surf.GetOutputPort() )
cc.Update()
c = cc.GetOutput()

print("Reading vol")
rd_vol = vtk.vtkXMLUnstructuredGridReader()
rd_vol.SetFileName(volmesh_filename)
rd_vol.Update()

sf = vtk.vtkDataSetSurfaceFilter()
sf.SetInputConnection( rd_vol.GetOutputPort() )
sf.Update()

T = vtk.vtkTransform()
T.Identity()
T.Scale(volmesh_rescale,volmesh_rescale,volmesh_rescale)
transformFilter = vtk.vtkTransformFilter()
transformFilter.SetInputConnection( sf.GetOutputPort() )
transformFilter.SetTransform(T)
transformFilter.Update()

volmesh = transformFilter.GetOutput()


v_array = volmesh.GetCellData().GetArray(cellarray)

s_array = vtk.vtkShortArray()
s_array.SetName(cellarray)
s_array.SetNumberOfComponents( v_array.GetNumberOfComponents() )
s_array.SetNumberOfTuples( surfmesh.GetNumberOfCells() )

loc = vtk.vtkCellLocator()
loc.SetDataSet(volmesh)
loc.BuildLocator()


closestPoint = [0,0,0]
cell = vtk.vtkGenericCell()
cellId = vtk.mutable(0)
subId = vtk.mutable(0)
dist2 = vtk.mutable(0)

for i in progressbar(range( surfmesh.GetNumberOfCells() )):
    loc.FindClosestPoint( c.GetPoint(i), closestPoint, cell, cellId, subId, dist2 ) 		
    s_array.SetTuple1( i, v_array.GetTuple1(cellId) )


surfmesh.GetCellData().AddArray( s_array )

print("writing")
wr = vtk.vtkXMLPolyDataWriter()
wr.SetInputData( surfmesh )
wr.SetFileName( outsurf_filename )
wr.Write()
