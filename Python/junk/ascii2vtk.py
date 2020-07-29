import vtk
import sys
import pandas as pd
from progressbar import progressbar

file_nodes = sys.argv[1]
file_tetra = sys.argv[2]
file_aha = sys.argv[3]
output_vtu = sys.argv[4]

print("reading")
nodes = pd.read_table(file_nodes, skiprows=1, skipfooter=1, header=None, names=["id","x","y","z"], delim_whitespace=True, engine='python')
tetra = pd.read_table(file_tetra, skiprows=1, skipfooter=1, header=None, names=["id","p1","p2","p3","p4"], delim_whitespace=True, engine='python')
aha   = pd.read_table(file_aha,  header=None, names=["id","aha"], delim_whitespace=True)

nodes_v = nodes.values
tetra_v = tetra.values
aha_v = aha.values

pts = vtk.vtkPoints()
pts.SetNumberOfPoints(nodes_v.shape[0])

print("Processing nodes")
for i in progressbar(range(nodes_v.shape[0])):
    pts.SetPoint(i, nodes_v[i,1:4])


print("Processing elements")
cells = vtk.vtkCellArray()
for i in progressbar(range(tetra_v.shape[0])):
    cells.InsertNextCell(4)
    cells.InsertCellPoint( tetra_v[i,1]-1 )
    cells.InsertCellPoint( tetra_v[i,2]-1 )
    cells.InsertCellPoint( tetra_v[i,3]-1 )
    cells.InsertCellPoint( tetra_v[i,4]-1 )


print("Processing aha")
array = vtk.vtkShortArray()
array.SetName("AHA")
array.SetNumberOfComponents(1)
array.SetNumberOfTuples(tetra_v.shape[0])

for i in progressbar(range(aha_v.shape[0])):
    array.SetTuple1(i, aha_v[i,1])

ug = vtk.vtkUnstructuredGrid()
ug.SetPoints(pts)
ug.SetCells(vtk.VTK_TETRA, cells)
ug.GetCellData().AddArray( array )

wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputData(ug)
wr.EncodeAppendedDataOff()
wr.SetFileName(output_vtu)
wr.Write()
