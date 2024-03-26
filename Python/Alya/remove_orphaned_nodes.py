import vtk
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import sys
from numpy_indexed import remap

filenamein = sys.argv[1]
filenameout = sys.argv[2]

print(f"Read {filenamein}")
rd = vtk.vtkXMLUnstructuredGridReader()
rd.SetFileName(filenamein)
rd.Update()
mesh = rd.GetOutput()

nodes = vtk_to_numpy(mesh.GetPoints().GetData())
offsets = mesh.GetCells().GetOffsetsArray()
cells = vtk_to_numpy(mesh.GetCells().GetConnectivityArray())

#find orphaned nodes
print(f"Find orphans")
node_mask = np.zeros(nodes.shape[0], dtype=np.uint8)
node_mask[cells] = 1
used_node_mask = node_mask==1
new_nodes = nodes[used_node_mask]

if( not any(node_mask==0) ):
   print("No orphans found. Nothing to output")
else:
   new_node_ids = np.arange(new_nodes.shape[0], dtype=int)
   used_old_node_ids = np.where(used_node_mask==True)[0]
   print(f"Found {np.sum(node_mask==0)} orphaned nodes. Remapping the nodes")
   new_cells = remap( cells, used_old_node_ids, new_node_ids )

   pts_vtk = vtk.vtkPoints()
   pts_vtk.SetData( numpy_to_vtk(new_nodes) )
   new_mesh = vtk.vtkUnstructuredGrid()
   new_mesh.SetPoints(pts_vtk)
   new_mesh.SetCells( mesh.GetCellTypesArray(),  mesh.GetCells() )
   new_mesh.GetCells().SetData( mesh.GetCells().GetOffsetsArray(), numpy_to_vtk(new_cells) )

   #copy cell arrays
   print("copy cell arrays")
   for i in range(mesh.GetCellData().GetNumberOfArrays()):
      new_mesh.GetCellData().AddArray(mesh.GetCellData().GetArray(i))

   #renumber point arrays
   print("renumber point arrays")
   for i in range(mesh.GetPointData().GetNumberOfArrays()):
      array = vtk_to_numpy(mesh.GetPointData().GetArray(i))
      if array.ndim == 1:
         new_array = array[used_node_mask]
      else:
         new_array = array[used_node_mask,:]

      array_vtk = numpy_to_vtk(new_array)
      name = mesh.GetPointData().GetArray(i).GetName()
      array_vtk.SetName(name)
      new_mesh.GetPointData().AddArray(array_vtk)

   print(f"Writing {filenameout}")
   wr = vtk.vtkXMLUnstructuredGridWriter()
   wr.SetInputData(new_mesh)
   wr.SetFileName(filenameout)
   wr.EncodeAppendedDataOff()
   wr.Write()

