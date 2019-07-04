import vtk
import numpy as np
from progressbar import progressbar

inmesh = 'volmesh+scar-clean.vtu'
outmesh = 'volmesh+scar-clean1.vtu'



print('Reading mesh')
#rd = vtk.vtkDataSetReader()
rd = vtk.vtkXMLUnstructuredGridReader()
rd.SetFileName(inmesh)
rd.Update()
mesh = rd.GetOutput()

nodemask = np.zeros(mesh.GetNumberOfPoints(), dtype=np.int32)

print('Checking elements')
for i in progressbar(range( mesh.GetNumberOfCells() ) ):
    ids = vtk.vtkIdList()
    mesh.GetCellPoints(i, ids) 

    for j in range( ids.GetNumberOfIds() ):
        nodemask[ ids.GetId(j) ] = 1



print('Number of orphan nodes: ', nodemask.shape[0]-nodemask.sum())


if nodemask.shape[0]-nodemask.sum()>0 :
    print('Removing orphan nodes')


    pts = vtk.vtkPoints()
    corresp = np.ones(mesh.GetNumberOfPoints(), dtype=np.int64)*(-1)

    for i in progressbar(range(mesh.GetNumberOfPoints())):
      if nodemask[i]==1:
        ptid = pts.InsertNextPoint(mesh.GetPoint(i))  
        corresp[i] = ptid


    celltypes = []
    cells = vtk.vtkCellArray()
    for i in progressbar(range(mesh.GetNumberOfCells())):
        celltypes.append( mesh.GetCellType(i) )
        cell = mesh.GetCell(i)
        cells.InsertNextCell( cell.GetPointIds().GetNumberOfIds() )

        for j in range( cell.GetPointIds().GetNumberOfIds() ):
          oldptid = cell.GetPointIds().GetId(j)
          cells.InsertCellPoint( corresp[oldptid] )

    newmesh = vtk.vtkUnstructuredGrid()
    newmesh.SetPoints(pts)
    newmesh.SetCells(celltypes, cells)
    newmesh.GetCellData().AddArray( mesh.GetCellData().GetArray('Scar') )


    print('Passing the point arrays')
    probe = vtk.vtkProbeFilter()
    probe.SetSourceData(mesh)
    probe.SetInputData(newmesh)
    probe.PassCellArraysOn ()
    probe.Update()




    wr = vtk.vtkXMLUnstructuredGridWriter()
    wr.SetFileName(outmesh)
    wr.SetInputData(probe.GetOutput())
    wr.SetDataModeToAppended()
    wr.EncodeAppendedDataOff()
    wr.Write()

