import vtk
import sys
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk

inputfilename = sys.argv[1]
outputfilename = sys.argv[2]
parts_to_exclude = ['fluid']

rd = vtk.vtkCGNSReader()
rd.SetFileName(inputfilename)
rd.EnableAllBases()
rd.LoadMeshOn()
rd.Update()
mbd = rd.GetOutput()

ug = vtk.vtkUnstructuredGrid()
merge = vtk.vtkMergeCells()
merge.MergeDuplicatePointsOn()
merge.SetUnstructuredGrid(ug)

npts = 0
ncells = 0
nparts = 0
for i in range(mbd.GetNumberOfBlocks()):
    blocks = mbd.GetBlock(i)
    for j in range(blocks.GetNumberOfBlocks()):
        block = blocks.GetBlock(j)
        npts   += block.GetNumberOfPoints()
        ncells += block.GetNumberOfCells()
        nparts += 1

merge.SetTotalNumberOfCells(ncells)
merge.SetTotalNumberOfPoints(npts)
merge.SetTotalNumberOfDataSets(nparts)

material = 1
materials = {}
for i in range(mbd.GetNumberOfBlocks()):
    blocks = mbd.GetBlock(i)

    for j in range(blocks.GetNumberOfBlocks()):
        block = blocks.GetBlock(j)
        metadata = { x.split(':')[0].strip():x.split(':')[1].strip() for x in blocks.GetMetaData(j).__str__().split('\n') if len(x.split(':'))==2 }
        name = metadata["NAME"]
        print(name)
        if name not in parts_to_exclude:
            array = numpy_to_vtk( material*np.ones(block.GetNumberOfCells(), dtype=int) )
            array.SetName('Material')
            block.GetCellData().AddArray(array)
            materials[name] = material
            material = material+1

            npts   += block.GetNumberOfPoints()
            ncells += block.GetNumberOfCells()
            nparts += 1
            merge.MergeDataSet(block)


merge.Finish()

print(materials)

wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputData(ug)
wr.SetFileName(outputfilename)
wr.Write()
