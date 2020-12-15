import vtk
import sys
import progressbar
import os

#!!!!! Assumes tetras !!!!!!
arrayname = 'id'


def write_ruben_polydata(pd, filename):
    with open(filename, 'w') as ff:
        ff.write("# vtk DataFile Version 4.2\n")
        ff.write("vtk output\n")
        ff.write("ASCII\n")
        ff.write("DATASET POLYDATA\n")
        ff.write(f"POINTS {pd.GetNumberOfPoints()} float\n")
        for i in progressbar.progressbar(range(pd.GetNumberOfPoints())):
            ff.write(" ".join([str(c) for c in pd.GetPoint(i)])+"\n")

        ff.write(f"POLYGONS {pd.GetNumberOfPolys()} {pd.GetNumberOfPolys()*4}\n")
        
        for i in progressbar.progressbar(range(pd.GetNumberOfCells())):
            pts = pd.GetCell(i).GetPointIds()
            ff.write("3 "+" ".join([str(pts.GetId(j)) for j in range(pts.GetNumberOfIds())])+"\n")

        ff.write(f"CELL_DATA {pd.GetNumberOfPolys()}\n")
        ff.write(f"SCALARS RegionId int 1\n")
        ff.write(f"LOOKUP_TABLE default\n")

        arr = pd.GetCellData().GetArray( arrayname )
        for i in progressbar.progressbar(range(pd.GetNumberOfCells())):
            ff.write(f"{int(arr.GetTuple1(i))}\n")



def write_ruben_tetmesh(pd, filename):
    with open(filename, 'w') as ff:
        ff.write("# vtk DataFile Version 4.2\n")
        ff.write("vtk output\n")
        ff.write("ASCII\n")
        ff.write("DATASET UNSTRUCTURED_GRID\n")
        ff.write(f"POINTS {pd.GetNumberOfPoints()} float\n")

        for i in progressbar.progressbar(range(pd.GetNumberOfPoints())):
            ff.write(" ".join([str(c) for c in pd.GetPoint(i)])+"\n")


        ff.write(f"CELLS {pd.GetNumberOfCells()} {pd.GetNumberOfCells()*5}\n")
        for i in progressbar.progressbar(range(pd.GetNumberOfCells())):
            pts = pd.GetCell(i).GetPointIds()
            ff.write("4 "+" ".join([str(pts.GetId(j)) for j in range(pts.GetNumberOfIds())])+"\n")

        ff.write(f"CELL_TYPES {pd.GetNumberOfCells()} \n")
        for i in progressbar.progressbar(range(pd.GetNumberOfCells())):
            ff.write("10\n")




rd = vtk.vtkPolyDataReader()
rd.SetFileName('surf-lbl.vtk')
rd.Update()

write_ruben_polydata(rd.GetOutput(),'surf-lbl-ruben.vtk')


rd = vtk.vtkDataSetReader()
rd.SetFileName('volmesh.vtk')
rd.Update()

write_ruben_tetmesh(rd.GetOutput(),'volmesh-ruben.vtk')

