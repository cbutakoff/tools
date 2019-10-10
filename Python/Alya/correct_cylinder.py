import vtk
import numpy as np
import mpmath
import sys
from progressbar import progressbar

volmesh_filename = sys.argv[1]
volmesh_filename_out = sys.argv[2]

mpmath.mp.dps = 30 #set decimal precision

R=mpmath.mpf(0.4) #cm
L=mpmath.mpf(4.0)  #cm

cylinder_axis = np.array([0,0,1], dtype=int)  #accepted only X,Y or Z  axes


print('Reading mesh')
rd = vtk.vtkXMLUnstructuredGridReader()
rd.SetFileName(volmesh_filename)
rd.Update()

volmesh = rd.GetOutput()

print('Generating point ids')
idf = vtk.vtkIdFilter()
idf.SetInputData(rd.GetOutput())
idf.PointIdsOn ()
idf.CellIdsOff ()
idf.Update()

ids_array = idf.GetIdsArrayName()

print('Extracting surface')
sf = vtk.vtkDataSetSurfaceFilter()
sf.SetInputData(idf.GetOutput())
sf.Update()


print('Generate cell normals to identify boundary')
nf = vtk.vtkPolyDataNormals()
nf.SetInputData(sf.GetOutput())
nf.SplittingOff ()
nf.ComputePointNormalsOff ()
nf.ComputeCellNormalsOn ()
nf.Update()


#normals = vtkArrayDownCast<vtkFloatArray>(output->GetCellData()->GetArray("Normals"))


print('Multiply face normal by the cylinder axis')
calc = vtk.vtkArrayCalculator()
calc.SetInputData(nf.GetOutput())
calc.SetAttributeTypeToCellData()
calc.AddScalarVariable ("nx", "Normals", 0)
calc.AddScalarVariable ("ny", "Normals", 1)
calc.AddScalarVariable ("nz", "Normals", 2)
calc.SetFunction(f"nx*{cylinder_axis[0]}+ny*{cylinder_axis[1]}+nz*{cylinder_axis[2]}")
calc.SetResultArrayName("dot")
calc.Update()

print('Extract tube')
th = vtk.vtkThreshold()
th.SetInputData(calc.GetOutput())
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "dot")
th.ThresholdBetween(-0.9999, 0.9999)
th.Update()

surf = th.GetOutput()
ptids = surf.GetPointData().GetArray( ids_array )


#loop over tube points
plane_axes = np.array([1,1,1],dtype=int) - cylinder_axis
for i in progressbar(range( surf.GetNumberOfPoints() )):
    pt = surf.GetPoint(i)
    ptid : int = int(ptids.GetTuple1(i))

    x = mpmath.mpf(pt[0])
    y = mpmath.mpf(pt[1])
    z = mpmath.mpf(pt[2])
    r = np.sqrt( x*x*plane_axes[0] + y*y*plane_axes[1] + z*z*plane_axes[2] )

    if plane_axes[0]!=0:
        x=x*R/r

    if plane_axes[1]!=0:
        y=y*R/r
    
    if plane_axes[2]!=0:
        z=z*R/r

    volmesh.GetPoints().SetPoint(ptid, x, y, z)        


print('Saving')
wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputData(volmesh)
wr.SetFileName(volmesh_filename_out)
wr.SetDataModeToAppended()
wr.EncodeAppendedDataOff()
wr.Write()


