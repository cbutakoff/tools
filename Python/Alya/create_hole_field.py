# Create the field on nodes of the holes, defined as edges of surfaces with defined codes, 
# This script is for alya cavity volume calculation using divergence or mixed product
# The output also generates a section to paste into sld.dat with the correct hole IDs for each cavity
import vtk
import numpy as np
import sys
from vtk.util.numpy_support import numpy_to_vtk
from progressbar import progressbar

surfmesh_filename = sys.argv[1]
volmesh_filename  = sys.argv[2]
output_filename   = sys.argv[3]
arrayname         = 'BoundaryId'
cavities          = [1,6,11,16] #codes that define the cavity surface
field_id          = 1           #number of the field to store in 
filename          = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format("xxx", field_id, 1)  #field filename

rd = vtk.vtkXMLPolyDataReader()
rd.SetFileName(surfmesh_filename)
rd.Update()
surfmesh = rd.GetOutput()
labels   = surfmesh.GetCellData().GetArray(arrayname)

rdv = vtk.vtkXMLUnstructuredGridReader()
rdv.SetFileName(volmesh_filename)
rdv.Update()
volmesh = rdv.GetOutput()

th = vtk.vtkThreshold()
th.SetInputData(surfmesh)
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, arrayname)
th.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)

loc = vtk.vtkStaticPointLocator()
loc.SetDataSet(volmesh)
loc.BuildLocator()


holes = np.zeros(volmesh.GetNumberOfPoints(), dtype=int)
hole_ids = {}
hole_id = 0

for cavity_label, cavity in enumerate(cavities):
    print("Extract cavity ",cavity)
    th.SetUpperThreshold(cavity+0.1)
    th.SetLowerThreshold(cavity-0.1)
    th.Update()
    
    sf = vtk.vtkDataSetSurfaceFilter()
    sf.SetInputConnection(th.GetOutputPort())
    sf.Update()

    fe = vtk.vtkFeatureEdges()
    fe.SetInputConnection( sf.GetOutputPort() )
    fe.BoundaryEdgesOn()
    fe.NonManifoldEdgesOff()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(fe.GetOutputPort())
    connect.Update()

    ncontours = connect.GetNumberOfExtractedRegions()
    print("Number of holes ",ncontours)


    if cavity not in hole_ids:
        hole_ids[cavity] = []

    for i in progressbar(range(ncontours)):
        connect.InitializeSpecifiedRegionList()
        connect.AddSpecifiedRegion(i)
        connect.SetExtractionModeToSpecifiedRegions()
        connect.Update()

        cl = vtk.vtkCleanPolyData()
        cl.SetInputConnection(connect.GetOutputPort())
        cl.Update()
        contour = cl.GetOutput()
    
        hole_id += 1
        hole_ids[cavity].append( hole_id )
    
        for i in range(contour.GetNumberOfPoints()):
            ptid = loc.FindClosestPoint( contour.GetPoint(i) ) 
            assert ptid>=0, 'Point locator failed'
            holes[ptid] = hole_id


holes_vtk = numpy_to_vtk(holes)
holes_vtk.SetName('holes')
volmesh.GetPointData().AddArray(holes_vtk)

wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputData(volmesh)
wr.SetFileName(output_filename)
wr.EncodeAppendedDataOff()
wr.Write()

for i, (a, b) in enumerate(hole_ids.items()):
    print(f"CAVITY: {i+1}")
    print(f"   BOUNDARY: {a}")
    print( "   VALSET: "," ".join([str(x) for x in b]))
    print(f"   VALFIELD: {field_id}")
    print( "END_CAVITY")

#=====================================================================
def write_vector_mpio(filename, vector, time):
    ndim = vector.shape[1]
    npts = vector.shape[0]

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(b'XFIEL00\0')

        if ndim > 1:
            f.write(b'VECTO00\0')
        else:
            f.write(b'SCALA00\0')

        f.write(b'NPOIN00\0')
        f.write(b'REAL000\0')
        f.write(b'8BYTE00\0')
        f.write(b'SEQUE00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array([0], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array([1], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array([0], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag
        np.array([time], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
        f.write(b'0000000\0')
        f.write(b'NONE000\0') #1
        f.write(b'NONE000\0') #2
        f.write(b'NONE000\0') #3
        f.write(b'NONE000\0') #4
        f.write(b'NONE000\0') #5
        f.write(b'NONE000\0') #6
        f.write(b'NONE000\0') #7
        f.write(b'NONE000\0') #8
        f.write(b'NONE000\0') #9
        f.write(b'NONE000\0') #10
        vector.astype('float64').tofile(f)
#=====================================================================


write_vector_mpio( filename, holes[:,np.newaxis], 0 )

