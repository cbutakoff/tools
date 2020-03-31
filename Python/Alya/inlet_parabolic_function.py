import vtk
import numpy as np
import sys

surface_filename = sys.argv[1]
boundary_array = sys.argv[2]
boundary_id = int(sys.argv[3])
volmesh_filename = sys.argv[4] #volumetric mesh that corresponds to the surface mesh 
volmesh_scale = 1 #if volmesh has a scale different to the surface mesh, put here the scaling factor for the transform
radius_decimals = 2 #how many decimals of the radius to keep in the output

print("Reading surface mesh")
rd = vtk.vtkPolyDataReader()
rd.SetFileName(surface_filename)
rd.Update()

print("Building normals")
normals = vtk.vtkPolyDataNormals()
normals.SetInputData(rd.GetOutput())
normals.SplittingOff ()
normals.ComputeCellNormalsOn ()
normals.ComputePointNormalsOff ()
normals.Update()


print("Extracting inlet")
th = vtk.vtkThreshold()
th.SetInputData(normals.GetOutput())
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, boundary_array)
th.ThresholdBetween(boundary_id-0.1, boundary_id+0.1)
th.Update()

sf = vtk.vtkDataSetSurfaceFilter()
sf.SetInputData(th.GetOutput())
sf.Update()

print("Extracting inlet boundary")
fe = vtk.vtkFeatureEdges()
fe.SetInputData(sf.GetOutput())
fe.BoundaryEdgesOn ()
fe.FeatureEdgesOff ()
fe.NonManifoldEdgesOff ()
fe.ManifoldEdgesOff ()
fe.Update() 

edges = fe.GetOutput()

#get centroid, radius
#extract points
pts = np.zeros((edges.GetNumberOfPoints(), 3))
for i in range(pts.shape[0]):
   pts[i,:] = edges.GetPoint(i)

c = pts.mean(axis=0)
r = pts-c
r = np.sqrt(np.sum(r**2, axis=1))

mean_r = np.round( np.mean(r), decimals=radius_decimals )

print( "=============================================================") 
print("Center: ",c)
print("R mean: ", mean_r)
print("R std: ", np.std(r))
cn = -c
print( "SPACE_&_TIME_FUNCTIONS")
print( "   FUNCTION=...")
print(f"      V * (1- ( (x {cn[0]:+})^2 + (y {cn[1]:+})^2 + (z {cn[2]:+})^2)/({mean_r}^2) )")  
print( "   END_FUNCTION")
print( "END_SPACE_&_TIME_FUNCTIONS")
print( "=============================================================") 

#get a normal
n = np.array(sf.GetOutput().GetCellData().GetNormals().GetTuple(0))
print(f"Inlet normal: {n}")
d = np.mean(r)/10
pt1 = c+n*d
pt2 = c-n*d
print('Candidate points for injector:',pt1, pt2)

#read volmesh
print("Reading vol mesh")
volrdr = vtk.vtkXMLUnstructuredGridReader()
volrdr.SetFileName(volmesh_filename)
volrdr.Update()

tff = vtk.vtkTransform()
tff.Scale(volmesh_scale,volmesh_scale,volmesh_scale)

print("Rescaling vol mesh")
tf = vtk.vtkTransformFilter()
tf.SetInputData(volrdr.GetOutput())
tf.SetTransform(tff)
tf.Update()

print("Building locator")
loc = vtk.vtkCellLocator()
loc.SetDataSet( tf.GetOutput() )
loc.BuildLocator()

print("Finding interior point")
GenCell = vtk.vtkGenericCell()
pcoords = [0,0,0]
weights = [0,0,0,0,0,0,0,0,0,0,0,0]
cellid1 = loc.FindCell ( pt1, 1e-10, GenCell, pcoords, weights) 
cellid2 = loc.FindCell ( pt2, 1e-10, GenCell, pcoords, weights) 
if ((cellid1<0) and (cellid2<0)) or ((cellid1>=0) and (cellid2>=0)):
    print("No location found for injection, something weird with the mesh")
else:
    if cellid1>0:
       print( "=============================================================") 
       print( "Particle injector:") 
       print( "INJECTOR=          XX")
       print(f"   GEOMETRY:     CIRCLE,PARAM= {pt1[0]}, {pt1[1]}, {pt1[2]}, {mean_r}, {n[0]}, {n[1]}, {n[2]}, XX");
       print( "   DISTRIBUTION: UNICA");
       print( "   NPARTICLES:   ASIS");
       print( "   TYPE=           XX");
       print( "   STOCA =         ON");
       print( "   INITIAL_TIME=   XX $s");
       print( "   PERIOD_TIME=    XX $s");
       print( "   FINAL_TIME=     XX $s");
       print( "   VELOCITY : ZERO");
       print( "END_INJECTOR");
       print( "=============================================================") 

    else:
       print( "=============================================================") 
       print( "Particle injector:") 
       print( "INJECTOR=          XX")
       print(f"   GEOMETRY:     CIRCLE,PARAM= {pt2[0]}, {pt2[1]}, {pt2[2]}, {mean_r}, {n[0]}, {n[1]}, {n[2]}, XX");
       print( "   DISTRIBUTION: UNICA");
       print( "   NPARTICLES:   ASIS");
       print( "   TYPE=           XX");
       print( "   STOCA =         ON");
       print( "   INITIAL_TIME=   XX $s");
       print( "   PERIOD_TIME=    XX $s");
       print( "   FINAL_TIME=     XX $s");
       print( "   VELOCITY : ZERO");
       print( "END_INJECTOR");
       print( "=============================================================") 



