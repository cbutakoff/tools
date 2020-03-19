import vtk
import numpy as np
import sys

surface_filename = sys.argv[1]
boundary_array = sys.argv[2]
boundary_id = int(sys.argv[3])


rd = vtk.vtkPolyDataReader()
rd.SetFileName(surface_filename)
rd.Update()


th = vtk.vtkThreshold()
th.SetInputData(rd.GetOutput())
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, boundary_array)
th.ThresholdBetween(boundary_id-0.1, boundary_id+0.1)
th.Update()

sf = vtk.vtkDataSetSurfaceFilter()
sf.SetInputData(th.GetOutput())
sf.Update()

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

print("Center: ",c)
print("R mean: ", np.mean(r))
print("R std: ", np.std(r))
cn = -c
print( "SPACE_&_TIME_FUNCTIONS")
print( "   FUNCTION=...")
print(f"      V * (1- ( (x {cn[0]:+})^2 + (y {cn[1]:+})^2 + (z {cn[2]:+})^2)/({np.mean(r)}^2) )")  
print( "   END_FUNCTION")
print( "END_SPACE_&_TIME_FUNCTIONS")


