import vtk
import numpy as np

line_pts = np.array([[1,2,3],[10,2,3]]).T #coords in cols
line_length = np.linalg.norm( line_pts[:,1]-line_pts[:,0] )
line_direction = (line_pts[:,1]-line_pts[:,0])/line_length
line_start = line_pts[:,0]
mesh_dist_from_ends = [0, 0]
N_div_len = 30
N_div_angle = 20
line_rad = 1

ll, tt = np.meshgrid( np.linspace(0+mesh_dist_from_ends[0], line_length-mesh_dist_from_ends[1], N_div_len), np.linspace(0,2*np.pi, N_div_angle) )
xx = line_rad*np.cos(tt) 
yy = line_rad*np.sin(tt) 
zz = np.zeros(yy.shape) + ll

#line parameters
v2 = line_direction

v_random = np.roll(v2,1)
v1 = np.cross(v2, v_random)
v1 = v1/np.linalg.norm(v1)
v0 = np.cross(v1, v2)
v0 = v0/np.linalg.norm(v0)
T = np.zeros([3,3])
T[:,0] = v0
T[:,1] = v1
T[:,2] = v2

xyz = T.dot( np.vstack( (xx.ravel(), yy.ravel(), zz.ravel()) ) ) + line_start[:,np.newaxis]

#create triangulation
pts = vtk.vtkPoints()
for t, l in zip( tt.ravel(), ll.ravel() ):
    pts.InsertNextPoint( t, l, 0 )

pd = vtk.vtkPolyData()
pd.SetPoints(pts)
dela = vtk.vtkDelaunay2D()
dela.SetInputData(pd)
dela.Update()


pts = vtk.vtkPoints()
for x,y,z in xyz.T:
    pts.InsertNextPoint(x,y,z)

pd = vtk.vtkPolyData()
pd.SetPoints(pts)
pd.SetPolys( dela.GetOutput().GetPolys() )

wr = vtk.vtkXMLPolyDataWriter()
wr.SetFileName('pts.vtp')
wr.SetInputData(pd)
wr.Write()

