import numpy as np
from scipy import sparse
import scipy.sparse.linalg as linalg_sp
import vtk

def CreatePolyData( pts, faces ):
    """
    Creates vtkPolyData from vertices and faces
    
    pts numpy.array: Nx3 array of vertices
    faces numpy.array: Mx3 array of faces

    Return vtkPolyData
    """
    (nv,mv) = pts.shape
    (nf,mf) = faces.shape
    cells = vtk.vtkCellArray()
    for j in range(nf):
        cell = vtk.vtkTriangle()
        cell.GetPointIds().SetNumberOfIds(3)
        cell.GetPointIds().SetId( 0, faces[j,0] )
        cell.GetPointIds().SetId( 1, faces[j,1] )
        cell.GetPointIds().SetId( 2, faces[j,2] )
        cells.InsertNextCell( cell )
    
    
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nv)
    for j in range(nv):
        points.SetPoint( j, pts[j,0], pts[j,1], pts[j,2] )
        
    new_mesh = vtk.vtkPolyData()
    new_mesh.SetPoints( points )
    new_mesh.SetPolys( cells )
    new_mesh.BuildCells()	
    
    return new_mesh

def ImportVTKPoints( pts, mesh ):
    """
    Iserts numpy array of points into vtkPolyData
    
    pts numpy.array: Nx3 array of vertices

    Return updated mesh
    """
    (nv,mv) = pts.shape
    
    
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nv)
    for j in range(nv):
        points.SetPoint( j, pts[j,0], pts[j,1], pts[j,2] )
        
    mesh.SetPoints( points )
    mesh.BuildCells()	
    
    #return mesh


def VTKPoints2PolyData( pts ):
    """
    Transforms numpy point array into  vtkPolyData
    
    pts numpy.array: Nx3 array of vertices

    Return vtkPolyData
    """
    (nv,mv) = pts.shape
    
    
    vertices = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nv)
    for j in range(nv):
        points.SetPoint( j, pts[j,0], pts[j,1], pts[j,2] )
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(j)    
 
    # Create a polydata object
    mesh = vtk.vtkPolyData()
 
    # Set the points and vertices we created as the geometry and topology of the polydata
    mesh.SetPoints(points)
    mesh.SetVerts(vertices)

    mesh.BuildCells()	
    
    return mesh




def ExtractVTKArray( mesh_array ):
    """
    Extracts scalars from the  mesh_array into numpy array

    mesh_data vtkPolyData: mesh
    filename string: filename

    Return array of scalars
    """
    result = np.zeros( mesh_array.GetNumberOfTuples() )
    for i in range( result.size ):
        result[i] = mesh_array.GetValue(i)
        
    return result
    
    
def SaveVTKPolyData( mesh, filename ):
    """
    Saves vtkPolyData

    mesh vtkPolyData: mesh
    filename string: filename

    Return Nothing
    """
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName( filename )
    writer.SetInputData( mesh )
    writer.SetFileTypeToBinary()
    writer.Update()






def ExtractVTKPoints( mesh ):
    """
    Extract points from vtk structures

    mesh vtk object: mesh

    Return the Nx3 numpy.array of the vertices.
    """
    n = mesh.GetNumberOfPoints();
    vertex = np.zeros( ( n, 3 ) )
    for i in range( n ):
        mesh.GetPoint( i, vertex[i,:] )
    
    return vertex    



def ExtractVTKTriFaces( mesh ):
    """
    Extract triangular faces from vtkPolyData

    mesh vtkPolyData: mesh

    Return the Nx3 numpy.array of the faces (make sure there are only triangles).
    """
    m = mesh.GetNumberOfCells()
    faces = np.zeros( (m, 3), dtype=int )
    for i in range( m ):
        ptIDs = vtk.vtkIdList()
        mesh.GetCellPoints( i, ptIDs )
        if ptIDs.GetNumberOfIds()!=3:
            raise Exception("Nontriangular cell!")
            
        faces[i,0] = ptIDs.GetId(0)
        faces[i,1] = ptIDs.GetId(1)
        faces[i,2] = ptIDs.GetId(2)
    
    return faces


def FlattenMesh( vertex, faces ):
    """
    Flattens a mesh (with one hole) by mapping the edge to a unit circle

    vertex 3xN numpy.array: vertices
    faces 3xM numpy.array: faces

    Return the vertices Nx2 numpy.array of the flattening (faces are the same).
    """
    n = vertex.shape[1]
    L = ComputeLaplacian( vertex, faces )
    boundary = np.array( GetSurfaceBoundary( faces ) )

    p = boundary.size
    t = np.linspace(0, 2*np.pi, p+1, endpoint=False)
    x0 = np.cos(t);
    y0 = np.sin(t);

    L = L.tolil()
    L[boundary,:] = 0;
    for i in range(p):
        L[boundary[i],boundary[i]] = 1;

    Rx = np.zeros(n); 
    Rx[boundary] = x0;
    Ry = np.zeros(n); 
    Ry[boundary] = y0;
    L = L.tocsr()

    result = np.zeros( (Rx.size, 2) )
    result[:,0] = linalg_sp.spsolve(L, Rx); #x
    result[:,1] = linalg_sp.spsolve(L, Ry); #y
    return result


def ComputeLaplacian( vertex, faces ):
    """
    Calculates the laplacian of a mesh

    vertex 3xN numpy.array: vertices
    faces 3xM numpy.array: faces

    Return the Laplacian (sparse matrix probably).
    """
    n = vertex.shape[1]
    m = faces.shape[1]
    
    #compute mesh weight matrix
    W = sparse.coo_matrix((n,n))
    for i in np.arange(1,4,1):
        i1 = np.mod(i-1,3)
        i2 = np.mod(i  ,3)
        i3 = np.mod(i+1,3)
        pp = vertex[:,faces[i2,:]] - vertex[:,faces[i1,:]]
        qq = vertex[:,faces[i3,:]] - vertex[:,faces[i1,:]]
        #% normalize the vectors
        pp = pp / np.sqrt(np.sum(pp**2, axis=0))
        qq = qq / np.sqrt(np.sum(qq**2, axis=0))

        #% compute angles
        ang = np.arccos(np.sum(pp*qq, axis=0))
        W = W + sparse.coo_matrix( (1 / np.tan(ang),(faces[i2,:],faces[i3,:])), shape=(n, n) )
        W = W + sparse.coo_matrix( (1 / np.tan(ang),(faces[i3,:],faces[i2,:])), shape=(n, n) )


    #compute laplacian
    d = W.sum(axis=0);
    D = sparse.dia_matrix((d, 0), shape=(n,n) );
    L = D - W;

    return L


def GetSurfaceBoundary( faces ):
    """
    Calculates the boundary of a surface

    faces 3xN numpy.array: faces of the mesh with only one hole

    Return List of point ids.
    """
    n = faces.max()+1
    m = faces.shape[1]

    #get mesh boundary
    A=sparse.lil_matrix((n,n), dtype=int);
    for i in range(m):
        f=faces[:,i];
        A[f[0],f[1]]=A[f[0],f[1]]+1;
        A[f[0],f[2]]=A[f[0],f[2]]+1;
        A[f[2],f[1]]=A[f[2],f[1]]+1;

    A=A+A.transpose();
    A=A.todense()

    for i in range(n):
        u = np.where( A[i,:]==1 )[1];
        if u.size > 0:
            boundary = [i, u[0]];
            break;

            
    s=boundary[1];
    i=1;
    while i<n:
        u=np.where(A[s,:]==1)[1];
        if u.size != 2:
            raise Exception("problem in boundary");
        if u[0]==boundary[i-1]:
            s=u[1];
        else:
            s=u[0];

        if s!=boundary[0]:
            boundary.append(s);
        else:
            break;
        i=i+1;
        
        
    return boundary



def CloseSurface( shape, output ):
    """
    Closes holes in a surface with a flat cover

    shape vtkPolyData: input shape
    output vtkPolyData: output shape

    Return nothing.
    """
    fe = vtk.vtkFeatureEdges()
    fe.SetInputData( shape )
    fe.BoundaryEdgesOn()
    fe.NonManifoldEdgesOff()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputData(fe.GetOutput())
    connect.Update()

    ncontours = connect.GetNumberOfExtractedRegions()

    append = vtk.vtkAppendPolyData()
    append.AddInput(shape)


    for i in range(ncontours):
        connect.AddSpecifiedRegion(i)
        connect.SetExtractionModeToSpecifiedRegions()
        connect.Update()
        edges = connect.GetOutput()
        cover = vtk.vtkPolyData()
        GenerateHoleCover(edges, cover)

        append.AddInput(cover)
        connect.DeleteSpecifiedRegion(i)


    append.Update()


    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(append.GetOutput())
    cleaner.Update()

    output.DeepCopy(cleaner.GetOutput())




def GenerateHoleCover(edges, cover):
    """
    Generates a flat cover for a convex hole defined by edges

    edges vtkPolyData: input edges
    cover vtkPolyData: output cover

    Return nothing.
    """
    #We'll create the building blocks of polydata including data attributes.
    polys = vtk.vtkCellArray()
    scalars = vtk.vtkFloatArray()
    points = vtk.vtkPoints()

    sur_filt = vtk.vtkCleanPolyData()
    sur_filt.SetInputData( edges )
    sur_filt.Update()

    points.DeepCopy(sur_filt.GetOutput().GetPoints())

    #add centroid
    centr = np.zeros(3)
    for i in range( points.GetNumberOfPoints() ):
        pt = np.zeros(3)
        points.GetPoint(i,pt)
        centr = centr + pt

    centr = centr/ points.GetNumberOfPoints()

    cnt_pt = points.InsertNextPoint(centr)

    #add cells
    for i in range( sur_filt.GetOutput().GetNumberOfCells() ):
        cell = sur_filt.GetOutput().GetCell(i)
        polys.InsertNextCell(3)
        polys.InsertCellPoint(cell.GetPointId(0))
        polys.InsertCellPoint(cell.GetPointId(1))
        polys.InsertCellPoint(cnt_pt)


    # We now assign the pieces to the vtkPolyData.
    cover.SetPoints(points)
    cover.SetPolys(polys)
