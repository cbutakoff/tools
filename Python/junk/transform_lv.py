import vtk
import numpy as np
import matplotlib.pyplot as plt
import mytools
import scipy.spatial as sci_spatial

def random_transform( meshin, meshout, trasformfile, Nsamples, sigma_eavm ):
	"""
	Apply a random transform to vtkPolyData and decimate to make it look like carto mesh

	meshin string: filename for the input mesh (make sure it has less than 15k points
	meshout string: filename for the output mesh (transformed and decimated)
	trasformfile string: filename for the output text file to save the transformation matrix
	Nsamples integer: number of vertices to sample from the mesh
	sigma_eavm float: std of gaussian noise to add to carto points (added along mesh normal)

	"""

	#------------------------------
	# parameters
	#------------------------------
	max_translation = 10 #mm
	max_rotation = 1 #rad
	#sigma_eavm = 1 #noise for carto points
	scar_sampling_proportion = 0.5 #proportion of points to be sampled from scar	
	PI_array = "PI" #array that stores pixel intensities of MRI
	




	#------------------------------------
	# read the mesh that will be transformed
	#---------------------------------------
	reader = vtk.vtkPolyDataReader()
	reader.SetFileName( meshin )
	reader.Update()

	trif = vtk.vtkTriangleFilter()
	trif.SetInput( reader.GetOutput() )
	trif.Update()

	mesh = trif.GetOutput()


	#---------------------------------------
	# create scar
	#--------------------------------------
	#read PI array, threshold (40% and 60%) and create scar labels
	
	#these are values for the scar labels, there is a hardcoded dependency below
	healthy_tissue = 0; #scalar value for healthy tissue
	border_zone = 1;#scalar value for borderzone tissue
	core_zone = 2;#scalar value for corezone tissue
		
	scar_array = mytools.ExtractVTKArray( mesh.GetPointData().GetArray( PI_array ) )
	imax = scar_array.max()
	imin = scar_array.min()
	t1 = (imax-imin)*0.4 + imin
	t2 = (imax-imin)*0.6 + imin
	scar_array[ scar_array <= t1 ] = healthy_tissue
	scar_array[ np.all( [scar_array > t1, scar_array <= t2 ], axis=0 ) ] = border_zone
	scar_array[ scar_array > t2 ] = core_zone


	#------------------------------------------------------
	# apply random similarity transform
	#--------------------------------------
	tform = vtk.vtkTransform()
	tform.Translate( (np.random.rand(3)-0.5)*2*max_translation )
	tform.RotateX( (np.random.rand(1)-0.5)*2*max_rotation )
	tform.RotateY( (np.random.rand(1)-0.5)*2*max_rotation )
	tform.RotateZ( (np.random.rand(1)-0.5)*2*max_rotation )
	tform.Update()

	tform_filter = vtk.vtkTransformPolyDataFilter()
	tform_filter.SetTransform(tform)
	tform_filter.SetInput( mesh )
	tform_filter.Update()

	mesh = tform_filter.GetOutput()

	#------------------------------------------------------
	# randomly sample the transformed mesh
	#------------------------------------------------------
	    
	scar_ids = np.where(scar_array>0)[0];
	healthy_ids = np.where(scar_array==0)[0];

	a = np.arange(scar_ids.size)
	np.random.shuffle(a)
	    
	npoints_scar = int(Nsamples * scar_sampling_proportion) #number of points to leave in the scar
	npoints_scar = min( a.size, npoints_scar ) #change number of points if there are not enough
	npoints_healthy = Nsamples-npoints_scar #number of points to leave in healthy area
	print 'Sampling {} in scar and {} everywhere else \n'.format( npoints_scar, npoints_healthy )

	#generate UNIQUE random vertex indices within scar
	scar_pts = np.zeros( ( npoints_scar, 3 ) )


	scar_pt_ids = scar_ids[a[:npoints_scar]]
	for i in range( npoints_scar ):
	    ptid = scar_pt_ids[i]
	    mesh.GetPoint( ptid, scar_pts[i,:] )

	healthy_pts = np.zeros( ( npoints_healthy, 3 ) )
	a = np.arange(healthy_ids.size)
	np.random.shuffle(a)

	npoints_healthy = min( a.size, npoints_healthy ) #change number of points if there are not enough

	healthy_pt_ids = healthy_ids[a[:npoints_healthy]]
	for i in range( npoints_healthy ):
	    ptid = healthy_pt_ids[i]
	    mesh.GetPoint( ptid, healthy_pts[i,:] )


	sampled_pts = np.vstack( (scar_pts, healthy_pts) )
	sampled_ids = np.concatenate( (scar_pt_ids, healthy_pt_ids) )
	sampled_ids = sampled_ids.astype(int)

	print 'Sampled {} points in total\n'.format(sampled_pts.shape)

	vertex = mytools.ExtractVTKPoints(mesh).T
	faces = mytools.ExtractVTKTriFaces(mesh).T

	#-----------------------------------------
	# calculate connectivity for the sampled points
	#-------------------------------------------
	print 'this is the slowest part, depends on the number of vertices \n'
	v = mytools.FlattenMesh( vertex, faces )
	print 'slowest part finished \n'

	faces_resampled = sci_spatial.Delaunay( v[sampled_ids,:] )

	#-------------------------------------------
	# add gauss. noise to vertices
	#---------------------------------------------
	normals = mesh.GetPointData().GetNormals()
	nn = np.zeros( ( 3, normals.GetNumberOfTuples()) )
	for i in range( normals.GetNumberOfTuples() ):
	    normals.GetTuple( i, nn[:,i] )

	vertex = vertex + nn*np.random.randn( normals.GetNumberOfTuples() ) * sigma_eavm

	vtkpoints = vtk.vtkPoints()
	vtkpoints.SetNumberOfPoints(sampled_pts.shape[0])

	flatpoints = vtk.vtkPoints()
	flatpoints.SetNumberOfPoints(sampled_pts.shape[0])

	for i in range( sampled_pts.shape[0] ):
	    vtkpoints.SetPoint(i, vertex[:,sampled_ids[i]])

	for i in range( sampled_pts.shape[0] ):
	    flatpoints.SetPoint(i, v[sampled_ids[i],0], v[sampled_ids[i],1], 0)

	    
	cells = vtk.vtkCellArray()
	for i in range( faces_resampled.vertices.shape[0] ):
	    cells.InsertNextCell(3)
	    cells.InsertCellPoint( faces_resampled.vertices[i,0] )
	    cells.InsertCellPoint( faces_resampled.vertices[i,1] )
	    cells.InsertCellPoint( faces_resampled.vertices[i,2] )
	    
	#pix_intens_resampled will contain pixel intensities for the subsampled mesh
	pix_intens = mesh.GetPointData().GetArray( PI_array )
	pix_intens_resampled = vtk.vtkFloatArray() #first store scar+boundary, then healthy
	pix_intens_resampled.SetName( PI_array )
	pix_intens_resampled.SetNumberOfValues( npoints_scar + npoints_healthy )
	for i in range( sampled_pts.shape[0] ):
	    pix_intens_resampled.SetValue(i, pix_intens.GetValue(sampled_ids[i]))

	    
	resampled_mesh = vtk.vtkPolyData()
	resampled_mesh.SetPoints( vtkpoints )
	resampled_mesh.SetPolys( cells )
	resampled_mesh.GetPointData().AddArray( pix_intens_resampled )
	resampled_mesh.Update()

	normalgen = vtk.vtkPolyDataNormals()
	normalgen.SetInput( resampled_mesh )
	normalgen.SplittingOff()
	normalgen.FlipNormalsOn()
	normalgen.Update()

	mytools.SaveVTKPolyData( normalgen.GetOutput(), meshout )

	#save the matrix
	matrix = tform.GetMatrix()
	np_matrix = np.zeros([4,4])
	for i in range(4):
	    for j in range(4):
		np_matrix[i,j] = matrix.GetElement(i,j)
	np.savetxt(trasformfile, np_matrix)

	#save indices (because of the decimation, these are not true)
	#np.savetxt("indices.txt", sampled_ids)

