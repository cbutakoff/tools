import transform_lv

meshes = ['mesh_decim.vtk']
npts_array = [500]
sigma_carto_array = [1]

for mesh in meshes:
	for npts in npts_array:
		for sigma_carto in sigma_carto_array:
			t = mesh.split('.')
			outfile_mesh = '{}_carto_n{:04d}_s{:+.2f}.vtk'.format( ''.join(t[:-1]), npts, sigma_carto )
			outfile_tform = '{}_tform_n{:04d}_s{:+.2f}.txt'.format( ''.join(t[:-1]), npts, sigma_carto )

			print 'Generating '+outfile_mesh+'\n'
			transform_lv.random_transform( mesh, outfile_mesh, outfile_tform, npts, sigma_carto )
