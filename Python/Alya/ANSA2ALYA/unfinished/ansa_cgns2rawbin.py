# PYTHON script
import os
import numpy as np
from os.path import join
import json
import h5py
import sys

#cython imports
import pyximport; pyximport.install()
import ansa_element_operations as ansaop


def ListRootEntities(hdf_file):
	group = hdf_file['Base']
	entities = []
	for key in group.keys():
		if type(group[key]) is h5py._hl.group.Group:
			g1 = group[key]
			if 'FamilyBC' in g1:
				entities.append({'Name':key, 'Type':'bound'})
			else:
				entities.append({'Name':key, 'Type':'solid'})

	return entities


def ExtractSolidWithBC(hdf_file, solid_name, shift_cellid):
	group = hdf_file['Base'][solid_name]

	print(group[' data'][()])
	tmpdata = np.squeeze( group[' data'][()] )
	npoints = tmpdata[0]
	ncells = tmpdata[1]

	coord = np.zeros( (npoints, 3), dtype=np.float64 )

	coord[:,0] = group['GridCoordinates']['CoordinateX'][' data'][()]
	coord[:,1] = group['GridCoordinates']['CoordinateY'][' data'][()]
	coord[:,2] = group['GridCoordinates']['CoordinateZ'][' data'][()]
	
	cells_old = group['GridElements']['ElementConnectivity'][' data'][()]
	print(cells_old)

	#in cells_raw replace the element type with the number of nodes
	assert ( cells_old.dtype == np.int64 ), f'Expecting nodes with int64 datatype. Datatype used: {cells_old.dtype.name}'

	cells, tmp = ansaop.element_list_CGNS2ANSABIN_int64( cells_old, shift_cellid, ncells )
	cells_old = None

	#iterate over the shells
	faceid_range = np.squeeze( group['GridShells']['ElementRange'][' data'][()] )
	nfaces = faceid_range[1] - faceid_range[0]+1

	#faceid_range[0] - first id, faceid_range[1] - last id (starting with 1). Usually cells are before the faces
	faces_old = group['GridShells']['ElementConnectivity'][' data'][()]
	faces, face_offsets = ansaop.element_list_CGNS2ANSABIN_int64( faces_old, faceid_range[0]-1+shift_cellid, nfaces, True )
	faces_old = None
	print(face_offsets)


def main(argv):
	cgns_filename = argv[1]	
	
	cgns_file = h5py.File(cgns_filename, 'r')

	props = ListRootEntities(cgns_file)
	print('Found properties: ', props)
	
	ExtractSolidWithBC(cgns_file, 'interior', 0)

	exit()

	#=================================================
	# 
	#
	#             MAIN CODE
	#
	info = {}
	
	
	file_flags = 'w'
	if writebinary:
		file_flags = 'wb'
		
	file_ext = '.csv'
	if writebinary:
		file_ext = '.bin'
	

	nodes_filename = name_prefix + '_' + nodes_name + file_ext

	
	#----------------------------------------------------------
	#
	#ge the main properties and split into surfaces and solids
	#
	#----------------------------------------------------------
	
	properties = base.CollectEntities(deck,None,'__PROPERTIES__')
	print('Found properties: ', properties)

	
	boundaries = []
	solids = []
	
	for i in range( len(properties) ):
		card = base.GetEntityCardValues( deck, properties[i], ['Name'])


		if base.GetEntityType(deck, properties[i]) == 'SHELL_PROPERTY' :
			if card['Name'] in properties2extract:
				boundaries.append({'Name':card['Name'],  'Index':i})
		elif  base.GetEntityType(deck, properties[i]) == 'SOLID_PROPERTY' :
			if card['Name'] in properties2extract:
				solids.append({'Name':card['Name'],  'Index':i})

	print('Boundaries: ', boundaries)
	print('Solids: ', solids)

	
	info['Solids'] = 	[]
	info['Boundaries'] = []

	#if len(solids)>1:
	#	print('More than 1 solid. Exitting')
	#	return -1
	utils.MainProgressBarSetVisible(1)

		
		
	#----------------------------------------------------------
	#
	#get all nodes
	#
	#----------------------------------------------------------
	ansa_nodes = base.CollectEntities(deck, None, "NODE", True)
	print('Number of nodes', len(ansa_nodes))
	
	utils.MainProgressBarSetValue(0)
	utils.MainProgressBarSetText( 'Extracting nodes' )
	progress_i =0
	progress_t = max( int(len(ansa_nodes)/100), 1 )
	
	

		
	with open( join(path2store, nodes_filename), file_flags) as file:
		if writebinary:
			file.write( b"NODES: ID (uint64), X, Y, Z (float64) (this line is 255 chars long)".ljust(255) )
		else:
			file.write( "NODES: ID (uint64), X, Y, Z\n" )
			
		for i, node	in enumerate(ansa_nodes):
			if (progress_i > progress_t):
				utils.MainProgressBarSetValue( 100*i/len(ansa_nodes) )
				progress_i = 0
			progress_i = progress_i+1

			vals = base.GetEntityCardValues( deck, node, ['ID','X','Y','Z'])
			if writebinary:
				np.array( [vals['ID']], dtype=np.uint64 ).tofile( file )
				np.array( [vals['X'], vals['Y'], vals['Z']], dtype=np.float64 ).tofile( file )
			else:
				file.write( "%d %e %e %e\n" % (vals['ID'], vals['X'], vals['Y'], vals['Z']) )

	#----------------------------------------------------------
	#
	#get solid cells
	#
	#----------------------------------------------------------
	number_of_nodes_max = 0;
	for solid_idx, solid in enumerate( solids ): #for each solid
		
		cells =  base.CollectEntities(deck, properties[solid['Index']], 'SOLID',True)
		print( solid['Name'], ": number of cells: ", len(cells) )
		info['Solids'].append( {'Name': solid['Name'], 'NCells' : len(cells) } )
		
		
		utils.MainProgressBarSetValue(0)
		utils.MainProgressBarSetText( 'Extracting %s' % (solid['Name']) )
		progress_i =0
		progress_t = max( int(len(cells)/100), 1 )
		
		with open( join( path2store, name_prefix+"_"+solid['Name']+file_ext ), file_flags) as file:
			
			if writebinary:
				file.write(b"SOLID: ID , NNODES, N1, N2, ... (all uint64) (this line is 255 chars long)".ljust(255) )
			else:
				file.write("SOLID: ID , NNODES, N1, N2, ... (all uint64) (this line is 255 chars long)\n" )

			#extract cells
			for cell_idx, cell in enumerate( cells ):
				if (progress_i > progress_t):
					utils.MainProgressBarSetValue( 100*cell_idx/len(cells) )
					progress_i = 0
				progress_i = progress_i+1
				
				nodelist = GetEntityNodes( deck, cell ) 
				if len(nodelist)>number_of_nodes_max:
					number_of_nodes_max = len(nodelist)
				
				if writebinary:
					np.array( [cell._id, len(nodelist)] + nodelist, dtype=np.uint64 ).tofile( file )
				else:
					file.write("%d %d %s\n" % (cell._id, len(nodelist), " ".join(map(str, nodelist))))  
		
	info['Max_nodes_per_solid_cell'] = number_of_nodes_max
	
	
	#----------------------------------------------------------
	#
	#process the shells and fing the solid cells that are connected ot the faces
	#
	#----------------------------------------------------------
	number_of_nodes_max = 0
	for boundary in boundaries: #boundary number in the list of boundaries
		#get all the faces
		faces =  base.CollectEntities(deck, properties[ boundary['Index'] ], 'SHELL', True)
		print( boundary['Name'], ": number of cells: ", len(faces) )
		info['Boundaries'].append( {'Name': boundary['Name'], 'NCells' : len(faces) } )
		
		
		
		with open( join( path2store, name_prefix+"_"+boundary['Name']+file_ext), file_flags ) as file:
			if writebinary:
				file.write(b"BOUNDARY: ID, BOUNDARY_SOLID_CELL_ID, NNODES, N1, N2, ... (all uint64) (this line is 255 chars long)".ljust(255))
			else:
				file.write("BOUNDARY: ID, BOUNDARY_SOLID_CELL_ID, NNODES, N1, N2, ... (all uint64) (this line is 255 chars long)\n")
			
			utils.MainProgressBarSetValue(0)
			utils.MainProgressBarSetText( 'Extracting %s' % (boundary['Name']) )
			progress_i =0
			progress_t = max( int(len(faces)/100), 1 )

			for face in faces:
				if (progress_i > progress_t):
					utils.MainProgressBarSetValue( 100*cell_idx/len(cells) )
					progress_i = 0
				progress_i = progress_i+1

				nodelist = GetEntityNodes( deck, face )
				#print('Boundary Face Id :', face._id,'  Nodes:',nodelist)
				if len(nodelist)>number_of_nodes_max:
					number_of_nodes_max = len(nodelist)

		
				#find the tetra using these nodes
				#tetra sharing at least one of the node ids
				boundary_tets = base.NodesToElements(nodelist,'entities')
				boundary_tet = -1
				for key, values in boundary_tets.items():
					for value in values:
						if base.GetEntityType( deck, value )=="SOLID":
							nodes = GetEntityNodes ( deck, value )
							if len( set( nodelist) - set(nodes ) ) == 0:
								boundary_tet = value
								break;  #there might be several with the same id
				
				#print('Boundary tets:',boundary_tet._id, ' touching face :', face._id)
				if writebinary:
					np.array( [int(face._id), boundary_tet._id, len(nodelist)] + nodelist, dtype=np.uint64 ).tofile( file )
				else:		
					file.write("%d %d %d %s\n" % (int(face._id), boundary_tet._id, len(nodelist), " ".join(map(str, nodelist))))  
				
	info['Max_nodes_per_boundary_face'] = number_of_nodes_max

	with open( join( path2store, name_prefix+'.json'), 'w' ) as file:
		json.dump(info, file, indent=3)

if __name__ == '__main__':
	main(sys.argv)


