# PYTHON script
import os
import ansa
import numpy as np
from os.path import join
import json
from ansa import *

deck = constants.CGNS

def GetEntityNodes(deck, entity):
	nodes = base.CollectEntities(deck, entity,'NODE')
	nodelist = []
	for node in nodes:
		nodelist.append( node._id )
	
	return nodelist

def appendPartListViewItems(listView):
	
	
	properties = base.CollectEntities(deck,None,'__PROPERTIES__')
	print('Found properties: ', properties)

	rn = [
	guitk.constants.BCRenameType_None,
	]


	cols = guitk.BCListViewColumns(listView)
	guitk.BCListViewSaveSortStateAndDisableSorting(listView)

	for i in range( len(properties) ):
		card = base.GetEntityCardValues( deck, properties[i], ['Name'])

		if (base.GetEntityType(deck, properties[i]) == 'SHELL_PROPERTY') or (base.GetEntityType(deck, properties[i]) == 'SOLID_PROPERTY') :
				guitk.BCListViewAddItem(listView, cols, [card['Name']], rn)
	

	guitk.BCListViewRestoreSortState(listView)



def main():
	#print(utils.SelectOpenDir(''))
	currentproject = base.DataBaseName()  
	print(os.path.dirname(currentproject))

	controls = {}

	window = guitk.BCWindowCreate("ANSA2BIN", guitk.constants.BCOnExitDestroy)
	hl = guitk.BCBoxLayoutCreate(window, guitk.constants.BCVertical)
	guitk.BCLabelCreate(hl, "Path to save to")
	lineedit = guitk.BCLineEditPathCreate(hl, guitk.constants.BCHistoryFolders, os.path.dirname(currentproject), guitk.constants.BCHistorySelect, "LineEditPath")

	controls['path'] = lineedit

	hbox = guitk.BCBoxLayoutCreate(hl, guitk.constants.BCHorizontal)
	bgShape = guitk.BCButtonGroupCreate(hbox, "Export format", guitk.constants.BCHorizontal)
	csvRadioButton = guitk.BCRadioButtonCreate(bgShape, "CSV", None, None)
	mpioRadioButton = guitk.BCRadioButtonCreate(bgShape, "BINARY", None, None)
	guitk.BCRadioButtonSetChecked(mpioRadioButton, True)


	guitk.BCLabelCreate(hl, "Alya problem name")
	intButtonLineEdit = guitk.BCButtonLineEditCreate(hl, os.path.basename(currentproject).split('.')[0] )

        	
	listView = guitk.BCListViewCreate(hl, 1, ['Properties to export'], 1)        	
	appendPartListViewItems(listView)
	guitk.BCListViewSetSelectionMode(listView, guitk.constants.BCMulti)
	guitk.BCListViewSelectAll(listView)

	controls['path'] = lineedit
	controls['format'] = bgShape
	controls['problem_name'] = intButtonLineEdit
	controls['prop_list'] = listView
	#BCButtonGroupGetSelectedId(bg)

	button = guitk.BCPushButtonCreate(hl, "Generate", generatePressed, controls)        	
            	
	guitk.BCLineEditPathSetFileDialogTitle(lineedit, "Select the path to store the meshes")
	guitk.BCWindowSetAcceptFunction(window, acceptFunction, None)
	guitk.BCWindowSetOnCloseFunction(window, closeFunction, None)
	guitk.BCShow( window )

def generatePressed(b, controls):
	path2store = guitk.BCLineEditGetText( guitk.BCComboBoxGetLineEdit ( guitk.BCLineEditPathGetCombo( controls['path'] )) ) 
	print( 'Path: ', path2store )

	format_id = guitk.BCButtonGroupGetSelectedId( controls['format'] )
	formats_dict = {0:"CSV", 1:"BINARY"}
	print( 'Format: ', formats_dict[format_id]  )  #0 - csv, 1 - binary
	writebinary = format_id==1

	name_prefix = guitk.BCLineEditGetText( controls['problem_name'] )
	print( 'Alya problem name: ', name_prefix)
        
	properties2extract = []
	guitk.BCListViewForEachItem( controls['prop_list'], guitk.constants.BCIterateSelected, _extractSelection, properties2extract)
	print('Properties to save: ', properties2extract)
	

	#writebinary = True
	#path2store = "/home/costa/Downloads/cyl/mesh"
	#name_prefix = 'cylinder'
	#properties2extract = ['ext_surf','inlet','outlet','layers','interior']
	nodes_name = 'nodes'
	
	
	
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
				utils.MainProgressBarSetValue( int(np.round(100*i/len(ansa_nodes))) )
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
					utils.MainProgressBarSetValue( int(np.round(100*cell_idx/len(cells))) )
					progress_i = 0
				progress_i = progress_i+1
				
				nodelist = GetEntityNodes( deck, cell ) 
				if len(nodelist)>number_of_nodes_max:
					number_of_nodes_max = len(nodelist)
				
				#vv = base.GetEntityCardValues( deck, cell, ['type'] )
				#print(vv)
				
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
			
			utils.MainProgressBarSetText( 'Extracting %s' % (boundary['Name']) )
			utils.MainProgressBarSetValue(0)
			progress_i =0
			progress_t = max( int(len(faces)/100), 1 )

			for face_idx, face in enumerate(faces):
				
				if (progress_i > progress_t):
					utils.MainProgressBarSetValue( int(np.round(100*face_idx/len(faces))) )
					progress_i = 0
				progress_i = progress_i+1


				nodelist = GetEntityNodes( deck, face )
				#print('Boundary Face Id :', face._id,'  Nodes:',nodelist)
				if len(nodelist)>number_of_nodes_max:
					number_of_nodes_max = len(nodelist)

		
				#find the tetra using these nodes
				#tetra sharing at least one of the node ids
				#this is too slow
#				boundary_tets = base.NodesToElements(nodelist,'entities')
#				boundary_tet = -1
#				for key, values in boundary_tets.items():
#					for value in values:
#						if base.GetEntityType( deck, value )=="SOLID":
#							nodes = GetEntityNodes ( deck, value )
#							if len( set( nodelist) - set(nodes ) ) == 0:
#								boundary_tet = value
#								break;  #there might be several with the same id
				
				
				#print('Boundary tets:',boundary_tet._id, ' touching face :', face._id)
				if writebinary:
#					np.array( [ face._id, boundary_tet._id, len(nodelist)] + nodelist, dtype=np.uint64 ).tofile( file )
					np.array( [ face._id, 0, len(nodelist)] + nodelist, dtype=np.uint64 ).tofile( file )
				else:		
#					file.write("%d %d %d %s\n" % (int(face._id), boundary_tet._id, len(nodelist), " ".join(map(str, nodelist))))  
					file.write("%d %d %d %s\n" % (int(face._id), 0, len(nodelist), " ".join(map(str, nodelist))))  
				
	info['Max_nodes_per_boundary_face'] = number_of_nodes_max

	with open( join( path2store, name_prefix+'.json'), 'w' ) as file:
		json.dump(info, file, indent=3)

	utils.MainProgressBarSetText( "Finished")
	utils.MainProgressBarSetValue(100)



	return 0

def _extractSelection(item_2, selected_props):
	selected_props.append ( guitk.BCListViewItemGetText(item_2,0) )
	return 0

def acceptFunction(window, data):	
	return 0

def closeFunction(window, data):
	utils.MainProgressBarSetText( "")
	utils.MainProgressBarSetValue(0)


if __name__ == '__main__':
		main()


