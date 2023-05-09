#!/usr/bin/env python
# coding: utf-8

# In[215]:


import vtk
import json
import numpy as np
import pandas as pd
import progressbar
import sys
import os
import argparse
import pathlib


#need python 3.6 or newer, mostly because ints are now 64bit
assert sys.version_info >= (3, 6)

supported_formats = ['mpio','vtk']

parser = argparse.ArgumentParser()
parser.add_argument("input_mesh_prefix", help='Prefix of the bin files exported from ANSA')
parser.add_argument("-op", "--output_mesh_prefix", help='Prefix of output files. For MPIO use problem name. Default = input_mesh_prefix')
parser.add_argument("-if", "--input_folder",  help='Folder with input bins (default: %(default)s)', default='.')
parser.add_argument("-of", "--output_folder", help='Folder for the output (default = input_folder)')
parser.add_argument("-s", "--scale",  metavar='N', type=float, help='Scale factor for the mesh (default: %(default)s)', default=1)
parser.add_argument('-c', '--infercode', action='store_true', required=False, help='Infer code from the bundary name if its like: 02_blah-blah (default: %(default)s)', default=False)
args = parser.parse_args()

input_path = args.input_folder
output_path = args.output_folder
if output_path is None:
    output_path = input_path

input_mesh_prefix = args.input_mesh_prefix

output_mesh_prefix = args.output_mesh_prefix
if output_mesh_prefix is None:
    output_mesh_prefix = input_mesh_prefix

redo_surface2vol_ids = False #recalculate vol elements that are connected to the surface
nodes_name = 'nodes'
ext = '.bin'



#print('input_path = ', input_path)
#print('output_path = ', output_path)
#print('input_mesh_prefix = ', input_mesh_prefix)
#print('output_mesh_prefix = ', output_mesh_prefix)
#print('outputformat = ', outputformat)



assert os.path.exists(input_path), f'{input_path} does not exist'

os.makedirs(output_path, exist_ok=True)


output_filename_vtk_boundary = os.path.join(output_path, f'{output_mesh_prefix}-boundary.vtp')
output_filename_info = os.path.join(output_path, f'{output_mesh_prefix}.info') #to save material and face names

input_jsonfile = os.path.join(input_path, input_mesh_prefix+'.json')
imput_nodesfilename = os.path.join(input_path, f"{input_mesh_prefix}_{nodes_name}{ext}")


assert os.path.exists(input_jsonfile), f"Can't find {input_mesh_prefix}.json in {input_path}"
assert os.path.exists(imput_nodesfilename), f"Can't find {imput_nodesfilename}.json in {input_path}"



with open( input_jsonfile, 'r') as file:
    info = json.load(file)


#read nodes
nodes_df = None
with open( imput_nodesfilename, 'rb') as f:
    header = f.read(255)
    data = np.fromfile(f, dtype=np.dtype([('id', np.uint64), ('x', np.float64), ('y', np.float64), ('z', np.float64)]))
    nodes_df = pd.DataFrame(data, columns=data.dtype.names)
    nodes_df = nodes_df.set_index('id')
    data = None
    
print('Nodes')
nodes_df.head()





def ReadBoundaryCellBIN( filename, ncells, nodemask ):
        
    print('Reading ', filename)
    with open(filename, 'rb') as f:
        header = f.read(255)
        data = np.fromfile(f, dtype=np.uint64 ) #id, nnodes

    cell_count : int = 0
    data_count : int = 0
    data_pp_ids = np.empty(shape=(ncells,), dtype=np.uint64)
    data_pp_volume_ids = np.empty(shape=(ncells,), dtype=np.uint64)
    data_pp_nodes = np.empty(shape=(ncells,), dtype=object)
    print('Processing ncells = ', ncells)
    while cell_count<ncells:
        data_pp_ids[cell_count] = data[data_count]
        data_pp_volume_ids [cell_count] = data[data_count+1] #reference to the volume mesh cell
        n :int = int( data[data_count+2] )

        data_pp_nodes[cell_count] = data[(data_count + 3):(data_count + 3 + n)] 

        nodemask[ data_pp_nodes[cell_count]  ] = 1

        data_count += 3+n
        cell_count += 1
    
    return pd.DataFrame({'nodes':data_pp_nodes, 'volumeCellId':data_pp_volume_ids}, index=data_pp_ids)



#mask to see which nodes are used
#all ids start with 1, ignore the nodemask[0]
nodemask = np.zeros( int(nodes_df.index.max())+1, dtype=np.int8 )



#read boundaries

bound_df = pd.DataFrame()

boundary_names = []

#sort the boundary names, so that they always appear in teh same order
info['Boundaries'] = sorted(info['Boundaries'], key = lambda i: i['Name']) 

for boundary_idx, boundary_info in enumerate( info['Boundaries'] ):
    filename = f"{input_mesh_prefix}_{boundary_info['Name']}{ext}"
    ncells = boundary_info['NCells']
    
    df = ReadBoundaryCellBIN( os.path.join(input_path, filename), ncells, nodemask )
    if 'Code' in boundary_info:
        boundary_id = boundary_info['Code']
    else:
        boundary_id = boundary_idx+1


    if args.infercode:
        boundary_id = int(boundary_info['Name'].split("_")[0])

    df['boundaryID'] = boundary_id 
    bound_df = pd.concat( (bound_df,df), axis=0, ignore_index=True )

    boundary_names.append( {'Name':boundary_info['Name'], 'Id':boundary_id } )

#free memory
df = None

print(boundary_names)



#leave only used nodes and create new ids
nodes_df = nodes_df.loc[ np.where(nodemask>0)[0] ]
nodemask = None
nodes_df['newPointId'] = np.arange(1,nodes_df.shape[0]+1,dtype=np.uint64)

print( nodes_df.head() )




bound_df['newFaceId'] = np.arange(1,bound_df.shape[0]+1,dtype=np.uint64)
bound_df.head()



#renumber the nodes in the cells
print("Renumbering boundary nodes")
bound_df['nodes'] = bound_df['nodes'].apply( lambda x: [ nodes_df.loc[i,'newPointId'] for i in x] )




with open(output_filename_info, 'w') as f:
    f.write('============================\n')
    f.write('Stats\n')
    f.write('============================\n')
    f.write(f'Number of nodes: {nodes_df.shape[0]}\n')
    f.write(f'Number of boundary faces: {bound_df.shape[0]}\n')


    f.write('\n')
    f.write('============================\n')
    f.write('Boundary codes\n')
    f.write('============================\n')
    for m in boundary_names:
        f.write( f"{m['Name']} - {m['Id']}\n" )


    


    # ====================================
    #
    #          save points
    #
    # ====================================
    nodes_df.sort_values(by='newPointId', ascending=True, inplace=True)

    pts = vtk.vtkPoints()
    pts.SetDataTypeToDouble()
    print('Extracting points', flush=True)
    with progressbar.ProgressBar(max_value=nodes_df.shape[0]) as bar:
        for idx, row in enumerate(nodes_df.itertuples(index=False, name='Pandas')):
            pts.InsertNextPoint( getattr(row, "x")*args.scale, getattr(row, "y")*args.scale, getattr(row, "z")*args.scale)
        


    # ====================================
    #
    #        process the boundary
    #
    # ================================
    bound_df.sort_values(by='newFaceId', ascending=True, inplace=True)

    boundary_id = vtk.vtkShortArray()
    boundary_id.SetName('BoundaryId')
    boundary_id.SetNumberOfComponents(1)
    boundary_id.SetNumberOfTuples(bound_df.shape[0])


    cells = vtk.vtkCellArray()
    #celltypes = np.zeros(bound_df.shape[0], dtype=np.uint64)

    print('Extracting boundary faces', flush=True)
    with progressbar.ProgressBar(max_value=bound_df.shape[0]) as bar:
        for idx, row in enumerate(bound_df.itertuples(index=False, name='Pandas')):
            nodes = getattr(row, "nodes")
            cells.InsertNextCell( len(nodes) )
            for node in nodes:
                cells.InsertCellPoint( int(node)-1 )

    #        celltypes[idx] = vtkfacetypes[len(nodes)]
            boundary_id.SetTuple1(idx,  int(getattr(row, "boundaryID")) )
            bar.update(idx)

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    pd.GetCellData().AddArray(boundary_id)

    print('Writing boundary vtk')
    wr = vtk.vtkXMLPolyDataWriter()
    wr.SetFileName(output_filename_vtk_boundary)
    wr.SetInputData(pd)
    wr.SetDataModeToBinary()
    wr.EncodeAppendedDataOff()
    wr.Write()
