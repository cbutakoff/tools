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
import mpio

#need python 3.6 or newer, mostly because ints are now 64bit
assert sys.version_info >= (3, 6)


# In[163]:

path = '.'
meshname = 'cylinder'
nodes_name = 'nodes'
ext = '.bin'
outputformat = 'alyampio' # 'vtk' or 'alyampio'
output_filename_vtk_boundary = 'boundary.vtk'
output_filename_vtk_volume = 'volume.vtk'
output_filename_info = 'info' #to save material and face names
mpio_problem_name = 'cylinder'

# In[164]:


with open( os.path.join(path, meshname+'.json'),'r') as file:
    info = json.load(file)


# In[165]:


#read nodes
nodesfilename = meshname+"_"+nodes_name+ext
nodes_df = None
with open( os.path.join(path, nodesfilename), 'rb') as f:
    header = f.read(255)
    data = np.fromfile(f, dtype=np.dtype([('id', np.uint64), ('x', np.float64), ('y', np.float64), ('z', np.float64)]))
    nodes_df = pd.DataFrame(data, columns=data.dtype.names)
    nodes_df = nodes_df.set_index('id')
    data = None
    
print('Nodes')
nodes_df.head()


# In[ ]:





# In[166]:





# In[167]:


def ReadVolumeCellBIN( filename, ncells, nodemask ):
        
    print('Reading ', filename)
    with open(filename, 'rb') as f:
        header = f.read(255)
        data = np.fromfile(f, dtype=np.uint64 ) #id, nnodes

        
    celltypes = set()
    
    cell_count : int = 0
    data_count : int = 0
    data_pp_ids = np.empty(shape=(ncells,), dtype=np.uint64)
    data_pp_nodes = np.empty(shape=(ncells,), dtype=object)
    print('Processing ncells = ', ncells)
    while cell_count<ncells:
        n :int = int( data[data_count+1] )
        celltypes.add(n)
        #cell_id : int = int( data[cell_count] )
        #print(cell_count + n)
        #node_ids = data[(cell_count + 2):(cell_count + 2 + n)]
        data_pp_ids[cell_count] = data[data_count]
        data_pp_nodes[cell_count] = data[(data_count + 2):(data_count + 2 + n)] 

        nodemask[ data_pp_nodes[cell_count]  ] = 1

        data_count += 2+n
        cell_count += 1
        
    
    return pd.DataFrame({'nodes':data_pp_nodes}, index=data_pp_ids), celltypes


def ReadBoundaryCellBIN( filename, ncells ):
        
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
        data_count += 3+n
        cell_count += 1
    
    return pd.DataFrame({'nodes':data_pp_nodes, 'volumeCellId':data_pp_volume_ids}, index=data_pp_ids)


# In[168]:


#read solid cells
cells_df = pd.DataFrame()

#mask to see which nodes are used
#all ids start with 1, ignore the nodemask[0]
nodemask = np.zeros( int(nodes_df.index.max())+1, dtype=np.int8 )

material_names = []

celltypes = set()
for solid_idx, solid_info in enumerate( info['Solids'] ):
    filename = meshname + '_' + solid_info['Name'] + ext
    ncells = solid_info['NCells']
    
    df, celltypes_local = ReadVolumeCellBIN( os.path.join(path, filename), ncells, nodemask )
    celltypes = celltypes.union(celltypes_local)

    df['solidID'] = solid_idx+1
    cells_df = cells_df.append(df)
    material_names.append( {'Name':solid_info['Name'], 'Id':solid_idx+1 } )

#free memory
df = None


print(material_names)

#create new id for solid cells
cells_df['newCellId'] = np.arange(1,cells_df.shape[0]+1,dtype=np.uint64)
print( cells_df.head() )

#leave only used nodes and create new ids
nodes_df = nodes_df.loc[ np.where(nodemask>0)[0] ]
nodemask = None
nodes_df['newPointId'] = np.arange(1,nodes_df.shape[0]+1,dtype=np.uint64)

print( nodes_df.head() )


# In[169]:


#read boundaries

bound_df = pd.DataFrame()

boundary_names = []

for boundary_idx, boundary_info in enumerate( info['Boundaries'] ):
    filename = meshname + '_' + boundary_info['Name'] + ext
    ncells = boundary_info['NCells']
    
    df = ReadBoundaryCellBIN( os.path.join(path, filename), ncells )
    df['boundaryID'] = boundary_idx+1
    bound_df = bound_df.append(df)

    boundary_names.append( {'Name':boundary_info['Name'], 'Id':boundary_idx+1 } )

#free memory
df = None

print(boundary_names)

bound_df['newFaceId'] = np.arange(1,bound_df.shape[0]+1,dtype=np.uint64)
bound_df.head()


# In[188]:


#renumber the nodes in the cells
cells_df['nodes'] = cells_df['nodes'].apply( lambda x: [ nodes_df.loc[i,'newPointId'] for i in x] )
bound_df['nodes'] = bound_df['nodes'].apply( lambda x: [ nodes_df.loc[i,'newPointId'] for i in x] )


# In[193]:


#renumber the volume_cell_ids in the boundary cells
bound_df['volumeCellId'] = bound_df['volumeCellId'].apply( lambda x: cells_df.loc[x,'newCellId'] ) 
bound_df.head()


with open(output_filename_info, 'w') as f:
    f.write('Material codes\n')
    f.write('============================\n')
    for m in material_names:
        f.write( f"{m['Name']} - {m['Id']}\n" )

    f.write('\n')
    f.write('============================\n')
    f.write('Volume cell types\n')
    f.write('============================\n')
    celltype_mames = {4:"TETRA04", 5:"PYRAMID05", 6:"WEDGE06", 8:"HEX08"}
    for m in celltypes:
        f.write( f"{celltype_mames[m]}\n" )



    f.write('\n')
    f.write('============================\n')
    f.write('Boundary codes\n')
    f.write('============================\n')
    for m in boundary_names:
        f.write( f"{m['Name']} - {m['Id']}\n" )


    celltypes


if outputformat == 'vtk':
    # ====================================
    #
    #          save points
    #
    # ====================================
    nodes_df.sort_values(by='newPointId', ascending=True, inplace=True)

    pts = vtk.vtkPoints()
    print('Extracting points', flush=True)
    with progressbar.ProgressBar(max_value=nodes_df.shape[0]) as bar:
        for idx, row in enumerate(nodes_df.itertuples(index=False, name='Pandas')):
            pts.InsertNextPoint( getattr(row, "x"), getattr(row, "y"), getattr(row, "z"))
        
    # ====================================
    #
    #          save volume cells
    #
    # ====================================
    cells_df.sort_values(by='newCellId', ascending=True, inplace=True)

    #number of nodes : type
    vtkcelltypes = {4:vtk.VTK_TETRA, 5:vtk.VTK_PYRAMID, 6:vtk.VTK_WEDGE, 8:vtk.VTK_HEXAHEDRON}
    #vtkfacetypes = {3:vtk.VTK_TRIANGLE, 4:vtk.VTK_QUAD}



    cells = vtk.vtkCellArray()
    celltypes = np.zeros(cells_df.shape[0], dtype=np.uint64)

    material_id = vtk.vtkShortArray()
    material_id.SetName('Material')
    material_id.SetNumberOfComponents(1)
    material_id.SetNumberOfTuples(cells_df.shape[0])

    print('Extracting volume cells', flush=True)
    with progressbar.ProgressBar(max_value=cells_df.shape[0]) as bar:
        for idx, row in enumerate(cells_df.itertuples(index=False, name='Pandas')):
            nodes = getattr(row, "nodes")
            cells.InsertNextCell( len(nodes) )
            for node in nodes:
                cells.InsertCellPoint( int(node)-1 )

            celltypes[idx] = vtkcelltypes[len(nodes)]
            material_id.SetTuple1(idx,  int(getattr(row, "solidID")) )
            bar.update(idx)

    ug = vtk.vtkUnstructuredGrid()
    ug.SetPoints(pts)
    ug.SetCells(celltypes, cells)
    ug.GetCellData().AddArray(material_id)

    print('Writing vtk unstructured grid')
    wr = vtk.vtkDataSetWriter()
    wr.SetFileName(output_filename_vtk_volume)
    wr.SetInputData(ug)
    wr.Write()

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

    vol_cell_id = vtk.vtkShortArray()
    vol_cell_id.SetName('VolCellId')
    vol_cell_id.SetNumberOfComponents(1)
    vol_cell_id.SetNumberOfTuples(bound_df.shape[0])

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
            vol_cell_id.SetTuple1(idx,  int(getattr(row, "volumeCellId"))-1 ) #correct for VTK numbering from 0
            bar.update(idx)

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    pd.GetCellData().AddArray(boundary_id)
    pd.GetCellData().AddArray(vol_cell_id)

    print('Writing boundary vtk')
    wr = vtk.vtkPolyDataWriter()
    wr.SetFileName(output_filename_vtk_boundary)
    wr.SetInputData(pd)
    wr.Write()

    # ====================================
    #
    #          ALYA
    #
    # ====================================

elif outputformat=='alyampio':
    # ====================================
    #
    #          save points
    #
    # ====================================
    nodes_df.sort_values(by='newPointId', ascending=True, inplace=True)

    mpio.MPIO_write_matrix( path, mpio_problem_name, nodes_df[['x','y','z']].values, 'COORD', 'NPOIN' )

    # ====================================
    #
    #          save LNODS, LTYPE, LMATER
    #
    # ====================================
    # TET04 30
    # PYR05 32
    # PEN06 34
    # HEX08 37
    alya_elem_types = {4:30, 5:32, 6:34, 8:37}
    cells_df.sort_values(by='newCellId', ascending=True, inplace=True)

    nodes = np.zeros( (cells_df.shape[0], info['Max_nodes_per_solid_cell']), dtype=np.uint64 )
    nodetypes = np.zeros( cells_df.shape[0], dtype=np.uint64 )
    materials = np.zeros( cells_df.shape[0], dtype=np.uint64 )
    print('Extracting volume cells', flush=True)
    with progressbar.ProgressBar(max_value=cells_df.shape[0]) as bar:
        for idx, row in enumerate(cells_df.itertuples(index=False, name='Pandas')):
            nodes_local = getattr(row, "nodes")
            n : int = len(nodes_local)
            nodes[idx, 0:n] = nodes_local
            nodetypes[idx] = alya_elem_types[n]
            materials[idx] = int(getattr(row, "solidID")) 
            bar.update(idx)

    #write LTYPE (Element types) on elems
    mpio.MPIO_write_matrix( path, mpio_problem_name, nodetypes, 'LTYPE', 'NELEM' )
    #write LNODS (Elements connectivity) on elems
    mpio.MPIO_write_matrix( path, mpio_problem_name, nodes, 'LNODS', 'NELEM' )            
    #write LMATE (material) on elems
    mpio.MPIO_write_matrix( path, mpio_problem_name, materials, 'LMATE', 'NELEM' )            



    # ====================================
    #
    #          save LNODB, LELBO
    #
    # ====================================
    bound_df.sort_values(by='newFaceId', ascending=True, inplace=True)

    nodes = np.zeros( (bound_df.shape[0], info['Max_nodes_per_boundary_face']), dtype=np.uint64 )
    boundary_ids = np.zeros( bound_df.shape[0], dtype=np.uint64 )
    boundary_volume_cell_ids = np.zeros( bound_df.shape[0], dtype=np.uint64 )

    print('Extracting boundary faces', flush=True)
    with progressbar.ProgressBar(max_value=bound_df.shape[0]) as bar:
        for idx, row in enumerate(bound_df.itertuples(index=False, name='Pandas')):
            nodes_local = getattr(row, "nodes")
            n : int = len(nodes_local)
            nodes[idx, 0:n] = nodes_local
               

            boundary_ids[idx] =  int(getattr(row, "boundaryID")) 
            boundary_volume_cell_ids[idx] = int(getattr(row, "volumeCellId")) 
            bar.update(idx)

    #write LNODB (boundary faces) vect
    mpio.MPIO_write_matrix( path, mpio_problem_name, nodes, 'LNODB', 'NBOUN' )            

    #write LELBO (Boundary elements) volume cell ids, NBOUN, LNODB, scal
    mpio.MPIO_write_matrix( path, mpio_problem_name, boundary_volume_cell_ids, 'LELBO', 'NBOUN' )            
    
    #write CODBO on NBOUN
    mpio.MPIO_write_matrix( path, mpio_problem_name, boundary_ids, 'CODBO', 'NBOUN' )            
