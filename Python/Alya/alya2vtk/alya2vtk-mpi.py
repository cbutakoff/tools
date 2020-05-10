import numpy as np   #needs install
import os
import pandas #needs install
import sys
from mpi4py import MPI #needs install
from queue import Queue
import vtk #needs install
from progressbar import progressbar

file_suffix = ".post.mpio.bin"
alya_id_type = np.int64

def read_header_mpio(f):
    magic = np.fromfile(f,count=1, dtype=np.int64)[0]
    if magic != 27093:
        print(f'File {filename} does not appear to be alya mpio file')
        
    format = str(f.read(8))
    if not ('MPIAL' in format):
        assert False,f'File {filename} does not appear to be alya mpio file'

    version = str(f.read(8))
    obj = str(f.read(8))
    dimension = str(f.read(8))
    association = str(f.read(8))
    datatype = str(f.read(8))
    datatypelen = str(f.read(8))    
    seq_par = str(f.read(8))
    filt = str(f.read(8))    
    sorting = str(f.read(8))    
    idd = str(f.read(8))    

    if not ('NOID' in idd):
        assert False, f'ID column in {filename} is not supported'

    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,   f'Lost alignment reding {filename}'
    

    columns = np.fromfile(f,count=1,dtype=np.int64)[0]
    lines = np.fromfile(f,count=1,dtype=np.int64)[0]
    timestep_no = np.fromfile(f,count=1,dtype=np.int64)[0]
    nsubdomains = np.fromfile(f,count=1,dtype=np.int64)[0]
    mesh_div = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag1 = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag2 = np.fromfile(f,count=1,dtype=np.int64)[0]
    time = np.fromfile(f,count=1,dtype=np.float64)[0]
    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,f'Lost alignment reding {filename}'

    junk = str(f.read(8))    #1
    junk = str(f.read(8))    #2
    junk = str(f.read(8))    #3
    junk = str(f.read(8))    #4
    junk = str(f.read(8))    #5
    junk = str(f.read(8))    #6
    junk = str(f.read(8))    #7
    junk = str(f.read(8))    #8
    junk = str(f.read(8))    #9
    junk = str(f.read(8))    #10
    
    if 'INT' in datatype:
        dt = 'int'
    elif 'REAL' in datatype:
        dt = 'float'
    else:
        assert False,f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt+'64'
    elif '4' in datatypelen:
        dt = dt+'32'
    else:
        assert False,f'Unsupported data type length {datatypelen}'


    header = {'DataType':dt, 'Lines':lines,'Columns':columns, 'TimeStepNo':timestep_no, 'Time':time, 'NSubdomains':nsubdomains}

    if 'ELEM' in association:
        header['Association'] = 'element'
    elif 'POIN' in association:
        header['Association'] = 'node'
    else:
        assert False,f'Unsupported association: {association}'


    if 'SCALA' in dimension:
        header['VariableType'] = 'scalar'
    elif( 'VECTO' in dimension  ):
        header['VariableType'] = 'vector'
    else:
        assert False, f"unsupported type of variable {variabletype}"

    assert ('NOFIL' in filt), "Filtered fields are not supported"


    return header



def generate_variables_table(inputfolder, project_name):
    #Parse the filelist of the fields
    with open(os.path.join(inputfolder, f'{project_name}.post.alyafil'),'r') as f:
        field_filelist = f.read().splitlines()

    #remove spaces and empty lines
    field_filelist = [x.strip() for x in field_filelist if x.strip()!='']
    new_field_filelist = []

    #extract array names and iteration numbers 
    fields = []
    iteration_numbers = []
    for filename in field_filelist:
        s1 = filename[len(project_name):].split('-');
        fields = fields + [s1[1]]
        iteration_numbers =  iteration_numbers + [ int(s1[2].split('.')[0]) ] #this will be long in python 3
        new_field_filelist = new_field_filelist + [filename]

    #add special files if present
    codnofile = f'{project_name}-CODNO{file_suffix}'
    if os.path.isfile( os.path.join(inputfolder, codnofile) ):
        new_field_filelist = new_field_filelist + [codnofile]
        iteration_numbers =  iteration_numbers + [ 0 ] #saved with the 0th mesh only
        fields = fields + ['CODNO']
        

    variable_info = pandas.DataFrame({'field':fields, 'iteration':iteration_numbers,'filename':new_field_filelist})

    return variable_info


def read_partitions(inputfolder, project_name):
    partitions = None
    #read the partitioning info
    partition_filename = os.path.join(inputfolder,f'{project_name}.post.alyapar')
    with open(partition_filename) as f:
        partitions = np.fromstring(f.read(), dtype=alya_id_type, sep=' ')

    #describes the mesh partitions
    #partition_id,  NumberOfElementsInPartition,  NumberOfPointsInPartition, NumberOfBoundariesInPartition
    partitions = np.reshape(partitions[1:],(partitions[0],4))
    number_of_blocks = partitions.shape[0]

    partitions = pandas.DataFrame(partitions, columns=['id','Elements','Points','Boundaries'])
    return partitions

def read_block_allvars(inputfolder, project_name, iteration, blockid, variable_table):
    vars = variable_table[variable_table[iteration]==iteration]
    if vars.shape[0]==0:
        return

    #unstructured_grid = read_block_mesh(inputfolder, project_name, blockid)
    #read_block_array(inputfolder, project_name, blockid, arrayid)
    #

    #cols:  field  iteration filename
    for index, row in vars.iterrows():
        varname = row['field']
        filename = row['filename']


#partition_id is the number from the table, not index
def read_alyampio_array(filename, partition_id, partition_table):
    if not os.path.isfile(filename):
        raise ValueError(f"File does not exist: {filename}")

    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        #get subtable to calculate cumulative sum of the number of rows per partition
        if partition_id > 1:
            subt = partition_table[partition_table['id']<partition_id]
            cumsum = subt.cumsum(axis=0)    
        else:
            cumsum = pandas.DataFrame({'Elements':[0],'Points':[0]})


        this_part = partition_table[partition_table['id']==partition_id]
        print("this_part = ", this_part.iloc[0,:]['Elements'])

        #skip to the block start depending association

        bytes_per_number = 8
        if '64' in header['DataType'] :
            bytes_per_number = 8
        elif '8' in header['DataType'] :
            bytes_per_number = 4
        else:
            raise ValueError(f"Unsupported number format in {filename}")

        ncols = header['Columns']
        rows2read = 0
        rows2skip = 0


        #Elements  Points 
        if header['Association'] == 'element':
            rows2skip = cumsum['Elements'].values[0]
            rows2read =  this_part.iloc[0,:]['Elements']
        elif header['Association'] == 'node':
            rows2skip = cumsum['Points'].values[0] 
            rows2read =  this_part.iloc[0,:]['Points']
        else:
            raise ValueError(f"Unsupported association in {filename}")


        f.seek(rows2skip*ncols*bytes_per_number, 1)                
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataType']), count=rows2read*ncols ), (rows2read, ncols) )

        return {'tuples':tuples, 'header':header};


def numpyarray2vtkarray(data, name, arraytype):
    #data is 1d (npoints) or 2d (npoints x n) numpy.ndarray
    #name to name the array

    if( data.shape[0] == 0 ):
        raise ValueError("numpyarray2vtkarray: empty data array")

    if arraytype == 'float':
        a = vtk.vtkFloatArray()
    elif arraytype == 'double':
        a = vtk.vtkDoubleArray()
    elif arraytype == 'short':
        a = vtk.vtkShortArray()
    elif arraytype == 'int64':
        a = vtk.vtkTypeInt64Array()
    else:
        raise ValueError("numpyarray2vtkarray: unknown array type")
    
    if len(data.shape)==1: #1d
        a.SetNumberOfComponents(1)
    else:
        a.SetNumberOfComponents(data.shape[1])

    a.SetNumberOfTuples(data.shape[0])
    a.SetName(name)

    print('data = ', data)

    for i, values in enumerate(data):
        a.SetTuple(i, values)
        

    return a


def read_mpio_partition_geometry(inputfolder, project_name, partition_id, partition_table):
    
    coords = read_alyampio_array(f'{project_name}-COORD{file_suffix}', partition_id, part)
    coords = coords['tuples']

    elem = read_alyampio_array(f'{project_name}-LNODS{file_suffix}', partition_id, part)
    elem = elem['tuples']

    eltype = read_alyampio_array(f'{project_name}-LTYPE{file_suffix}', partition_id, part)
    eltype = eltype['tuples']


    eltype_alya2vtk = {37:{'Name': vtk.VTK_HEXAHEDRON,'Vertices':8}, 30:{'Name':vtk.VTK_TETRA,'Vertices':4}, \
                       32:{'Name': vtk.VTK_PYRAMID   ,'Vertices':5}, 34:{'Name':vtk.VTK_WEDGE,'Vertices':6},\
                       10:{'Name': vtk.VTK_TRIANGLE  ,'Vertices':3}, 12:{'Name':vtk.VTK_QUAD, 'Vertices':4}}


    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(coords.shape[0])
    for ptid, xyz in enumerate(coords):
        pts.SetPoint(ptid, xyz)


    element_types = [eltype_alya2vtk[t[0]]['Name'] for t in eltype]  #eltype is like [[1],[2],[3],..]
    element_nvert = [eltype_alya2vtk[t[0]]['Vertices'] for t in eltype]  #eltype is like [[1],[2],[3],..]

    cells = vtk.vtkCellArray()
    for cellid, cellpts in enumerate(elem):
        nnodes = element_nvert[cellid]
        cells.InsertNextCell(nnodes)
        for i in range(nnodes):
            cells.InsertCellPoint(cellpts[i]-1)



    ug = vtk.vtkUnstructuredGrid()
    ug.SetPoints(pts)
    ug.SetCells(element_types, cells)
    
    return ug


def write_vtk_all_arrays(input_folder, project_name, output_folder, partition_id, partition_table, iteration, variable_table):

    ug = read_mpio_partition_geometry(input_folder, project_name, partition_id, partition_table)

    #extarct all variables for this timestep
    vars_this_step = variable_table[variable_table['iteration']==iteration]
    print(vars_this_step)

    if( vars_this_step.shape[0]==0 ):
        return  #nothing to save

    for index, row in vars_this_step.iterrows():
        filename = row['filename']
        varname = row['field']
        
        array = read_alyampio_array(os.path.join(input_folder, filename), partition_id, partition_table)
        
        if( 'int' in array['header']['DataType'] ):
            vtkarray = numpyarray2vtkarray(array['tuples'], varname, 'short')
        else:
            vtkarray = numpyarray2vtkarray(array['tuples'], varname, 'double')

        if array['header']['Association'] == 'element':
            ug.GetCellData().AddArray( vtkarray )
        elif array['header']['Association'] == 'node':
            ug.GetPointData().AddArray( vtkarray )



    wr = vtk.vtkXMLUnstructuredGridWriter()
    wr.SetFileName( os.path.join(output_folder, f"{project_name}_{partition_id}.vtu" ) )
    wr.SetInputData( ug )
    wr.SetDataModeToAppended()
    wr.EncodeAppendedDataOff()
    wr.Write()  



path = sys.argv[1]
name = sys.argv[2]
vars = generate_variables_table(path, name)
part = read_partitions(path, name)

write_vtk_all_arrays(path, name, './', 1, part, 0, vars)



    
