import numpy as np   #needs install
import os
import pandas #needs install
import sys
import vtk #needs install
from mpi4py import MPI

file_suffix = ".post.mpio.bin"
ncpus = 10
time_roundoff = 12
points_as_double = True


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
        partitions = np.fromstring(f.read(), dtype=np.int64, sep=' ')

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
    if partition_id<1:
        raise ValueError("partition_id must be > 0")

    if not os.path.isfile(filename):
        raise ValueError(f"File does not exist: {filename}")

    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        #get subtable to calculate cumulative sum of the number of rows per partition
        if partition_id > 1:
            subt = partition_table[partition_table['id']<partition_id]
            cumsum = subt.sum(axis=0)    
        else:
            cumsum = {'Elements':0,'Points':0}


        this_part = partition_table[partition_table['id']==partition_id]

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
            rows2skip = cumsum['Elements']
            rows2read =  this_part.iloc[0,:]['Elements']
        elif header['Association'] == 'node':
            rows2skip = cumsum['Points'] 
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

    for i, values in enumerate(data):
        a.SetTuple(i, values)
        

    return a


def read_mpio_partition_geometry(inputfolder, project_name, partition_id, partition_table):
    
    coords = read_alyampio_array( os.path.join(inputfolder, f'{project_name}-COORD{file_suffix}'), partition_id, part)
    coords = coords['tuples']

    elem = read_alyampio_array( os.path.join(inputfolder, f'{project_name}-LNODS{file_suffix}'), partition_id, part)
    elem = elem['tuples']

    eltype = read_alyampio_array( os.path.join(inputfolder, f'{project_name}-LTYPE{file_suffix}'), partition_id, part)
    eltype = eltype['tuples']


    eltype_alya2vtk = {37:{'Name': vtk.VTK_HEXAHEDRON,'Vertices':8}, 30:{'Name':vtk.VTK_TETRA,'Vertices':4}, \
                       32:{'Name': vtk.VTK_PYRAMID   ,'Vertices':5}, 34:{'Name':vtk.VTK_WEDGE,'Vertices':6},\
                       10:{'Name': vtk.VTK_TRIANGLE  ,'Vertices':3}, 12:{'Name':vtk.VTK_QUAD, 'Vertices':4}}


    pts = vtk.vtkPoints()

    if points_as_double:
        pts.SetDataTypeToDouble()

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


def write_vtk_all_arrays(input_folder, project_name, output_folder, partition_id_list, partition_table, iteration, variable_table):

    af = vtk.vtkAppendFilter()

    variable_pairs = {}
    for partition_id in partition_id_list:
        ug = read_mpio_partition_geometry(input_folder, project_name, partition_id, partition_table)

        #extarct all variables for this timestep
        vars_this_step = variable_table[variable_table['iteration']==iteration]

        if( vars_this_step.shape[0]==0 ):
            return  #nothing to save

        for index, row in vars_this_step.iterrows():
            filename = row['filename']
            varname = row['field']
            
            array = read_alyampio_array(os.path.join(input_folder, filename), partition_id, partition_table)
            
            #if( 'int' in array['header']['DataType'] ):
            #    vtkarray = numpyarray2vtkarray(array['tuples'], varname, 'short')
            #else:
            vtkarray = numpyarray2vtkarray(array['tuples'], varname, 'double')

            if array['header']['Association'] == 'element':
                ug.GetCellData().AddArray( vtkarray )
            elif array['header']['Association'] == 'node':
                ug.GetPointData().AddArray( vtkarray )

            variable_pairs[varname] = array['header']['Association']
    

        af.AddInputData(ug)


    af.Update()

    wr = vtk.vtkXMLUnstructuredGridWriter()
    wr.SetFileName( os.path.join(output_folder, f"{project_name}_{partition_id_list[0]}.vtu" ) )
    wr.SetInputData( af.GetOutput() )
    wr.SetDataModeToAppended()
    wr.EncodeAppendedDataOff()
    wr.Write()  

    return array['header']['Time']


def extract_vars_vtu(vtufilename):
    
    vars = []
    extract = False
    pd = False
    cd = False
    with open(vtufilename, 'rb') as f:
        while 1>0:
            line = f.readline()

            if b"<PointData>" in line:
                extract = True
            elif b"<CellData>" in line:
                extract = True


            if pd and cd:
                break

            if extract:
                vars.append(line)


            if b"</PointData>" in line:
                pd = True
                extract = False
            elif b"</CellData>" in line:
                cd = True
                extract = False


    return vars

#
#
#
#   MAIN CODE
#
#
#
#

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
MPIsize = comm.Get_size()
output_partitions = 2


path = None
name = None
outpath = None
vars = None
joined_partitions = None
iterations = None
part = None

if rank == 0:
    path = sys.argv[1]
    name = sys.argv[2]
    outpath = sys.argv[3]

    vars = generate_variables_table(path, name)
    part = read_partitions(path, name)

    iterations = vars['iteration'].unique()
    partition_ids = part['id'].unique()

    #python's int is platform dependent and does not have int64. To avoid dealing with numpy scatter, I will use string just for scattering
    iterations_str = [str(vv) for vv in iterations]

    joined_partitions = [  partition_ids[i::output_partitions] for i in range(output_partitions) ]
    #print(joined_partitions)
    #print(len(joined_partitions))

    #print("iterations:",iterations_str)
    iterations_split = [  iterations_str[i::MPIsize] for i in range(MPIsize) ] 
    #print("Master iterations_split=",iterations_split)
    iterations = iterations_split
    

path = comm.bcast(path, root=0)
name = comm.bcast(name, root=0)
outpath = comm.bcast(outpath, root=0)
vars = comm.bcast(vars, root=0)
joined_partitions = comm.bcast(joined_partitions, root=0)
part = comm.bcast(part, root=0)


iterations = comm.scatter( iterations, root=0)

iterations = np.array(iterations).astype(np.int64)

#print("Rank =", rank," vars=",iterations)

#
#
#
#   Slaves 
#
#
#

iter_step = {}
for iter in iterations:
    subpath = os.path.join(outpath, f"vtk/{iter:08d}") 
    os.makedirs( subpath, exist_ok=True )
    for part_list in joined_partitions:
        print("Mpi =", rank, ' Iter = ',iter)
        #print('Iter = ',iter,', part = ',part_list)
        time = write_vtk_all_arrays(path, name, subpath, part_list, part, iter, vars)

        if not iter in iter_step:
            iter_step[iter] = {'Time':np.round(time,time_roundoff)} 

        #print(iter_step)

#
#
#
#
#
#
#
iter_step = comm.gather(iter_step, root=0)


if rank == 0:
    #at this point iter_step is a list of dicts. Conver to a single dict
    iter_step = {key: value for d in iter_step for (key, value) in d.items() }

    #print("iter_step= ",iter_step)

    with open(os.path.join(outpath,'vtk',f'{name}.pvd'),'w') as ff:
        ff.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        ff.write('<Collection>\n')

        for iter, value in iter_step.items():

            iter_filename = f"{name}_{iter:08d}.pvtu"
            ff.write(f'   <DataSet timestep="{value["Time"]}" file="{iter_filename}"/>\n')

            with open(os.path.join(outpath,'vtk',iter_filename),'w') as ff_pvtu:
                ff_pvtu.write( '<?xml version="1.0"?>\n' )
                ff_pvtu.write( '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n') 
                ff_pvtu.write( '<PUnstructuredGrid GhostLevel="0">\n' )

                #these vars must be in te same order as in vtu
                allvars = extract_vars_vtu( os.path.join(outpath,'vtk',f"{iter:08d}/{name}_{joined_partitions[0][0]}.vtu") )
                for line in allvars:
                    ff_pvtu.write(line.decode('utf-8').replace('PointData','PPointData').replace('CellData','PCellData'))
                    
                ff_pvtu.write( '<PPoints>\n')

                if points_as_double:
                    ff_pvtu.write( '    <PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>\n')
                else:
                    ff_pvtu.write( '    <PDataArray type="Float32" Name="Points" NumberOfComponents="3"/>\n')

                ff_pvtu.write( '</PPoints>\n')
                for pp in joined_partitions:
                    ff_pvtu.write( f'<Piece Source="{iter:08d}/{name}_{pp[0]}.vtu"/>\n')
                ff_pvtu.write( '</PUnstructuredGrid>\n')
                ff_pvtu.write( '</VTKFile>\n')

        ff.write('</Collection>\n')
        ff.write('</VTKFile>\n')
        

    



