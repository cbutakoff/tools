
# coding: utf-8

# In[17]:


#writing ensight Gold binary (comaptible with vtk 8)
#parallelization borowed here https://gist.github.com/fspaolo/51eaf5a20d6d418bd4d0
import numpy as np   #needs install
import os
import pandas #needs install
import sys
from mpi4py import MPI #needs install
import json


vtk_installed = True
try:
    import vtk #needs install
except:
    vtk_installed = False
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("Failed to load VTK, will not save VTK with double precision for point coordinates")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

#Generally you don't need VTK double precision for point coordinates



WORKTAG = 1
DIETAG = 0

    



comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
my_name = MPI.Get_processor_name()


# In[9]:
#-------------------------------------------------
#
# Parse arguments
#
#-------------------------------------------------

inputfolder = None
project_name = None
outputfolder = None


#If true - read mpio, if false, read alyabin
MPIO = None


try:
    if my_rank==0:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("task_name", help='Name of the alya task')
        parser.add_argument("input_folder", help='Folder with input alyabins')
        parser.add_argument("output_folder", help='Folder for the output ensight case')
        parser.add_argument("-f","--format", help='Format of the data to expect: mpio, alyabin, auto(default)', default = 'auto')
        parser.add_argument("-v","--vtu", action='store_true', required=False, help='Generate VTU with double precision of coordinates?')
        args = parser.parse_args()

        vtk_installed = vtk_installed & args.vtu
        if vtk_installed:
            print("VTU will be generated")
        else:
            print("VTU will NOT be generated")


        inputfolder = args.input_folder
        project_name = args.task_name
        outputfolder = args.output_folder
        
        if args.format == 'alyabin':
            MPIO = False
        elif args.format == 'mpio':
            MPIO = True
        elif args.format == 'auto':
            #check what kind of files are present
            filename_mpio = os.path.join(inputfolder,f'{project_name}-LNODS.post.mpio.bin');
            filename_alyabin = os.path.join(inputfolder,f'{project_name}-LNODS.post.alyabin');
            mpio_file_present = os.path.isfile(filename_mpio) 
            alyabin_file_present = os.path.isfile(filename_alyabin) 
            assert mpio_file_present!=alyabin_file_present, "Found both alyabin and mpio files. Specify format to use in the --format argument"
            
            if mpio_file_present:
                MPIO = True
            else:
                MPIO = False

        else:
            assert False, f'Unsupported format: {args.format}'

	#check if input path exists
        import pathlib
        path = pathlib.Path(inputfolder)
        assert path.exists(), f'{inputfolder} does not exist'

        #create output folders
        path = pathlib.Path(outputfolder)
        path.mkdir(parents=True, exist_ok=True)
        

        print(f'-------------------------------------------------------');
        print(f'Alya task: {project_name}');
        print(f'Input path: {inputfolder}');
        print(f'Output path: {outputfolder}');
        if MPIO:
            print(f'Format: MPIO');
        else:
            print(f'Format: ALYABIN');
        print(f'-------------------------------------------------------');
  
        sys.stdout.flush()
finally:
    inputfolder = comm.bcast(inputfolder, root=0)
    project_name = comm.bcast(project_name, root=0)
    outputfolder = comm.bcast(outputfolder, root=0)
    MPIO = comm.bcast(MPIO, root=0)
    vtk_installed = comm.bcast(vtk_installed, root=0)



if MPIO:
    file_suffix = '.post.mpio.bin'
else:
    file_suffix = '.post.alyabin'


#-------------------------------------------------
#
#  Important variable to adjust if needed
#
#-------------------------------------------------


ensight_id_type = np.int32
ensight_float_type = np.float32
iterationid_number_of_digits = 6  #how many digits to use for the iteration id in variable files

#-------------------------------------------------
#
#  Important variable to adjust if needed
#
#-------------------------------------------------



# # Functions

# In[10]:


def identify_alya_id_type(project_name):
    #read he header where element ids are stored and see if it's int8 or int4
    filename = os.path.join(inputfolder,f'{project_name}-LNODS{file_suffix}');


    with open(filename, 'rb') as f:
        header = read_header(f)
        if '32' in header['DataType']:
            alya_id_type = np.int32
        elif '64' in header['DataType']:
            alya_id_type = np.int64
        else:
            assert False, f'Alya id type {header[6]} is not supported'
    return alya_id_type        



def read_one_fp90_record(file_object, number_of_elements, datatype):
    #fortran stores length with every block, in the beginning and the end
    count_read = 0
    record = []
    while count_read<number_of_elements:
        #in case the record is stored as several blocks join them
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)
        #block_len is in bytes
        block = np.fromfile(file_object, dtype=datatype, count=block_len[0]//np.dtype(datatype).itemsize)
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)    
        count_read = count_read+block_len
        record = record + [block]
        
    return np.concatenate(record)


# In[11]:
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



def read_header(file_object):
    if MPIO:
        return read_header_mpio(file_object)
    else:
        return read_header_alyabin(file_object)


def read_header_alyabin(file_object):
    #a sanity check
    assert hasattr(file_object, 'read'), "read_header: argument is not a file object"
    
    ihead = read_one_fp90_record(file_object, 1, np.int32)    #! Header: 1234, int32
    assert ihead[0] ==1234, "Header is not 1234"
    strings = [];
    integers = [];
    for i in range(0,9):
        strings = strings +[read_one_fp90_record(file_object, 8, np.uint8).tostring().decode().strip()]   
        #read(ii) strings(1) ! AlyaPost, char8bytes
        #read(ii) strings(2) ! Version, char8bytes
        #read(ii) strings(3) ! NAME, char8bytes
        #read(ii) strings(4) ! SCALA/VECTO, char8bytes
        #read(ii) strings(5) ! NELEM/NPOIN/NBOUN, char8bytes
        #read(ii) strings(6) ! INTEG/REAL,char8bytes
        #read(ii) strings(7) ! 4BYTE/8BYTE, char8bytes -- 4/8 byte integers used subsequently for ids/element type
        #read(ii) strings(8) ! SEQUE/PARAL, char8bytes      
        #read(ii) strings(9) ! NOFIL/FILTE, char8bytes
        
    #for i in range(len(strings)):
    #    print(strings[i])
        
        
    for i in range(0,5):
        integers = integers +[read_one_fp90_record(file_object, 1, np.int32)]   



    #for i in range(len(integers)):
    #    print(integers[i])

    
        
    #read(ii) integers(1) ! ??? int32
    #read(ii) integers(2) ! nelem_total, int32
    #read(ii) integers(3) ! # Subdomains, number of parts?
    if( strings[1][0:5] != 'V0001' ):
        integers = integers +[read_one_fp90_record(file_object, 1, np.int32)]   
        #read(ii) integers(4) ! Time step, int32
        
    reals = read_one_fp90_record(file_object, 1, np.float64)
    #read(ii) reals(1)    ! Time, float64

    if( strings[1][0:5] == 'V0001' ):
        integers[3] = int(reals)  #! floor()?

    #return {'strings':strings, 'integers':integers, 'reals':reals}

    number_of_dimensions = integers[0][0]
    number_of_tuples_total = integers[1][0]
    time_instant_int = integers[3][0]
    time_instant_real = reals[0]

    if( strings[5] == 'REAL' ):
        field_dtype = 'float'
    if( strings[5] == 'INTEG' ):
        field_dtype = 'int'

    if( strings[4] == 'NPOIN'):
        association = 'node'
    else:
        association = 'element'

    if strings[6] == '4BYTE':
        field_dtype = field_dtype + '32'
    elif strings[6] == '8BYTE':
        field_dtype = field_dtype + '64'
    else:
        assert False, f'Alya id type {header[6]} is not supported'


    if( strings[3] == 'SCALA' ):
        variabletype = 'scalar'
    elif( strings[3] == 'VECTO' ):
        variabletype = 'vector'
    else:
        assert False, "unsupported type of variable"


    assert ( strings[8] == 'NOFIL'), "Filtered types not supported"

    header = {'DataType':field_dtype, 'Lines':number_of_tuples_total,'Columns':number_of_dimensions, 'TimeStepNo':time_instant_int, \
                'Time':time_instant_real, 'Association': association, 'VariableType': variabletype}

    return header



# In[12]:
def read_alya_array(filename, number_of_blocks):
    if MPIO:
        return read_alyampio_array(filename, number_of_blocks)
    else:
        return read_alyabin_array(filename, number_of_blocks)

def read_alyampio_array(filename, number_of_blocks):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataType']) ), (header['Lines'], header['Columns']) )
        
        return {'tuples':tuples, 'header':header, 'tuples_per_block':[]};


def read_alyabin_array(filename, number_of_blocks):
    with open(filename,'rb') as f:
        header = read_header_alyabin(f)
        number_of_dimensions = header['Columns']
        number_of_tuples_total = header['Lines']
        time_instant_int = header['TimeStepNo']
        time_instant_real = header['Time']
        #print(f'Reading array: {number_of_dimensions} dim, {number_of_tuples_total} tuples\n')

        datatype = np.dtype( header['DataType'] )
        
        tuples = np.zeros((number_of_tuples_total,number_of_dimensions), dtype=datatype)

        c = 0;
        tuples_per_block = np.zeros(number_of_blocks, dtype=np.int32)
        
        for i in range(number_of_blocks):
            number_of_tuples_in_block = read_one_fp90_record(f, 1, alya_id_type)[0] #stored by alya
            tuples_per_block[i] = number_of_tuples_in_block

            #print(f'Block {i}/{number_of_blocks}: {(number_of_tuples_in_block)} tuples\n')
            tuples_temp = read_one_fp90_record(f, number_of_dimensions*number_of_tuples_in_block, datatype)
            
            tuples[c:c+number_of_tuples_in_block, :] =                 np.reshape(tuples_temp, (number_of_tuples_in_block,number_of_dimensions))
            c = c+number_of_tuples_in_block

    return {'tuples':tuples, 'header':header, 'tuples_per_block':tuples_per_block};


# In[13]:


def read_alya_variable(variable_name, iteration, number_of_blocks):
    if iteration>=0:
        field_filename = os.path.join(inputfolder, '%s-%s-%08d%s'% (project_name, variable_name, iteration, file_suffix)) 
    else:
        field_filename = os.path.join(inputfolder, '%s-%s%s'% (project_name, variable_name, file_suffix)) 
    #print(field_filename)
    

    field_data = read_alya_array(field_filename, number_of_blocks)
        
        
    return field_data
 


# In[14]:

def write_geometry_vtk_double(point_coordinates, elements, element_types):

    print('Saving VTK mesh with double precision')

    element_alya2vtk = {37:{'vtkid':vtk.VTK_HEXAHEDRON,'Vertices':8}, 30:{'vtkid':vtk.VTK_TETRA,'Vertices':4}, \
                         32:{'vtkid':vtk.VTK_PYRAMID,'Vertices':5}, 34:{'vtkid':vtk.VTK_WEDGE,'Vertices':6},\
                         10:{'vtkid':vtk.VTK_TRIANGLE ,'Vertices':3}, 12:{'vtkid':vtk.VTK_QUAD,'Vertices':4}}



    pts = vtk.vtkPoints()
    pts.SetDataTypeToDouble()
    pts.SetNumberOfPoints( point_coordinates.shape[0] )
    
    if point_coordinates.shape[1]==2:
        for i in range(point_coordinates.shape[0]):
            pts.SetPoint(i, point_coordinates[i,0], point_coordinates[i,1], 0)
    else:
        for i in range(point_coordinates.shape[0]):
            pts.SetPoint(i, point_coordinates[i,:])

    elems = vtk.vtkCellArray()
    celltypes = []
    for i in range(elements.shape[0]):
        eltype = element_types[i]
        el = elements[ i, 0:element_alya2vtk[ eltype ]['Vertices'] ]
        celltypes.append( element_alya2vtk[ eltype ]['vtkid'] )
        elems.InsertNextCell( el.shape[0] )
        for j in range( el.shape[0] ):
            elems.InsertCellPoint( el[j]-1 )

    ug = vtk.vtkUnstructuredGrid()
    ug.SetPoints(pts)
    ug.SetCells( celltypes, elems )
    
    wr = vtk.vtkXMLUnstructuredGridWriter()
    wr.SetFileName( os.path.join(outputfolder, f'{project_name}.vtu') )
    wr.SetInputData(ug)
    wr.SetDataModeToAppended()
    wr.EncodeAppendedDataOff()	
    wr.Write()



def write_geometry(number_of_blocks):
    point_coordinates = read_alya_array(os.path.join(inputfolder,f'{project_name}-COORD{file_suffix}'), \
                                        number_of_blocks)
    element_types = read_alya_array(os.path.join(inputfolder,f'{project_name}-LTYPE{file_suffix}'),  \
                                    number_of_blocks)
    #Read connectivity (indices inside start with 1)
    connectivity = read_alya_array(os.path.join(inputfolder,f'{project_name}-LNODS{file_suffix}'),    \
                                   number_of_blocks)


    #np.savetxt( 'connectivity.txt', connectivity['tuples'].astype(np.int32), fmt='%d' )
    #np.savetxt( 'inverse.txt', inverse_pt_correspondence.astype(np.int32), fmt='%d' )

    #elements have ids local to each block, tranform them to global ids
    #a = connectivity['tuples_per_block'][0]
    a = partitions['Elements'][0]
    npts =  partitions['Points'][0]
    for i in range(1, partitions['Elements'].shape[0]): #for each block, skip 0
        b = a + partitions['Elements'][i]
        connectivity['tuples'][a:b,:] = connectivity['tuples'][a:b,:] + npts
        a = b
        npts = npts + partitions['Points'][i]

    #print("Connectivty dimensions:", connectivity['tuples'].shape)
        
    
    #Ensight groups elements by type. Create grouping
    element_alya2ensi = {37:{'Name':b'hexa8','Vertices':8}, 30:{'Name':b'tetra4','Vertices':4}, \
                         32:{'Name':b'pyramid5','Vertices':5}, 34:{'Name':b'penta6','Vertices':6},\
                         10:{'Name':b'tria3','Vertices':3}, 12:{'Name':b'quad4','Vertices':4}}
    

    
    point_coordinates1 = point_coordinates['tuples']

    ndim = point_coordinates1.shape[1]
    point_coordinates2 = np.zeros( (inverse_pt_correspondence.max()+1, ndim ), dtype=ensight_float_type)
    point_coordinates2[inverse_pt_correspondence,:] = point_coordinates1 
    #print('Point ',inverse_pt_correspondence[0], ', ', point_coordinates2[inverse_pt_correspondence[0],:])
    
    point_coordinates = point_coordinates2 

    connectivity = connectivity['tuples'];
    for i in range(connectivity.shape[1]):
        #here -1 to transform to python array, and +1 to ensight array indexing
        connectivity[:,i] = inverse_pt_correspondence[connectivity[:,i]-1]+1                 
    

    #dump vtk mesh witoh double precision, since ENSIGHT supports only single precision
    if(vtk_installed):
        write_geometry_vtk_double( point_coordinates, connectivity, element_types['tuples'].ravel() )

    
    #geometry ensight
    with open(os.path.join(outputfolder,f'{project_name}.ensi.geo'),'wb') as f:
        f.write(b'C Binary'.ljust(80))
        f.write(b'description line 1'.ljust(80))
        f.write(b'description line 2'.ljust(80))
        f.write(b'node id given'.ljust(80))
        f.write(b'element id given'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=ensight_id_type))   #int
        f.write(b'description line 1'.ljust(80))
        f.write(b'coordinates'.ljust(80))

        number_of_points = point_coordinates.shape[0]
        f.write(np.array([number_of_points], dtype=ensight_id_type))   #int
        f.write(np.arange(1,number_of_points+1, dtype=ensight_id_type))

        ndim = point_coordinates.shape[1]

        iii = 0
        while iii<ndim:
            f.write( point_coordinates[:,iii].ravel().astype(ensight_float_type) )  #x coord
            iii = iii+1

        #fill the rest with 0
        while iii<3:
            f.write( 0*point_coordinates[:,0].ravel().astype(ensight_float_type) )  #x coord
            iii = iii+1


        for elem_alya_id, elem_ensi_id in element_alya2ensi.items():        
            #print("Saving elements ", elem_alya_id, " as ", elem_ensi_id)

            element_locations = np.where( element_types['tuples']==elem_alya_id )[0] #returns 2 sets, take first
        
        
            elements = connectivity[element_locations,0:elem_ensi_id['Vertices']]
            #print("Locations: ", element_locations)

        
            number_of_elements = elements.shape[0]
            #print("Number of elements = ",number_of_elements)


            f.write(elem_ensi_id['Name'].ljust(80))  #tetra4 or hexa8
            f.write( np.array([number_of_elements], dtype=ensight_id_type) )   #int
            f.write( np.array( element_locations, dtype=ensight_id_type)+1 )
            f.write( elements.ravel().astype(ensight_id_type) )

        



# In[15]:



def write_variable_percell(varname, iteration, number_of_blocks):
    try:
        data = read_alya_variable(varname, iteration, number_of_blocks)
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('!!! An error occured reading variable ',varname,' iteration ', iteration)
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        sys.stdout.flush()
        return {'time_real':-1, 'time_int':-1,             'variable_type':'FAILED', 'variable_association':'FAILED'}
    
    #variable ensight
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';
    with open( os.path.join(outputfolder, fmt % (project_name, varname, iteration)),'wb') as f:
        f.write(b'description line 1'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=ensight_id_type))   #int

        if 'MATER' in varname:  #this is very particular          
            #data2write = np.zeros(inverse_el_correspondence.max()+1, dtype = ensight_float_type)
            #data2write[inverse_el_correspondence] = data['values']['tuples'].ravel()
            data2write = data['tuples'].ravel()

            for elem_alya_id, elem_ensi_id in element_alya2ensi.items():        
                #print("Saving elements ", elem_alya_id, " as ", elem_ensi_id)

                element_locations = np.where( element_types==elem_alya_id )[0] #returns 2 sets, take first
            
        
                values = data2write[element_locations]
                number_of_values = values.shape[0]


                f.write(elem_ensi_id['Name'].ljust(80))  #tetra4 or hexa8
                f.write( values.ravel().astype(ensight_float_type) )


        else:
            assert False, f'For now only saving MATER is implmeneted, not {varname}'       
        
    return {'time_real':data['header']['Time'], 'time_int':data['header']['TimeStepNo'], \
            'variable_type':data['header']['VariableType'], 'variable_association':data['header']['Association']}


def write_variable_pernode(varname, iteration, number_of_blocks):
    #print("Writing variable: ",varname)
    
    try:
        data = read_alya_variable(varname, iteration, number_of_blocks)
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('!!! An error occured reading variable ',varname,' iteration ', iteration)
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        sys.stdout.flush()
        return {'time_real':-1, 'time_int':-1,             'variable_type':'FAILED', 'variable_association':'FAILED'}
    
    #variable ensight
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';

    if iteration>=0:
        filename = os.path.join(outputfolder, fmt % (project_name, varname, iteration))
    else:
        filename = os.path.join(outputfolder, fmt % (project_name, varname, 0))
        data['header']['Time'] = 0
        data['header']['TimeStepNo'] = 0

    with open( filename,'wb') as f:
        f.write(b'description line 1'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=ensight_id_type))   #int
        f.write(b'coordinates'.ljust(80))


        if data['header']['VariableType']=='scalar':            
            data2write = np.zeros(inverse_pt_correspondence.max()+1, dtype = ensight_float_type)
            #print('Data22write: ',data2write.shape)
            #print('Ravel: ',data['values']['tuples'].ravel().shape)
            #print('Corresp:',inverse_pt_correspondence.shape )
            #print("Writing variable: ",varname)
            data2write[inverse_pt_correspondence] = data['tuples'].ravel()
            f.write( data2write )  #z coord    
        elif data['header']['VariableType']=='vector':
            #data has coordinates in the order [[x,y,z],[x,y,z],...]
            #expected order of coordinates
            #vx_n1 vx_n2 ... vx_nn nn floats
            #vy_n1 vy_n2 ... vy_nn nn floats
            #vz_n1 vz_n2 ... vz_nn nn floats
            #Rearrange the  matrix
            data2write = np.zeros( [inverse_pt_correspondence.max()+1, 3], dtype = ensight_float_type)
            ncomponents = data['tuples'].shape[1] #for 1d, 2d and 3d problems
            data2write[inverse_pt_correspondence,0:ncomponents] = data['tuples']
            f.write( data2write.ravel(order='F').astype(ensight_float_type) )  #z coord    
        else:
            assert False, f"Unknown varibale type: {data['variabletype']}"
        
        
        
        
    return {'time_real':data['header']['Time'], 'time_int':data['header']['TimeStepNo'], \
            'variable_type':data['header']['VariableType'], 'variable_association':data['header']['Association']}


# # Main program

#Ensight groups elements by type. Create grouping
element_alya2ensi = {37:{'Name':b'hexa8','Vertices':8}, 30:{'Name':b'tetra4','Vertices':4}, \
                     32:{'Name':b'pyramid5','Vertices':5}, 34:{'Name':b'penta6','Vertices':6},\
                     10:{'Name':b'tria3','Vertices':3}, 12:{'Name':b'quad4','Vertices':4}}


# # Read the partitioning info

# In[16]:

#print(f'My rank:{my_rank}')
alya_id_type = None

if my_rank == 0:
    #identify id type
    alya_id_type = identify_alya_id_type(project_name)
    #print(f'Node 0: Using Alya id type {np.dtype(np.int32).name}')
    

#broadcast ALYA id type to all the nodes
alya_id_type = comm.bcast(alya_id_type, root=0);


#if my_rank != 0:
#    print(f'Node {my_rank}: Using Alya id type {alya_id_type}')


partitions = None
if my_rank == 0:
    #read the partitioning info
    partition_filename = os.path.join(inputfolder,f'{project_name}.post.alyapar')
    with open(partition_filename) as f:
        partitions = np.fromstring(f.read(), dtype=alya_id_type, sep=' ')

    #describes the mesh partitions
    #partition_id,  NumberOfElementsInPartition,  NumberOfPointsInPartition, NumberOfBoundariesInPartition
    partitions = np.reshape(partitions[1:],(partitions[0],4))
    number_of_blocks = partitions.shape[0]

    partitions = pandas.DataFrame(partitions, columns=['id','Elements','Points','Boundaries'])

partitions = comm.bcast(partitions, root=0)


# In[64]:


#if my_rank == 0:
    #this script does not handle boundaries yet
    #assert (partitions[:,3]==0).all(),  'this script does not handle boundaries yet'


# # Identify variables

# In[65]:

if my_rank == 0:
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
        iteration_numbers =  iteration_numbers + [ -1 ] #this will be long in python 3
        fields = fields + ['CODNO']
        

    variable_info = pandas.DataFrame({'field':fields, 'iteration':iteration_numbers,'filename':new_field_filelist})
    variable_info['time_int'] = 0
    variable_info['time_real'] = 0
    variable_info['variabletype']=''
    variable_info['association']=''


inverse_pt_correspondence = None
element_types = None
#inverse_el_correspondence = None

if my_rank == 0:    
    #read correct element arrangement
    LNINV = read_alya_array(os.path.join(inputfolder,f'{project_name}-LNINV{file_suffix}'), \
                                number_of_blocks)
    inverse_pt_correspondence = (LNINV['tuples']-1).ravel(); #convert ids to python
    
    #verify the point correspondence
    max_id = inverse_pt_correspondence.max();
    pt_ids = np.zeros(max_id+1)
    pt_ids[inverse_pt_correspondence] = 1    
    assert not (LNINV['tuples']<0).any(), "Negative elements in LNINV, wrong mesh"    
    assert (pt_ids>0).all(), "Some points in the mesh do not have a correspondence in the parittions"
    pt_ids = None #free memeory


#    LEINV = read_alya_array(os.path.join(inputfolder,f'{project_name}-LEINV.post.alyabin'), \
#                          number_of_blocks, alya_id_type)
#    inverse_el_correspondence = (LEINV['tuples']-1).ravel(); #convert ids to python

    LTYPE = read_alya_array(os.path.join(inputfolder,f'{project_name}-LTYPE{file_suffix}'),  \
                                    number_of_blocks)
    element_types =    LTYPE['tuples'].ravel()



inverse_pt_correspondence = comm.bcast(inverse_pt_correspondence, root=0)
element_types = comm.bcast(element_types, root=0)
#inverse_el_correspondence = comm.bcast(inverse_el_correspondence, root=0)


# # Unknown stuff

# In[67]:


#god knows what are these
#LNINV = read_alya_array(os.path.join(inputfolder,f'{project_name}-LNINV.post.alyabin'), \
#                                number_of_blocks, alya_id_type)
#LELCH = read_alya_array(os.path.join(inputfolder,f'{project_name}-LELCH.post.alyabin'), \
#                                number_of_blocks, alya_id_type)
#LEINV = read_alya_array(os.path.join(inputfolder,f'{project_name}-LEINV.post.alyabin'), \
#                                number_of_blocks, alya_id_type)


# # Write geometry and variables

# In[7]:


#blocks are mesh partitions
#
#
#
#for index, row in variable_info.iterrows():
#    info = write_variable_pernode(row.field, row.iteration)
#    variable_info.loc[index, 'time_real'] = info['time_real']
#    variable_info.loc[index, 'time_int']= info['time_int']
#    variable_info.loc[index, 'variabletype'] = info['variable_type']
#    variable_info.loc[index, 'association'] = info['variable_association']
#


# # MPI stuff

# In[ ]:


#
#
#   Functions in the middle of the code. The queue is filled before
#
#


def CreateWork(variable_info, number_of_blocks):
    q = []
    q.append( json.dumps({'filetype':'geometry','number_of_blocks':number_of_blocks}) )

    for index, row in variable_info.iterrows():
        q.append( json.dumps({'filetype':'variable','name':row.field, \
               'iteration':row.iteration,'table_index':index,\
               'number_of_blocks':number_of_blocks, \
               'table_index':index}) )

    return q


# In[ ]:


def do_work(work):
    #print(f'Node: Using Alya id type {alya_id_type}')

    info = {}
    if work['filetype'] == 'geometry':
        write_geometry(work['number_of_blocks'])
    elif work['filetype'] == 'variable':
        if 'MATER' in work['name']:
           info = write_variable_percell(work['name'], work['iteration'], work['number_of_blocks'])
        else:
           info = write_variable_pernode(work['name'], work['iteration'], work['number_of_blocks'])

        info['table_index'] = work['table_index'];       
    else:
        assert False, f'Unsupported file type {work["filetype"]}'
    
    info['filetype'] = work['filetype']
    return info


def process_result(result):
    if result['filetype'] == 'variable':
        index = result['table_index']
        variable_info.loc[index, 'time_real'] = result['time_real']
        variable_info.loc[index, 'time_int']= result['time_int']
        variable_info.loc[index, 'variabletype'] = result['variable_type']
        variable_info.loc[index, 'association'] = result['variable_association']

    return 


# In[ ]:

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def chunk(xs,n):
    return [xs[index::n] for index in range(n)]

#
#
#   End of functions in the middle of the code. The queue is filled before
#
#



#=================================================================
#
#  Scatter the queue and process
# In[ ]:

comm.Barrier()

queue2scatter = None
can_continue = 1

if my_rank == 0:
    # generate work queue
    print('Variables')
    print(variable_info)
    print('Number of blocks ', number_of_blocks)
    queue = CreateWork(variable_info, number_of_blocks)

    num_procs = comm.Get_size()
    queue2scatter = chunk(queue, num_procs)

    if len(queue)<num_procs:
        print(f'Number of conversion tasks {len(queue)} smaller than number of CPUs {comm.Get_size()}. Reduce number of CPUs in mpirun.')
        can_continue = 0

    sys.stdout.flush()
    
can_continue = comm.bcast(can_continue, root=0)

if can_continue>0:
    if my_rank == 0:
        print("Scattering work")

    queue2scatter = comm.scatter(queue2scatter, root=0)


    results = []
    for iwork, work in enumerate(queue2scatter):
        if my_rank == 0:
            print( f"Master processing {iwork}/{len(queue2scatter)}" )

        results.append( json.dumps( do_work( json.loads(work) ), cls=NumpyEncoder ) )


    if my_rank == 0:
        print("Gathering results")

    results = comm.gather(results, root=0)


    if my_rank == 0:
        for resultset in results:
            for result in resultset:
                process_result( json.loads(result) )





    #=================================================================
    #
    # # Write ensight case file

    # In[ ]:

    comm.Barrier()
    if my_rank == 0:
        print(variable_info)
        print('Saving case')
        sys.stdout.flush()
        
        case_file = f'{project_name}.ensi.case'
        with open(os.path.join(outputfolder, case_file), 'w') as f:
            f.write('# Converted from Alya\n')
            f.write('# Ensight Gold Format\n')
            f.write('#\n')
            f.write(f'# Problem name: {project_name}\n')
            f.write('FORMAT\n')
            f.write('type: ensight gold\n')
            f.write('\n')
            f.write('GEOMETRY\n')
            f.write(f'model: 1 {project_name}.ensi.geo\n')
            f.write('\n')
            f.write('VARIABLE\n')

            variables = variable_info.field.unique();

            #define dataframe for each variable
            #define timeline for each variable base on the number of timesteps
            variable_info_per_variable = []
            timeline_id = 1
            timelines = {}; #each timeline will be identified only by the number of elements
            c = 0
            for varname in variables:
                df = variable_info[variable_info.field==varname]
                ntimesteps = df.shape[0]

                variable_info_per_variable = variable_info_per_variable + \
                    [{'varname':varname, \
                      'data':variable_info[variable_info.field==varname].sort_values(by='iteration'),\
                      'timeline':timeline_id}]
                c=c+1        
                timeline_id += 1
      
            for var_data in variable_info_per_variable:
                varname = var_data['varname']
                timeline_id = var_data['timeline']
                print(f'Varname {varname}\n')
                df = var_data['data'].iloc[0] #get one record with this varibale
                line = f'{df.variabletype} per {df.association}: {timeline_id} {varname} {project_name}.ensi.{varname}-'+                '*'*iterationid_number_of_digits+'\n'       
                f.write(line)



            #to make sure the whole array gets printed
            #otherwise numpy converts to string a summary, e.g. (1,2,3,...,5,6,7)
            np.set_printoptions(threshold=np.inf)

            f.write('\n')
            f.write('TIME\n')
            for var_data in variable_info_per_variable:
                f.write(f'time set: {var_data["timeline"]}\n')
                number_of_timesteps = var_data['data'].shape[0]
                f.write(f'number of steps: {number_of_timesteps}\n')
                f.write(f'filename numbers: \n')
                f.write(str(var_data['data'].time_int.values)[1:-1]+'\n')
                f.write('time values:\n')
                f.write(str(var_data['data'].time_real.values)[1:-1]+'\n')


# In[33]:




