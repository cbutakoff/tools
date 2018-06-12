
# coding: utf-8

# In[17]:


#writing ensight Gold binary (comaptible with vtk 8)
#parallelization borowed here https://gist.github.com/fspaolo/51eaf5a20d6d418bd4d0
import numpy as np   #needs install
import os
import pandas #needs install

from mpi4py import MPI #needs install
from queue import Queue

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

try:
    if my_rank==0:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("task_name", help='Name of the alya task')
        parser.add_argument("input_folder", help='Folder with input alyabins')
        parser.add_argument("output_folder", help='Folder for the output ensight case')
        args = parser.parse_args()

        inputfolder = args.input_folder
        project_name = args.task_name
        outputfolder = args.output_folder

	#check if input path exists
        import pathlib
        path = pathlib.Path(inputfolder)
        assert path.exists(), f'{inputfolder} does not exist'

        #create output folders
        path = pathlib.Path(outputfolder)
        path.mkdir(parents=True, exist_ok=True)

        import sys
        sys.stdout.flush()
        print(f'-------------------------------------------------------');
        print(f'Alya task: {project_name}');
        print(f'Input path: {inputfolder}');
        print(f'Output path: {outputfolder}');
        print(f'-------------------------------------------------------');
finally:
    inputfolder = comm.bcast(inputfolder, root=0)
    project_name = comm.bcast(project_name, root=0)
    outputfolder = comm.bcast(outputfolder, root=0)

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
    filename = os.path.join(inputfolder,f'{project_name}-LNODS.post.alyabin');
    with open(filename, 'rb') as f:
        header = read_header(f)
        if header['strings'][6] == '4BYTE':
            alya_id_type = np.int32
        elif header['strings'][6] == '8BYTE':
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


def read_header(file_object):
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

    return {'strings':strings, 'integers':integers, 'reals':reals}


# In[12]:


def read_alya_array(filename, number_of_blocks, datatype):
    with open(filename,'rb') as f:
        header = read_header(f)
        number_of_dimensions = header['integers'][0][0]
        number_of_tuples_total = header['integers'][1][0]
        time_instant_int = header['integers'][3][0]
        time_instant_real = header['reals'][0]
        print(f'Reading array: {number_of_dimensions} dim, {number_of_tuples_total} tuples\n')

        tuples = np.zeros((number_of_tuples_total,number_of_dimensions), dtype=datatype)

        c = 0;
        tuples_per_block = np.zeros(number_of_blocks, dtype=np.int32)
        
        for i in range(number_of_blocks):
            number_of_tuples_in_block = read_one_fp90_record(f, 1, alya_id_type)[0] #stored by alya
            tuples_per_block[i] = number_of_tuples_in_block

            print(f'Block {i}/{number_of_blocks}: {(number_of_tuples_in_block)} tuples\n')
            tuples_temp = read_one_fp90_record(f, number_of_dimensions*number_of_tuples_in_block, datatype)
            
            tuples[c:c+number_of_tuples_in_block, :] =                 np.reshape(tuples_temp, (number_of_tuples_in_block,number_of_dimensions))
            c = c+number_of_tuples_in_block

    return {'tuples':tuples, 'time_real':time_instant_real, 'time_int':time_instant_int, 'tuples_per_block':tuples_per_block};


# In[13]:


def read_alya_variable(variable_name, iteration, number_of_blocks):
    field_filename = os.path.join(inputfolder, '%s-%s-%08d.post.alyabin'% (project_name, variable_name, iteration)) 
    print(field_filename)
    
    field_dtype = alya_id_type
    association = '';
    variabletype = ''; #scalar, vector, ...

    with open(field_filename,'rb') as f:
        header = read_header(f)


    
    if( header['strings'][5] == 'REAL' ):
        field_dtype = np.float64
    if( header['strings'][5] == 'INTEG' ):
        field_dtype = alya_id_type

    if( header['strings'][4] == 'NPOIN' ):
        association = 'node'
    else:
        association = 'element'


    if( header['strings'][3] == 'SCALA' ):
        variabletype = 'scalar'
    elif( header['strings'][3] == 'VECTO' ):
        variabletype = 'vector'
        print('Reading vectors, this has not been tested yet')
    else:
        assert False, "unsupported type of variable"


    if( header['strings'][8] == 'NOFIL' ):
        field_data = read_alya_array(field_filename, number_of_blocks, field_dtype)
    else: 
        assert False, "Filtered types not supported"
        
        
    return {'values':field_data, 'association':association, 'variabletype':variabletype}
 


# In[14]:


def write_geometry(number_of_blocks):
    point_coordinates = read_alya_array(os.path.join(inputfolder,f'{project_name}-COORD.post.alyabin'), \
                                        number_of_blocks, np.float64)
    element_types = read_alya_array(os.path.join(inputfolder,f'{project_name}-LTYPE.post.alyabin'),  \
                                    number_of_blocks, alya_id_type)
    #Read connectivity (indices inside start with 1)
    connectivity = read_alya_array(os.path.join(inputfolder,f'{project_name}-LNODS.post.alyabin'),    \
                                   number_of_blocks, alya_id_type)


    #np.savetxt( 'connectivity.txt', connectivity['tuples'].astype(np.int32), fmt='%d' )
    #np.savetxt( 'inverse.txt', inverse_pt_correspondence.astype(np.int32), fmt='%d' )

    #elements have ids local to each block, tranform them to global ids
    a = connectivity['tuples_per_block'][0]
    npts =  point_coordinates['tuples_per_block'][0]
    for i in range(1,connectivity['tuples_per_block'].shape[0]): #for each block, skip 0
        b = a + connectivity['tuples_per_block'][i]
        connectivity['tuples'][a:b,:] = connectivity['tuples'][a:b,:] + npts
        a = b
        npts = npts + point_coordinates['tuples_per_block'][i]

    print("Connectivty dimensions:", connectivity['tuples'].shape)
        

    #assume all elements are the same
#    element_type = b'hexa8';
#    if element_types['tuples'][0] == 37: #alya hex08 element
#        element_type = b'hexa8';
#    elif element_types['tuples'][0] == 30: #alya tet04 element
#        element_type = b'tetra4';
#    elif element_types['tuples'][0] == 32: #alya pyr05 element
#        element_type = b'pyramid5';
#    elif element_types['tuples'][0] == 34: #alya pen06 element
#        element_type = b'penta6';
    
    #Ensight groups elements by type. Create grouping
    element_alya2ensi = {37:{'Name':b'hexa8','Vertices':8}, 30:{'Name':b'tetra4','Vertices':4}, \
                         32:{'Name':b'pyramid5','Vertices':5}, 34:{'Name':b'penta6','Vertices':6}}
    

    #print(f'Id {connectivity["tuples"][0,0]} transforms to {inverse_pt_correspondence[connectivity["tuples"][0,:]]}'  )
    #print('Element 0, original connectivity: ', connectivity['tuples'][0,:])
    #print('Element 0, original points: ', point_coordinates['tuples'][connectivity['tuples'][0,:],:])
    
    point_coordinates1 = point_coordinates['tuples']

    #print('Point 0: ', point_coordinates1[0,:])
    #print('0 maps to ',inverse_pt_correspondence[0])

    point_coordinates2 = np.zeros( (inverse_pt_correspondence.max()+1,3), dtype=ensight_float_type)
    point_coordinates2[inverse_pt_correspondence,:] = point_coordinates1 
    #print('Point ',inverse_pt_correspondence[0], ', ', point_coordinates2[inverse_pt_correspondence[0],:])
    
    point_coordinates = point_coordinates2 

    connectivity = connectivity['tuples'];
    for i in range(connectivity.shape[1]):
        #here -1 to transform to python array, and +1 to ensight array indexing
        connectivity[:,i] = inverse_pt_correspondence[connectivity[:,i]-1]+1                 
    
    print('Element 0, transformed connectivity: ', connectivity[0,:])
    print('Element 0, transformed points: ', point_coordinates[connectivity[0,:],:])

    
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
        f.write( point_coordinates[:,0].ravel().astype(ensight_float_type) )  #x coord
        f.write( point_coordinates[:,1].ravel().astype(ensight_float_type) )  #y coord
        f.write( point_coordinates[:,2].ravel().astype(ensight_float_type) )  #z coord

        for elem_alya_id, elem_ensi_id in element_alya2ensi.items():        
            print("Saving elements ", elem_alya_id, " as ", elem_ensi_id)

            element_locations = np.where( element_types['tuples']==elem_alya_id )[0] #returns 2 sets, take first
        
        
            elements = connectivity[element_locations,0:elem_ensi_id['Vertices']]
            #print("Locations: ", element_locations)

        
            number_of_elements = elements.shape[0]
            print("Number of elements = ",number_of_elements)


            f.write(elem_ensi_id['Name'].ljust(80))  #tetra4 or hexa8
            f.write( np.array([number_of_elements], dtype=ensight_id_type) )   #int
            f.write( np.array( element_locations, dtype=ensight_id_type)+1 )
            f.write( elements.ravel().astype(ensight_id_type) )

        



# In[15]:

def write_material(number_of_blocks):
    #this is per element variable
    materials_file = os.path.join(inputfolder,f'{project_name}-LMATE.post.alyabin')

    if not os.path.isfile( materials_file ) :
        return                
        
    materials = read_alya_array(materials_file, number_of_blocks, alya_id_type)

    print("Writing variable: LMATE")

    
    #variable ensight
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';
    with open( os.path.join(outputfolder, fmt % (project_name, varname, iteration)),'wb') as f:
        f.write(b'description line 1'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=ensight_id_type))   #int
        f.write(b'coordinates'.ljust(80))


        if data['variabletype']=='scalar':            
            data2write = np.zeros(inverse_pt_correspondence.max()+1, dtype = ensight_float_type)
            print('Data22write: ',data2write.shape)
            print('Ravel: ',data['values']['tuples'].ravel().shape)
            print('Corresp:',inverse_pt_correspondence.shape )
            print("Writing variable: ",varname)
            data2write[inverse_pt_correspondence] = data['values']['tuples'].ravel()
            f.write( data2write )  #z coord    
        elif data['variabletype']=='vector':
            #data has coordinates in the order [[x,y,z],[x,y,z],...]
            #expected order of coordinates
            #vx_n1 vx_n2 ... vx_nn nn floats
            #vy_n1 vy_n2 ... vy_nn nn floats
            #vz_n1 vz_n2 ... vz_nn nn floats
            #Rearrange the  matrix
            data2write = np.zeros( [inverse_pt_correspondence.max()+1, 3], dtype = ensight_float_type)
            data2write[inverse_pt_correspondence,:] = data['values']['tuples']
            f.write( data2write.ravel(order='F').astype(ensight_float_type) )  #z coord    
        else:
            assert False, f"Unknown varibale type: {data['variabletype']}"
        
        
        
        
    return {'time_real':data['values']['time_real'], 'time_int':data['values']['time_int'],             'variable_type':data['variabletype'], 'variable_association':data['association']}





def write_variable_pernode(varname, iteration, number_of_blocks):
    print("Writing variable: ",varname)
    data = read_alya_variable(varname, iteration, number_of_blocks)

    
    #variable ensight
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';
    with open( os.path.join(outputfolder, fmt % (project_name, varname, iteration)),'wb') as f:
        f.write(b'description line 1'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=ensight_id_type))   #int
        f.write(b'coordinates'.ljust(80))


        if data['variabletype']=='scalar':            
            data2write = np.zeros(inverse_pt_correspondence.max()+1, dtype = ensight_float_type)
            print('Data22write: ',data2write.shape)
            print('Ravel: ',data['values']['tuples'].ravel().shape)
            print('Corresp:',inverse_pt_correspondence.shape )
            print("Writing variable: ",varname)
            data2write[inverse_pt_correspondence] = data['values']['tuples'].ravel()
            f.write( data2write )  #z coord    
        elif data['variabletype']=='vector':
            #data has coordinates in the order [[x,y,z],[x,y,z],...]
            #expected order of coordinates
            #vx_n1 vx_n2 ... vx_nn nn floats
            #vy_n1 vy_n2 ... vy_nn nn floats
            #vz_n1 vz_n2 ... vz_nn nn floats
            #Rearrange the  matrix
            data2write = np.zeros( [inverse_pt_correspondence.max()+1, 3], dtype = ensight_float_type)
            data2write[inverse_pt_correspondence,:] = data['values']['tuples']
            f.write( data2write.ravel(order='F').astype(ensight_float_type) )  #z coord    
        else:
            assert False, f"Unknown varibale type: {data['variabletype']}"
        
        
        
        
    return {'time_real':data['values']['time_real'], 'time_int':data['values']['time_int'],             'variable_type':data['variabletype'], 'variable_association':data['association']}


# # Main program



# # Read the partitioning info

# In[16]:

print(f'My rank:{my_rank}')
alya_id_type = None

if my_rank == 0:
    #identify id type
    alya_id_type = identify_alya_id_type(project_name)
    print(f'Node 0: Using Alya id type {np.dtype(np.int32).name}')
    

#broadcast ALYA id type to all the nodes
alya_id_type = comm.bcast(alya_id_type, root=0);


if my_rank != 0:
    print(f'Node {my_rank}: Using Alya id type {alya_id_type}')



if my_rank == 0:
    #read the partitioning info
    partition_filename = os.path.join(inputfolder,f'{project_name}.post.alyapar')
    with open(partition_filename) as f:
        partitions = np.fromstring(f.read(), dtype=alya_id_type, sep=' ')

    #describes the mesh partitions
    #partition_id,  NumberOfElementsInPartition,  NumberOfPointsInPartition, NumberOfBoundariesInPartition
    partitions = np.reshape(partitions[1:],(partitions[0],4))
    number_of_blocks = partitions.shape[0]
    

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

    #extract array names and iteration numbers 
    fields = []
    iteration_numbers = []
    for filename in field_filelist:
        s1 = filename.split('-');
        fields = fields + [s1[1]]
        iteration_numbers =  iteration_numbers + [ int(s1[2].split('.')[0]) ] #this will be long in python 3

    variable_info = pandas.DataFrame({'field':fields, 'iteration':iteration_numbers,'filename':field_filelist})
    variable_info['time_int'] = 0
    variable_info['time_real'] = 0
    variable_info['variabletype']=''
    variable_info['association']=''


inverse_pt_correspondence = None
if my_rank == 0:    
    #read correct element arrangement
    LNINV = read_alya_array(os.path.join(inputfolder,f'{project_name}-LNINV.post.alyabin'), \
                                number_of_blocks, alya_id_type)
    inverse_pt_correspondence = (LNINV['tuples']-1).ravel(); #convert ids to python
    
    #verify the point correspondence
    max_id = inverse_pt_correspondence.max();
    pt_ids = np.zeros(max_id+1)
    pt_ids[inverse_pt_correspondence] = 1    
    assert not (LNINV['tuples']<0).any(), "Negative elements in LNINV, wrong mesh"    
    assert (pt_ids>0).all(), "Some points in the mesh do not have a correspondence in the parittions"
    pt_ids = None #free memeory

inverse_pt_correspondence = comm.bcast(inverse_pt_correspondence, root=0)


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


class Work(object):
    def __init__(self, variable_info, number_of_blocks):
        #files are tuples of filename_input, filetype
        #filenames are the
        q = Queue()
        q.put({'filetype':'geometry','number_of_blocks':number_of_blocks})
        for index, row in variable_info.iterrows():
            q.put({'filetype':'variable','name':row.field, \
                   'iteration':row.iteration,'table_index':index,\
                   'number_of_blocks':number_of_blocks, \
                   'table_index':index})
        self.work = q

    def get_next(self):
        if self.work.empty():
            return None
        return self.work.get()

    def get_size(self):
        return self.work.qsize()


# In[ ]:


def do_work(work):
    print(f'Node: Using Alya id type {alya_id_type}')

    info = {}
    if work['filetype'] == 'geometry':
        write_geometry(work['number_of_blocks'])
    elif work['filetype'] == 'variable':
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


def master(comm):
    status = MPI.Status()
    
    # generate work queue
    print('Variables')
    print(variable_info)
    print('number_of_blocks')
    print(number_of_blocks)
    wq = Work(variable_info, number_of_blocks)

    num_procs = min(wq.get_size(), comm.Get_size())
    print(f'Number of processes: {num_procs}')

    # Seed the slaves, send one unit of work to each slave (rank)
    sent_jobs = 0
    for rank in range(1, num_procs):
        work = wq.get_next()
        comm.send(work, dest=rank, tag=WORKTAG)
        sent_jobs = sent_jobs+1
    
    # Loop over getting new work requests until there is no more work to be done
    while True:
        work = wq.get_next()
        if not work: break

        # Receive results from a slave
        result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        sent_jobs = sent_jobs-1
        process_result(result)

        # Send the slave a new work unit
        comm.send(work, dest=status.Get_source(), tag=WORKTAG)
        sent_jobs = sent_jobs+1
    
    print(f'Jobs still running: {sent_jobs}')

    # No more work to be done, receive all outstanding results from slaves
    for rank in range(sent_jobs): 
        print(f'Receiving {rank}')
        result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        print(f'Received {rank}')
        process_result(result)

    print(f'Shutting down')
    # Tell all the slaves to exit by sending an empty message with DIETAG
    for rank in range(1, comm.Get_size()):
        print(f'Shutting down {rank}')
        comm.send(0, dest=rank, tag=DIETAG)


# In[ ]:


def slave(comm):
    my_rank = comm.Get_rank()
    status = MPI.Status()

    print(f'running processor: {my_rank}')

    while True:
        # Receive a message from the master
        print(f'Proc {my_rank}: receiving')
        work = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        print(f'Proc {my_rank}: received')


        print(f'Process {my_rank}, kill signal: {status.Get_tag()==DIETAG}')
        # Check the tag of the received message
        if status.Get_tag() == DIETAG: break 

        # Do the work
        result = do_work(work)

        # Send the result back
        comm.send(result, dest=0, tag=0)
        print(f'Proc {my_rank}: sent')


# In[ ]:

comm.Barrier()


if my_rank == 0:
    master(comm)
else:
    slave(comm)


# # Write ensight case file

# In[ ]:

comm.Barrier()
if my_rank == 0:
    print(variable_info)
    
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
        max_timeline_id = 1
        timelines = {}; #each timeline will be identified only by the number of elements
        c = 0
        for varname in variables:
            df = variable_info[variable_info.field==varname]
            ntimesteps = df.shape[0]
            if ntimesteps in timelines:
                timeline_id = timelines[ntimesteps]['timeline']
            else:
                timeline_id = max_timeline_id
                timelines[ntimesteps] = {'timeline':timeline_id, 'array_loc':c}
                max_timeline_id = max_timeline_id+1

            variable_info_per_variable = variable_info_per_variable + \
                [{'varname':varname, \
                  'data':variable_info[variable_info.field==varname].sort_values(by='iteration'),\
                  'timeline':timeline_id}]
            c=c+1        
    
    
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
        for timeset in timelines:
            var_data = variable_info_per_variable[ timelines[timeset]['array_loc'] ]
            f.write(f'time set: {var_data["timeline"]}\n')
            number_of_timesteps = var_data['data'].shape[0]
            f.write(f'number of steps: {number_of_timesteps}\n')
            f.write(f'filename numbers: \n')
            f.write(str(var_data['data'].time_int.values)[1:-1]+'\n')
            f.write('time values:\n')
            f.write(str(var_data['data'].time_real.values)[1:-1]+'\n')


# In[33]:




