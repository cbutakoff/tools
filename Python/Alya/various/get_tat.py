import numpy as np  # needs install
import os
import pandas  # needs install
import sys
import json


path = sys.argv[1]
json_file = 'tat.json'


file_suffix_mpio = '.post.mpio.bin'
file_suffix_alya = '.post.alyabin'
file_isoch = 'ISOCH'  # filename must contain this


MPIO = True


def read_one_fp90_record(file_object, number_of_elements, datatype):
    # fortran stores length with every block, in the beginning and the end
    count_read = 0
    record = []
    while count_read < number_of_elements:
        # in case the record is stored as several blocks join them
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)
        # block_len is in bytes
        block = np.fromfile(file_object, dtype=datatype,
                            count=block_len[0]//np.dtype(datatype).itemsize)
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)
        count_read = count_read+block_len
        record = record + [block]

    return np.concatenate(record)

def read_header_mpio(f):
    magic = np.fromfile(f, count=1, dtype=np.int64)[0]
    if magic != 27093:
        print(f'File does not appear to be alya mpio file')

    format = str(f.read(8))
    if not ('MPIAL' in format):
        assert False, f'File does not appear to be alya mpio file'

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

    columns = np.fromfile(f, count=1, dtype=np.int64)[0]
    lines = np.fromfile(f, count=1, dtype=np.int64)[0]
    timestep_no = np.fromfile(f, count=1, dtype=np.int64)[0]
    nsubdomains = np.fromfile(f, count=1, dtype=np.int64)[0]
    mesh_div = np.fromfile(f, count=1, dtype=np.int64)[0]
    tag1 = np.fromfile(f, count=1, dtype=np.int64)[0]
    tag2 = np.fromfile(f, count=1, dtype=np.int64)[0]
    time = np.fromfile(f, count=1, dtype=np.float64)[0]

    junk = str(f.read(8))
    if not ('0000000' in junk):
        assert False, f'Lost alignment reding {filename}'

    junk = str(f.read(8))  # 1
    junk = str(f.read(8))  # 2
    junk = str(f.read(8))  # 3
    junk = str(f.read(8))  # 4
    junk = str(f.read(8))  # 5
    junk = str(f.read(8))  # 6
    junk = str(f.read(8))  # 7
    junk = str(f.read(8))  # 8
    junk = str(f.read(8))  # 9
    junk = str(f.read(8))  # 10

    if 'INT' in datatype:
        dt = 'int'
    elif 'REAL' in datatype:
        dt = 'float'
    else:
        assert False, f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt+'64'
    elif '4' in datatypelen:
        dt = dt+'32'
    else:
        assert False, f'Unsupported data type length {datatypelen}'

    header = {'DataType': dt, 'Lines': lines, 'Columns': columns,
              'TimeStepNo': timestep_no, 'Time': time, 'NSubdomains': nsubdomains}

    if 'ELEM' in association:
        header['Association'] = 'element'
    elif 'POIN' in association:
        header['Association'] = 'node'
    else:
        assert False, f'Unsupported association: {association}'

    if 'SCALA' in dimension:
        header['VariableType'] = 'scalar'
    elif('VECTO' in dimension):
        header['VariableType'] = 'vector'
    else:
        assert False, f"unsupported type of variable {header['VariableType']}"

    assert ('NOFIL' in filt), "Filtered fields are not supported"

    return header


def read_header(file_object):
    if MPIO:
        return read_header_mpio(file_object)
    else:
        return read_header_alyabin(file_object)


def read_header_alyabin(file_object):
    # a sanity check
    assert hasattr(
        file_object, 'read'), "read_header: argument is not a file object"

    ihead = read_one_fp90_record(
        file_object, 1, np.int32)  # ! Header: 1234, int32
    assert ihead[0] == 1234, "Header is not 1234"
    strings = []
    integers = []
    for i in range(0, 9):
        strings = strings + \
            [read_one_fp90_record(
                file_object, 8, np.uint8).tostring().decode().strip()]
        # read(ii) strings(1) ! AlyaPost, char8bytes
        # read(ii) strings(2) ! Version, char8bytes
        # read(ii) strings(3) ! NAME, char8bytes
        # read(ii) strings(4) ! SCALA/VECTO, char8bytes
        # read(ii) strings(5) ! NELEM/NPOIN/NBOUN, char8bytes
        # read(ii) strings(6) ! INTEG/REAL,char8bytes
        # read(ii) strings(7) ! 4BYTE/8BYTE, char8bytes -- 4/8 byte integers used subsequently for ids/element type
        # read(ii) strings(8) ! SEQUE/PARAL, char8bytes
        # read(ii) strings(9) ! NOFIL/FILTE, char8bytes

    # for i in range(len(strings)):
    #    print(strings[i])

    for i in range(0, 5):
        integers = integers + [read_one_fp90_record(file_object, 1, np.int32)]

    # for i in range(len(integers)):
    #    print(integers[i])

    # read(ii) integers(1) ! ??? int32
    # read(ii) integers(2) ! nelem_total, int32
    # read(ii) integers(3) ! # Subdomains, number of parts?
    if(strings[1][0:5] != 'V0001'):
        integers = integers + [read_one_fp90_record(file_object, 1, np.int32)]
        # read(ii) integers(4) ! Time step, int32

    reals = read_one_fp90_record(file_object, 1, np.float64)
    # read(ii) reals(1)    ! Time, float64

    if(strings[1][0:5] == 'V0001'):
        integers[3] = int(reals)  # ! floor()?

    # return {'strings':strings, 'integers':integers, 'reals':reals}

    number_of_dimensions = integers[0][0]
    number_of_tuples_total = integers[1][0]
    time_instant_int = integers[3][0]
    time_instant_real = reals[0]

    if(strings[5] == 'REAL'):
        field_dtype = 'float'
    if(strings[5] == 'INTEG'):
        field_dtype = 'int'

    if(strings[4] == 'NPOIN'):
        association = 'node'
    else:
        association = 'element'

    if strings[6] == '4BYTE':
        field_dtype = field_dtype + '32'
    elif strings[6] == '8BYTE':
        field_dtype = field_dtype + '64'
    else:
        assert False, f'Alya id type {header[6]} is not supported'

    if(strings[3] == 'SCALA'):
        variabletype = 'scalar'
    elif(strings[3] == 'VECTO'):
        variabletype = 'vector'
    else:
        assert False, "unsupported type of variable"

    assert (strings[8] == 'NOFIL'), "Filtered types not supported"

    header = {'DataType': field_dtype, 'Lines': number_of_tuples_total, 'Columns': number_of_dimensions, 'TimeStepNo': time_instant_int,
              'Time': time_instant_real, 'Association': association, 'VariableType': variabletype}

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

        tuples = np.reshape(np.fromfile(f, dtype=np.dtype(
            header['DataType'])), (header['Lines'], header['Columns']))

        return {'tuples': tuples, 'header': header, 'tuples_per_block': []}


def read_alyabin_array(filename, number_of_blocks):
    alya_int_type = identify_alya_id_type(path)


    with open(filename, 'rb') as f:
        header = read_header_alyabin(f)
        number_of_dimensions = header['Columns']
        number_of_tuples_total = header['Lines']
        time_instant_int = header['TimeStepNo']
        time_instant_real = header['Time']
        #print(f'Reading array: {number_of_dimensions} dim, {number_of_tuples_total} tuples\n')

        datatype = np.dtype(header['DataType'])

        tuples = np.zeros(
            (number_of_tuples_total, number_of_dimensions), dtype=datatype)

        c :int= 0
        tuples_per_block = np.zeros(number_of_blocks, dtype=np.int32)


        for i in range(number_of_blocks):
            number_of_tuples_in_block = read_one_fp90_record(
                f, 1, alya_int_type)[0]  # stored by alya
            tuples_per_block[i] = number_of_tuples_in_block

            print(f'Block {i}/{number_of_blocks}: {(number_of_tuples_in_block)} tuples. Ndimensions: {number_of_dimensions}\n')
            tuples_temp = read_one_fp90_record(
                f, number_of_dimensions*number_of_tuples_in_block, datatype)

            tuples[c:c+number_of_tuples_in_block, :] = np.reshape( tuples_temp, (number_of_tuples_in_block, number_of_dimensions) )
            c = c+number_of_tuples_in_block

    return {'tuples': tuples, 'header': header, 'tuples_per_block': tuples_per_block}


def read_alya_variable(path, filename, number_of_blocks):
    field_filename = os.path.join( path, filename )
    return read_alya_array(field_filename, number_of_blocks)




def find_isoch( path ):
    #find the ishoch file with the largest timestep number
    global MPIO
    files = os.listdir( path )
    isoch_names = []
    isoch_numbers = []
    for file in files:
        if file_isoch in file.upper():
            data = file.split('.')
            if file_suffix_mpio in file.lower():
                filename = data[-4]
                MPIO = True
            elif file_suffix_alya in file.lower():
                filename = data[-3]
                MPIO = False

            number = int(filename.split('-')[-1])
            isoch_numbers.append( number )
            isoch_names.append( file )


    ind = np.argmax( np.array(isoch_numbers) )
    print('Isoch: ', isoch_names[ind] )

    return isoch_names[ind]


def find_file(path, pattern):
    files = os.listdir( path )
    for file in files:
        if pattern in file.lower():
            return file

    return None


def identify_alya_id_type(path):
    #read he header where element ids are stored and see if it's int8 or int4
    lninv_file = find_file(path, 'lninv')
    filename = os.path.join( path, lninv_file )


    with open(filename, 'rb') as f:
        header = read_header(f)
        if '32' in header['DataType']:
            alya_id_type = np.int32
        elif '64' in header['DataType']:
            alya_id_type = np.int64
        else:
            assert False, f'Alya id type {header[6]} is not supported'
    return alya_id_type        




iso_file = find_isoch(path)
part_file = find_file(path, ".post.alyapar")




if part_file is None:
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('!! Partition file *.post.alyapar not found ')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    raise 

# read the partitioning info
partition_filename = os.path.join( path, part_file )
with open(partition_filename) as f:
    partitions = np.fromstring(f.read(), dtype=np.int64, sep=' ')

# describes the mesh partitions
#partition_id,  NumberOfElementsInPartition,  NumberOfPointsInPartition, NumberOfBoundariesInPartition
partitions = np.reshape(partitions[1:], (partitions[0], 4))
number_of_blocks = partitions.shape[0]

partitions = pandas.DataFrame(
    partitions, columns=['id', 'Elements', 'Points', 'Boundaries'])


try:
    data = read_alya_variable(path, iso_file, number_of_blocks)
except:
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('!!! An error occured reading variable ', iso_file)
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    sys.stdout.flush()
    raise
    
iso_data = data['tuples'].ravel()
max_t = iso_data.max()
iso_data[iso_data<=1e-14] = np.finfo('double').max

activated_perc = ( iso_data<np.finfo('double').max ).sum()*100/iso_data.shape[0]

perc = np.percentile(iso_data, [95, 98, 99])
string = f'TAT: {max_t:.5f}, Activated: {activated_perc:.2f}%'
json_data = {}

json_data['TAT'] = max_t
json_data['Activated'] = activated_perc

for name, p in zip( ['Perc95','Perc98','Perc99'], perc):
    if p<np.finfo('double').max:
        string += f', {name}: {p:.5f}'
        json_data[name] = p

print( string )

with open(json_file, 'w') as f:
    json.dump( json_data, f )


#import matplotlib.pyplot as plt
#plt.hist(iso_data, bins=100)
#plt.show()

