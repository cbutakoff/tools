import numpy as np  # needs install
import os
import pandas  # needs install
import sys
import vtk
#from progressbar import progressbar
import multiprocessing as mp
import argparse
import json



#replacement for progressbar for mutltithreading
def progressbar(a):
    return a

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
def read_alya_array(filename, number_of_blocks, first_line, last_line):
    if MPIO:
        return read_alyampio_array(filename, number_of_blocks, first_line, last_line)
    else:
        return read_alyabin_array(filename, number_of_blocks, first_line, last_line)


def read_alyampio_array(filename, number_of_blocks, first_line, last_line):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)

        f.seek(first_line*header['Columns'],1)
        tuples = np.reshape(np.fromfile(f, dtype=np.dtype(
            header['DataType'])), (last_line-first_line+1, header['Columns']))

        return {'tuples': tuples, 'header': header }


def read_alyabin_array(filename, number_of_blocks, first_line, last_line):
    alya_int_type = identify_alya_id_type(path)


    with open(filename, 'rb') as f:
        header = read_header_alyabin(f)
        number_of_dimensions : int = header['Columns']
        number_of_tuples_total : int = header['Lines']
        time_instant_int : int = header['TimeStepNo']
        time_instant_real = header['Time']
        #print('Reading ', filename)
        #print(f'Reading array: {number_of_dimensions} dim, {number_of_tuples_total} tuples\n')

        datatype = np.dtype(header['DataType'])

        tuples = np.zeros(
            (last_line-first_line+1, number_of_dimensions), dtype=datatype)

        c = 0
        d = 0

        #print('First line: ', first_line)
        #print('Last line: ', last_line)

        for i in range(number_of_blocks):
            number_of_tuples_in_block = read_one_fp90_record(
                f, 1, alya_int_type)[0]  # stored by alya

            tuples_temp = read_one_fp90_record(f, number_of_dimensions*number_of_tuples_in_block, datatype)

            #print(f'Block {i}/{number_of_blocks}: {(number_of_tuples_in_block)} tuples\n')
            #print('c+number_of_tuples_in_block=',c+number_of_tuples_in_block)
            #print('c=',c)
            if ( c+number_of_tuples_in_block > first_line) and (c<=last_line):
                #print('number_of_tuples_in_block:', number_of_tuples_in_block)
                #print('number_of_dimensions:', number_of_dimensions)
                tuples_temp = np.reshape( tuples_temp, (number_of_tuples_in_block, number_of_dimensions) )

                #print('max(first_line-c,0)=',max(first_line-c,0))
                #print('min(last_line-c,tuples_temp1.shape[0])=',min(last_line-c+1,tuples_temp.shape[0]))
                tuples_temp1 = tuples_temp[ max(first_line-c,0):min(last_line-c+1,tuples_temp.shape[0]) ,: ]
                #print('tuples_temp1: ',tuples_temp1)
                #print('tuples_temp1.shape: ',tuples_temp1.shape)

                tuples[ d : d+tuples_temp.shape[0] , :] = tuples_temp1
                            
                d +=  tuples_temp1.shape[0]
                #print('d= ',d)
                #print('tuples.shape[0]=',tuples.shape[0])

            if d>=tuples.shape[0]:
                break

            c += number_of_tuples_in_block

    return {'tuples': tuples, 'header': header }



def read_alya_variable(path, filename, number_of_blocks, block_number, association):
    field_filename = os.path.join( path, filename )

    if block_number>0:
        first_line = partitions[association].values[0:block_number].sum()
    else:
        first_line = 0

    last_line = (first_line + partitions[association].values[block_number] - 1)
    #print('first_line=',first_line)
    #print('last_line=',last_line)


    return read_alya_array(field_filename, number_of_blocks, first_line, last_line)


def find_isoch( path ):
    #find the ishoch file with the largest timestep number
    global MPIO
    files = os.listdir( path )
    isoch_names = []
    isoch_numbers = []
    for file in files:
        if file_isoch in file.upper():
            if problem_name in file:
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

    if(len(isoch_numbers)!=0):
        ind = np.argmax( np.array(isoch_numbers) )
        #print('Isoch: ', isoch_names[ind] )

        return isoch_names[ind]
    else:
        return None



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



def ProcessPartition(part_id):    
    coords = read_alya_variable(path, 'wedge_scars-COORD.post.alyabin', number_of_blocks, part_id, 'Points')['tuples']
    connectivity = read_alya_variable(path, 'wedge_scars-LNODB.post.alyabin', number_of_blocks, part_id, 'Boundaries')['tuples']
    eltype = read_alya_variable(path, 'wedge_scars-LTYPB.post.alyabin', number_of_blocks, part_id, 'Boundaries')['tuples'].ravel()
    bnd_codes = read_alya_variable(path, 'wedge_scars-CODBO.post.alyabin', number_of_blocks, part_id, 'Boundaries')['tuples'].ravel()
    iso_data = read_alya_variable(path, iso_file, number_of_blocks, part_id, 'Points')['tuples'].ravel()
        

    element_alya2vtk = {37:{'Name':b'hexa8','Vertices':8,'VTK':vtk.VTK_HEXAHEDRON}, 30:{'Name':b'tetra4','Vertices':4,'VTK':vtk.VTK_TETRA}, \
                            32:{'Name':b'pyramid5','Vertices':5,'VTK':vtk.VTK_PYRAMID}, 34:{'Name':b'penta6','Vertices':6,'VTK':vtk.VTK_WEDGE},\
                            10:{'Name':b'tria3','Vertices':3,'VTK':vtk.VTK_TRIANGLE}, 12:{'Name':b'quad4','Vertices':4,'VTK':vtk.VTK_QUAD}}




    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(coords.shape[0])

    iso = vtk.vtkFloatArray()
    iso.SetNumberOfComponents(1)
    iso.SetNumberOfTuples(pts.GetNumberOfPoints())
    iso.SetName('ISOCH')

    #print('Creating points')
    for i in progressbar(range(coords.shape[0])):
        pts.SetPoint(i, coords[i,:])
        iso.SetTuple1(i, iso_data[i])

    codes = vtk.vtkShortArray()
    codes.SetNumberOfComponents(1)
    codes.SetNumberOfTuples(connectivity.shape[0])
    codes.SetName('BC')


    cells = vtk.vtkCellArray()
    #print('Creating booundary elements')
    for i in progressbar( range(connectivity.shape[0]) ):
        ptids = connectivity[ i, 0:element_alya2vtk[eltype[i]]['Vertices'] ]
        cells.InsertNextCell(ptids.shape[0])
        for pt in ptids:
            cells.InsertCellPoint(pt-1)

        codes.SetTuple1(i, bnd_codes[i])

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    pd.GetCellData().AddArray(codes)
    pd.GetPointData().AddArray(iso)

    th = vtk.vtkThreshold()
    th.SetInputData(pd)
    th.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,'BC')
    th.ThresholdBetween(epicardium_boundary_code-0.1, epicardium_boundary_code+0.1)
    th.Update()

    vel = np.array([])
    if th.GetOutput().GetNumberOfPoints()>0:
        grad = vtk.vtkGradientFilter()
        grad.SetInputConnection(th.GetOutputPort())
        grad.ComputeGradientOn ()
        grad.ComputeDivergenceOff ()
        grad.ComputeVorticityOff ()
        grad.ComputeQCriterionOff ()
        grad.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,'ISOCH')
        grad.Update()

        calc = vtk.vtkArrayCalculator()
        calc.SetInputConnection(grad.GetOutputPort())
        calc.SetAttributeTypeToPointData()
        calc.AddVectorArrayName('Gradients')
        calc.AddVectorVariable("V", 'Gradients')
        calc.SetFunction("1/mag(V)")
        calc.ReplaceInvalidValuesOn()
        calc.SetReplacementValue(0)
        calc.SetResultArrayName("VelocityMag")
        calc.Update()

        vel_vtk = calc.GetOutput().GetPointData().GetArray("VelocityMag")

        vel = np.zeros(calc.GetOutput().GetNumberOfPoints())
        for i in range(vel.shape[0]):
            vel[i] =  vel_vtk.GetTuple1(i)

        #wr = vtk.vtkDataSetWriter()
        #wr.SetFileName( f'{output_name}_{part_id:03d}.vtk' )
        #wr.SetInputConnection(calc.GetOutputPort())
        #wr.Write()
    
    return vel


def convert(o):
    if isinstance(o, np.int64): return int(o)  
    raise TypeError


#----------------------------------------------------------------
#
#                  Main program
#
#--------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("project_name", help='Name of the alya project')
    parser.add_argument("code", type=int, help='Boundary code for epicardium')
    parser.add_argument("--ncpus", help='Number of CPUs to use', type=int, default=1)
    parser.add_argument("--input_folder","-i", help='Folder with alyabins/mpios', default = '.')
    parser.add_argument("--csv", help='Filename of the csv to save velocities')
    parser.add_argument("--min", help='Min velocty', type=float)
    parser.add_argument("--max", help='Max velocty', type=float)
    parser.add_argument("--bins", help='Number of bins for mode calculation', type=int, default=100)
    parser.add_argument("--json", help='JSON file to store the historgram and the numbers', type=str)

    args = parser.parse_args()

    problem_name = args.project_name
    path = args.input_folder

    epicardium_boundary_code = args.code
    ncpus = args.ncpus

    print(f'Using {ncpus} cpus. Historgram will use {args.bins} bins.')

    file_suffix_mpio = '.post.mpio.bin'
    file_suffix_alya = '.post.alyabin'
    file_isoch = 'ISOCH'  # filename must contain this

    MPIO = True

    iso_file = find_isoch(path)

    if iso_file is None:
        print(f'No isochrones found in {path}')
        raise FileNotFoundError

    part_file = find_file(path, ".post.alyapar")




    if part_file is None:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('!! Partition file *.post.alyapar not found ')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        raise FileNotFoundError

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


    p = mp.Pool(ncpus)
    vel = np.concatenate(p.map(ProcessPartition, range(partitions.shape[0]))).ravel()
    #for part_id in range(partitions.shape[0]):
    #    print("===================================",part_id)
    #    ProcessPartition(part_id)


    if args.min is not None:
        min_veloc = args.min
    else:
        min_veloc = vel.min()

    if args.max is not None:
        max_veloc = args.max
    else:
        max_veloc = vel.max()


    if args.csv is not None:
        np.savetxt(args.csv, vel)
    if args.min is not None:
        vel = vel[ vel >= args.min ]
    if args.max is not None:
        vel = vel[ vel <= args.max ]


    #mode
    hh = np.histogram(vel, bins=args.bins)   
    bin_centers = (hh[1][1:]-hh[1][0:-1])/2+hh[1][1:]
    max_ind = np.argmax(hh[0])  
    mode = bin_centers[max_ind]  

    #median
    median = np.median(vel)
    mean = vel.mean()
    sd = vel.std()


    print(f'Conduction velocity mean={mean}, median={median}, mode={mode}, in range=[{min_veloc}, {max_veloc}]' )
    
    if args.json:
      json_data = {}
      json_data['Histogram'] = {'Bins': list(bin_centers), 'Counts':list(hh[0])}
      json_data['Mode'] = mode
      json_data['Mean'] = mean
      json_data['Median'] = median
      json_data['STD'] = sd

      
      with open(args.json, 'w') as f:
         json.dump(json_data, f, default=convert)
