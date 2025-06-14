import numpy as np
import os
import progressbar
import sys


def read_alyampio_array(filename):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataTypePython']) ), (header['Lines'], header['Columns']) )
        
        return {'tuples':tuples, 'header':header};


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


    header = {
        'Version':version, 
        'Object':obj,
        'Dimension':dimension,
        'Columns':columns,
        'Lines':lines,
        'Association':association,
        'DataType':datatype,
        'DataTypeLength':datatypelen,
        'TimeStepNo':timestep_no, 
        'Time':time, 
        'NSubdomains':nsubdomains,              
        'Div':mesh_div,
        'DataTypePython':dt}


    assert ('NOFIL' in filt), "Filtered fields are not supported"


    return header



def write_vector_mpio(filename, vector, header):
    ndim = vector.shape[1]
    npts = vector.shape[0]

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(eval(header['Version']))
        f.write(eval(header['Object']))

        if ndim==1:
            f.write(b'SCALA00\0')
        else:    
            f.write(b'VECTO00\0')

        f.write(eval(header['Association']))
        f.write(eval(header['DataType']))
        f.write(eval(header['DataTypeLength']))
        f.write(b'SEQUE00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array(header['TimeStepNo'], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array(header['NSubdomains'], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array(header['Div'], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag2
        np.array(header['Time'], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
        f.write(b'0000000\0')
        f.write(b'NONE000\0') #1
        f.write(b'NONE000\0') #2
        f.write(b'NONE000\0') #3
        f.write(b'NONE000\0') #4
        f.write(b'NONE000\0') #5
        f.write(b'NONE000\0') #6
        f.write(b'NONE000\0') #7
        f.write(b'NONE000\0') #8
        f.write(b'NONE000\0') #9
        f.write(b'NONE000\0') #10

        vector.astype(header['DataTypePython']).tofile(f)




ref_t = float(sys.argv[1])
files = sorted([f for f in os.listdir(".") if "INTRA" in f])

data1 = read_alyampio_array(files[0])

zero_cross = -1*np.ones( data1['tuples'].shape )
for f1, f2 in zip( files[:-1], files[1:] ):
    data1 = read_alyampio_array(f1)
    t1    = data1['header']['Time']

    if t1>=ref_t:
        print(f1)

        data2 = read_alyampio_array(f2)
        t2    = data2['header']['Time']

        t = t1 - (t2-t1) * data1['tuples'] / (data2['tuples']-data1['tuples']) #0 crossing of INTRA
        valid_t = (t>=t1) & (t<=t2) & ( data2['tuples']>data1['tuples']  ) & (zero_cross<0)

        if any(valid_t): 
            zero_cross[ valid_t ] = t[ valid_t ]

        if all(zero_cross>=0):
            break; #all isochrones set

data1['header']['Time'] = 0.0
data1['header']['TimeStepNo'] = 0
write_vector_mpio( "VT_4mat_ohara_nomodel_ab-ISOCH-00000000.post.mpio.bin", zero_cross, data1['header'] )

