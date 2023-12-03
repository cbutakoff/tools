import numpy as np  # needs install
import os
import pandas  # needs install
import sys
import json
from pathlib import Path
import progressbar as pb


path = Path(sys.argv[1])
alya_problem = sys.argv[2] #e.g. nostenosis_nolayers
varname = sys.argv[3] #e.g YPLUS

file_suffix_mpio = '.post.mpio.bin'


#-----------------------------------


def read_alyampio_header(filename):
    with open(filename, 'rb') as f:
        return read_header_mpio(f)


def read_alyampio_array(filename):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)

        tuples = np.reshape(np.fromfile(f, dtype=np.dtype(header['DataTypePython'])),
                            (header['Lines'], header['Columns']))

        return {'tuples': tuples, 'header': header};


def read_header_mpio(f):
    magic = np.fromfile(f, count=1, dtype=np.int64)[0]
    if magic != 27093:
        print(f'File {f} does not appear to be alya mpio file')

    format = str(f.read(8))
    if not ('MPIAL' in format):
        assert False, f'File {f} does not appear to be alya mpio file'

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
        assert False, f'ID column in {f} is not supported'

    junk = str(f.read(8))
    if not ('0000000' in junk):
        assert False, f'Lost alignment reding {f}'

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

    options = []
    for i in range(10):
        options.append(f.read(8))

    if 'INT' in datatype:
        dt = 'int'
    elif 'REAL' in datatype:
        dt = 'float'
    else:
        assert False, f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt + '64'
    elif '4' in datatypelen:
        dt = dt + '32'
    else:
        assert False, f'Unsupported data type length {datatypelen}'

    header = {
        'Version': version,
        'Object': obj,
        'Dimension': dimension,
        'Columns': columns,
        'Lines': lines,
        'Association': association,
        'DataType': datatype,
        'DataTypeLength': datatypelen,
        'TimeStepNo': timestep_no,
        'Time': time,
        'NSubdomains': nsubdomains,
        'Div': mesh_div,
        'DataTypePython': dt,
        'Options': options}

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

        for i in range(10):
            option = header["Options"][i] if header["Options"][i] else b'NONE000\0'
            if isinstance( option, str ):
                option = option.encode()
            assert len(option)==8, f"{i}-th option \'{header['Options'][i]}\' has to be 8 symbols long"
            f.write(option)

        vector.astype(header['DataTypePython']).tofile(f)

#----------------------------------------------------------

with open(path/f"{alya_problem}.post.alyafil", 'r') as ff:
    files = [x.strip() for x in ff if varname in x]

if len(files)==0:
    print("Nothing to postprocess")
    exit(1)


data = read_alyampio_array(files[0])
maxvalue = data["tuples"].copy()


progress = pb.ProgressBar()
for filename in progress(files):
    data = read_alyampio_array(filename)
    maxvalue = np.fmax( maxvalue, data["tuples"] )

varname_new = list(varname)
varname_new[-1] = "M"
varname_new = "".join(varname_new)
output_filename = path/f"{alya_problem}-{varname_new}-00000000.post.mpio.bin"
data['header']['Time']=0.0
data['header']['TimeStepNo']=0
write_vector_mpio(output_filename, maxvalue, data['header'])
