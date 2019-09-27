import numpy as np
import os


def GenerateMPIOFieldFilename(problemname, field_number):
    return f'{problemname}-XFIEL.{field_number:08d}.00000001.mpio.bin'


def MPIO_write_field(path, problem_name, field_number, vector, time):
    #vector is Npoints x ndim
    #association is points
    ndim = vector.shape[1]
    npts = vector.shape[0]

    filename = os.path.join( path, GenerateMPIOFieldFilename(problem_name, field_number) )

    realtype =  ( vector.dtype == np.dtype('float64')) | ( vector.dtype == np.dtype('float32') )

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(b'XFIEL00\0')

        if ndim == 1:
            f.write(b'SCALA00\0')
        else:
            f.write(b'VECTO00\0')


        f.write(b'NPOIN00\0')

        if realtype:
            f.write(b'REAL000\0')
        else:
            f.write(b'INTEG00\0')
        
        f.write(b'8BYTE00\0')
        f.write(b'SEQUE00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array([0], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array([1], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array([0], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag2
        np.array([time], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
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

        if realtype:
            vector.astype('float64').tofile(f)
        else:
            vector.astype('int64').tofile(f)



def MPIO_write_matrix(path, problem_name, vector, variable, association):
    #vector is Npoints x ndim 
    #variable = COORD,  LNODS, LTYPE, LTYPE, LNODS, LMATE, LNODB, LELBO
    #association = NPOIN, NELEM, NBOUN

    if len(vector.shape)==1:
        ndim = 1
    else:
        ndim = vector.shape[1]

    npts = vector.shape[0]

    var = variable.upper().ljust(7,'0').encode()+b'\0'

    filename = os.path.join( path, f'{problem_name}-{variable}.mpio.bin' ) 

    realtype =  ( vector.dtype == np.dtype('float64')) | ( vector.dtype == np.dtype('float32') )

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(var)

        if ndim == 1:
            f.write(b'SCALA00\0')
        else:
            f.write(b'VECTO00\0')


        f.write( association.upper().ljust(7,'0').encode()+b'\0' )

        if realtype:
            f.write(b'REAL000\0')
        else:
            f.write(b'INTEG00\0')
        
        f.write(b'8BYTE00\0')
        f.write(b'SEQUE00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array([0], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array([1], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array([0], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag2
        np.array([0], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
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

        if realtype:
            vector.astype('float64').tofile(f)
        else:
            vector.astype('int64').tofile(f)

