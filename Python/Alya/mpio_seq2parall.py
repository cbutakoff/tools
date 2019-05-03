#!/usr/bin/env python
# coding: utf-8

# In[39]:


import numpy as np
import os
import multiprocessing as mp
import progressbar

project_name = '2d'
inputfolder = './'
outputfolder = './paral_fields'
ncpus = 10


# In[51]:




def read_alyampio_array(filename):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataType']) ), (header['Lines'], header['Columns']) )
        
        return {'tuples':tuples, 'header':header};

    xfiels[0]
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
    else:xfiels[0]
        assert False,f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt+'64'
    elif '4' in datatypelen:
        dt = dt+'32'
    else:
        assert False,f'Unsupported data type length {datatypelen}'


    header = {'DataType':dt, 'Lines':lines,'Columns':columns,               'TimeStepNo':timestep_no, 'Time':time, 'NSubdomains':nsubdomains,              'Div':mesh_div}

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



def write_vector_mpio(filename, vector, header):
    ndim = vector.shape[1]
    npts = vector.shape[0]

    
    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(b'XFIEL00\0')
        f.write(b'VECTO00\0')
        f.write(b'NPOIN00\0')
        f.write(b'REAL000\0')
        f.write(b'8BYTE00\0')
        f.write(b'PARAL00\0')
        f.write(b'NOFIL00\0')
        f.write(b'ASCEN00\0')
        f.write(b'NOID000\0')
        f.write(b'0000000\0')
        np.array([ndim], dtype=np.int64).tofile(f)
        np.array([npts], dtype=np.int64).tofile(f)
        np.array([ header['TimeStepNo'] ], dtype=np.int64).tofile(f)  # Time step number (signed int 64 bits):    ittim
        np.array([ header['NSubdomains'] ], dtype=np.int64).tofile(f)  # n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
        np.array([ header['Div'] ], dtype=np.int64).tofile(f)  #  Mesh division (signed int 64 bits):       divi
        np.array([0], dtype=np.int64).tofile(f)  #Tag 1 (signed int 64 bits):               tag1
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag2
        np.array([ header['Time'] ], dtype=np.float64).tofile(f)  #Time (real 64 bits):                      time
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
        vector.astype('float64').tofile(f)


def list_xfiel(path):
    allfiles = os.listdir(path=path)
    return [s for s in allfiles if ('XFIEL' in s.upper()) and not ('POST' in s.upper()) ]


def save_one_file_by_name(file_name):
  infilename = os.path.join( inputfolder, file_name )
  outfilename = os.path.join( outputfolder, file_name )
  save_one_file( infilename, outfilename )

def save_one_file(infilename, outfilename):
    data = read_alyampio_array(infilename)

    ncomponents = data['tuples'].shape[1] #for 1d, 2d and 3d problems
    field2write = data['tuples'][inverse_pt_correspondence,:]
    
    write_vector_mpio(outfilename, field2write, data['header'])


# In[53]:


os.makedirs(outputfolder, exist_ok=True)

lninv_filename = os.path.join(inputfolder,f'{project_name}-LNINV.post.mpio.bin')
print('Reading paritioning info from ',lninv_filename)
LNINV = read_alyampio_array( lninv_filename )
inverse_pt_correspondence = (LNINV['tuples']-1).ravel(); #convert ids to python

xfiels = list_xfiel(inputfolder)

print("main process")
bar = progressbar.ProgressBar(max_value=len(xfiels))
with mp.Pool(processes = ncpus) as p:
  for i, _ in enumerate( p.imap_unordered(save_one_file_by_name, xfiels), 1 ):
    bar.update(i)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




