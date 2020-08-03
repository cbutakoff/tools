#!/usr/bin/env python
# coding: utf-8


import numpy as np
import os
import progressbar
import sys
import pandas as pd

input_filename = sys.argv[1]
output_filename = sys.argv[2]



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



data = read_alyampio_array(input_filename)

with open(output_filename+'.hdr','w') as fff:
    for key, value in data['header'].items():
        fff.write(f'{key} : {value}\n')


df = pd.DataFrame( data['tuples'] )
df.index = list(range(1,df.shape[0]+1))
with open(output_filename+'.csv', 'w') as fff:
    df.to_csv(fff, header=False, sep=' ')

    
