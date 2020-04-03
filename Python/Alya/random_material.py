#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import progressbar
import io
import os
import sys

# In[3]:


problem_name = sys.argv[1]
field_id = int(sys.argv[2])
npoints = int(sys.argv[3])
nmaterials = int(sys.argv[4])



def write_vector_mpio(filename, vector, time):
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
        vector.astype('float64').tofile(f)



filename = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format(problem_name, field_id, 1)
v = np.random.randint(1,nmaterials+1,[npoints,1])
write_vector_mpio( filename, v, 0 )


