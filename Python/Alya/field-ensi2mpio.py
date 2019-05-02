#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import progressbar
import io
import os


# In[3]:


problem_name = '2d'
path = f'ensi/{problem_name}.ensi.case'
outpath = 'fields/'
varname = 'VELOC'
ndimensions = 2
starttime = 5
endtime = 6
field_id = 1



#figure out the timeline for the variable
timeset_id = -1
with open(path,'r') as f:
    for line in f:
        if varname in line:
            digits = [int(s) for s in line.split() if s.isdigit()]
            timeset_id = digits[0]
            break

if timeset_id<0:
    print('Timeset of the variable ',varname, ' not found in the case. Check the spelling of the variable name')
    exit()

timeset_id = timeset_id-1 #vtk makes them start from 0

print(f'Found {varname} in timeset {timeset_id}')

assert( os.path.isfile(path) )

os.makedirs(outpath,exist_ok=True)


rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(path)
rdr.ReadAllVariablesOn()
rdr.Update()



ntimesteps = rdr.GetTimeSets().GetItem(timeset_id).GetNumberOfTuples()
print(f'Number of timesteps: {ntimesteps}')


timesteps = []
for time_idx in range(ntimesteps):
    time_instant  = rdr.GetTimeSets().GetItem(timeset_id).GetTuple1(time_idx)
    timesteps.append(np.round(time_instant,8))

timesteps = np.array(timesteps)
t20 = timesteps[ (timesteps>=starttime) & (timesteps<=endtime) ]

if t20.shape[0]==0:
    print('No timesteps in the specified range')
    exit()



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



corrected_times = np.round(t20-t20[0], 15)

for i in progressbar.progressbar(range(t20.shape[0])):
    t = t20[i]
    rdr.SetTimeValue(t)
    rdr.Update()
    mesh = rdr.GetOutput().GetBlock(0)
    
    veloc = mesh.GetPointData().GetArray(varname)

    v = np.zeros((mesh.GetNumberOfPoints(),ndimensions), dtype = np.float64)
    for j in range(veloc.GetNumberOfTuples()):
        vv = veloc.GetTuple3(j)
        v[j,:] = vv[:ndimensions]

    filename = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format(problem_name, field_id, i+1)

    write_vector_mpio( os.path.join( outpath, filename), v, corrected_times[i] )

#save times
with open( os.path.join( outpath,'times.in' ) ,'w') as ftimes:
    for t in progressbar.progressbar(corrected_times):
        ftimes.write(f'{t}\n')



