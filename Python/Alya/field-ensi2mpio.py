#

import vtk
import numpy as np
import progressbar
import io
import os
import argparse

# In[3]:

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("problem_name", help='Name of the alya problem')
parser.add_argument("input_ensicase", help='Path to the ensi.case file')
parser.add_argument("output_folder", help='Folder for the output fields case')
parser.add_argument("varname", help='Variable associated to the NODES to save')
parser.add_argument("-f", "--field",type=int, default=1, help = "Field id for the generated field (>0)")
parser.add_argument("-t", "--timesteps", help='Comvert only timesteps in the specified range, use time_step=-1 to save all timesteps', nargs=3, type=float, metavar=('time_start', 'time_end','time_step'))
args = parser.parse_args()

assert args.field > 0, "Field ID must be > 0"

problem_name = args.problem_name   #'geo200216_c6c7_straight'
path =         args.input_ensicase #f'/media/gpfs_bsc/gpfs/projects/bsc21/WORK-BEATRIZ/spinal/catheters/C6-C7/straight/ensi/geo200216_c6c7_straight.ensi.case'
outpath =      args.output_folder  #'./'
varname =      args.varname        #'VELOC'
field_id =     args.field



assert( os.path.isfile(path) )

os.makedirs(outpath,exist_ok=True)


rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(path)
rdr.ReadAllVariablesOn()
rdr.Update()



#figure out the timeline for the variable
timeset_id = -1

for i in range(rdr.GetTimeSets().GetNumberOfItems()):
    if varname == rdr.GetPointArrayName(i):
        timeset_id = i

assert timeset_id>=0, f'Timeset of the variable {varname} not found in the case. Check the spelling of the variable name'
print(f'Found {varname} in timeset {timeset_id}')

ntimesteps = rdr.GetTimeSets().GetItem(timeset_id).GetNumberOfTuples()
print(f'Total number of timesteps in the mesh for the variable {varname}: {ntimesteps}')


#get all the timesteps
timesteps = []
for time_idx in range(ntimesteps):
    time_instant  = rdr.GetTimeSets().GetItem(timeset_id).GetTuple1(time_idx)
    timesteps.append(np.round(time_instant,8))


timesteps = np.array(timesteps)
t20 = timesteps
if not args.timesteps is None:
    if args.timesteps[2]<0: #all timesteps in the range   
        t20 = timesteps[ (timesteps>=args.timesteps[0]) & (timesteps<=args.timesteps[1]) ]
    else
        t20 = np.arange(args.timesteps[0], args.timesteps[1]+args.timesteps[2],args.timesteps[2]) 

    assert t20.shape[0]>0 f'No timesteps in the specified range {args.timesteps}'
else:
    assert t20.shape[0]>0 f'There are no timesteps for the specified variable'

print(f"Extracting {t20.shape[0]} timesteps: {t20}")

def write_vector_mpio(filename, vector, time):
    ndim = vector.shape[1]
    npts = vector.shape[0]

    with open(filename, 'wb') as f:
        np.array([27093], dtype=np.int64).tofile(f)
        f.write(b'MPIAL00\0')
        f.write(b'V000400\0')
        f.write(b'XFIEL00\0')

        if ndim > 1:
            f.write(b'VECTO00\0')
        else:
            f.write(b'SCALA00\0')

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
        np.array([0], dtype=np.int64).tofile(f)  #Tag 2 (signed int 64 bits):               tag
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
    assert not veloc is None, f'{varname} array not found in the mesh'

    v = np.zeros((mesh.GetNumberOfPoints(), veloc.GetNumberOfComponents()), dtype = np.float64)
    for j in range(veloc.GetNumberOfTuples()):
        v[j,:] = veloc.GetTuple(j)

    filename = '{:s}-XFIEL.{:08d}.{:08d}.mpio.bin'.format(problem_name, field_id, i+1)

    write_vector_mpio( os.path.join( outpath, filename), v, corrected_times[i] )

#save times
with open( os.path.join( outpath,'times.in' ) ,'w') as ftimes:
    for t in progressbar.progressbar(corrected_times):
        ftimes.write(f'{t}\n')



