import vtk
import numpy as np
import sys
from progressbar import progressbar
from scipy.stats import entropy
import pandas as pd

casefile = sys.argv[1]
variable_name = sys.argv[2]
timeline = int(sys.argv[3])-1
nbins = int(sys.argv[4])
outputfile = sys.argv[5]

rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(casefile)
rdr.ReadAllVariablesOn()

print("Verify if the file is valid: ", casefile)
if (not rdr.CanReadFile(casefile)):
    print("Invalid file: ")
    raise

rdr.DebugOff()
rdr.Update()

print("Number of timesets: ", rdr.GetTimeSets().GetNumberOfItems() )
print("Time values for timeset ",timeline)
    
n_timesteps = rdr.GetTimeSets().GetItem(timeline).GetNumberOfTuples()
times = []
for i in range(n_timesteps):
    if rdr.GetTimeSets().GetItem(timeline).GetTuple1(i)<0:
      continue

    times.append( rdr.GetTimeSets().GetItem(timeline).GetTuple1(i));

print("Times: ",times)

t=[]
e=[]
for i in progressbar(range(n_timesteps)):
    rdr.SetTimeValue(times[i])
    rdr.Update()

    volmesh = rdr.GetOutput().GetBlock(0)
    scalars = volmesh.GetPointData().GetArray( variable_name ) 
    np_scalars = np.zeros( scalars.GetNumberOfTuples() )
    #extract scalars
    for k in range(np_scalars.shape[0]):
        np_scalars[k] = scalars.GetTuple1(k)
    
    hist,edges = np.histogram(np_scalars, bins=nbins, density=True)
    entr = entropy(hist)

    t.append(np.round(times[i],13))
    e.append(entr)

df = pd.DataFrame( {"Time":t, "Entropy":e})

df.to_csv(outputfile)
