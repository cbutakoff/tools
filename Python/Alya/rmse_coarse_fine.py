import vtk
import numpy as np
import sys

# case2 and case1 are compared at the nodes from case1  
case_coarse_filename = sys.argv[1]  #1.ensi.case -- coarse mesh
case_fine_filename = sys.argv[2]  #2.ensi.case -- fine mesh
variable = sys.argv[3]        #VELOC

rounding_digits = 7

def FindTimeSet(reader, var):
    for i in range( reader.GetNumberOfVariables() ):
        if reader.GetDescription(i) == var:
            break

    return i

def CreateReader(case1_filename):
    ensi_rdr1 = vtk.vtkEnSightGoldBinaryReader()
    ensi_rdr1.SetCaseFileName(case1_filename)
    ensi_rdr1.ReadAllVariablesOn()
    ensi_rdr1.Update()
    
    return ensi_rdr1


def GetTimesteps(reader, timeset):
    ntimesteps=reader.GetTimeSets().GetItem(timeset).GetNumberOfTuples()

    times1 = []
    for time_idx in range(ntimesteps):
        times1.append( np.round(reader.GetTimeSets().GetItem(timeset).GetTuple1(time_idx), decimals=rounding_digits))

    return times1

def GetVariable(reader, variable, t, ptids = None):
    reader.SetTimeValue(t)
    reader.Update()

    mesh = reader.GetOutput().GetBlock(0)
    var  = mesh.GetPointData().GetArray(variable)
    ncomp = var.GetNumberOfComponents()

    v = None

    if ptids is None:
        v = np.zeros( (var.GetNumberOfTuples(), ncomp) )
        for i in range( v.shape[0] ):
            v[i,:] = var.GetTuple(i)
    else:
        v = np.zeros( (ptids.shape[0], ncomp) )
        for i, ptid in enumerate(ptids):
            v[i,:] = var.GetTuple(ptid)

    return v

def GetCoarsePointOnFine(reader_coarse, reader_fine):
    mesh_coarse = reader_coarse.GetOutput().GetBlock(0)
    mesh_fine   = reader_fine.GetOutput().GetBlock(0)

    loc = vtk.vtkStaticPointLocator()
    loc.SetDataSet( mesh_fine )


    ids = np.zeros( mesh_coarse.GetNumberOfPoints(), dtype=np.int64 )
    for i in range( ids.shape[0] ):
        ptid = loc.FindClosestPoint( mesh_coarse.GetPoint(i) )
        ids[i] = ptid
    
    return ids


####################################################################3
# READING
print("Reading coarse mesh", case_coarse_filename)
print("Reading fine mesh ", case_fine_filename)
rdr_coarse = CreateReader(case_coarse_filename)
rdr_fine   = CreateReader(case_fine_filename)


timeset_coarse = FindTimeSet(rdr_coarse, variable)
timeset_fine   = FindTimeSet(rdr_fine, variable)

print(f"Timelines {timeset_coarse}, {timeset_fine}")


times_coarse = GetTimesteps(rdr_coarse, timeset_coarse)
times_fine   = GetTimesteps(rdr_fine,   timeset_fine)

time_intersection=np.sort(list(set(times_coarse).intersection(times_fine)))

coarse_ids_on_fine = GetCoarsePointOnFine(rdr_coarse, rdr_fine)

print("Time intersection (n=",len(time_intersection),"):", time_intersection)
print("--------------------------")    

####################################################################3
# compute
for t in time_intersection:
    v1 = GetVariable( rdr_coarse, variable, t )
    v2 = GetVariable( rdr_fine  , variable, t, coarse_ids_on_fine )

    for j in range(v1.shape[1]):
        print (f"time: {t} | RMSE({j}) = {np.sqrt(np.mean((v1[:,j]-v2[:,j])**2))}" )
        print (f"time: {t} | MIN({j}) = {np.min(v1[:,j]-v2[:,j])}" )
        print (f"time: {t} | MAX({j}) = {np.max(v1[:,j]-v2[:,j])}" )
