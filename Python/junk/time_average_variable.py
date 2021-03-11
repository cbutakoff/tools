import vtk
import numpy as np
import progressbar
import io
import os
import argparse
import sys


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_ensicase", help='Path to the ensi.case file')
parser.add_argument("output_folder", help='Folder for the output fields case')
parser.add_argument("output_prefix", help='Prefix for the output files', type=str )
parser.add_argument("-t", help='Time interval to average (s)', nargs=3, type=float, metavar="timestart timeend dt")
args = parser.parse_args()

assert (np.array(args.t) > 0).all(), "Time must be  > 0"
assert args.t[1] > args.t[0], "End time > start time"
assert args.output_prefix != '', "Output prefix cannot be empty string"

path =         args.input_ensicase #f'/media/gpfs_bsc/gpfs/projects/bsc21/WORK-BEATRIZ/spinal/catheters/C6-C7/straight/ensi/geo200216_c6c7_straight.ensi.case'
outpath =      args.output_folder  #'./'
varname_intra = 'INTRA'
dt =           args.t

assert( os.path.isfile(path) )

os.makedirs(outpath,exist_ok=True)


rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(path)
rdr.ReadAllVariablesOn()
rdr.Update()
mesh = rdr.GetOutput().GetBlock(0)


#figure out the timeline for the variable
timeset_id = -1

for i in range(rdr.GetTimeSets().GetNumberOfItems()):
    if varname_intra == rdr.GetPointArrayName(i): timeset_id = i

assert timeset_id>=0, f'Timeset of the variable {varname_intra} not found in the case. Check the spelling of the variable name'
print(f'Found {varname_intra} in timeset {timeset_id}')

ntimesteps = rdr.GetTimeSets().GetItem(timeset_id).GetNumberOfTuples()
print(f'Total number of timesteps in the mesh for the variable {varname_intra}: {ntimesteps}')


#get all the timesteps
timesteps = []
for time_idx in range(ntimesteps):
    time_instant  = rdr.GetTimeSets().GetItem(timeset_id).GetTuple1(time_idx)
    timesteps.append(np.round(time_instant,8))

timesteps = np.array(timesteps)

#lump timesteps according to the dt
time_intervals = []
for t in np.arange( dt[0], dt[1], dt[2] ):
    time_intervals.append( timesteps[ (timesteps>=t) & (timesteps<(t+dt[2])) ].tolist() )

#print (time_intervals)
assert len(time_intervals)>0, 'No time intervals identified'

file_times = [] #these are the times for the generated files
with progressbar.ProgressBar(max_value=len(time_intervals)) as bar:
    for time_id, time_list in enumerate(time_intervals):

        intra = np.zeros( (mesh.GetNumberOfPoints(),len(time_list)) ) #npoints x nttimesteps
        mesh = None
        for ti, t in enumerate(time_list):
            rdr.SetTimeValue(t)
            rdr.Update()
            mesh = rdr.GetOutput().GetBlock(0)

            vtk_var = mesh.GetPointData().GetArray(varname_intra)
            intra[:,ti] = [ vtk_var.GetTuple1(i) for i in range(vtk_var.GetNumberOfTuples()) ]
           

        ar1 = vtk.vtkFloatArray()        
        ar1.SetNumberOfComponents(1)
        ar1.SetNumberOfTuples(vtk_var.GetNumberOfTuples())
        ar1.SetName("INTRA_MAX")
        [ ar1.SetTuple1(i, x) for i, x in enumerate( intra.max(axis=1) ) ] #get max INTRA for each node

        mesh.GetPointData().AddArray(ar1)

        wr = vtk.vtkXMLUnstructuredGridWriter()
        wr.SetFileName(f'{args.output_prefix}_{time_id:010d}.vtu')
        wr.SetInputData(mesh)
        wr.Write()

        file_times.extend( [time_list[-1]] ) #add the last time to the times

        bar.update(time_id)


#wrte the PVTU to open in paraview


print('Generating main pvd file')
with open(f'{args.output_prefix}.pvd','w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write(' <Collection>\n')

    for time_id in range(len(time_intervals)):
        filename = f"{args.output_prefix}_{time_id:010d}.vtu"
        if os.path.isfile(filename):
            f.write(f'<DataSet timestep="{np.round(file_times[time_id],10)}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')



