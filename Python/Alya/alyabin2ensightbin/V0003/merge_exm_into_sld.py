"""
pass all exmedi arrays to solidz mesh
"""

import vtk
import sys
import os
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from progressbar import progressbar, ProgressBar
import shutil
import glob

casename_exm = sys.argv[1]
casename_sld = sys.argv[2]
out_path     = sys.argv[3] #will copy solidz case and add new variables

apply_displacement = False
displ_varname = 'DISPL' #variable name for displacement
time_check_tolerance = 1e-7 #tolerance when time values are compared

#==============================================================
#
#   FUNCTIONS
#

def find_displ_timeset(casefileaname):
  variables_section = False
  timeline = None
  with open(casefileaname) as ff:
    for line in ff:
      if "VARIABLE" in line.upper():
        variables_section = True
      if variables_section:
        print (line)
        if " DISPL " in line.upper():
          data = np.array(line.upper().split())
          print(data)
          idx = np.where( data=='DISPL' )[0]
          timeline = int( data[idx-1] )-1
          break

  return timeline

def calculate_interpolation_weights( mesh_exm, mesh_sld ):
    """ Calculate intepolation matrices, !!!!!!ONLY TRIANGLES!!!!!! AND ONLY POINT DATA !!!!
    
        For each point of mesh_sld use vtkCellLocator to find closest point on mesh_exm, interpolation weights and point ids for interpolation

        Return interp_weights - npts x 3 matrix of weights, interp_ptids - npts x 3 matrix of point ids 
        so that interpolated value of array X is:
        X[i] = interp_weights[i,0] * a[interp_ptids[i,0]] + interp_weights[i,1] * a[interp_ptids[i,1]] + interp_weights[i,2] * a[interp_ptids[i,2]]
        or 
        X = np.sum(interp_weights * a[interp_ptids], axis=1)

    """

    npts = mesh_sld.GetNumberOfPoints()
    interp_weights = np.zeros( (npts, 3), dtype=np.float32 ) 
    interp_ptids   = np.zeros( (npts, 3), dtype=np.int64   ) 

    loc = vtk.vtkCellLocator()
    loc.SetDataSet( mesh_exm ) 
    loc.BuildLocator()

    for i in progressbar( range(npts) ):
        closest_point = [0,0,0]
        cellId = vtk.mutable(0)
        subId  = vtk.mutable(0)
        dist2  = vtk.mutable(0)
        loc.FindClosestPoint( mesh_sld.GetPoint(i), closest_point, cellId, subId, dist2 )
            
        #get interpolation weights
        cell = mesh_exm.GetCell( cellId )
        closestPoint1 = [0,0,0]
        pcoords       = [0,0,0]
        weights       = [0]*cell.GetNumberOfPoints() 
        reval = cell.EvaluatePosition( closest_point, closestPoint1, subId, pcoords, dist2, weights )
        assert reval>=0, 'mesh_exm.GetCell() failed for unknown reason'

        interp_weights[i,:] = weights

        #get point ids
        ptids = cell.GetPointIds()
        for j in range( ptids.GetNumberOfIds() ):
            interp_ptids[i, j] = ptids.GetId(j)

    return interp_weights, interp_ptids

def interpolate_array( array, interp_weights, interp_ptids ):
    """ Interpolate array "array"
        
        array          : np.ndarray(K) or np.ndarray((K,N))  vector or columns of values to interpolate from. Must contain at least max(interp_ptids.ravel()) values
        interp_weights : np.ndarray((npts,3)) matrix of weights
        interp_ptids   : np.ndarray((npts,3)) matrix of point ids 

        returns np.ndarray( npts ) of interpolate values
    """

    result = None

    if array.ndim==1:
        result = np.sum(interp_weights * array[interp_ptids], axis=1)
    elif array.ndim==2:
        result = np.zeros( interp_weights.shape[0], array.shape[1] )

        for i in range( a.shape[1] ):
            aaa = array[:,i]
            result[:,i] = np.sum(interp_weights * aaa[interp_ptids], axis=1)

    assert not result is None, f'interpolate_array: array dimension has to be 1 or 2, passed {array.ndim}, shape: {array.shape}'
    
    return result

#
#  END OF FUNTIONS
#
#===============================================================

#==============================================================
#
#   ENSI Writer
#
def write_variable_pernode(data, outputfolder, project_name, varname, iteration):
    #variable ensight
    iterationid_number_of_digits = 6
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';
    filename = os.path.join(outputfolder, fmt % (project_name, varname, iteration))

    with open( filename,'wb') as f:
        f.write(b'description line 1'.ljust(80))
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=np.int32))   #int
        f.write(b'coordinates'.ljust(80))


        if data.ndim==1:            
            f.write( data )  
        elif data.ndim==2:
            #data has coordinates in the order [[x,y,z],[x,y,z],...]
            #expected order of coordinates
            #vx_n1 vx_n2 ... vx_nn nn floats
            #vy_n1 vy_n2 ... vy_nn nn floats
            #vz_n1 vz_n2 ... vz_nn nn floats
            #Rearrange the  matrix
            f.write( data.ravel(order='F').astype(np.float32) )  
        else:
            assert False, f"write_variable_pernode: array must have dimension 1 or 2, got {data.ndim}, shape {data.shape}"      

def get_case_name( filename ):
    assert filename.endswith('.ensi.case'), f'get_case_name was passed not ensi.case, passed: {filename}'
    return filename.replace('.ensi.case','')

def get_DISPL_filenumbers( casename ) -> int:
    """ returns filename_numbers, time_values, timeset_number for DISPL """
    variable_section    = False
    time_section        = False
    in_displ_timeset    = False
    in_filename_numbers = False
    in_time_values      = False
    displ_timeset = None   
    nsteps        = None

    filename_numbers = []
    time_values = []

    with open( casename, 'r' ) as ff:
        for line in ff:
            if len(line.split()) == 0:
                continue        
            
            if line.split()[0][0] == '#':
                continue

            if line.split()[0]=='VARIABLE':
                variable_section = True

            if line.split()[0]=='TIME':
                time_section = True

            if variable_section:
                if displ_varname in line:
                    data = line.split()
                    for i, v in enumerate(data):
                        if v==displ_varname:
                            displ_timeset = int( data[i-1] )
                            variable_section = False
        
            if time_section:
                if "time set" in line.lower():
                    if int( line.split()[-1] ) == displ_timeset:
                        in_displ_timeset = True
                        time_section = False

            if in_displ_timeset:
                if "number of steps" in line.lower():
                    nsteps = int( line.split()[-1] )
                elif "filename numbers" in line.lower():
                    in_filename_numbers = True
                    in_time_values = False
                elif "time values" in line.lower():
                    in_filename_numbers = False
                    in_time_values = True
                else:
                    if in_filename_numbers:
                        filename_numbers.extend(line.split())
                        in_filename_numbers = len(filename_numbers)<nsteps                             
                    elif in_time_values:
                        time_values.extend(line.split())
                        in_time_values = len(time_values)<nsteps

    return np.array( [int(x) for x in filename_numbers], dtype=int), np.array( [float(x) for x in time_values] ), displ_timeset



def update_casefile_newvariables( casefilename, vars, displ_filename_numbers, displ_timeset_no  ):
    """ vars = { 'INTRA': {'file_numbers' : [1,2,3..], 'time_values' : [0,0.1,0.2,...], 'ndim':1, 'ntimeset':1 }, ... } 

        displ_filename_numbers, displ_timeset_no  are needed to decide if I can use the DISPL timeset or do I need to generate a new timeset
    """

    with open( casefilename, 'r' ) as fff:
        filedata = fff.read()        

    #get min timset number
    min_timeset = 100000000 
    for varname, vardata in vars.items():
        if vardata['ntimeset'] < min_timeset:
            min_timeset = vardata['ntimeset'] 

    #see which variables can reuse DISPL timeset
    #renumber other timesets
    timset_no = min_timeset
    for varname, vardata in vars.items():
        if np.array_equal( vardata['file_numbers'], displ_filename_numbers ):
            vars[varname]['ntimeset'] = displ_timeset_no
        else:
            vars[varname]['ntimeset'] = timset_no
            timset_no += 1



    vars_added_lines = []
    times_added_lines = []
    for varname, vardata in vars.items():
        line = ""
        if vardata['ndim']==1:
            line += "scalar"
        else:
            line += "vector"
       
        line += f" per node: {vardata['ntimeset']} {varname} {os.path.basename(casefilename).replace('.case','.')}{varname}-******" 
        
        vars_added_lines.append(line)
    
        if vardata['ntimeset'] != displ_timeset_no:
            line = []
            line.append( f"time set: {vardata['ntimeset']}" )
            line.append( f"number of steps: {len(vardata['file_numbers'])}" )
            line.append(  "filename numbers:" )
            line.append( " ".join( [ str(x)+"\n" if i%12==0 else str(x) for i,x in enumerate(vardata['file_numbers']) ] ) )
            line.append( "time values:" )
            line.append( " ".join( [ str(x)+"\n" if i%12==0 else str(x) for i,x in enumerate(vardata['time_values']) ] ) )

            times_added_lines.extend(line)

    #print(vars_added_lines)
    #print(times_added_lines)

    filedata = filedata.replace( 'VARIABLE', 'VARIABLE\n'+"\n".join( vars_added_lines ) )
    filedata = filedata.replace( 'TIME',     'TIME\n'+"\n".join( times_added_lines ) )

    with open( casefilename, 'w' ) as fff:
        fff.write(filedata)        



#
#   ENSI Writer
#
#==============================================================




print(f"Reading {casename_exm}")
assert os.path.isfile(casename_exm)
rdr_exm = vtk.vtkEnSightGoldBinaryReader()
rdr_exm.SetCaseFileName(casename_exm)
rdr_exm.ReadAllVariablesOn()
rdr_exm.Update()

print(f"Reading {casename_sld}")
assert os.path.isfile(casename_sld)
rdr_sld = vtk.vtkEnSightGoldBinaryReader()
rdr_sld.SetCaseFileName(casename_sld)
rdr_sld.ReadAllVariablesOn()
rdr_sld.Update()


#copy the sld case to the output path
print("Copying solidz case")
os.makedirs( out_path, exist_ok = True )
#shutil.copytree( os.path.dirname(casename_sld), out_path )
for filename in glob.glob(os.path.join(os.path.dirname(casename_sld), f'{get_case_name(os.path.basename(casename_sld))}*'), recursive=False):
    if not os.path.isdir(filename):
        shutil.copy(filename, out_path)

casefile_out = os.path.join( out_path, os.path.basename(casename_sld) )

print("Getting filenumbers") 
filename_numbers, time_values_tmp, displ_timeset_no = get_DISPL_filenumbers( casefile_out )

#print( filename_numbers )
#print( time_values_tmp ) 

#list exm variables
varnames = []
for i in range( rdr_exm.GetNumberOfVariables() ):
    varnames.append( rdr_exm.GetDescription(i) )

assert len(varnames)>0, 'Exmedi mesh is missing variables'

#find index of "DISPL" in sld mesh
displ_index = None
displ_index = find_displ_timeset(casename_sld)
assert not displ_index is None, 'Solidz does not have DISPL'

#get the timesteps of solidz
solidz_times = vtk_to_numpy( rdr_sld.GetTimeSets().GetItem(displ_index) )
#assert (np.abs(solidz_times - time_values_tmp)<time_check_tolerance).all(), 'Time values read by VTK for DISPL are different than read by this script. This should not have happened'

#check if by any chance solidz is one timestep further than exmedi and skip that timestep
#can happen with alya
if solidz_times.max() != time_values_tmp.max():
  print( "Max time in solidz ({solidz_times.max()}) and exmedi ({time_values_tmp.max()}) are different")
solidz_times = solidz_times[solidz_times<=time_values_tmp.max()]


#---------------------------------------------------------
#
#   precalculate interpolation weights
#

# extract the points of the exmedi mesh for interpolation later on
print("Loading geometry")
rdr_exm.SetTimeValue( rdr_exm.GetTimeSets().GetItem(0).GetTuple1(0) )
rdr_exm.Update()

rdr_sld.SetTimeValue( rdr_sld.GetTimeSets().GetItem(0).GetTuple1(0) )
rdr_sld.Update()

print("Precalcualting interpolation weights")
interp_weights, interp_ptids = calculate_interpolation_weights( rdr_exm.GetOutput().GetBlock(0), rdr_sld.GetOutput().GetBlock(0) )

#---------------------------------------------------------
#
#   Main loop:

# For each timestep:
written_vars = {}
for ivar, varname in enumerate(varnames):
    written_vars[varname]={}
    written_vars[varname]['file_numbers'] = []
    written_vars[varname]['time_values']  = []
    written_vars[varname]['ndim']         = None
    written_vars[varname]['ntimeset']     = ivar + rdr_sld.GetTimeSets().GetNumberOfItems() + 1


print("Saving variables")
bar = ProgressBar(max_value=len(solidz_times))
for ti, t in enumerate(solidz_times):
    rdr_exm.SetTimeValue(t)
    rdr_exm.Update()
    mesh_exm = rdr_exm.GetOutput().GetBlock(0) 

    rdr_sld.SetTimeValue(t)
    rdr_sld.Update()
    #mesh_sld = vtk.vtkUnstructuredGrid()
    #mesh_sld.DeepCopy( rdr_sld.GetOutput().GetBlock(0) )

    for varname in varnames:
        exm_array = mesh_exm.GetPointData().GetArray(varname)
        if exm_array is None:
            #print(f'Exmedi timestep {t} does not have array {varname}')
            pass
        else:
            #arr = numpy_to_vtk( interpolate_array(  vtk_to_numpy( mesh_exm.GetPointData().GetArray(varname) ), interp_weights, interp_ptids ) )
            #arr.SetName( varname )
            #mesh_sld.GetPointData().AddArray(arr)   
            arr = interpolate_array(  vtk_to_numpy( mesh_exm.GetPointData().GetArray(varname) ), interp_weights, interp_ptids )
            write_variable_pernode (  arr, out_path, get_case_name(os.path.basename(casename_sld)), varname, filename_numbers[ti] )
            written_vars[varname]['file_numbers'].append( filename_numbers[ti] )
            written_vars[varname]['time_values'].append( t )
            written_vars[varname]['ndim'] = arr.ndim

        #wr = vtk.vtkXMLUnstructuredGridWriter()
        #wr.SetInputData( mesh_sld )
        #wr.SetFileName(out_prefix)
        #wr.EncodeAppendedDataOff()
        #wr.Write()

    bar.update(ti)

print("")

#print(written_vars)
print(f"Updating casefile {casefile_out}")
update_casefile_newvariables( casefile_out, written_vars, filename_numbers, displ_timeset_no )

