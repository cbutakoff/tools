#!/home/costa/anaconda3/envs/py36/bin/python3.6
# -------- parameters ------------

input_filename_mesh = 'mesh_x4_hexfibers.vtk'
output_filename_mesh = 'PL6_hex.mesh.in'
output_filename_fibers = 'PL6_hex.fibers.in'
fiber_array_name = 'Fibers'
scale_factor = 0.00104 #cm
#----------------------------------


#---- main code ------------------------

import vtk
from TextProgressBar import TextProgressBar 


print('Reading mesh')
rd = vtk.vtkDataSetReader()
rd.SetFileName(input_filename_mesh)
rd.Update()
mesh = rd.GetOutput()

file = open(output_filename_mesh,"w")

print('Writing elements');
file.write( "ELEMENTS\n" )

pb = TextProgressBar( mesh.GetNumberOfCells(), prefix = 'Progress:', suffix = 'Complete', length = 50)

if mesh.GetCell(0).GetNumberOfPoints() == 8:
	id_order = [0,1,3,2,4,5,7,6]  #to reorder the points
	for i in range(mesh.GetNumberOfCells()):
	    ids = vtk.vtkIdList()
	    mesh.GetCellPoints(i, ids)
	    ids_list = [str(i+1)] + [ str(ids.GetId(id_order[k])+1) for k in range(ids.GetNumberOfIds())]
	    file.write( " ".join(ids_list) + "\n" )
	    pb.UpdateProgress(i+1)

elif  mesh.GetCell(0).GetNumberOfPoints() == 4:
	id_order = [0,1,2,3]  #to reorder the points
	for i in range(mesh.GetNumberOfCells()):
	    ids = vtk.vtkIdList()
	    mesh.GetCellPoints(i, ids)
	    ids_list = [str(i+1)] + [ str(ids.GetId(id_order[k])+1) for k in range(ids.GetNumberOfIds())]
	    file.write( " ".join(ids_list) + "\n" )	
	    pb.UpdateProgress(i+1)
	

file.write( "END_ELEMENTS\n" )


print('Writing points');
file.write( "COORDINATES\n" )

pb = TextProgressBar( mesh.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 50)
for i in range(mesh.GetNumberOfPoints()):
    pt = mesh.GetPoint(i)
    file.write("%d %f %f %f\n"%(i+1, pt[0]*scale_factor, pt[1]*scale_factor, pt[2]*scale_factor))
    pb.UpdateProgress(i+1)

file.write( "END_COORDINATES\n" )
file.close()



print('Writing fibers');
file = open(output_filename_fibers,"w")

fibers = mesh.GetPointData().GetArray(fiber_array_name)

pb = TextProgressBar( fibers.GetNumberOfTuples(), prefix = 'Progress:', suffix = 'Complete', length = 50)

for i in range(fibers.GetNumberOfTuples()):
    data = list(fibers.GetTuple(i))
    for k in range(3): #discard large numbers
        if abs(data[k])>1:
            data[k] = 0
    file.write( "%d %f %f %f\n"%(i+1, data[0], data[1], data[2]) )

    pb.UpdateProgress(i+1)
    
file.close()
