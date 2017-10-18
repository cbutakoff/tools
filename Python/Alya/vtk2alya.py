#!/home/costa/anaconda3/envs/py36/bin/python3.6
# -------- parameters ------------

input_filename_mesh = 'mesh_x4_hexfibers.vtk'
output_filename_mesh = '19w_hex.mesh.dat'
output_filename_fibers = '19w_hex.point.fibers.dat'
fiber_array_name = 'Fibers'
#----------------------------------


#---- main code ------------------------

import vtk
from TextProgressBar import TextProgressBar 



rd = vtk.vtkDataSetReader()
rd.SetFileName(input_filename_mesh)
rd.Update()
mesh = rd.GetOutput()

file = open(output_filename_mesh,"w")

print('Writing elements');
file.write( "ELEMENTS\n" )

pb = TextProgressBar( mesh.GetNumberOfCells(), prefix = 'Progress:', suffix = 'Complete', length = 50)

for i in range(mesh.GetNumberOfCells()):
    ids = vtk.vtkIdList()
    mesh.GetCellPoints(i, ids)
    ids_list = [str(i+1)] + [ str(ids.GetId(k)+1) for k in range(ids.GetNumberOfIds())]
    file.write( " ".join(ids_list) + "\n" )
    pb.UpdateProgress(i+1)
	

file.write( "END_ELEMENTS\n" )


print('Writing points');
file.write( "COORDINATES\n" )

pb = TextProgressBar( mesh.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 50)
for i in range(mesh.GetNumberOfPoints()):
    pt = mesh.GetPoint(i)
    file.write("%d %f %f %f\n"%(i+1, pt[0], pt[1], pt[2]))
    pb.UpdateProgress(i+1)

file.write( "END_COORDINATES\n" )
file.close()



print('Writing fibers');
file = open(output_filename_fibers,"w")

fibers = mesh.GetPointData().GetArray(fiber_array_name)

pb = TextProgressBar( fibers.GetNumberOfTuples(), prefix = 'Progress:', suffix = 'Complete', length = 50)

for i in range(fibers.GetNumberOfTuples()):
    data = fibers.GetTuple(i)
    file.write( "%d %f %f %f\n"%(i+1, data[0], data[1], data[2]) )

    pb.UpdateProgress(i+1)
    
file.close()
