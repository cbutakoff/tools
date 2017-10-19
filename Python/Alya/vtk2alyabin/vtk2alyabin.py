#!/home/costa/anaconda3/envs/py36/bin/python3.6
# -------- parameters ------------

input_filename_mesh = 'hex.vtk'
output_filename_mesh = 'hex1.dom.bin'
output_filename_fibers = 'hex1.fibers'
fiber_array_name = 'Fibers'
scale_factor = 1 #cm
#----------------------------------


#---- main code ------------------------

import vtk
import numpy as np
from TextProgressBar import TextProgressBar 

alya_int = np.int64;


def write_f90_number(file, data):
    dtype = data.dtype
    count = data.size
    file.write(np.array([count*np.dtype(dtype).itemsize], dtype=np.int32)) #length
    file.write(data.ravel())
    file.write(np.array([count*np.dtype(dtype).itemsize], dtype=np.int32)) #length

def write_f90_string(file, data):
    data1 = data+b'12345678'
    count = len(data1)
    file.write(np.array([count], dtype=np.int32)) #length
    file.write(data1)
    file.write(np.array([count], dtype=np.int32)) #length


print('Reading mesh')
rd = vtk.vtkDataSetReader()
rd.SetFileName(input_filename_mesh)
rd.Update()
mesh = rd.GetOutput()

with open(output_filename_mesh,"wb") as file:
    print('Writing elements');

    one = np.array([1], dtype=alya_int)
    write_f90_number( file, np.array([123456], dtype=np.int32) )
    write_f90_string( file, b'DIMEN'.ljust(20) )
    write_f90_number( file, np.array([mesh.GetNumberOfCells(), mesh.GetNumberOfPoints(), 0, 0, 0, 0, 0], dtype=alya_int) )

    element_type_id = 30;
    number_of_pts_per_cell = mesh.GetCell(0).GetPointIds().GetNumberOfIds();
    if number_of_pts_per_cell ==4:
        element_type_id = 30 #tet04
    elif number_of_pts_per_cell ==8:
        element_type_id = 37 #hex08


    write_f90_string( file, b'LTYPE'.ljust(20) )
    element_types = np.ones(mesh.GetNumberOfCells(), dtype=alya_int)*element_type_id;
    print(f'Number of cells {element_types.shape}')
    write_f90_number( file, element_types )

    write_f90_string( file, b'LNNOD'.ljust(20) )  #elemnt ids
    write_f90_number( file, np.arange(1, mesh.GetNumberOfCells()+1, dtype=alya_int) )

    pb = TextProgressBar( mesh.GetNumberOfCells(), prefix = 'Progress:', suffix = 'Complete', length = 50)

    cells = np.zeros((mesh.GetNumberOfCells(), number_of_pts_per_cell), dtype=alya_int)

    if  number_of_pts_per_cell == 4:
        id_order = [0,1,2,3]  #to reorder the points
    elif  number_of_pts_per_cell == 8:
        id_order = [0,1,3,2,4,5,7,6]  #to reorder the points

    print('Writing cells');

    print(cells.shape)
    for i in range(mesh.GetNumberOfCells()):
        ids = vtk.vtkIdList()
        mesh.GetCellPoints(i, ids)
        ids_list = [ ids.GetId(id_order[k])+1 for k in range(ids.GetNumberOfIds())]

        for k in range(len(ids_list)):
            cells[i,k] = ids_list[k]
            
        pb.UpdateProgress(i+1)

    write_f90_string( file, b'LNODS'.ljust(20) )  #elemnt ids
    write_f90_number( file, cells )


    print('Writing points');

    pb = TextProgressBar( mesh.GetNumberOfPoints(), prefix = 'Progress:', suffix = 'Complete', length = 50)
    pts = np.zeros(mesh.GetNumberOfPoints()*3, dtype=np.float64)
    c= 0 
    for i in range(mesh.GetNumberOfPoints()):
        pt = mesh.GetPoint(i)
        for k in range(3):
            pts[c] = pt[k]*scale_factor
            c=c+1

        pb.UpdateProgress(i+1)

    write_f90_string( file, b'COORD'.ljust(20) )  
    write_f90_number( file, pts )

    write_f90_string( file, b'END_FILE'.ljust(20) )  
    


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

