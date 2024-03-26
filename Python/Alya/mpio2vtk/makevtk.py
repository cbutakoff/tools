""" Script to reconstruct a mesh from the mpios, for debugging material generation and other things on the mesh """

import vtk
import numpy as np
import sys
import os
from vtk.util.numpy_support import numpy_to_vtk

name = sys.argv[1]
coords = f"{name}-COORD.mpio.bin"
elems  = f"{name}-LNODS.mpio.bin"
lmate  = f"{name}-LMATE.mpio.bin"
leset  = f"{name}-LESET.mpio.bin"
cellt  = f"{name}-CELLT.mpio.bin"
out    = "mesh.vtu"



def read_alyampio_array(filename):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataTypePython']) ), (header['Lines'], header['Columns']) )
        
        return {'tuples':tuples, 'header':header};


def read_header_mpio(f):
    magic = np.fromfile(f,count=1, dtype=np.int64)[0]
    if magic != 27093:
        print(f'File {filename} does not appear to be alya mpio file')
        
    format = str(f.read(8))
    if not ('MPIAL' in format):
        assert False,f'File {filename} does not appear to be alya mpio file'

    version = str(f.read(8))
    obj = str(f.read(8))
    dimension = str(f.read(8))
    association = str(f.read(8))
    datatype = str(f.read(8))
    datatypelen = str(f.read(8))    
    seq_par = str(f.read(8))
    filt = str(f.read(8))    
    sorting = str(f.read(8))    
    idd = str(f.read(8))    

    if not ('NOID' in idd):
        assert False, f'ID column in {filename} is not supported'

    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,   f'Lost alignment reding {filename}'
    

    columns = np.fromfile(f,count=1,dtype=np.int64)[0]
    lines = np.fromfile(f,count=1,dtype=np.int64)[0]
    timestep_no = np.fromfile(f,count=1,dtype=np.int64)[0]
    nsubdomains = np.fromfile(f,count=1,dtype=np.int64)[0]
    mesh_div = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag1 = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag2 = np.fromfile(f,count=1,dtype=np.int64)[0]
    time = np.fromfile(f,count=1,dtype=np.float64)[0]
    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,f'Lost alignment reding {filename}'

    junk = str(f.read(8))    #1
    junk = str(f.read(8))    #2
    junk = str(f.read(8))    #3
    junk = str(f.read(8))    #4
    junk = str(f.read(8))    #5
    junk = str(f.read(8))    #6
    junk = str(f.read(8))    #7
    junk = str(f.read(8))    #8
    junk = str(f.read(8))    #9
    junk = str(f.read(8))    #10
    
    if 'INT' in datatype:
        dt = 'int'
    elif 'REAL' in datatype:
        dt = 'float'
    else:
        assert False,f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt+'64'
    elif '4' in datatypelen:
        dt = dt+'32'
    else:
        assert False,f'Unsupported data type length {datatypelen}'


    header = {
        'Version':version, 
        'Object':obj,
        'Dimension':dimension,
        'Columns':columns,
        'Lines':lines,
        'Association':association,
        'DataType':datatype,
        'DataTypeLength':datatypelen,
        'TimeStepNo':timestep_no, 
        'Time':time, 
        'NSubdomains':nsubdomains,              
        'Div':mesh_div,
        'DataTypePython':dt}


    assert ('NOFIL' in filt), "Filtered fields are not supported"


    return header



c = read_alyampio_array(coords)['tuples']
v = read_alyampio_array(elems)['tuples']

pts = vtk.vtkPoints()
pts.SetData( numpy_to_vtk(c) )

eles = vtk.vtkCellArray()
for pt in v:
    eles.InsertNextCell(pt.shape[0])
    [eles.InsertCellPoint(ptc-1) for ptc in pt]

ug = vtk.vtkUnstructuredGrid()
ug.SetPoints( pts )
ug.SetCells( vtk.VTK_TETRA, eles )

for filename, arrayname in zip([lmate, leset, cellt],['lmate', 'leset', 'cellt']):  
    if os.path.isfile(filename):
        m = read_alyampio_array(filename)['tuples'].astype(int)
        mvtk = numpy_to_vtk(m) 
        mvtk.SetName(arrayname)
        ug.GetCellData().AddArray( mvtk )

for i in range(10):
    filename = f"{name}-XFIEL.{i:08d}.00000001.mpio.bin"
    if os.path.isfile(filename):
        print(f"reading {filename}")
        m = read_alyampio_array(filename)['tuples']
        mvtk = numpy_to_vtk(m)
        mvtk.SetName(f"XFIEL{i:02d}")
        ug.GetPointData().AddArray( mvtk )


wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetFileName(out)
wr.SetInputData(ug)
wr.Write()



