import vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np
from os.path import isfile
from sys import exit

main_filename = 'fluidda.pvd'
pts_filename0 = '../fluidda.pts.res'
pts_filename = '../fluidda.pts.{:08d}.res'
ncpus = 10

timestep = 5e-2

filenumbers = range(10,40010,10)


def read_file( filename ):
    dt = np.dtype([('T', 'f8'), ('ILAGR', 'i8'), ('ITYPE', 'i8'), ('EXIST', 'i8'),\
     ('COORX', 'f8'), ('COORY', 'f8'), ('COORZ', 'f8'),\
     ('VELOX', 'f8'), ('VELOY', 'f8'), ('VELOZ', 'f8'),\
     ('DTK', 'f8'), ('CD', 'f8') ])

    with open(filename, mode='rb') as file:
        hdr = file.read(255)
        data = np.fromfile(file, dtype=dt)
    
    df = pd.DataFrame.from_records(data, exclude=['VELOX','VELOY','VELOZ','DTK', 'CD'])

    return df
    

def save_one_file_by_number(file_number):
  infilename = pts_filename.format(file_number) 
  outfilename = "pts_{:010d}.vtp".format(file_number)
  save_one_file( infilename, outfilename )


def save_one_file(infilename, outfilename):

    df1 = read_file( infilename )

    if df1.shape[0]==0:
      return

    pts = vtk.vtkPoints()
    x = df1['COORX']
    y = df1['COORY']
    z = df1['COORZ']
    tp = df1['ITYPE']
    il = df1['ILAGR']
    ex = df1['EXIST']
    tt = df1['T']
    pts.SetNumberOfPoints(df1.shape[0])

    for j in range(df1.shape[0]):         
        pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))

    types = vtk.vtkShortArray()
    types.SetNumberOfComponents(1)
    types.SetNumberOfTuples(pts.GetNumberOfPoints())
    types.SetName('ITYPE')
  
    ilagr = vtk.vtkIntArray()
    ilagr.SetNumberOfComponents(1)
    ilagr.SetNumberOfTuples(pts.GetNumberOfPoints())
    ilagr.SetName('ILAGR')

    exist = vtk.vtkShortArray()
    exist.SetNumberOfComponents(1)
    exist.SetNumberOfTuples(pts.GetNumberOfPoints())
    exist.SetName('EXIST')

    T = vtk.vtkFloatArray()
    T.SetNumberOfComponents(1)
    T.SetNumberOfTuples(pts.GetNumberOfPoints())
    T.SetName('T')

    for j in range(df1.shape[0]):         
        types.SetTuple1(j, tp.iloc[j])
        ilagr.SetTuple1(j, il.iloc[j])
        exist.SetTuple1(j, ex.iloc[j])
        T.SetTuple1(j, tt.iloc[j])

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.GetPointData().AddArray(types)
    pd.GetPointData().AddArray(ilagr)
    pd.GetPointData().AddArray(exist)
    pd.GetPointData().AddArray(T)

    wr= vtk.vtkXMLPolyDataWriter()
    wr.SetInputData(pd)
    wr.SetDataModeToBinary()
    wr.SetFileName(outfilename)
    wr.Write()


#save 0-th timestep
save_one_file( pts_filename0, "pts_{:010d}.vtp".format(0) )


bar = progressbar.ProgressBar(max_value=len(filenumbers))
with mp.Pool(processes = ncpus) as p:
  for i, _ in enumerate( p.imap_unordered(save_one_file_by_number, filenumbers), 1 ):
    bar.update(i)



print('Generating main pvd file')
with open(main_filename,'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write(' <Collection>\n')

    filename = "pts_{:010d}.vtp".format(0)
    if isfile(filename):
        f.write(f'<DataSet timestep="0" group="" part="0"\n')
        f.write(f'     file="{filename}"/>\n')


    for filenumber in progressbar.progressbar(filenumbers):
        filename = "pts_{:010d}.vtp".format(filenumber)
        if isfile(filename):
            f.write(f'<DataSet timestep="{np.round(filenumber*timestep,10)}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')




