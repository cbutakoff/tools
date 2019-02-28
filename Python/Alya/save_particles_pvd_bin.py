import vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np
from os.path import isfile

main_filename = 'fluidda.pvd'
pts_filename = '../fluidda.pts.res'
ncpus = 10

print('Parsing ',pts_filename)


dt = np.dtype([('T', 'f8'), ('ILAGR', 'i8'), ('ITYPE', 'i8'), ('EXIST', 'i8'),\
 ('COORX', 'f8'), ('COORY', 'f8'), ('COORZ', 'f8'),\
 ('VELOX', 'f8'), ('VELOY', 'f8'), ('VELOZ', 'f8'),\
 ('DTK', 'f8'), ('CD', 'f8') ])

with open(pts_filename, mode='rb') as file: # b is important -> binary
    hdr = file.read(255)
    data = np.fromfile(file, dtype=dt)

df = pd.DataFrame.from_records(data)


#print(df)
#
#df = pd.read_csv(pts_filename, delim_whitespace=True,\
# names=['T', 'ILAGR', 'ITYPE', 'EXIST', 'COORX', 'COORY', 'COORZ', 'VELOX', 'VELOY', 'VELOZ', 'DTK', 'CD'],\
# header=0, comment='#')
#df = pd.DataFrame({'T':T,'ILAGR':ILAGR, 'ITYPE':ITYPE, 'EXIST':EXIST, 'COORX':COORX, 'COORY':COORY, 'COORZ':COORZ})

df['T']=  df['T'].apply(lambda x: np.round(x,10))

#print(df['T'])

#unique_times = df['T'].unique()
unique_times = list(np.round(np.arange(0,100+(5e-2),5e-2),10))

#print(unique_times)


print('Saving the particles')


def save_one_file(i):
    pts = vtk.vtkPoints()
    df1 = df[df['T']==unique_times[i]]
    x = df1['COORX']
    y = df1['COORY']
    z = df1['COORZ']
    tp = df1['ITYPE']
    il = df1['ILAGR']
    pts.SetNumberOfPoints(df1.shape[0])

    if df1.shape[0]!=0:
        for j in range(df1.shape[0]):         
            pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))

        types = vtk.vtkShortArray()
        types.SetNumberOfComponents(1)
        types.SetNumberOfTuples(pts.GetNumberOfPoints())
        types.SetName('Type')
      
        ilagr = vtk.vtkIntArray()
        ilagr.SetNumberOfComponents(1)
        ilagr.SetNumberOfTuples(pts.GetNumberOfPoints())
        ilagr.SetName('ILAGR')


        for j in range(df1.shape[0]):         
            types.SetTuple1(j, tp.iloc[j])
            ilagr.SetTuple1(j, il.iloc[j])

        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.GetPointData().AddArray(types)
        pd.GetPointData().AddArray(ilagr)

        wr= vtk.vtkXMLPolyDataWriter()
        wr.SetInputData(pd)
        wr.SetDataModeToBinary()
        filename = "pts_{:010d}.vtp".format(i)
        wr.SetFileName(filename)
        wr.Write()

bar = progressbar.ProgressBar(max_value=len(unique_times))
with mp.Pool(processes = ncpus) as p:
    for i, _ in enumerate(p.imap_unordered(save_one_file, range(len(unique_times))), 1):
        bar.update(i)




print('Generating main pvd file')
with open(main_filename,'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write(' <Collection>\n')

    for i in progressbar.progressbar(range(len(unique_times))):
        filename = "pts_{:010d}.vtp".format(i)
        if isfile(filename):
            f.write(f'<DataSet timestep="{unique_times[i]}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')


