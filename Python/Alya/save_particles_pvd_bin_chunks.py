import vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np
from os.path import isfile
from sys import exit

main_filename = 'fluidda.pvd'
pts_filename = '../fluidda.pts.res'
ncpus = 10
lines_per_chunk = 100000

starttime = 0
endtime = 200
timestep = 5e-2

halfstep = timestep/2.0

times = np.array( np.round(np.arange(starttime,endtime+timestep,timestep),10) )

df_global = pd.DataFrame([])

def read_chunk( file, nlines ):
    dt = np.dtype([('T', 'f8'), ('ILAGR', 'i8'), ('ITYPE', 'i8'), ('EXIST', 'i8'),\
     ('COORX', 'f8'), ('COORY', 'f8'), ('COORZ', 'f8'),\
     ('VELOX', 'f8'), ('VELOY', 'f8'), ('VELOZ', 'f8'),\
     ('DTK', 'f8'), ('CD', 'f8') ])
    data = np.fromfile(file, count=nlines, dtype=dt)
    print('Read ',len(data),' lines')
    df = pd.DataFrame.from_records(data, exclude=['VELOX','VELOY','VELOZ','DTK', 'CD'])
#    df['T']=  df['T'].apply(lambda x: np.round(x,10))
#    print(df['T'])

    return df


def read_vtk_file(filename):
    rd = vtk.vtkXMLPolyDataReader()
    rd.SetFileName(filename)
    rd.Update()
    data = rd.GetOutput()

    N = data.GetNumberOfPoints()
    x = [];
    y = [];
    z = [];
    itype = [];
    ilagr = [];
    exist = [];
    T = [];

    for i in range(N):
        pt = data.GetPoint(i)
        x.append(pt[0])
        y.append(pt[1])
        z.append(pt[2])
        itype.append( data.GetPointData().GetArray('ITYPE').GetTuple1(i) )
        ilagr.append( data.GetPointData().GetArray('ILAGR').GetTuple1(i) )
        exist.append( data.GetPointData().GetArray('EXIST').GetTuple1(i) )
        T.append( data.GetPointData().GetArray('T').GetTuple1(i) )

    df = pd.DataFrame({'COORX':x, 'COORY':y, 'COORZ':z, 'ITYPE':itype, 'ILAGR':ilagr, 'EXIST':exist, 'T':T })
    return df
    

def save_one_file(i):
    #df1 = df[df['T']==unique_times[i]]
    df1 = df_global[ (df['T']>=times[i]-halfstep) & (df['T']<times[i]+halfstep) ]

    if df1.shape[0]==0:
      return

    filename = "pts_{:010d}.vtp".format(i)

    if isfile(filename):
        dft = read_vtk_file(filename)
        df1 = dft.append(df1, ignore_index=True, sort=False) #add df to the end of the read data

    pts = vtk.vtkPoints()
    x = df1['COORX']
    y = df1['COORY']
    z = df1['COORZ']
    tp = df1['ITYPE']
    il = df1['ILAGR']
    ex = df1['EXIST']
    tt = df1['T']
    pts.SetNumberOfPoints(df1.shape[0])

    if df1.shape[0]!=0:
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
        wr.SetFileName(filename)
        wr.Write()




print('Parsing ',pts_filename)

filesize=0


with open(pts_filename, mode='rb') as file: # b is important -> binary

    file.seek(0,2)
    filesize = file.tell()
    file.seek(0,0)

    hdr = file.read(255)


    while True:
        df = read_chunk(file,lines_per_chunk)
        if df.shape[0]==0:
            break

        df_global = df

        print( 'Read ', file.tell()*100/filesize,'%' )

        local_time_global_ids = []   

        df_times_min = df['T'].min()
        df_times_max = df['T'].max()

        #find the ids if time values that correcpond to this chunk of data
        for i, t in enumerate(times):
            if( ( t >= df_times_min-halfstep ) & ( t< df_times_max+halfstep ) ):
              local_time_global_ids.append( i )
            



        bar = progressbar.ProgressBar(max_value=len(local_time_global_ids))
        with mp.Pool(processes = ncpus) as p:
            for i, _ in enumerate( p.imap_unordered(save_one_file, local_time_global_ids), 1 ):
                bar.update(i)



print('Generating main pvd file')
with open(main_filename,'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write(' <Collection>\n')

    for i in progressbar.progressbar(range(len(times))):
        filename = "pts_{:010d}.vtp".format(i)
        if isfile(filename):
            f.write(f'<DataSet timestep="{times[i]}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')




