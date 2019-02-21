import vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np

main_filename = '2d.pvd'
pts_filename = '../2d.pts.res'
ncpus = 10

print('Parsing ',pts_filename)

T =[]
ILAGR =[]
ITYPE = []
EXIST = []
COORX = []
COORY = []
COORZ = []
with open(pts_filename, mode='rb') as file: # b is important -> binary
    hdr = file.read(255)
    print(hdr)
    while True:
        T1 = np.fromfile(file, dtype=np.float64, count=1)
        if T1.shape[0]==0: break #end of file
        T.append(T1)

        ILAGR.append( np.fromfile(file, dtype=np.int64, count=1) )
        ITYPE.append( np.fromfile(file, dtype=np.int64, count=1) )
        EXIST.append( np.fromfile(file, dtype=np.int64, count=1) )
        COORX.append( np.fromfile(file, dtype=np.float64, count=1) )
        COORY.append( np.fromfile(file, dtype=np.float64, count=1) )
        COORZ.append( np.fromfile(file, dtype=np.float64, count=1) )
        VELOX = np.fromfile(file, dtype=np.float64, count=1)
        VELOY = np.fromfile(file, dtype=np.float64, count=1)
        VELOZ = np.fromfile(file, dtype=np.float64, count=1)
        DTK = np.fromfile(file, dtype=np.float64, count=1)
        CD = np.fromfile(file, dtype=np.float64, count=1)

#
#df = pd.read_csv(pts_filename, delim_whitespace=True,\
# names=['T', 'ILAGR', 'ITYPE', 'EXIST', 'COORX', 'COORY', 'COORZ', 'VELOX', 'VELOY', 'VELOZ', 'DTK', 'CD'],\
# header=0, comment='#')
df = pd.DataFrame({'T':T,'ILAGR':ILAGR, 'ITYPE':ITYPE, 'EXIST':EXIST, 'COORX':COORX, 'COORY':COORY, 'COORZ':COORZ})

df['T']=  df['T'].apply(lambda x: np.round(x,10))

print(df)

#unique_times = df['T'].unique()
unique_times = list(np.round(np.arange(0,1,0.01),10))

#print(unique_times)

empty_file = np.zeros(len(unique_times))

print('Saving the particles')


def save_one_file(i):
    pts = vtk.vtkPoints()
    df1 = df[df['T']==unique_times[i]]
    x = df1['COORX']
    y = df1['COORY']
    z = df1['COORZ']
    tp = df1['ITYPE']
    pts.SetNumberOfPoints(df1.shape[0])

    if df1.shape[0]==0:
        empty_file[i]=1
    else:
        for j in range(df1.shape[0]):         
            pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))

        types = vtk.vtkShortArray()
        types.SetNumberOfComponents(1)
        types.SetNumberOfTuples(pts.GetNumberOfPoints())
        types.SetName('Type')
        
        for j in range(df1.shape[0]):         
            types.SetTuple1(j, tp.iloc[j])

        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.GetPointData().AddArray(types)

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
        if empty_file[i]==0:
            filename = "pts_{:010d}.vtp".format(i)
            f.write(f'<DataSet timestep="{unique_times[i]}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')


