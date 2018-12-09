import vtk
import pandas as pd
import  progressbar 

main_filename = 'fluidda.vtk'
pts_filename = 'fluidda.pts.res'
ncpus = 10

print('Parsing ',pts_filename)
df = pd.read_csv(pts_filename, delim_whitespace=True)


pts = vtk.vtkPoints()
T = df['T']
x = df['COORX']
y = df['COORY']
z = df['COORZ']
tp = df['ITYPE']
ilagr = df['ILAGR']
exist = df['EXIST']

N=df.shape[0]
pts.SetNumberOfPoints(N)

types = vtk.vtkShortArray()
types.SetNumberOfComponents(1)
types.SetNumberOfTuples(pts.GetNumberOfPoints())
types.SetName('Type')

times = vtk.vtkFloatArray()
times.SetNumberOfComponents(1)
times.SetNumberOfTuples(pts.GetNumberOfPoints())
times.SetName('Time')

vtkilagr = vtk.vtkIdTypeArray()
vtkilagr.SetNumberOfComponents(1)
vtkilagr.SetNumberOfTuples(pts.GetNumberOfPoints())
vtkilagr.SetName('ILAGR')

vtkexist = vtk.vtkShortArray()
vtkexist.SetNumberOfComponents(1)
vtkexist.SetNumberOfTuples(pts.GetNumberOfPoints())
vtkexist.SetName('EXIST')


#verts = vtk.vtkCellArray()

for j in progressbar.progressbar(range(N)):         
    pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))
    types.SetTuple1(j, tp.iloc[j])
    times.SetTuple1(j, T.iloc[j])
    vtkilagr.SetTuple1(j, ilagr.iloc[j])    
    vtkexist.SetTuple1(j, exist.iloc[j])    
    
#    verts.InsertNextCell(1)
#    verts.InsertCellPoint(j)
    


#create connectivity
#get unique IDs
uIds = pd.unique(ilagr)

lines = vtk.vtkCellArray()

for pid in progressbar.progressbar(uIds):
    df_id = df[df['ILAGR']==pid]
    if df_id.shape[0]>1:
        pt0id = df_id.index[0] #use the fact that points in VTK have the same order as index
        for j in range(1,df_id.shape[0]):
            pt1id = df_id.index[j]
            
            lines.InsertNextCell(2)
            lines.InsertCellPoint(pt0id)
            lines.InsertCellPoint(pt1id)
            
            pt0id = pt1id


pd = vtk.vtkPolyData()
pd.SetPoints(pts)
pd.GetPointData().AddArray(types)
pd.GetPointData().AddArray(times)
pd.GetPointData().AddArray(vtkilagr)
pd.GetPointData().AddArray(vtkexist)
pd.SetLines(lines)




print('Writing ', main_filename)
wr= vtk.vtkPolyDataWriter()
wr.SetInputData(pd)
wr.SetFileTypeToBinary()
wr.SetFileName(main_filename)
wr.Write()
