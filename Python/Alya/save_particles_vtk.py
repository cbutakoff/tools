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

pts.SetNumberOfPoints(df.shape[0])

types = vtk.vtkShortArray()
types.SetNumberOfComponents(1)
types.SetNumberOfTuples(pts.GetNumberOfPoints())
types.SetName('Type')

times = vtk.vtkFloatArray()
times.SetNumberOfComponents(1)
times.SetNumberOfTuples(pts.GetNumberOfPoints())
times.SetName('Time')



for j in progressbar.progressbar(range(df.shape[0])):         
    pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))
    types.SetTuple1(j, tp.iloc[j])
    times.SetTuple1(j, T.iloc[j])
    
pd = vtk.vtkPolyData()
pd.SetPoints(pts)
pd.GetPointData().AddArray(types)
pd.GetPointData().AddArray(times)


print('Writing ', main_filename)
wr= vtk.vtkPolyDataWriter()
wr.SetInputData(pd)
wr.SetDataModeToBinary()
wr.SetFileName(main_filename)
wr.Write()


