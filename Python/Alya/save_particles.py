import vtk
import pandas as pd

main_filename = 'fluidda.pvd'

df = pd.read_csv('fluidda.pts.res', delim_whitespace=True)

unique_times = df['T'].unique()
print(unique_times)

with open(main_filename,'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write(' <Collection>\n')
    


    for i in range(len(unique_times)):
        pts = vtk.vtkPoints()
        df1 = df[df['T']==unique_times[i]]
        x = df1['COORX']
        y = df1['COORY']
        z = df1['COORZ']
        pts.SetNumberOfPoints(df1.shape[0])
    
        for j in range(df1.shape[0]):         
            pts.SetPoint(j, (x.iloc[j],y.iloc[j],z.iloc[j]))
    
        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
    
        wr= vtk.vtkXMLPolyDataWriter()
        wr.SetInputData(pd)
        wr.SetDataModeToBinary()
        filename = "pts_{:04d}.vtp".format(i)
        wr.SetFileName(filename)
        wr.Write()

        f.write(f'<DataSet timestep="{unique_times[i]}" group="" part="0"\n')
        f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')
