
# coding: utf-8

# In[1]:


import vtk
import numpy as np
from progressbar import progressbar, ProgressBar


# In[2]:


grdfile = 'fensap.grd'
vtkfile = 'wedge-vol.vtk'


# In[3]:


with open(grdfile,'r') as f:
    line = f.readline().split()
    npoints = int(line[0])
    
    line = f.readline().split()
    nelement_types = int(line[0])
    
    nelements = 0
    for k in range(nelement_types):
        line = f.readline().split()
        nelements = nelements + int(line[1])
        
    f.readline()  #1.0 1.0
    f.readline()  #fensap input file
    

    
    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(npoints)
    for i in progressbar(range(npoints)):
        line = f.readline().split()
        pts.SetPoint(i, float(line[0]), float(line[1]), float(line[2]))
        
    cells = vtk.vtkCellArray()
    celltypes = []
    c = 0
    with ProgressBar(max_value=nelements) as bar:    
        while True:
            line = f.readline()
            if line=='':
                break

            line = line.split()

            nverts = len(line)
            
            cells.InsertNextCell( nverts )
            for k in line:
                cells.InsertCellPoint( int(k)-1 )

            if nverts == 4:
                celltypes.append(vtk.VTK_TETRA)
            elif nverts == 8:
                celltypes.append(vtk.VTK_HEXAHEDRON)
            elif nverts == 6:
                celltypes.append(vtk.VTK_WEDGE)
            elif nverts == 5:
                celltypes.append(vtk.VTK_PYRAMID)
            elif nverts == 10:
                celltypes.append(vtk.VTK_PENTAGONAL_PRISM)
            elif nverts == 12:
                celltypes.append(vtk.VTK_HEXAGONAL_PRISM)
            else:
                print('Unknown cell type:', line)
            
                
            c=c+1
            bar.update(c)
        


# In[4]:


ug = vtk.vtkUnstructuredGrid()
ug.SetPoints(pts)
ug.SetCells(celltypes, cells)


# In[5]:


wr = vtk.vtkDataSetWriter()
wr.SetFileName(vtkfile)
wr.SetInputData(ug)
wr.SetFileTypeToBinary()
wr.Write()


# In[20]:




