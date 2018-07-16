
# coding: utf-8

# In[1]:


import vtk
import numpy as np
import os
import progressbar
import pandas


# In[14]:


single_file = 'fensap.dom.geo'
elemtype_file = 'fluidda.elem.type'
coords_file = 'fluidda.coord'
elem_file = 'fluidda.elem'
boundary_file = 'fluidda.boundary'
boundary_label_file = 'fensap.fix.dat'

vtk_boundary_file = 'boundary.vtk'
vtk_volmesh_file = 'volmesh.vtk'


# In[18]:


split_geometryfile = True


# In[22]:


def ExtractData(f, filename, pb):
    with open(filename,'w') as fout:
        while True:
            line=f.readline()
            if not line: break
            if 'END_' in line: break
            fout.write(' '.join(line.split())+'\n')
            pb.update(bar.value+len(line))
    

if split_geometryfile:   
    #for the progressbar
    statinfo = os.stat(single_file)
    filesize = statinfo.st_size   
    
    with progressbar.ProgressBar(max_value=filesize) as bar:
        with open(single_file, 'r') as f:
            while True:
                line=f.readline()
                if not line: break
                bar.update(bar.value+len(line))

                line_start = line.strip()[0:5].upper()
                if line_start == 'TYPES':
                    print('Found types')
                    ExtractData(f, elemtype_file, bar)
                elif  line_start == 'ELEME':
                    print('Found elements')
                    ExtractData(f, elem_file, bar)
                elif  line_start == 'COORD': 
                    print('Found coordinates')
                    ExtractData(f, coords_file, bar)
                elif  line_start == 'BOUND':
                    print('Found boundary')
                    ExtractData(f, boundary_file, bar)
            
            
    


# In[8]:


#save boundary to vtk
print('Reading coordinates')
coords = pandas.read_csv(coords_file, header=None, delim_whitespace=True)
print('Reading boundary codes')
codes = pandas.read_csv(boundary_label_file, header=None, delim_whitespace=True)


# In[9]:


print('Passing coordinates to vtk')

pts = vtk.vtkPoints()
pts.SetNumberOfPoints(coords.shape[0])

vv = coords.values

for i in progressbar.progressbar(range(coords.shape[0])):
    pts.SetPoint(i, vv[i,1:])
    


# In[25]:


print('Passing boundaries to vtk')

boundaries = vtk.vtkCellArray()

with open(boundary_file,'r') as f:
    for line in f:
        data = line.strip().split()
        
        boundaries.InsertNextCell(len(data)-2)
        for i in range(1,len(data)-1):
             boundaries.InsertCellPoint( int(data[i])-1 )           


# In[26]:


print('Creating vtk array of boundary codes')

arr = vtk.vtkShortArray()
arr.SetName('Code')
arr.SetNumberOfComponents(1)
arr.SetNumberOfTuples(boundaries.GetNumberOfCells())

vv  = codes.values

for i in range(arr.GetNumberOfTuples()):
    arr.SetTuple1(i, vv[i,1])


# In[28]:


print('Saving  vtk boundary to ',vtk_boundary_file)

pd = vtk.vtkPolyData()
pd.SetPoints(pts)
pd.SetPolys(boundaries)
pd.GetCellData().AddArray(arr)

wr = vtk.vtkPolyDataWriter()
wr.SetFileName(vtk_boundary_file)
wr.SetFileTypeToBinary()
wr.SetInputData(pd)
wr.Write()


# In[19]:


print('Reading volumetric elments')

#save the volumetric mesh to vtk
#read the elments
elements = vtk.vtkCellArray()
types = []

eltypes = {4:vtk.VTK_TETRA, 5:vtk.VTK_PYRAMID, 6: vtk.VTK_WEDGE, 8:vtk.VTK_HEXAHEDRON}

statinfo = os.stat(elem_file)
filesize = statinfo.st_size   
    
with progressbar.ProgressBar(max_value=filesize) as bar:
    with open(elem_file,'r') as f:
        for line in f:
            pp = line.strip().split()
            bar.update(bar.value+len(line))

            n = len(pp)-1 #1st is the id

            types += [eltypes[n]]

            elements.InsertNextCell(n)
            for ptid in pp[1:]:
                elements.InsertCellPoint(int(ptid)-1)
            
            


# In[20]:


ug = vtk.vtkUnstructuredGrid()
ug.SetPoints(pts)
ug.SetCells(types, elements)


# In[22]:


print('Saving  vtk vol mesh ',vtk_volmesh_file)

wr = vtk.vtkDataSetWriter()
wr.SetFileName(vtk_volmesh_file)
wr.SetInputData(ug)
wr.SetFileTypeToBinary()
wr.Write()

