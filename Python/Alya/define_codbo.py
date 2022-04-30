""" Read all boundary labels and identify label intersections for alya on nodes """
import vtk
import numpy as np
import sys
from progressbar import progressbar
import pprint


mesh_filename    = sys.argv[1]
bc_info_filename = sys.argv[2]
arrayname        = sys.argv[3]

#extract BC codes and names
with open(bc_info_filename, 'r') as ff:
    lines = ff.readlines()

bc_line = None
for i in range(len(lines)):
    if 'Boundary codes' in lines[i]:
        bc_line = i+2
        break

bcs = {}
for i in range(bc_line,len(lines)):
    dd = lines[i].split('-')
    bcs[int(dd[-1])] = lines[i].strip()

print('Boundaries :')
pprint.pprint(bcs)

print('Read mesh ',mesh_filename)
rd = vtk.vtkXMLPolyDataReader()
rd.SetFileName(mesh_filename)
rd.Update()

mesh = rd.GetOutput()
ss   = mesh.GetCellData().GetArray(arrayname)

alya_boundaries = set()
cellids = vtk.vtkIdList()
for i in progressbar(range(mesh.GetNumberOfPoints())):
    mesh.GetPointCells( i, cellids	)
    labels = set() 		
    for j in range(cellids.GetNumberOfIds()):
        labels.add( int(ss.GetTuple1(cellids.GetId(j)) ) )
  
    alya_boundaries.add( tuple(sorted(list(labels))) )

print('Alya codes:')
alya_boundaries = sorted(list(alya_boundaries))
#pprint.pprint(alya_boundaries)
for cc in alya_boundaries:
    name = ''
    if len(cc)==1:
        if cc[0] in bcs.keys():
            name = bcs[cc[0]]

    print(' & '.join([str(x) for x in cc]).ljust(15,' ') + " 000 0.0 0.0 0.0    $ "+name)



