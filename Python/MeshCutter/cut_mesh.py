
# coding: utf-8

# In[1]:


import vtk
import numpy  as np

input_mesh = 'patch.vtk'
output_mesh = 'patch_regions.vtk'
line_textfile = 'line3.txt' 
array_name = 'Regions'

#The line is specified by ordered line points
#for example 0 2 3 4 0 -- connected line passing through the points 0 2 3 4 0 in that order 


# In[2]:


rd = vtk.vtkPolyDataReader()
rd.SetFileName(input_mesh)
rd.Update()
mesh = rd.GetOutput()


# In[3]:


#line points have to be ordered
lines = []
#np.loadtxt(line_textfile,dtype=int)
#line
with open(line_textfile, 'r') as f:
    for line in f:
        l = line.replace('\n','').strip()
        lines.append( np.array([int(x) for x in l.split(' ')]) )


# In[4]:


lines


# In[5]:


#extract connectivity
tri = np.zeros([mesh.GetNumberOfCells(),3], dtype=np.int64)

for i in range(tri.shape[0]):
    ids = mesh.GetCell(i).GetPointIds()
    for j in range(3): 
        tri[i,j] = ids.GetId(j)


# In[9]:


def find_triangles(p1_id, p2_id, tri):
    tt = (tri-p1_id) * (tri-p2_id)
    
    return np.where((tt==0).sum(axis=1)==2)[0]

def find_triangle_point_loc(p1_id, p2_id, tri):
    p1_loc = np.where((tri-p1_id)==0)
    p2_loc = np.where((tri-p2_id)==0)
    return [p1_loc[0][0], p2_loc[0][0]]

def find_celledge_neighbors(tri_id, tri):
    (p1_id, p2_id, p3_id) = tri[tri_id,:]
    t1 = find_triangles(p1_id, p2_id, tri)
    t2 = find_triangles(p1_id, p3_id, tri)
    t3 = find_triangles(p2_id, p3_id, tri)
    t = ( set(t1).union(set(t2)).union(set(t3)) ) -{tri_id}
    
    return list(t)

def triangle_common_edge(tri1, tri2):
    common_pts = set(tri1).intersection(set(tri2))
    if len(common_pts)<2:
        return {}
    else:
        return common_pts

    
    
def triangles_on_one_line(t1, t2, tri, line):
    edge = triangle_common_edge(tri[t1,:], tri[t2,:]);
    
    on_line = False;
    
    for i in range(line.shape[0]-1):
        segm = {line[i], line[i+1]}
        if len(segm-edge)==0:
            on_line = True
            break
            
    return on_line


def triangles_on_any_line(t1, t2, tri, lines):
    on_line = False;
    for line in lines:
        on_line = triangles_on_one_line(t1, t2, tri, line)
        if on_line:
            break;
            
    return on_line


# In[10]:


trilabel = np.zeros(mesh.GetNumberOfCells(), dtype=np.int64)


region_id = 0
for i in range(mesh.GetNumberOfCells()):
    if trilabel[i]==0:
        tri_stack = [i]  #traingles to process
        region_id = region_id+1
    
        while tri_stack: #whle not empty
            tri_id = tri_stack.pop()
            #print('Triangle ', tri_id)

            if(trilabel[tri_id]==0): #if not labeled yet
                trilabel[tri_id]=region_id
                neighb = find_celledge_neighbors(tri_id, tri)

                #print('Neighbors ', neighb)

                for j in range( len(neighb) ):
                    if trilabel[neighb[j]]==0:
                        #see if the triangles tri_id and neighb[j] are on the different sides of the line
                        #i.e. if they share any pair of points of the line
                        if not triangles_on_any_line(tri_id, neighb[j], tri, lines):
                            tri_stack.append(neighb[j])


# In[11]:


array = vtk.vtkIdTypeArray()
array.SetName(array_name)
array.SetNumberOfComponents(1)
array.SetNumberOfTuples(trilabel.shape[0])

for i in range(trilabel.shape[0]):
    array.SetTuple1(i, trilabel[i])
    
mesh.GetCellData().AddArray(array)

wr = vtk.vtkPolyDataWriter()
wr.SetFileName(output_mesh)
wr.SetInputData(mesh)
wr.Write()

