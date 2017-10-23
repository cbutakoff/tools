#!/home/costa/anaconda3/envs/py36/bin/python3.6

# coding: utf-8

# In[31]:
from mpi4py import MPI
import numpy as np
import vtk


mesh_filename = 'mesh_x4_hexfibers.vtk'
points_filename = 'x4_activation_points.ply'

out_table = 'activation.hex.in'
outmesh_filename = mesh_filename.split('.')[0]+'_act.vtk'



#alya params for activation points
radius = 10
CURDENSITY= -80
LAPSE=    0.002 #[s]



save_mesh = True


#MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print(f'Running process {rank}/{size}\n')

Ncpu = size;



# In[32]:


rdm = vtk.vtkDataSetReader()
rdm.SetFileName(mesh_filename)
rdm.Update()
mesh = rdm.GetOutput()


# In[33]:


rdp = vtk.vtkPLYReader()
rdp.SetFileName(points_filename)
rdp.Update()
pts = rdp.GetOutput()


# In[34]:


loc = vtk.vtkPointLocator()
loc.SetDataSet(mesh)
loc.BuildLocator()


#subdivide the index array (not continuous)
indices = list(range(pts.GetNumberOfPoints()))
ranges = [indices[i::Ncpu] for i in range(Ncpu)]

print(f'Number of inidces:{len(indices)}, number of ranges:{len(ranges)}\n')

index_range = ranges[rank]  #this range will be used on this node

activation_ptids = []  #local storage for activation point ids

#extract the neighbors on each nodes
print(f'Running on {len(index_range)} elements: {index_range[0]}, {index_range[1]}, ..., {index_range[-1]}\n')
for ptid in index_range:
    cp = vtk.vtkIdList()
    print(f'{ptid}/{pts.GetNumberOfPoints()}\n')
    #sys.stdout.flush()
    loc.FindPointsWithinRadius(radius, pts.GetPoint(ptid), cp)
    if cp.GetNumberOfIds()>0:
        activation_ptids = activation_ptids + [cp.GetId(i)+1 for i in range(cp.GetNumberOfIds())]

#in the root node 
if rank==0:
    #initialize global storage
    activation_ptids_global = activation_ptids
    
    #receive activation ids from all the nodes (suboptimal as it waits for every process)
    for i in range(1, Ncpu):
        buffer_size = np.zeros(1, dtype=np.uint64) #get buffer length
        comm.Recv(buffer_size, source=i)
        buffer = np.zeros(buffer_size[0], dtype=np.uint64)
        comm.Recv(buffer, source=i) #get the buffer
        activation_ptids_global = activation_ptids_global+list(buffer)
        print(f'CPU{rank}: received {buffer_size} elements')

else: #any other node
    length = np.zeros(1, dtype = np.uint64)
    length[0] = len(activation_ptids)
    comm.Send(length, dest=0) #send length

    buffer = np.array(activation_ptids, dtype = np.uint64)
    comm.Send(buffer, dest=0) #send buffer
    print(f'CPU{rank}: sent {length} elements')


#finalize
if rank == 0:
    #remove duplicate nodes
    activation_ptids_global = list(set(activation_ptids_global))


    # In[38]:


    f = open(out_table,'w')
    for ptid in activation_ptids:
        f.write(f'{ptid} {CURDENSITY} 0 {LAPSE}\n')

    f.close()


    # In[37]:


    #save the points with the oriinal mesh for validation
    if save_mesh:
        a = vtk.vtkShortArray()
        a.SetName('activation_points')
        a.SetNumberOfComponents(1)
        a.SetNumberOfTuples(mesh.GetNumberOfPoints())
        a.Fill(0);

        for ptid in activation_ptids:
            a.SetTuple1(ptid-1,1)

        mesh.GetPointData().AddArray(a)
        
        wr = vtk.vtkUnstructuredGridWriter()
        wr.SetFileName(outmesh_filename)
        wr.SetInputData(mesh)
        wr.Write()

