import vtk
import sys
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from progressbar import progressbar

input_image_filename = sys.argv[1]
output_vtu_filename  = sys.argv[2]
labels_to_extract = [3]


rd = vtk.vtkMetaImageReader()
rd.SetFileName(sys.argv[1])
rd.Update()

half_voxel_size = -np.array(rd.GetDataSpacing())/2.0

#transform the image to explicit structured grid, voxes are nodes
grid = vtk.vtkImageToStructuredGrid()
grid.SetInputConnection( rd.GetOutputPort() )
grid.Update()

t = vtk.vtkTransform()
t.Identity()
t.Translate(half_voxel_size.tolist())

#shift the grid by half voxel, os that the voxels are in the centers of each element
tf = vtk.vtkTransformFilter()
tf.SetInputConnection(grid.GetOutputPort())
tf.SetTransform(t)
tf.Update()

#sample the image at the centers of the elements
cc = vtk.vtkCellCenters()
cc.SetInputConnection(tf.GetOutputPort())
cc.Update()


interpolator = vtk.vtkImageInterpolator()
interpolator.SetInterpolationModeToNearest()
probe = vtk.vtkImageProbeFilter()
probe.SetSourceConnection(rd.GetOutputPort())
probe.SetInputConnection(cc.GetOutputPort())
probe.SetInterpolator(interpolator)
probe.Update()

#pass the sampled intensitis to array "Label"
#and associate it to the structured grid
probe_array = vtk_to_numpy( probe.GetOutput().GetPointData().GetScalars() ).astype(np.uint8)
pixel_data = numpy_to_vtk( probe_array )
pixel_data.SetName('Label')


mesh = tf.GetOutput()
mesh.GetCellData().AddArray(pixel_data)

#threshold and combine the elements corresponding to each label
th = vtk.vtkThreshold()
th.SetInputData(mesh)
th.SetInputArrayToProcess(0,0,0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, 'Label')
th.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)

append = vtk.vtkAppendDataSets()
append.MergePointsOn()
for label in progressbar(labels_to_extract):
    th.SetUpperThreshold(label+0.1)
    th.SetLowerThreshold(label-0.1)
    th.Update()
    subpart = vtk.vtkUnstructuredGrid()
    subpart.DeepCopy( th.GetOutput() )
    append.AddInputData(subpart)


append.Update()

wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputConnection(append.GetOutputPort())
wr.SetFileName(sys.argv[2])
wr.EncodeAppendedDataOff()
wr.Write()
