import vtk
import sys
import numpy as np
import os
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from progressbar import progressbar

input_image_filename = sys.argv[1]
output_vtu_filename  = sys.argv[2]
labels_to_extract = [3]


rd = vtk.vtkMetaImageReader()
rd.SetFileName(sys.argv[1])
rd.Update()

bounds = np.array( rd.GetOutput().GetBounds() )

#pad image by 1 voxel on all sides
extent = np.array( rd.GetOutput().GetExtent() )
extent[0::2] = extent[0::2]-1
extent[1::2] = extent[1::2]+1

pad = vtk.vtkImageConstantPad()
pad.SetInputConnection( rd.GetOutputPort() )
pad.SetConstant(0)
pad.SetOutputWholeExtent(extent)
pad.Update()

#Running vtkImageToStructuredGrid on vtkImageConstantPad output directaly fails to reconstrcut the grid correctly. vtkImageToStructuredGrid does not like negative extent
image = pad.GetOutput()
spacing = image.GetSpacing()
origin = np.array(image.GetOrigin())
image.SetOrigin( origin - spacing ) #shift by 1 voxel
extent[1::2] = extent[1::2]-extent[0::2]
extent[0::2] = 0
image.SetExtent( extent )

half_voxel_size = -np.array(image.GetSpacing())/2.0

#transform the image to explicit structured grid, voxes are nodes
grid = vtk.vtkImageToStructuredGrid()
grid.SetInputData( image )
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
probe.SetSourceData(image)
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

#clip the result to image bounds
box = vtk.vtkBox()
box.SetBounds(bounds)

clip = vtk.vtkClipDataSet()
clip.SetInputConnection(append.GetOutputPort())
clip.SetClipFunction (box)
clip.InsideOutOn()
clip.Update()


wr = vtk.vtkXMLUnstructuredGridWriter()
wr.SetInputConnection(clip.GetOutputPort())
wr.SetFileName(sys.argv[2])
wr.EncodeAppendedDataOff()
wr.Write()
