import os
import shutil
import pydicom
import glob
import numpy as np
import pandas as pd
import SimpleITK as sitk
from progressbar import progressbar

slice_files = glob.glob("*.dcm")

#order slices
times = []
spacing = np.zeros(3)
slice_files = sorted(slice_files)
inst = []
pos =[]
slice_shape = None
pixeltype = None

n_time_instants = 24
n_slices = int( len(slice_files)/n_time_instants )



for slicefile in slice_files:
    ds = pydicom.read_file(slicefile)

    pixels = np.array(ds.pixel_array)
    slice_shape = pixels.shape
    pixeltype = pixels.dtype
#    print(slicefile)
#    print("ds.SliceThickness = ",ds.SliceThickness)
#    print("ds.ImageOrientationPatient = ",ds.ImageOrientationPatient)
#    print("ds.ImagePositionPatient = ",   ds.ImagePositionPatient)
#    print("ds.PixelSpacing = ",ds.PixelSpacing)
    spacing[0:2] = [float(x) for x in ds.PixelSpacing]
    spacing[2] = float(ds.SpacingBetweenSlices)
    times.append( ds.TriggerTime ) #strings, intentionally
    inst.append(int(ds.InstanceNumber))
    pos.append(float(ds.SliceLocation))

    
df = pd.DataFrame({'Time':times, 'InstanceNumber':inst, 'File':slice_files, 'Position':pos})
df1 = df.sort_values( by = 'InstanceNumber', ascending=True )
df1['TimeInstant'] = np.tile( np.arange(n_time_instants), n_slices )

for i in progressbar( range(n_time_instants) ): 
    subdf = df[df1['TimeInstant']==i].sort_values( by = 'Position', ascending=True )

    stack = np.zeros( list(slice_shape)+[subdf.shape[0]] , dtype = pixeltype)

    time = subdf['Time'].iloc[0]
    for j in range(subdf.shape[0]):
        ds = pydicom.read_file( subdf['File'].iloc[j] )
        stack[:,:,j] = ds.pixel_array

    #filename = f"im_{time:06.0f}.nrrd"
    filename = f"im_{i:03d}.nrrd"
    img = sitk.GetImageFromArray(stack, isVector=False)
    img.SetSpacing(spacing[::-1])
    sitk.WriteImage(img, filename)
