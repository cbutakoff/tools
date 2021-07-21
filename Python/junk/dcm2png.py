import pydicom
import numpy as np
from PIL import Image
import sys

ds = pydicom.dcmread(sys.argv[1])

new_image = ds.pixel_array.astype(float)
scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0
scaled_image = np.uint8(scaled_image)
final_image = Image.fromarray(scaled_image)
final_image.save(sys.argv[2])
