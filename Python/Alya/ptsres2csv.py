import vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np
from os.path import isfile, join
from os import listdir
from sys import exit
import re

import sys

filenamein = sys.argv[1]
filenameout = sys.argv[2]


def read_file( filename ):
    dt = np.dtype([('T', 'f8'), ('ILAGR', 'i8'), ('ITYPE', 'i8'), ('EXIST', 'i8'),\
     ('COORX', 'f8'), ('COORY', 'f8'), ('COORZ', 'f8'),\
     ('VELOX', 'f8'), ('VELOY', 'f8'), ('VELOZ', 'f8'),\
     ('DTK', 'f8'), ('CD', 'f8') ])

    with open(filename, mode='rb') as file:
        hdr = file.read(255)
        data = np.fromfile(file, dtype=dt)
    
    df = pd.DataFrame.from_records(data)

    return df
    

read_file( filenamein ).to_csv(filenameout, index=False)




