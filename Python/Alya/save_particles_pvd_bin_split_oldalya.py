import vtk
from vtk.util.numpy_support import numpy_to_vtk
import pandas as pd
import  progressbar 
import multiprocessing as mp
import numpy as np
from os.path import isfile, join
from os import listdir, SEEK_END
from sys import exit
import re
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("case_name", help='Name of the alya task')
parser.add_argument("input_folder", help='Folder with input pts.res')
parser.add_argument("-o","--output_folder", default="./", help='Folder for the output vtps')
parser.add_argument("-n","--ncpus", help='Number of cpus to parallelize (1 -- run sequential)', default = 10, type=int)
parser.add_argument("-f","--force", action='store_true', required=False, help='Force overwriting vtps')
args = parser.parse_args()


case_name       = args.case_name
case_path       = args.input_folder
output_path     = args.output_folder 
ncpus           = args.ncpus
force_overwrite = args.force
main_filename   = f'{case_name}.pvd'
dat_file        = join(case_path, f'{case_name}.dat')
pts_filename0   = join(case_path, f'{case_name}.pts.res')
pts_filename    = join(case_path, f'{case_name}'+'.pts.{:08d}.res')


#extract the timestep from the dat
timestep = -1
with open(dat_file,'r') as f:
    for line in f:
        if "time_step_size" in line.lower():
            numbers = [float(s) for s in line.split() if s[0].isdigit()]
            timestep = numbers[-1]

if timestep<0:
    print ('Error extracting timestep')
    exit            

#extract filenumbers
filenumbers = []
for f in listdir(case_path):
    m = re.match(f'{case_name}.pts.(\d+).res',f) 
    if m:
        filenumbers.append( int(m.group(1)) )

filenumbers = sorted(filenumbers) 


def read_file( filename ):
    dt = np.dtype([('T', 'f8'), ('ILAGR', 'i8'), ('ITYPE', 'i8'), ('EXIST', 'i8'),\
     ('COORX', 'f8'), ('COORY', 'f8'), ('COORZ', 'f8'),\
     ('VELOX', 'f8'), ('VELOY', 'f8'), ('VELOZ', 'f8'),\
     ('DTK', 'f8'), ('CD', 'f8') ])

    with open(filename, mode='rb') as file:
        hdr = file.read(255)
        data = np.fromfile(file, dtype=dt)
    
    df = pd.DataFrame.from_records(data, exclude=['VELOX','VELOY','VELOZ','DTK', 'CD'])

    return df
    

def save_one_file_by_number(file_number):
    infilename = pts_filename.format(file_number) 
    outfilename = join( output_path, "pts_{:010d}.vtp".format(file_number) )

    if not is_vtp_complete(outfilename):
        save_one_file( infilename, outfilename )


def save_one_file(infilename, outfilename):

    df1 = read_file( infilename )

    if df1.shape[0]==0:
        return

    pts = vtk.vtkPoints()
    pts.SetData( numpy_to_vtk(df1[["COORX","COORY","COORZ"]].to_numpy(), array_type = vtk.VTK_FLOAT) )  

    types = numpy_to_vtk(df1["ITYPE"].to_numpy(), array_type = vtk.VTK_SHORT)
    types.SetName('ITYPE')


    ilagr = numpy_to_vtk(df1["ILAGR"].to_numpy(), array_type = vtk.VTK_INT) 
    ilagr.SetName('ILAGR')


    exist = numpy_to_vtk(df1["EXIST"].to_numpy(), array_type = vtk.VTK_SHORT)
    exist.SetName('EXIST')


    T = numpy_to_vtk(df1["T"].to_numpy(), array_type = vtk.VTK_FLOAT)
    T.SetName('T')

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.GetPointData().AddArray(types)
    pd.GetPointData().AddArray(ilagr)
    pd.GetPointData().AddArray(exist)
    pd.GetPointData().AddArray(T)

    wr= vtk.vtkXMLPolyDataWriter()
    wr.SetInputData(pd)
    wr.SetDataModeToBinary()
    wr.SetFileName(outfilename)
    wr.Write()


def read_last_block(in_file, block_size=1024):
    in_file.seek(-block_size, SEEK_END)
    return in_file.read(block_size).decode("utf-8") 


def is_vtp_complete(filename):
    result = False

    if not force_overwrite:
        pattern = "</VTKFile>"
        if isfile(filename):
            with open(filename,'rb') as ff:
                try:
                    result = pattern in read_last_block(ff, block_size=len(pattern)+100)
                except:
                    result = False
                    pass

    return result

#save 0-th timestep
vtp_filename = join( output_path, "pts_{:010d}.vtp".format(0) )
if not is_vtp_complete(vtp_filename):
    save_one_file( pts_filename0, vtp_filename )


if ncpus>1:
    bar = progressbar.ProgressBar(max_value=len(filenumbers))
    with mp.Pool(processes = ncpus) as p:
        for i, _ in enumerate( p.imap_unordered(save_one_file_by_number, filenumbers), 1 ):
            bar.update(i)
else:
    for filenumber in progressbar.progressbar(filenumbers):
        save_one_file_by_number(filenumber)


print('Generating main pvd file')
with open( join( output_path, main_filename), 'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1"\n')
    f.write('         byte_order="LittleEndian"\n')
    f.write('         compressor="vtkZLibDataCompressor">\n')
    f.write('  <Collection>\n')

    filename = "pts_{:010d}.vtp".format(0)
    if isfile( join( output_path, filename) ):
        f.write(f'<DataSet timestep="0" group="" part="0"\n')
        f.write(f'     file="{filename}"/>\n')


    for filenumber in progressbar.progressbar(filenumbers):
        filename = "pts_{:010d}.vtp".format(filenumber)
        if isfile( join( output_path, filename) ):
            f.write(f'<DataSet timestep="{np.round(filenumber*timestep,10)}" group="" part="0"\n')
            f.write(f'     file="{filename}"/>\n')
    

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')




