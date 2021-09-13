import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import vtk
from scipy import interpolate
from io import StringIO
from vtk.util.numpy_support import vtk_to_numpy
import os
import subprocess

ecg_filename = sys.argv[1]
ensicase     = sys.argv[2]
output_csv   = sys.argv[3]
ndeci        = 5
timeline_var = 'INTRA'
LA_label     = 'LA'
LL_label     = 'LL'
RA_label     = 'RA'

def CleanECGFile( filename_in : str, filename_out : str ):
    """ Remove ALya rubbish from the vin file, leave only header with column names"""
    #subprocess.call( (f"grep -vwE '(\-\-\|)' {filename_in} | sed 's/#/ /' > {filename_out}").split() )
    with open(filename_out, 'w') as fout:
        ps = subprocess.Popen( ("grep", "-vwE", "(\-\-\|)", filename_in), stdout=subprocess.PIPE)
        ps1 = subprocess.Popen( ("sed", "s/#/ /"), stdin=ps.stdout, stdout=fout)
        ps1.wait()


def find_timeset(casefileaname, timeline_var):
  variables_section = False
  timeline = None
  with open(casefileaname) as ff:
    for line in ff:
      if "VARIABLE" in line.upper():
        variables_section = True
      if variables_section:
        print (line)
        if f" {timeline_var} " in line.upper():
          data = np.array(line.upper().split())
          print(data)
          idx = np.where( data==timeline_var )[0]
          timeline = int( data[idx-1] )-1
          break

  return timeline
             

def resample_dataframe( df : pd.DataFrame, timeline : np.ndarray, timecolumnname : str ) -> pd.DataFrame:
    x = df[timecolumnname].to_numpy()
    df1 = df.loc[:, df.columns != timecolumnname]
    y = df1.to_numpy()
    f = interpolate.interp1d( x, y, axis=0, kind='linear',bounds_error=False,fill_value="extrapolate")
    ensi_vals1 = pd.DataFrame( f(timeline), columns = df1.columns )
    ensi_vals1.insert(0, 'Time', timeline)
    return ensi_vals1


ecgfile_clean = f"{ecg_filename}_clean"
CleanECGFile( ecg_filename, ecgfile_clean )
df_ecg = pd.read_csv( ecgfile_clean, delim_whitespace=True)
df_ecg['Time'] = np.round( df_ecg['Time'], ndeci )
os.remove(ecgfile_clean)
print(df_ecg.head())


print('Reading ',ensicase)
rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(ensicase)
rdr.ReadAllVariablesOn()
rdr.Update()

displ_timeline = find_timeset( ensicase, timeline_var )
assert not displ_timeline is None, 'INTRA not found in ensi case'

ensi_times = np.round( vtk_to_numpy(rdr.GetTimeSets().GetItem(displ_timeline)),  ndeci)


df_ecg = resample_dataframe( df_ecg, ensi_times, 'Time')
df_ecg["L1"] = df_ecg[LA_label] -  df_ecg[RA_label]
df_ecg["L2"] = df_ecg[LL_label] -  df_ecg[RA_label]
df_ecg["L3"] = df_ecg[LL_label] -  df_ecg[LA_label]

df_ecg.to_csv(output_csv,index=False)

