import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import vtk
from scipy import interpolate
from io import StringIO

pv_filename  = sys.argv[1]
ecg_filename = sys.argv[2]
ensicase     = sys.argv[3]
output_csv   = sys.argv[4]
ndeci        = 5

print('Reading ',pv_filename)
df = pd.read_csv( pv_filename, delim_whitespace=True, header=None, skiprows=17, names=['I','Time','cav1','cycle1','LVphase','LVV','LVP','LVWK','cav2','cycle2','RVphase','RVV','RVP','RVWK'] )
df['Time'] = np.round( df['Time'], ndeci )
print(df.head())


print('Reading ',ecg_filename)
with open( ecg_filename, 'r' ) as fff:
  lines = fff.readlines()

header_line = None
for i in range(len(lines)):
  if ('Step' in lines[i]) & (not '--|' in lines[i]):
    header_line = i
    break

print(header_line)
lines = lines[header_line:]
lines[0] = lines[0].replace('#',' ')
print(lines[:5])

df_ecg = pd.read_csv( StringIO("\n".join(lines)), delim_whitespace=True )
df_ecg['Time'] = np.round( df_ecg['Time'], ndeci )
print(df_ecg.head())
del lines


print('Reading ',ensicase)
rdr = vtk.vtkEnSightGoldBinaryReader()
rdr.SetCaseFileName(ensicase)
rdr.ReadAllVariablesOn()
rdr.Update()

displ_timeline = None
for i in range(rdr.GetNumberOfVariables()):
    if rdr.GetDescription(i)=='DISPL':
        displ_timeline = i
        break

assert not displ_timeline is None, 'DISPL not found in ensi case'

ensi_times = [np.round(rdr.GetTimeSets().GetItem(displ_timeline).GetTuple1(i),ndeci) for i in range(rdr.GetTimeSets().GetItem(displ_timeline).GetNumberOfTuples())]

x = df['Time'].to_numpy()
y = df[['LVV','LVP','LVWK','RVV','RVP','RVWK']].to_numpy()
f = interpolate.interp1d( x, y, axis=0, kind='linear',bounds_error=False,fill_value="extrapolate")
ensi_vals1 = pd.DataFrame( f(ensi_times), columns = ['LVV','LVP','LVWK','RVV','RVP','RVWK'] )
ensi_vals1['Time'] = ensi_times



x = df['Time']
y = df[['LVphase','RVphase']].to_numpy()
f = interpolate.interp1d( x, y, axis=0 , kind='nearest',bounds_error=False,fill_value="extrapolate")
ensi_vals2 = pd.DataFrame( f(ensi_times), columns = ['LVphase','RVphase'] )
ensi_vals2['Time'] = ensi_times


x = df_ecg['Time'].to_numpy()
y = df_ecg.drop(['Step','Time'], axis=1)
f = interpolate.interp1d( x, y.to_numpy(), axis=0, kind='linear',bounds_error=False,fill_value="extrapolate")
ensi_vals3 = pd.DataFrame( f(ensi_times), columns = y.columns )
ensi_vals3['Time'] = ensi_times




merged = pd.merge( pd.merge(ensi_vals1,ensi_vals2,on='Time'), ensi_vals3,on='Time')
time_column = merged.pop('Time')
merged.insert(0, 'Time', time_column)
merged.to_csv(output_csv,index=False)

