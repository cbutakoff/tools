#!/usr/bin/env python
# coding: utf-8

# In[66]:

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy.interpolate import interp1d

# In[2]:

casefile = sys.argv[1]
ecgfile = sys.argv[2]
output_csv = sys.argv[3]
imagefilename = sys.argv[4]



#read the timesteps of the INTRA
class BreakIt(Exception): pass
timelineid_found = False
timeline_started = False
reading_time_values = False;
times = []
nsteps=0
with open(casefile) as f:

    try:
        for line in f:
            #print ("Line :",line)
            if reading_time_values:
                for v in line.split():
                    #print( len(times), nsteps )
                    if len(times)>=nsteps:
                        raise BreakIt #break the loops, all timesteps extracted

                    times.append(float(v))


            if timelineid_found:
                line_low = line.lower()                
                if "time" in line_low:
                    if ( "set" in line_low ):
                        if int(line_low.split()[-1])==timeline_id:
                            timeline_started = True;

            if timeline_started:
                ll = line.lower()
                if "number" in ll:
                    if "steps" in ll:
                        nsteps = int(ll.split()[-1])
                        
                if "time" in ll:
                    if "values" in ll:
                        reading_time_values = True;

            if "INTRA" in line:
                timeline_id = int(line.split()[3])
                print('INTRA timeline id: ', timeline_id)
                timelineid_found = True;
        
    except BreakIt:
        pass
                

times = np.array(times)
print("Times: ", times)



#calculate ECG
df = pd.read_csv(ecgfile, header=None, delim_whitespace=True, names=['iter','t','0','LL','RA','LA'])
df1 = df.loc[0].to_frame().transpose()
if df1.loc[0,'iter'] != 0:
    df1.loc[0,'iter'] =0
    df1.loc[0,'t'] =0
    df = df1.append(df, ignore_index=True)

print(df.head())
df['L1'] = df['LA']-df['RA']
df['L2'] = df['LL']-df['RA']
df['L3'] = df['LL']-df['LA']

l1_new = np.interp(times, df['t'], df['L1'])   # use interpolation function returned by `interp1d`
l2_new = np.interp(times, df['t'], df['L2'])   # use interpolation function returned by `interp1d`
l3_new = l_interp(times, df['t'], df['L3'])   # use interpolation function returned by `interp1d`

df_interp = pd.DataFrame( {'Time':times, 'L1':l1_new, 'L2':l2_new, 'L3':l3_new} )
df_interp.to_csv(output_csv, index=False)


df_interp.plot(x='Time',y=['L1','L2','L3'])
plt.ylabel('mV')
plt.axhline(y=0, xmin=0, xmax=1, linestyle=':',c='r')


plt.savefig(imagefilename)

