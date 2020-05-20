
# coding: utf-8

# In[1]:


import progressbar
import pandas as pd
import glob
import matplotlib.pyplot as plt
from math import ceil
import sys
import numpy
import os


density = 1.06  #g/cm3
flowrate_cm2s = False
set2plot = [2,3]
meanp_column = 2
flowrate_column = 3
columns_per_table = 4 #deault number of cols in .set



def extract_flowrates(path):
    #find the results file
    nsi_set_file = glob.glob(os.path.join(path, '*.nsi.set'), recursive=False)

    nsi_set_file = nsi_set_file[0]
    print(f'Using input from {nsi_set_file}')    



    #read header
    numsets = 0
    with open(nsi_set_file, 'r') as f:
        data = f.readlines()    

    time = 0
    reading_set_table = False
    flowrates0 = []
    flowrates1 = []
    meanp0 = []
    meanp1 = []
    times = []
    for line in data:
        if numsets==0:
            if 'NUMSETS' in line:
                numsets = int(line.split(':')[1])

        if 'Time' in line:
            time = float(line.split(' ')[-1])
            times.append(time)
        elif line[0]=='#':
            pass
        else:
            try:
                row = line.strip().split()
                boundary = int(row[0])
                if boundary == set2plot[0]:   
                    flowrates0.append( row[flowrate_column] )
                    meanp0.append( row[meanp_column] )
                elif boundary == set2plot[1]:
                    flowrates1.append( row[flowrate_column] )
                    meanp1.append( row[meanp_column] )
            except:
                pass
            


    df = pd.DataFrame({'Time':times, f'Flowrate_{set2plot[0]}':flowrates0, f'P_{set2plot[0]}':meanp0, f'Flowrate_{set2plot[1]}':flowrates1, f'P_{set2plot[1]}':meanp1}, dtype=float)


    return df



df1 = extract_flowrates(sys.argv[1])
df2 = extract_flowrates(sys.argv[2])
print("===================================================================")
print(f"{sys.argv[1]} stats: ", df1.describe())
print("===================================================================")
print(f"{sys.argv[2]} stats: ", df2.describe())
print("===================================================================")
print(f"{sys.argv[2]}-{sys.argv[1]} stats: ", df2.describe()-df1.describe())


#df['Time'] = df['Time']*1000

# create figure and axis objects with subplots()
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)

axs[0,0].plot(df1['Time'],df1[f'Flowrate_{set2plot[0]}'],label=f'case:{sys.argv[1]}')
axs[0,0].plot(df2['Time'],df2[f'Flowrate_{set2plot[0]}'],label=f'case:{sys.argv[2]}')
axs[0,0].set_title(f'FR_{set2plot[0]}')

axs[1,0].plot(df1['Time'],df1[f'Flowrate_{set2plot[1]}'],label=f'case:{sys.argv[1]}')
axs[1,0].plot(df2['Time'],df2[f'Flowrate_{set2plot[1]}'],label=f'case:{sys.argv[2]}')
axs[1,0].set_title(f'FR_{set2plot[1]}')

axs[0,1].plot(df1['Time'],df1[f'P_{set2plot[0]}'],label=f'case:{sys.argv[1]}')
axs[0,1].plot(df2['Time'],df2[f'P_{set2plot[0]}'],label=f'case:{sys.argv[2]}')
axs[0,1].set_title(f'P_{set2plot[0]}')

axs[1,1].plot(df1['Time'],df1[f'P_{set2plot[1]}'],label=f'case:{sys.argv[1]}')
axs[1,1].plot(df2['Time'],df2[f'P_{set2plot[1]}'],label=f'case:{sys.argv[2]}')
axs[1,1].set_title(f'P_{set2plot[1]}')

# make a plot
plt.legend()
fig.tight_layout()
plt.show()



    

