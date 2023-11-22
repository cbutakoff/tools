
# coding: utf-8

# In[1]:


import progressbar
import pandas as pd
import glob
import matplotlib.pyplot as plt
from math import ceil
import sys
import numpy


density = 1.06  #g/cm3
flowrate_cm2s = False
set2plot = int(sys.argv[1])
meanp_column = 2
flowrate_column = 3
columns_per_table = 4 #deault number of cols in .set

# In[2]:


#find the results file
nsi_set_file = glob.glob('*.nsi.set', recursive=False)

nsi_set_file = nsi_set_file[0]
print(f'Using input from {nsi_set_file}')    



#read header
numsets = 0
with open(nsi_set_file, 'r') as f:
    data = f.readlines()    

time = 0
reading_set_table = False
flowrates = []
meanp = []
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
            if boundary == set2plot:   
                flowrates.append( row[flowrate_column] )
                meanp.append( row[meanp_column] )
        except:
            pass
        


df = pd.DataFrame({'Time':times, 'Flowrate':flowrates, 'Mean P':meanp}, dtype=float)


if flowrate_cm2s:
    df['Flowrate'] = df['Flowrate']/density

#df['Time'] = df['Time']*1000

# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(df['Time'],-df['Flowrate'], color="red")
ax.set_xlabel('Time[s]')

# set y-axis label
if flowrate_cm2s:
    ax.set_ylabel("Flowrate[cm3/s]",color="red")
else:
    ax.set_ylabel("Flowrate[g/s]",color="red")

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(df['Time'],df['Mean P'],color="blue")
ax2.set_ylabel("Mean Pressure [Baryes]",color="blue")

plt.title(f'Set {set2plot}')
fig.tight_layout()
plt.show()

# save the plot as a file
#fig.savefig('two_different_y_axis_for_single_python_plot_with_twinx.jpg',
#            format='jpeg',
#            dpi=100,
#            bbox_inches='tight')



    

