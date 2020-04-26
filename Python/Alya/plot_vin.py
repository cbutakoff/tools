import numpy
import matplotlib.pyplot as plt
import pandas as pd
import sys

print('Parameters: inputfile.exm.vin L1|L2|L3 [image_filename.ext]')
print('Optional: image_filename.ext -- ext = png, jpg, svg,...')


columns = ['Iter','Time','rubbish','Ca','LL','RA','LA' ]

lead2plot = sys.argv[2]
saveimage = ''
if len(sys.argv)>3:
    saveimage = sys.argv[3]

df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=columns, comment='#')

print(df.head())
df['L1'] = df['LA']-df['RA']
df['L2'] = df['LL']-df['RA']
df['L3'] = df['LL']-df['LA']

fig, ax = plt.subplots()

ax.plot( df['Time'],df[lead2plot],'r',label=lead2plot)
#ax.plot( df['Time'],df['L2'],'g',label='L2')
#ax.plot( df['Time'],df['L3'],'b',label='L3')
ax.set_ylabel('Voltage(mV)')
ax.set_xlabel('Time(s)')

ax2 = ax.twinx()
ax2.plot(df['Time'],df['Ca'],'b',label='Ca')
ax2.plot([], [], 'r', label = lead2plot)
ax2.set_ylabel('Calcium transient (umol)')

plt.legend()


if len(saveimage)>0:
    plt.savefig(saveimage)
else:
    plt.show()

