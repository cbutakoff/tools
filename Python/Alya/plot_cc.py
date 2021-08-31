import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

print(sys.argv[1])
df = pd.read_csv( sys.argv[1], delim_whitespace=True, header=None, skiprows=17 )
print(df.head())


fig, ((ax0),(ax1), (ax2),(ax10)) = plt.subplots(4, 1)

ax0.plot ( df[1], df[4] )
ax10.plot( df[1], df[10] )

lns1 = ax1.plot( df[1], df[5],  '-', color='red', label = 'LV volume')
ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
lns2 = ax3.plot( df[1], df[6],  '-', color='blue', label = 'LV pressure')
lns3 = ax3.plot( df[1], df[7],  '-', color='green', label = 'LV WK pressure')

lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)



lns4 = ax2.plot( df[1], df[11], '-', color='red', label = 'RV volume')
ax4  = ax2.twinx()  # instantiate a second axes that shares the same x-axis
lns5 = ax4.plot( df[1], df[12], '--', color='blue', label = 'RV pressure')
lns6 = ax4.plot( df[1], df[13], '--', color='green', label = 'RV WK pressure')

# added these three lines
lns = lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs, loc=0)

plt.show()
