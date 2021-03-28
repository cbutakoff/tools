from svgpathtools import svg2paths 
#0,0 -- top left
paths, attributes = svg2paths('drawing.svg')         



import matplotlib.pyplot as plt
import numpy as np

def ecg_from_path( path ):
    xy = []
    for line in path:
        xy.append( [np.real(line.start),  -np.imag(line.start)] )
    
    line = paths[-1]
    xy.append( [np.real(line.end),  -np.imag(line.end)] )

    xy = np.array(xy)
    
    baseline = np.median(xy[:,1])
    xy[:,1] = xy[:,1] - baseline
    xy[:,0] = (xy[:,0]-xy[:,0].min())*0.2/5
    xy[:,1] = xy[:,1]*0.5/5
    
    return xy



ecg = ecg_from_path( paths[-1] )
plt.plot( ecg[:,0], ecg[:,1])
plt.axis('equal')
plt.xlabel('Time[s]')
plt.ylabel('V[mV]')

plt.show()
