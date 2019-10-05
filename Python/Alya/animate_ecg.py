import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

fig, ax = plt.subplots(figsize=(5,2))

ecg = pd.read_csv('dm-400-290.csv')
maxtime = 1.2
ecg = ecg[ecg.Time<=1.2]

line, = ax.plot(ecg.Time, ecg.L1)
point, = ax.plot(0,0,'or')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (ms)')
ax.set_title('L1')

def init():  # only required for blitting to give a clean slate.
    line.set_ydata(np.zeros(ecg.Time.shape[0]))
    return line,


def animate(t):
    mask = ecg.Time<=t
    line.set_xdata( ecg.Time[mask] )  # update the data.
    line.set_ydata( ecg.L1[mask] )  # update the data.
    point.set_xdata( ecg.Time[mask].values[-1] )
    point.set_ydata( ecg.L1[mask].values[-1] )
    return line,


ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=0.01, frames = ecg.Time )

# To save the animation, use e.g.
#
# ani.save("ecg.mp4")
#
# or
#
from matplotlib.animation import FFMpegWriter
writer = FFMpegWriter(fps=10, metadata=dict(artist='Me'), bitrate=1800)
ani.save("ecg.mp4", writer=writer)

plt.show()