import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("vin_file", metavar="problem.exm.vin", help='File with volume integrals')
parser.add_argument("lead", metavar="L1|L2|L3", choices=['L1', 'L2', 'L3'], help='Lead to plot')
parser.add_argument("-ca","--calcium", action='store_true', required=False, help='Add to plot calcium')
parser.add_argument("-cl","--cyclelength", help='Cycle length in seconds', type=float, default = -1, required=False)
parser.add_argument("-c","--cycles", metavar="N", nargs='*', type=int, required=False, help='Cycle numbers to plot from 0. Default all')
parser.add_argument("-i","--image", metavar="image.png", required=False, help='Filename to save image (pg, jpg, svg,...)', default=None)
args = parser.parse_args()


columns = ['Iter','Time','rubbish','Ca','LL','RA','LA' ]


df = pd.read_csv(args.vin_file, header=None, delim_whitespace=True, comment='#')
cols = list(df.columns)
cols[0:2] = ['Iter', 'Time']
cols[-4:] = ['Ca','LL','RA','LA']

df.columns = cols

df['L1'] = df['LA']-df['RA']
df['L2'] = df['LL']-df['RA']
df['L3'] = df['LL']-df['LA']

if args.cycles is None:


    fig, ax = plt.subplots()

    ax.plot( df['Time'],df[args.lead],'r',label=args.lead)
    ax.set_ylabel('Voltage(mV)')
    ax.set_xlabel('Time(s)')

    if args.calcium:
        ax2 = ax.twinx()
        ax2.plot(df['Time'],df['Ca'],'b',label='Ca')
        ax2.plot([], [], 'r', label = args.lead)
        ax2.set_ylabel('Calcium transient (umol)')

else:
    df1 = pd.DataFrame()
    df1["Time"] = df[df['Time']<=args.cyclelength]['Time']

    if args.calcium:
        for i in args.cycles:
            ta = np.round(i*args.cyclelength,10)
            tb = np.round((i+1)*args.cyclelength,10)
            dd = df[(df['Time']>ta) & (df['Time']<=tb)]['Ca'].values
            df1[f'Ca{i}']=np.NaN
            df1[f'Ca{i}'].iloc[0:dd.shape[0]]=dd

        #print(df1.head())
        df1.plot(x='Time')
        plt.ylabel('Calcium transient (umol)')
    else:
        for i in args.cycles:
            ta = np.round(i*args.cyclelength,10)
            tb = np.round((i+1)*args.cyclelength,10)
            dd = df[(df['Time']>ta) & (df['Time']<=tb)][args.lead].values
            #print(df[(df['Time']>i*args.cyclelength) & (df['Time']<=(i+1)*args.cyclelength)]['Time'])
            df1[f'{args.lead}{i}']=np.NaN
            df1[f'{args.lead}{i}'].iloc[0:dd.shape[0]]=dd

        df1.plot(x='Time')
        plt.ylabel('Voltage(mV)')
  

plt.legend()


if not args.image is None:
    plt.savefig(args.image)
else:
    plt.show()

