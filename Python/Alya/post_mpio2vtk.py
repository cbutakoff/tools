import sys
import re
import shutil
import os
import numpy as np
import pandas as pd

main_pvd_file = sys.argv[1]
output_pvd_file = sys.argv[2]
ecg_file = ''
ecg_file_out = ''
if len(sys.argv)>3:
    ecg_file = sys.argv[3]
    ecg_file_out = sys.argv[4]


variables2keep = ['CALCI','QNETO']
timesteps = []

with open(main_pvd_file, 'r') as fin:
    with open(output_pvd_file, 'w') as fout:
        for line in fin:
            if "DataSet" in line:
                pattern = r"[\'\"]([A-Za-z0-9_\./\\-]*)[\'\"]"                
                m = re.findall(pattern, line)    
                filename = ''
#                for text in m:
#                    if 'pvtu' in text.lower():
#                        filename = text
#                        break
                #m should have only 2 values: time and filename
                timesteps.append(round(float(m[0]),10))
                filename = m[1]
                

                keepfile = False
                with open(filename) as pvtu:
                    for pvtu_line in pvtu:
                        if("PDataArray" in pvtu_line):
                            for var in variables2keep:
                                if var in pvtu_line:
                                    keepfile = True
                                    break
                
                if keepfile:
                    fout.write(line)
#                else:
#                    os.remove(filename)
#                    shutil.rmtree(filename[:-5], ignore_errors=True)
                
            else:
                fout.write(line)                



#calculate ECG
if len( ecg_file )>0:
    print('reading')
    df = pd.read_csv(ecg_file, header=None, delim_whitespace=True, names=['iter','t','0','LL','RA','LA'])
    df1 = df.loc[0].to_frame().transpose()
    if df1.loc[0,'iter'] != 0:
        df1.loc[0,'iter'] =0
        df1.loc[0,'t'] =0
        df = df1.append(df, ignore_index=True)

    print(df.head())
    df['L1'] = df['LA']-df['RA']
    df['L2'] = df['LL']-df['RA']
    df['L3'] = df['LL']-df['LA']

    l1_new = np.interp(timesteps, df['t'], df['L1'] )   # use interpolation function returned by `interp1d`
    l2_new = np.interp(timesteps, df['t'], df['L2'])   # use interpolation function returned by `interp1d`
    l3_new = np.interp(timesteps, df['t'], df['L3'])   # use interpolation function returned by `interp1d`

    df_interp = pd.DataFrame( {'Time':timesteps, 'L1':l1_new, 'L2':l2_new, 'L3':l3_new} )
    df_interp.to_csv(ecg_file_out, index=False)
