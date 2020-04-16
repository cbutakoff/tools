import sys
import re
import shutil
import os

main_pvd_file = sys.argv[1]
output_pvd_file = sys.argv[2]
variables2keep = ['CALCI','QNETO']

with open(main_pvd_file, 'r') as fin:
    with open(output_pvd_file, 'w') as fout:
        for line in fin:
            if "DataSet" in line:
                pattern = r"[\'\"]([A-Za-z0-9_\./\\-]*)[\'\"]"                
                m = re.findall(pattern, line)    
                filename = ''
                for text in m:
                    if 'pvtu' in text.lower():
                        filename = text
                        break

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
                else:
                    os.remove(filename)
                    shutil.rmtree(filename[:-5], ignore_errors=True)
                
            else:
                fout.write(line)                
