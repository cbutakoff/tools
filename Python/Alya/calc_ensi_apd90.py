import numpy as np
import os
from progressbar import progressbar, ProgressBar
import sys
import matplotlib.pyplot as plt

iterationid_number_of_digits = 6  #how many digits to use for the iteration id in ensi variable files


def get_percentile(x, y, perc):
    """perc in [0,100]

        return position, x[position], height of y[position], peak hight, time of the peak hight 
    """
    v_min  = y.min()
    v_max  = y.max()
    v_perc = v_min + (v_max-v_min)*(1-perc/100)
    perc_pos = np.where( np.diff((y<v_perc)>0) )[0][-1]

    v_max_time = x[ np.argmax(y) ]

    return [perc_pos, x[perc_pos], y[perc_pos]-v_min, v_max-v_min, v_max_time]


def find_zerocross(y):
    return np.where( ((y[1:]-y[:-1])>0) & (y[:-1]<=0) & (y[1:]>0) )[0]  


def write_scalar_variable_pernode(data, path, project_name, varname, iteration, header_string):
    #variable ensight
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';

    filename = os.path.join(path, fmt % (project_name, varname, iteration))


    with open( filename,'wb') as f:
        f.write( header_string.ljust(80) )
        f.write(b'part'.ljust(80))
        f.write(np.array([1], dtype=np.int32))   #int
        f.write(b'coordinates'.ljust(80))

        f.write( data.ravel().astype(np.float32) )  


def read_scalar_variable_pernode(path, project_name, varname, iteration):
    fmt = '%s.ensi.%s-'+f'%0{iterationid_number_of_digits}d';

    filename = os.path.join(path, fmt % (project_name, varname, iteration))
    header = ""    
    data = None

    with open( filename,'rb') as f:
        header = f.read( 80 )
        _      = f.read( 80 ) #part
        _      = f.read( 4 )  #dtype=np.int32
        _      = f.read( 80 ) #coordinates
        data   = np.fromfile( f, dtype=np.float32 )

    return header, np.array(data).ravel()

def extract_intra_data( path, project_name ):
    variables_section = False
    timeline = None
    casefileaname = os.path.join( path, f"{project_name}.ensi.case" )

    with open(casefileaname) as ff:
        for line in ff:
            if "VARIABLE" in line.upper():
                variables_section = True
            if variables_section:
                if " INTRA " in line.upper():
                    data = np.array(line.upper().split())
                    idx = np.where( data=='INTRA' )[0]
                    timeline = int( data[idx-1] )
                    break

    #extract file numbers and times
    with open(casefileaname) as ff:
        lines = ff.readlines()

    numbers_line       = None
    times_line         = None
    time_set_line      = None
    next_time_set_line = len(lines)+1
    for i, line in enumerate(lines):
        if "time set" in line:
            if int(line.split(":")[1]) == timeline:
                time_set_line = i

            if int(line.split(":")[1]) == timeline+1:
                next_time_set_line = i
                break

        if not time_set_line is None:
            if "filename numbers" in line:
                numbers_line = i+1
            elif "time values" in line:
                times_line = i+1


    #extract file numbers and time values
    file_numbers = [int(number) for number in (" ".join( lines[numbers_line:times_line-1] )).split()]
    time_values  = [float(t)    for t      in (" ".join( lines[times_line  :next_time_set_line] )).split()]
    
    return np.array(sorted(file_numbers)), np.array(sorted(time_values))


path     = sys.argv[1]
casename = sys.argv[2]
#ref_t    = -1 #float(sys.argv[3])
#ref_t1   = 100 #float(sys.argv[4])
fileno   = int(sys.argv[3])

numbers, times = extract_intra_data(path, casename) 

#idx     = (times >= ref_t) & (times <= ref_t1)
#numbers = numbers[ idx ] 
#times   = np.array(times  [ idx ])  
times   = np.array(times)  


#read first 
header, data1 = read_scalar_variable_pernode( path, casename, "INTRA", numbers[0] )
data = np.zeros( [data1.shape[0], len(numbers)], dtype=np.single )
data[:,0] = data1.astype(np.single)

apd90 = -1*np.ones( data.shape[0] )

with ProgressBar(max_value=len(numbers)-1 ) as bar:
    for i, nn in enumerate(numbers[1:]):
        _, data1 = read_scalar_variable_pernode( path, casename, "INTRA", nn )
        data[:,i+1] = data1.astype(np.single)
        bar.update(i)

#calculate apds
for i in progressbar( range(data.shape[0]) ):
    zc = find_zerocross(data[i,:])
    
    #if len(zc)==0:
    #    print(zc)
    #    plt.plot(times,data[i,:])
        #plt.vlines([apd90_time, times[zc[-1]]], data[i,:].min(), data[i,:].max(), colors=['r','g'])
    #    plt.show()

    if len(zc)>0:
        _, apd90_time , _, _, _ = get_percentile(times[zc[-1]:], data[i,zc[-1]:], 90) 
        apd90[i]     = apd90_time - times[zc[-1]]   

    #plot = True
    #if plot:
    #    import matplotlib.pyplot as plt
    #    plt.plot(times,data[i,:])
    #    plt.vlines([apd90_time, times[zc[-1]]], data[i,:].min(), data[i,:].max(), colors=['r','g'])
    #    plt.show()


write_scalar_variable_pernode(apd90[:,np.newaxis], path, casename, "REPOL", fileno, header)

