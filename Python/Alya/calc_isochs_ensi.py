import numpy as np
import os
import progressbar
import sys

iterationid_number_of_digits = 6  #how many digits to use for the iteration id in ensi variable files


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

    return header, data

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
ref_t    = float(sys.argv[3])
ref_t1   = float(sys.argv[4])

numbers, times = extract_intra_data(path, casename) 

idx     = (times >= ref_t) & (times <= ref_t1)
numbers = numbers[ idx ] 
times   = times  [ idx ] 

print(times)

header, data1 = read_scalar_variable_pernode( path, casename, "INTRA", numbers[0] )

zero_cross = -1*np.ones( data1.shape )
for f1, f2, t1, t2 in zip( numbers[:-1], numbers[1:], times[:-1], times[1:] ):
    print(t1)
    _, data1 = read_scalar_variable_pernode( path, casename, "INTRA", f1 )
    _, data2 = read_scalar_variable_pernode( path, casename, "INTRA", f2 )

    t = t1 - (t2-t1) * data1 / (data2 - data1) #0 crossing of INTRA
    valid_t = (t>=t1) & (t<=t2) & ( data2 > data1  ) & (zero_cross<0)

    if any(valid_t): 
        zero_cross[ valid_t ] = t[ valid_t ]

    if all(zero_cross>=0):
        break; #all isochrones set

write_scalar_variable_pernode(zero_cross, path, casename, "ISOCH", 0, header)

