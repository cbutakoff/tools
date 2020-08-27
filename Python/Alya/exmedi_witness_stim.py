import numpy as np
import re
import sys
from io import StringIO

problem_name = sys.argv[1]
added_distance = float(sys.argv[2])
nrandom_samples = int(sys.argv[3])


def read_alyampio_array(filename):
    with open(filename, 'rb') as f:
        header = read_header_mpio(f)
    
        tuples = np.reshape( np.fromfile(f, dtype=np.dtype(header['DataTypePython']) ), (header['Lines'], header['Columns']) )
        
        return {'tuples':tuples, 'header':header};


def read_header_mpio(f):
    magic = np.fromfile(f,count=1, dtype=np.int64)[0]
    if magic != 27093:
        print(f'File {filename} does not appear to be alya mpio file')
        
    format = str(f.read(8))
    if not ('MPIAL' in format):
        assert False,f'File {filename} does not appear to be alya mpio file'

    version = str(f.read(8))
    obj = str(f.read(8))
    dimension = str(f.read(8))
    association = str(f.read(8))
    datatype = str(f.read(8))
    datatypelen = str(f.read(8))    
    seq_par = str(f.read(8))
    filt = str(f.read(8))    
    sorting = str(f.read(8))    
    idd = str(f.read(8))    

    if not ('NOID' in idd):
        assert False, f'ID column in {filename} is not supported'

    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,   f'Lost alignment reding {filename}'
    

    columns = np.fromfile(f,count=1,dtype=np.int64)[0]
    lines = np.fromfile(f,count=1,dtype=np.int64)[0]
    timestep_no = np.fromfile(f,count=1,dtype=np.int64)[0]
    nsubdomains = np.fromfile(f,count=1,dtype=np.int64)[0]
    mesh_div = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag1 = np.fromfile(f,count=1,dtype=np.int64)[0]
    tag2 = np.fromfile(f,count=1,dtype=np.int64)[0]
    time = np.fromfile(f,count=1,dtype=np.float64)[0]
    
    junk = str(f.read(8))    
    if not ('0000000' in junk):
        assert False,f'Lost alignment reding {filename}'

    junk = str(f.read(8))    #1
    junk = str(f.read(8))    #2
    junk = str(f.read(8))    #3
    junk = str(f.read(8))    #4
    junk = str(f.read(8))    #5
    junk = str(f.read(8))    #6
    junk = str(f.read(8))    #7
    junk = str(f.read(8))    #8
    junk = str(f.read(8))    #9
    junk = str(f.read(8))    #10
    
    if 'INT' in datatype:
        dt = 'int'
    elif 'REAL' in datatype:
        dt = 'float'
    else:
        assert False,f'Unsupported data type {datatype}'

    if '8' in datatypelen:
        dt = dt+'64'
    elif '4' in datatypelen:
        dt = dt+'32'
    else:
        assert False,f'Unsupported data type length {datatypelen}'


    header = {
        'Version':version, 
        'Object':obj,
        'Dimension':dimension,
        'Columns':columns,
        'Lines':lines,
        'Association':association,
        'DataType':datatype,
        'DataTypeLength':datatypelen,
        'TimeStepNo':timestep_no, 
        'Time':time, 
        'NSubdomains':nsubdomains,              
        'Div':mesh_div,
        'DataTypePython':dt}


    assert ('NOFIL' in filt), "Filtered fields are not supported"


    return header


def ExtractPotentialTable(problem_name):
    #process the stimulation points
    with open(f'{problem_name}.exm.dat','r') as f:
        exm_data = f.read().upper()

    m = exm_data.partition("STARTING_POTENTIAL")[2].partition("END_STARTING_POTENTIAL")[0]
    if len(m)==0:
        raise Exception(f"Starting potential table not found in {problem_name}.exm.dat")

    lines = []
    start_line_extraction = False
    for line in m.split('\n'):
        if len(line.strip())>0:
            if line.strip()[0]!="$":
                if start_line_extraction:
                    lines += [line.split('$')[0].strip()]

                if 'NSTIM' in line:
                    start_line_extraction = True

    c = StringIO("\n".join(lines))

    matrix = np.loadtxt(c)    

    #extract unique activation points
    pts = []
    r = []
    for i in range(matrix.shape[0]):
        pt = matrix[i,1:4]
        if i>1:
            pts_np = np.array(pts)
            d = np.linalg.norm(pts_np - pt, axis=1)
            if not ( d<1e-14 ).any():
                pts += [pt.tolist()]
                r += [ matrix[i,4] ]
        else:
            pts += [pt.tolist()]
            r += [ matrix[i,4] ]
        

    return pts, r


def GetRandomPointsAtR( points, center, r1, r2, npts ):
    #gets point at distance >r1 and <r2
    d = np.linalg.norm( points-center, axis=1)
    pts_in = points[ (d<r2) & (d>r1) ]
    
    if pts_in.shape[0] == 0:
        raise Exception(f"No points found around {center} between radii {r1} and {r2}")

    idx = list(range(pts_in.shape[0]))
    if pts_in.shape[0] > npts:
        idx = np.random.randint(0, high=pts_in.shape[0]-1, size=npts )

    return pts_in[idx,:]


def GetScale( line, keyword ):
    xscale_part = line.partition(keyword)[2].split()
    xscale = None
    for v in xscale_part:
        vv = v.strip()
        try:
            xscale = float(vv)
        except:
            pass

        if not xscale is None:
            break

    return xscale




def GetScaleFactors( problem_name ):
    scale_line = ""
    with open( f"{problem_name}.dom.dat" ,'r') as ff:
        for line in ff:           
            if "XSCALE" in line.partition('$')[0]:
                scale_line = line.partition('$')[0].upper()
                break

    scales = None
    if len(scale_line) > 0:
        scales = [0,0,0]
        scales[0] = GetScale( scale_line, 'XSCALE' )
        scales[1] = GetScale( scale_line, 'YSCALE' )
        scales[2] = GetScale( scale_line, 'ZSCALE' )

        if ( (scales[0] is None) or (scales[1] is None ) or (scales[2] is None) ):
            raise Exception('Failed to find scales')

    return scales


#================================================
#
#  Main program
#
stim_pts, stim_r = ExtractPotentialTable(problem_name)
print ( "Unique stimulation locations (x,y,z,r):" )
print ( np.hstack([np.array(stim_pts), np.array(stim_r)[:,np.newaxis] ]) ) 

print ("Read mesh points")
mesh_point_data = read_alyampio_array(f"{problem_name}-COORD.mpio.bin")
coords = mesh_point_data['tuples']

print ("Check if there is scaling factors")
scales = GetScaleFactors( problem_name )
if not scales is None:
    print(f'Scaling factors found: {scales}')
    coords = coords*np.array(scales)

ppts = []
ref_pts = []
for i in range(len(stim_pts)):
    pts = GetRandomPointsAtR( coords, np.array(stim_pts[i]), stim_r[i]+added_distance-added_distance/2, stim_r[i]+added_distance+added_distance/2, nrandom_samples )
    ppts += pts.tolist()
    ref_pts += [stim_pts[i]]

witness_lines = [f"WITNESS_POINTS, NUMBER={len(ppts)}"]
for i in range(len(ppts)):
    pt = ppts[i]
    ref_pt = ref_pts[i]

    witness_lines += [f"   {pt[0]} {pt[1]} {pt[2]} $$$ for stim point {ref_pt[0]} {ref_pt[1]} {ref_pt[2]}"]

witness_lines += ["END_WITNESS_POINTS"]

print ( "\n".join(witness_lines) )



