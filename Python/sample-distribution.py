import numpy as np
import matplotlib.pyplot as plt 
from scipy import interpolate
from progressbar import progressbar
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('Radius', type=float,
                    help='Radius of the disc within which to sample')
parser.add_argument('Amplitude', type=float,
                    help='Amplitude of the parabolic density')
parser.add_argument('NSamples', type=int,
                    help='Number of samples to sample')
parser.add_argument('Output_Filename', type=str,
                    help='CSV file to save the particles.')
parser.add_argument('--grid_nsamples', '-gs', metavar='N', type=int, default=100,
                    help='Number of samples in along each axis to take in the grid discretizing the distribution, default 100.')
parser.add_argument('--normal','-n', metavar='N', nargs=3, type=float, default = [0,0,1],
                    help='Normal to the plane of the disc, default plane is XY.')
parser.add_argument('--center','-c', metavar='N', nargs=3, type=float, default = [0,0,0],
                    help='Center of the particles, default (0,0,0).')
parser.add_argument('--plot','-p', action='store_true',
                    help='Show plot of the distribution of the samples')
parser.add_argument('--extra_samples', '-e', metavar='N', type=int, default=10000,
                    help='Extra samples, to account for samples removed for being outside of the disc, default 10000.')

args = parser.parse_args()

print (args)




ExtraSamples = args.extra_samples #sample more, since some samples will be outside of the disc, to compensate for the removed particles
R = args.Radius
N = args.NSamples
M = args.Amplitude
output_filename = args.Output_Filename
grid_nsamples = args.grid_nsamples



########################################################################################
#
#
#  Functions
#
#
########################################################################################

def ParabolicProfile(X, Y, M, R):
    #returns
    Z = M * (1 - np.sqrt(X**2 + Y**2)/R )
    Z[Z<0]=0
    return  Z




def sample_conditional(X,F,N): 
    #following https://en.wikipedia.org/wiki/Inverse_transform_sampling
    
    #normalizes F to have integral =1, for continuous distributions
    area_under_F = np.trapz(F,X) #use trapezoidal rule for integration
    if area_under_F>=0.0000001:
        F=F/area_under_F
    cdf = F.cumsum() #cumulative distribution function
    samples_uniform = np.random.rand(N)*cdf.max() #uniform samples of the Y axis of the CDF
    
    #inverse CDF
    cdf_inverse = interpolate.interp1d(cdf, X, kind='linear', fill_value=0, bounds_error=False)
    
    #sample and return the values of the inverse CDF
    return cdf_inverse(samples_uniform)


def gibbs(X,Y,Z,N):
    xmin = X.min()
    xmax = X.max()
    ymin = Y.min()
    ymax = Y.max()
    
    f = interpolate.interp2d(X, Y, Z, kind='linear', fill_value=0, bounds_error=False)
    XY = np.hstack([X.ravel(),Y.ravel()])
    xy = np.zeros((N,2))
    
    #initial Y
    y = np.random.rand(1)*(ymax-ymin)+ymin
    
    for i in progressbar(range(N)):
        x = sample_conditional( X[0,:].ravel(), f(X[0,:], y).ravel(), 1 )
        y = sample_conditional( Y[:,0].ravel(), f(x,Y[:,0]).ravel(), 1 )
        xy[i,0] = x
        xy[i,1] = y
        
    return xy


########################################################################################
#
#
#  End of Functions
#
#
########################################################################################



#create a grid of X,Y points for plotting and interpolation
ttt = np.linspace(-R, R, num=grid_nsamples,  endpoint=True)
X,Y = np.meshgrid( ttt, ttt )

#define the function
Z = ParabolicProfile(X, Y, M, R)






#transfer points to the global coordinate system
C = np.array([args.center]).transpose()
n = np.array(args.normal)
n = n / np.linalg.norm(n,2)

#find one of the axes perpendicular to n by solving
# axis1.n = 0
#axis1[0]*n[0]+axis[1]*n[1]+axis[2]*n[2] = 0
if n[0]>0:
    axis1 = [ 0, n[0], n[1] ]
    axis1[0] = ( -axis1[1]*n[1] - axis1[2]*n[2] )/n[0]
elif n[1]>0:
    axis1 = [ n[1], 0, n[2] ]
    axis1[1] = ( -axis1[0]*n[0] - axis1[2]*n[2] ) / n[1]
elif n[2]>0:
    axis1 = [ n[1], n[2], 0 ]
    axis1[2] = (-axis1[0]*n[0] - axis1[1]*n[1])/n[2]
else:
    print('Normal cannot be (0,0,0)')
    exit(1)

axis1 = axis1 / np.linalg.norm(axis1, 2);
axis1 = axis1[:,np.newaxis]


axis2 = np.cross( n.ravel(), axis1.ravel() )
axis2 = axis2[:, np.newaxis]

if( sum( (axis1*axis2).ravel() )>1e-10 ):
    print('Axes generated incorrectly. Not orthogonal. Contact developer')
    exit(1)


#sample the particles
xy = gibbs(X, Y, Z, N+ExtraSamples)  #Number of samples to sample

if args.plot:
    plt.hist2d(xy[:,0],xy[:,1], bins=[50,50])
    plt.axis('equal')
    plt.colorbar()
    plt.show()


#remove the particles that are outside the disc
rr = np.sqrt(xy[:,0]**2+xy[:,1]**2)
xy_inside = xy[rr<=R,:]

if ( xy_inside.shape[0]<N ):
    print(f'The number of samples {xy_inside.shape[0]} inside the disc is smaller than requested {N}. Increase the extra_samples')


xy_inside = xy_inside[np.random.randint(xy_inside.shape[0], size=N),:]

coords = xy_inside[:,0]*axis1 + xy_inside[:,1]*axis2 + C



np.savetxt(output_filename, coords.transpose(), delimiter=",")

