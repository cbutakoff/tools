#see https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

nSample = 900
n = np.array([1,1,1])
c = np.array([15,15,15])
R = 10

n = n/np.linalg.norm(n)
v3 = np.random.random(3)
v3 = v3/np.linalg.norm(v3)
v2 = np.cross(n, v3)
v2 = v2/np.linalg.norm(v2)
v1 = np.cross(v2,n)

T = np.column_stack( (v2,n,v1) )

xyz = np.zeros((3,nSample))

theta = 0
phi = 0


#sphere
#theta = 2*np.pi*np.random.random(nSample);
#hemisphere
theta = np.pi*np.random.random(nSample);
phi = np.arccos(2*np.random.random(nSample)-1.0);
xyz[0,:] = np.cos(theta)*np.sin(phi);
xyz[1,:] = np.sin(theta)*np.sin(phi);
xyz[2,:] = np.cos(phi);

xyz = R*T.dot(xyz)+c[:,np.newaxis]

np.savetxt("points.csv", xyz.transpose(), delimiter=" ")


x = xyz[0,:]
y = xyz[1,:]
z = xyz[2,:]




fig = plt.figure(figsize=plt.figaspect(1)*1.5)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0

mid_x = (x.max()+x.min()) * 0.5
mid_y = (y.max()+y.min()) * 0.5
mid_z = (z.max()+z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()

