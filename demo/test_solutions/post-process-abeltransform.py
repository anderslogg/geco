
'''
USE:
    
    pip install --user PyAbel

Perform forward abel transform on RHO data
Called in vpsolver.py
outputs image of the RHO_array and forward transformation
'''
import abel
import numpy as np
import matplotlib.pyplot as plt
import sys, os, csv
from dolfin import *
from timeit import default_timer as timer


# Load metric functions
cwd = os.getcwd()
print(cwd)
filename = [f for f in os.listdir(cwd) if f.startswith('RHO_')][0]
print(filename)
R = filename.split('.')[0].split('_')[1]

print(os.listdir(cwd))

# Read mesh and create function space
mesh = Mesh(cwd + '/mesh.xml.gz')
V = FunctionSpace(mesh, 'P', 1)

# Read and save RHO field

RHO  = Function(V)
try:
    File('RHO_{:}.xml.gz'.format(R)) >> RHO
except:
    print('RHO_{:}.xml.gz file not found.'.format(R) )


def forward_abel_transform(RHO):	
    #r_max,z_max - dimensions of quarter-image
    r_max = 2
    #resolution of images
    res = 500
    rvals = np.linspace(0,r_max,res)
    
    z_max = 2
    zvals = np.linspace(0,z_max,res)
	
    #build 
    RHO_array_A = np.zeros((len(rvals), len(zvals)))

    for i in range(len(rvals)):
        for j in range(len(zvals)):
            r = rvals[i]
            z = zvals[j]
            RHO_array_A[j,i] = RHO(r,z)
	
	
	#These two lines create a reflection of the array across horiz. axis
	RHO_array_B = np.flipud(RHO_array_A)
    RHO_array_C = np.concatenate((RHO_array_B,RHO_array_A), axis=0)
	
	#this line mirrors over vert. axis
    RHO_array = np.concatenate((np.fliplr(RHO_array_C),RHO_array_C), axis=1)
	
	
	#Using 'hansenlaw' is much faster than 'direct' without cython implementation
    forward_abel = abel.Transform(RHO_array, direction='forward', method='hansenlaw').transform
    
	#Output saved in "demo/abel_out" directory
    #Constant multiple applied to second paramater alters contrast
    
    #fig, axs = plt.subplots(1, 2, figsize=(6, 4))
    #axs[0].imshow(RHO_array, clim=(0, np.max(RHO_array)*1.15), origin='lower')
    #axs[1].imshow(forward_abel, clim=(0, np.max(forward_abel)*1.15), origin='lower')
    #axs[0].set_title('Original')
    #axs[1].set_title('Forward Transform')
    #plt.tight_layout()
    #plt.savefig("abel_out/out.png")
    #plt.close()
    
    plt.title("Forward Transform, resolution = %d" % res)
    plt.imshow(forward_abel)
    plt.savefig("abel_out/out.png")
    plt.close()
    

start = timer()
forward_abel_transform(RHO)
end = timer()
print(end - start)

# res   |   time (s)
# 500       3.14830684662
# 1000      11.8434147835
# 1500      27.00735116
# 3000      173.100042105
# 5000      Killed