#!/usr/bin/env python

'''
geco-postprocess-rotation_curve

A script to compute velocity of particles as a function of radius

Required file: U_R.xml.gz, mesh.xml.gz

Usage:

with setup.py installed in a docker container run

geco-postprocess-rotation_curve

from the solution directory.

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys, os, csv
from dolfin import *

cwd = os.getcwd()

dir_list = filter(os.path.isdir, os.listdir(os.getcwd()))
num_dirs = len(dir_list)

#This sets the figure size. It is static currently
# something like 20/3 * num_components
fig = plt.figure(figsize=(20,num_dirs*5))
gs = gridspec.GridSpec(num_dirs, 4)


def rotation_curve(cur_dir):
    subdir = os.getcwd()
    potential = [f for f in os.listdir(subdir) if (f.startswith('U_') and not (f.startswith('U_R')))]
    R = potential[0].split('.')[0].split('_')[1]

    # Read mesh and create function space
    mesh = Mesh(subdir + '/mesh.xml.gz')
    V = FunctionSpace(mesh, 'P', 1)
    x = SpatialCoordinate(mesh)

    # Read and save U field
    U = Function(V)
    U.set_allow_extrapolation(True)

    try:
        File('U_{:}.xml.gz'.format(R)) >> U
    except:
        print('U_{:}.xml.gz file not found.'.format(R))

    components = [f for f in os.listdir(subdir) if (f.startswith('RHO_comp_') and (f.endswith('.xml.gz')))]
    print(components)

    comp_densities = []

    for i in range(len(components)):
        RHO_COMP = Function(V)
        RHO_COMP.set_allow_extrapolation(True)
        C = components[i].split('.')[0].split('_')[1]
        try:
            #str = 'RHO_comp_%d_{}.xml.gz' %i
            str = 'RHO_comp_%d_50.xml.gz' %i
            File(str.format(C)) >> RHO_COMP
            comp_densities.append(RHO_COMP)
        except:
            print(str + ' file not found.'.format(C))

    #v = sqrt(rU'(r))
    #z value is constant
    #r_max,z_max - dimensions of quarter-image
    r_max = 35
    #resolution of images
    res = 500
    rvals = np.linspace(0,r_max,res)

    z_max = 35
    zvals = np.linspace(0,z_max,res)

    #build RHO components
    for c in range(len(comp_densities)):
        RHO_array = np.zeros((len(rvals), len(zvals)))
        for i in range(len(zvals)):
            for j in range(len(rvals)):
                r = rvals[j]
                z = zvals[i]
                RHO_array[len(zvals)-1-j,i] = comp_densities[c](z,r)

    #Add RHO components to figure
        f = fig.add_subplot(gs[cur_dir,c])
        f.imshow(RHO_array, cmap='gist_heat')
        f.set_title("RHO component %d" % c)

    #From here, calculate rotation curve
    z = 0
    dr = 10 ** -10
    r_max = 50
    res = 1000
    rvals = np.linspace(0,r_max,res)
    v = np.zeros(len(rvals))

    inv_r = np.zeros(len(rvals))

    for i in range(1, len(rvals)):
            dUdr = (U(z,rvals[i] + dr) - U(z,rvals[i] - dr))/(2*dr)
            v[i] = np.sqrt(rvals[i]*dUdr)
            inv_r[i] = 1/np.sqrt(rvals[i])

    f = fig.add_subplot(gs[cur_dir,2:4])
    f.plot(rvals, v)
    f.plot(rvals[150:len(rvals)], inv_r[150:len(rvals)])
    f.set_title('Velocity vs. Radius')
    f.set_xlabel('radius')
    f.set_ylabel('velocity')
###########################################################3
for cur_dir in range(num_dirs):
    os.chdir(cwd+"/"+dir_list[cur_dir])
    rotation_curve(cur_dir)

os.chdir(cwd)
fig.savefig("rotation_curve.png")
