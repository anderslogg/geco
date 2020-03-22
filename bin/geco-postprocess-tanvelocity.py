#!/usr/bin/env python
'''
geco-postprocess-tanvelocity

A script to calculate tangential velocity

Required files: mesh.xml.gz, U_##.xdmf, RHO_comp_#_##.xdmf,
                data.csv, parameters_#.csv

Usage:

with setup.py installed in a docker container run

geco-postprocess-tanvelocity [-t] [-z]

-t  Optional user defined tolerance cutoff
    Default is relative to max value of RHO Components

-z  Optional user defined z axis value about which
    to calculate the rotation curve
    Default z = 0

from the solution directory.

'''
import argparse
from geco import *
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import sys, os, csv
from dolfin import *
from itertools import izip



def vel_tangential(tol=0,zed=0):

    dir = os.getcwd()

    potential = [f for f in os.listdir(dir) if (f.startswith('U_') and not (f.startswith('U_R')) and (f.endswith('.xml.gz')))]
    _U = potential[0].split('.')[0].split('_')[1]

    # Read mesh and create function space
    mesh = Mesh(dir + '/mesh.xml.gz')
    V = FunctionSpace(mesh, 'P', 1)
    x = SpatialCoordinate(mesh)

    # Read and save U field
    U = Function(V)
    U.set_allow_extrapolation(True)

    try:
        File(potential[0]) >> U
    except:
        print('U_{:}.xml.gz file not found.'.format(_U))


    components = [f for f in os.listdir(dir) if (f.startswith('RHO_comp_') and (f.endswith('.xml.gz')))]
    comp_densities = []
    C = components[0].split('.')[0].split('_')[1]


    for i in range(len(components)):
        rad_sup = components[i].split('_')[3].split('.')[0]
        RHO_COMP = Function(V)
        RHO_COMP.set_allow_extrapolation(True)
        try:
            str = 'RHO_comp_{}_{}.xml.gz'.format(i, rad_sup)
            File(str.format(C)) >> RHO_COMP
            comp_densities.append(RHO_COMP)
        except:
            print(str + ' file not found.'.format(C))

    data_dict = dict()
    param_dict_list = []
    parameters = [f for f in os.listdir(dir) if (f.endswith('.csv'))]
    parameters.sort()

    for f in parameters:
        if(f.startswith("data")):
            with open(f) as f:
                data_in=filter(None,csv.reader(f))
                data_dict=dict(izip(data_in[0],data_in[1]))
            f.close()
        elif(f.startswith("parameters")):
            with open(f) as f:
                param_dict_list.append(dict(filter(None, csv.reader(f))))
            f.close()
        else:
            continue

    r_max = np.ceil(float(data_dict["radius_of_support"]))

    #resolution of images
    res = 500
    rvals = np.linspace(0,r_max,res)

    z_max = np.ceil(float(data_dict["radius_of_support"]))

    zvals = np.linspace(0,z_max,res)
    for param_dict,rho in izip(param_dict_list, comp_densities):

        model = None

        model,weight,name = build_model(param_dict)
        num_steps = 150
        model.set_fields(U)
        model.set_integration_parameters(num_steps)
        model.read_parameters()

        vel_integral = project(model, V)
        vel_integral.set_allow_extrapolation(True)

        #If the tolerance is not user defined, set it
        #based on the maximum value of the current RHO
        reset_tol=False
        if tol == 0:
            reset_tol = True
            nodal_values = rho.vector().get_local()
            tol = np.amax(nodal_values) * 0.00001


        #Perform the integration
        #creates both the tanvelocity 2D numpy array
        #AND outputs the original RHO as
        RHO_array = np.zeros((len(rvals), len(zvals)))
        VEL_array = np.zeros((len(rvals), len(zvals)))
        for i in range(len(zvals)):
            for j in range(len(rvals)):
                r = rvals[j]
                z = zvals[i]
                RHO_array[j,i] = rho(r,z)
                if(RHO_array[j,i] > tol):
                    VEL_array[j,i] = vel_integral(z,r)/rho(r,z)
                else:
                    VEL_array[j,i] = 0

        VEL_array = VEL_array*(float(weight))
        RHO_array = RHO_array*(float(weight))

        plt.imshow(RHO_array, cmap='gist_heat', extent=(0,r_max,0,z_max), origin='lower')
        plt.title("RHO - " + name)
        plt.colorbar(format='%.0e')
        plt.savefig("RHO " + name)
        plt.close()


        plt.imshow(VEL_array, cmap='gist_heat', extent=(0,r_max,0,z_max), origin='lower')
        plt.title("Velocity - " + name)
        plt.colorbar(format='%.0e')
        plt.savefig("Velocity " + name)
        plt.close()


        res = 500
        v = np.zeros(len(rvals))
        inv_r = np.zeros(len(rvals))
        rvals = np.linspace(0,r_max,res)

        for i in range(1,len(rvals)):
                v[i] = VEL_array[i,zed]
                inv_r[i-1] = 1/np.sqrt(rvals[i])

        alpha = inv_r[1]*v[1]

        plt.plot(rvals, v)
        plt.plot(rvals[20:len(rvals)], alpha*inv_r[20:len(rvals)])
        plt.title('Velocity vs. Radius ' + name + '\n tolerance = {}'.format(tol))
        plt.xlabel('radius')
        plt.ylabel('velocity')
        plt.savefig("Rotation Curve" + name + ".png")
        plt.close()

        if reset_tol==True:
            tol=0
            reset_tol=False

    return ""


def build_model(param_dict):
#Extract weight and model name parameters, then remove them
    weight = param_dict.get('weight')
    del param_dict['weight']
    name = param_dict.get('model')
    print(name)
    del param_dict['model']

    curModel = TangentialVelocityModel(name)
    params = "curModel.parameters."
    for key in param_dict:
        print(params + key + " = " + param_dict[key])
        exec(params + key + " = " + param_dict[key])

    return curModel, weight, name

def TangentialVelocityModel(model):
    "Create given material model"

    # Fixes subtle issue when model comes in as unicode and should be a string.
    # This pops up when running from the GUI and then calling Parameters.rename.
    model = str(model)

    # Get model data
    if not model in model_data:
        error("Unknown material model: \"%s\"." % str(model))
    model_filename = model + ".h"
    template_filename = "tan_velocity.h"

    library_dir = os.path.join(os.path.dirname(__file__), '..','geco/')
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code
    template_code = open(os.path.join(cppcode_dir, template_filename)).read()
    model_code = open(os.path.join(cppcode_dir, model_filename)).read()

    # Extract relevant model code
    member_functions = model_code.split("// Member functions")[1].split("// Member variables")[0]
    member_variables = model_code.split("// Member variables")[1]

    # Stick specialized code into template and return
    # FIXME: python3 version:
    # template_code.format(member_functions=member_functions, member_functions=member_functions),
    # and in XXAnsatz.h as {member_functions} and {member_variables}
    cppcode = template_code % {"member_functions": member_functions,
                               "member_variables": member_variables}
    # Create expression
    rho = Expression(cppcode=cppcode, degree=1)

    # Set name of ansatz
    rho.parameters.rename(model)

    return rho


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tolerance', help="enter a tolerance cutoff",
                        type=float, required=False)
parser.add_argument('-z', '--zaxis', help="enter a z-axis ",
                        type=int, required=False)
args=vars(parser.parse_args())

tol = 0 if args['tolerance'] is None else args['tolerance']
zed = 0 if args['zaxis'] is None else args['zaxis']

vel_tangential(tol,zed)
