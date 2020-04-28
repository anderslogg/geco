'''
A collection of useful post processing functions
'''
from geco import *
import csv
import numpy as np
from itertools import izip
import sys, os, csv
from dolfin import *


def GetTitle(parameters, variables=["model", "weight"], sep="\n"):
    param_dict=GetParametersDicts(parameters)
    title=sep.join(map(str, [param_dict[x] for x in variables]))
    return title

def GetParametersStrings(parameters, sep="\n"):
    param_dict=GetParametersDicts(parameters)
    del param_dict['model']
    param_str=""
    for key, value in param_dict.iteritems():
        temp = key + ': ' + value +'\n'
        param_str+= temp
    param_str = os.getcwd().split('/')[-2] + sep +param_str
    return param_str

def GetParametersDicts(parameters):
    param_dict=None
    with open(parameters) as p:
        param_dict = dict(filter(None, csv.reader(p)))
    p.close()
    return param_dict

def GatherFiles(dir):
    ''' Returns often used postprocessing objects. All "RHO...xml.gz" are
    evaluated and returned in the list "components". The list is sorted first
    by length of filename then alphabetically so for example:
    components[0] = RHO_100.xml.gz
    components[1] = RHO_comp_0_100.xml.gz
    components[2] = RHO_comp_1_100.xml.gz

    parameters is returned as a list of the "parameters_#.csv" files. Also sorted:
    parameters[0] = ../parameters_0.csv
    parameters[1] = ../parameters_1.csv

    Note: the first model in component does not have a corresponding parameter.csv
    file so the component-parameter indices are off by one.
    '''
    try:
        mesh = Mesh(dir + '/mesh.xml.gz')
    except:
        print('mesh.xml.gz not found in: ' + dir)

    V=FunctionSpace(mesh,'P',1)

    try:
        potential = [f for f in os.listdir(dir) if (f.startswith('U_')
         and not (f.startswith('U_R')) and (f.endswith('.xml.gz')))]
        _U = potential[0].split('.')[0].split('_')[1]
        U = Function(V)
        U.set_allow_extrapolation(True)
        File(potential[0]) >> U
    except:
        print('U_{:}.xml.gz file not found.'.format(_U))

    try:
        parameters = [f for f in os.listdir(dir) if (f.startswith('param'))]
    except:
        print('component CSV files not found in found in: ' + dir)
    parameters.sort()

    try:
        c = [f for f in os.listdir(dir) if (f.startswith('RHO_') and (f.endswith('.xml.gz')))]
    except:
        print('RHO component files not found in found in: ' + dir)
    c.sort(key=len)
    c.sort()

    components = []
    for i in range(len(c)):
        rho = Function(V)
        rho.set_allow_extrapolation(True)
        try:
            File(c[i]) >> rho
            components.append(rho)
        except:
            print(c[i] + '{} file not found.'.format(i))

    return mesh, parameters, components, U, V


def ToNumpyArray(comp_density, r_max=20, res=250):
    rvals = np.linspace(0,r_max,res)
    zvals = np.linspace(0,r_max,res)
    rho_array = np.zeros((len(rvals), len(zvals)))
    for i in range(len(zvals)):
        for j in range(len(rvals)):
            r = rvals[j]
            z = zvals[i]
            rho_array[j,i] = comp_density(z,r)

    return rho_array

def TangentialVelocityArray(model, V, rho, tol=1e-3, r_max=20, res=250):
     vel_integral = project(model, V)
     vel_integral.set_allow_extrapolation(True)
     rvals = np.linspace(0,r_max,res)
     zvals = np.linspace(0,r_max,res)
     vel_array = np.zeros((len(rvals), len(zvals)))
     for i in range(len(zvals)):
         for j in range(len(rvals)):
            r = rvals[j]
            z = zvals[i]
            if(rho[j,i] > tol):
                vel_array[j,i] = vel_integral(z,r)/rho[j,i]
            else:
                vel_array[j,i] = 0

     return vel_array

def TangentialVelocityCurve(vel_array, r_max, z=0):
    r_length = np.shape(vel_array)[1]
    vel=np.zeros(r_length)
    inv_r = np.zeros(r_length)
    r_vals=np.linspace(0,r_max,r_length)
    for i in range(1,r_length):
        vel[i]=vel_array[i,z]
        inv_r[i-1] = 1/np.sqrt(r_vals[i])

    return vel, inv_r, r_vals

def RotationCurve(U, r_max, res=500, z=0):

    #v = sqrt(rU'(r, z))
    #z value is constant
    dr = 10 ** -10
    rvals = np.linspace(0,r_max,res)
    v = np.zeros(len(rvals))
    inv_r = np.zeros(len(rvals))

    for i in range(1,len(rvals)):
        dUdr = (U(rvals[i] + dr,z) - U(rvals[i] - dr,z))/(2*dr)
        v[i] = np.sqrt(rvals[i]*dUdr)
        inv_r[i] = 1/np.sqrt(rvals[i])

    return rvals, inv_r, v



def FilesDirsByName(filename):
    '''returns sorted list of files in all subdirectories
    matching "filename" as the first list and a list of the
    subdirectories themselves as the second list'''
    file_list = []
    dir_list = []
    for root, dirs, files in os.walk(os.getcwd(), topdown=False):
       for name in files:
          if filename in name:
              file_list.append(os.path.join(root, name))
              dir_list.append(root)

    dir_list.sort()
    file_list.sort()
    return file_list, dir_list

def MaxSupport(data_list = None):
    '''Returns the largest radius of support found
     after searching all subdirectory data.csv files'''
    if data_list == None:
        data_list, dir_list = FilesDirsByName("data.csv")

    cur_support = 1
    for d in data_list:
        try:
            header = np.genfromtxt(d, max_rows=1, delimiter=',', dtype=str)
            data_lines = np.genfromtxt(d, delimiter=',', dtype=None)
            data_lines = data_lines[-1]
            data = dict( zip(header, data_lines) )
        except:
            print('Unable to find ' + d)

        if float(data['radius_of_support']) > cur_support:
            cur_support = np.ceil(float(data['radius_of_support']))

    return cur_support




###########################################################
#FROM TAN VELOCITY

def BuildTanVelModel(param_dict, U, num_steps=150):
    '''This reconstructs a model by evaluating passed parameter file
    It calls '''

#in case a string to the parameter.csv file is passed instead of a dictionary
    if type(param_dict) is not dict:
        param_dict=GetParametersDicts(param_dict)
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

    curModel.set_fields(U)
    curModel.set_integration_parameters(num_steps)
    curModel.read_parameters()
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
    cppcode = template_code % {"member_functions": member_functions,
                               "member_variables": member_variables}
    rho = Expression(cppcode=cppcode, degree=1)
    rho.parameters.rename(model)

    return rho


#####################################################
''' From ROTATION CURVE '''
