#!/usr/bin/env python

"""
GECo Data Plot

A collection of tools for plotting data from GECo and adaptivesolver.py
solutions or a sequence of such solutions

Last Modified:  May 12 2017

"""
import numpy as np
import sys, os, glob
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import rcParams


def normalized_redshift_vs_radius_ratio(sol_dirs):
    "Plots the normalized redshift vs the ri/ro, as in Ansorg et al."
    
    e0_index       = get_data_index('E0', sol_dirs[0])
    rinner_index   = get_data_index('r_inner', sol_dirs[0])
    router_index   = get_data_index('r_outer', sol_dirs[0])
    redshift_index = get_data_index('central_redshift', sol_dirs[0])    

    for sol_dir in sol_dirs:
        e0_data, ri_data, ro_data, rs_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, rinner_index,router_index, redshift_index), \
          unpack=True)

        r_ratio = ri_data/ro_data
        zc_norm = rs_data/(1. + rs_data)
        plt.plot(r_ratio, zc_norm,':ro')
        
    plt.xlabel('Inner radius/Outer radius') # $\rho_i / \rho_o$')
    plt.ylabel('Normalized central redshift, $Z_C/(1 + Z_C)$')

    plt.show()

# TODO: Get labeling by E0-values working.
"""
        e0_iter = np.array([e0_data])
        ri_iter = np.array([ri_data])
        ro_iter = np.array([ro_data])
        rs_iter = np.array([rs_data])
        r_ratio = [rid/rod for rid, rod in zip(ri_iter, ro_iter)]
        zc_norm = [rsd/(1. + rsd) for rsd in rs_iter]
        plt.plot(r_ratio, zc_norm,':ro')
        for label, x, y in zip(e0_iter, r_ratio, zc_norm):
            plt.annotate(str(label), xy=(x, y), xytext=(-20, 20),
                     textcoords='offset points', ha='right', va='bottom')
"""        


def get_data_index(data, sol_dir):
    "Returns the column number corresponding to data in the solution sol_dir"
    
    header = np.loadtxt(sol_dir + '/data.csv', dtype=str, delimiter=',')
    names = header[0,:].tolist()
    data_index = names.index(data)
    
    return data_index

def angular_momentum_vs_E0(sol_dirs):
    "Plots the angular momentum vs E0"
    
    e0_index  = get_data_index('E0', sol_dirs[0])
    J_index   = get_data_index('total_angular_momentum', sol_dirs[0])

    for sol_dir in sol_dirs:
        e0_data, j_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, J_index), \
          unpack=True)
          
        plt.plot(e0_data, j_data,':ro')
        
    plt.xlabel('$E_0$') # $\rho_i / \rho_o$')
    plt.ylabel('Total Angular Momentum, $J$')

    plt.show()

def M_over_Rcirc_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"
    
    e0_index      = get_data_index('E0', sol_dirs[0])
    Rcirc_index   = get_data_index('Rcirc', sol_dirs[0])
    m_index       = get_data_index('mass', sol_dirs[0])    

    for sol_dir in sol_dirs:
        e0_data, rc_data, m_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, Rcirc_index, m_index), \
          unpack=True)

        m_over_rc = m_data/rc_data  
        plt.plot(e0_data, m_over_rc,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('$M/R_{circ}$')

    plt.show()    
