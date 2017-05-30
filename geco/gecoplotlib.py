#!/usr/bin/env python

"""
GECo Data Plot

A collection of tools for plotting data from GECo and adaptivesolver.py
solutions or a sequence of such solutions. 

TODO: This would be great as a GUI.

Last Modified:  May 25 2017

"""
import numpy as np
import sys, os, glob
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import rcParams


def normalized_redshift_vs_radius_ratio(sol_dirs, labels):
    "Plots the normalized redshift vs the ri/ro, as in Ansorg et al."

    for sol_dir in sol_dirs:

        try:
            e0_index       = get_data_index('E0', sol_dir)
            rinner_index   = get_data_index('r_inner', sol_dir)
            router_index   = get_data_index('r_outer', sol_dir)
            redshift_index = get_data_index('central_redshift', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue            
            
        e0_data, ri_data, ro_data, rs_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, rinner_index,router_index, redshift_index), \
          unpack=True)

        e0_iter = np.array([e0_data]).reshape(-1)
        ri_iter = np.array([ri_data]).reshape(-1)
        ro_iter = np.array([ro_data]).reshape(-1)
        rs_iter = np.array([rs_data]).reshape(-1)
        r_ratio = [rid/rod for rid, rod in zip(ri_iter, ro_iter)]
        zc_norm = [rsd/(1. + rsd) for rsd in rs_iter]
        plt.plot(r_ratio, zc_norm,':ro')
        if labels:
            for label, x, y in zip(e0_iter, r_ratio, zc_norm):
                plt.annotate(str(label), xy=(x, y), xytext=(-2, 2),
                            textcoords='offset points', ha='left', va='bottom')        
        
    plt.xlabel('Inner radius/Outer radius') # $\rho_i / \rho_o$')
    plt.ylabel('Normalized central redshift, $Z_C/(1 + Z_C)$')

    plt.show()
       
def get_data_index(data, sol_dir):
    "Returns the column number corresponding to data in the solution sol_dir"
    
    header = np.loadtxt(sol_dir + '/data.csv', dtype=str, delimiter=',')
    names = header[0,:].tolist()
    try:
        data_index = names.index(data)
        return data_index
    except ValueError:
        raise ValueError

def angular_momentum_vs_E0(sol_dirs):
    "Plots the angular momentum vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index  = get_data_index('E0', sol_dir)
            J_index   = get_data_index('total_angular_momentum', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue
        
        e0_data, j_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, J_index), \
          unpack=True)
          
        plt.plot(e0_data, j_data,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('Total Angular Momentum, $J$')

    plt.show()

def frac_binding_energy_vs_E0(sol_dirs):
    "Plots the fractional binding energy vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index('E0', sol_dir)
            Eb_index = get_data_index('frac_binding_energy', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue
        
        e0_data, eb_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, Eb_index), \
          unpack=True)
          
        plt.plot(e0_data, eb_data,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('Fractional Binding Energy, $E_b$')

    plt.show()

def normalized_central_redshift_vs_E0(sol_dirs):
    "Plots the normalized central redshift vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index('E0', sol_dir)
            Zc_index = get_data_index('central_redshift', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue             
    
        e0_data, zc_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, Zc_index), \
          unpack=True)

        zc_norm = zc_data/(1. + zc_data)
          
        plt.plot(e0_data, zc_norm,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('Normalized Central Redshift, $Z_c/(1 + Z_c)$')

    plt.show()    

def normalized_central_redshift_vs_frac_binding_energy(sol_dirs, labels):
    "Plots the normalized redshift vs fractional binding energy"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index('E0', sol_dir)            
            Zc_index = get_data_index('central_redshift', sol_dir)    
            Eb_index = get_data_index('frac_binding_energy', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue             
            
        e0_data, zc_data, eb_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, Zc_index, Eb_index), \
          unpack=True)

        e0_iter = np.array([e0_data]).reshape(-1)
        eb_iter = np.array([eb_data]).reshape(-1)
        zc_iter = np.array([zc_data]).reshape(-1)           
        zc_norm = [rsd/(1. + rsd) for rsd in zc_iter]
          
        plt.plot(eb_iter, zc_norm,':ro')

        if labels:
            for label, x, y in zip(e0_iter, eb_iter, zc_norm):
                plt.annotate(str(label), xy=(x, y), xytext=(-2, 2),
                            textcoords='offset points', ha='left', va='bottom')         
        
    plt.xlabel('Fractional Binding Energy $E_b$')
    plt.ylabel('Normalized Central Redshift, $Z_c/(1 + Z_c)$')

    plt.show()

def normalized_central_redshift_vs_M_squared_over_J(sol_dirs):
    "Plots the normalized redshift vs M^2/J"

    for sol_dir in sol_dirs:

        try:
            Zc_index = get_data_index('central_redshift', sol_dir)    
            Jt_index = get_data_index('total_angular_momentum', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue            
        
        zc_data, jt_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(Zc_index, Jt_index), \
          unpack=True)

        zc_norm = zc_data/(1. + zc_data)
        m_over_j = 1/(jt_data)
          
        plt.plot(m_over_j, zc_norm,':ro')
        
    plt.xlabel('Mass square over Total Angular Momentum  $M^2/J$')
    plt.ylabel('Normalized Central Redshift, $Z_c/(1 + Z_c)$')

    plt.show()      

def M_over_Rcirc_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"     

    for sol_dir in sol_dirs:

        try:
            e0_index      = get_data_index('E0', sol_dir)
            Rcirc_index   = get_data_index('Rcirc', sol_dir)
            m_index       = get_data_index('mass', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue
        
        e0_data, rc_data, m_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, Rcirc_index, m_index), \
          unpack=True)

        m_over_rc = m_data/rc_data  
        plt.plot(e0_data, m_over_rc,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('$M/R_{circ}$')

    plt.show()


def r_peak_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"     

    for sol_dir in sol_dirs:

        try:
            e0_index      = get_data_index('E0', sol_dir)
            rpeak_index   = get_data_index('r_peak', sol_dir)
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue
        
        e0_data, rp_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, rpeak_index), \
          unpack=True)
 
        plt.plot(e0_data, rp_data,':ro')
        
    plt.xlabel('$E_0$')
    plt.ylabel('Radius peak density')

    plt.show()


def r_values_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"     

    for sol_dir in sol_dirs:

        try:
            e0_index      = get_data_index('E0', sol_dir)
            rp_index   = get_data_index('r_peak', sol_dir)
            ri_index   = get_data_index('r_inner', sol_dir)
            ro_index   = get_data_index('r_outer', sol_dir)            
        except ValueError:
            print "Desired data not found in '%s'" %(sol_dir)
            continue
        
        e0_data, rp_data, ri_data, ro_data = \
          np.loadtxt(sol_dir + '/data.csv', delimiter=',', skiprows=1, \
          usecols=(e0_index, rp_index, ri_index, ro_index), \
          unpack=True)
 
        plt.plot(e0_data, rp_data,':ro')
        plt.plot(e0_data, ri_data,':bo')
        plt.plot(e0_data, ro_data,':go')                
        
    plt.xlabel('$E_0$')
    plt.ylabel('Peak, Outer, and Inner radius of matter support')
    plt.legend(['peak rad.', 'inner rad.', 'outer rad.'], loc='upper left')

    plt.show()    


