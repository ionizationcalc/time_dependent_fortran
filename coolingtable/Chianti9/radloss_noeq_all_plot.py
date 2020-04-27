#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name
    radloss_noeq_all_plot
Purpose:
    Read saved files and plot total radative loss rates including all elements.
Author:
    Chengcai
Update:
    Created on 2020-04-13
    2020-04-14: Plot ChiantiPy results.
    2020-04-24.
"""
import ChiantiPy.core as ch
import numpy as np
from scipy import linalg

import os
import pickle
import glob
import time as systime
import matplotlib.pyplot as plt

# ==============================================================================
# function: func_readtable_list
# ==============================================================================
# Get all elements in the folder
def func_readtable_list(root_path):
    # a. FreeFree
    ff_path = root_path + 'freefree_rate/'
    files = sorted(glob.glob(ff_path+'*.p'))
    nfile = len(files)
    for ifile in range(nfile):
        print('a. reading ', files[ifile])
        ffrate_dict = pickle.load(open(files[ifile], "rb"))
        ni = ffrate_dict['abundance']*np.array(ffrate_dict['eqfraction'])
        if (ifile == 0):
            temperature_table = np.array(ffrate_dict['temperature'])
            nte_table = len(temperature_table)
            ffrate = np.zeros(nte_table)
            fbrate = np.zeros(nte_table)
            bbrate = np.zeros(nte_table)
            twophorate = np.zeros(nte_table)
        ffrate += np.array(ffrate_dict['rate'])*ni
        
    # b. FreeBound
    fb_path = root_path + 'freebound_rate/'
    files = sorted(glob.glob(fb_path+'*.p'))
    nfile = len(files)
    for ifile in range(nfile):
        print('b. reading ', files[ifile])
        fbrate_dict = pickle.load(open(files[ifile], "rb"))
        ni = fbrate_dict['abundance']*np.array(fbrate_dict['eqfraction'])
        fbrate += fbrate_dict['rate']*ni

    # c. BondBond
    bb_path = root_path + 'boundbound_rate/'
    files = sorted(glob.glob(bb_path+'*.p'))
    nfile = len(files)
    for ifile in range(nfile):
        print('c. reading ', files[ifile])
        bbrate_dict = pickle.load(open(files[ifile], "rb"))
        ni = bbrate_dict['abundance']*np.array(bbrate_dict['eqfraction'])  
        bbrate += bbrate_dict['rate']*ni

    # d. TwoPhoton
    twopho_path = root_path + 'twophoton_rate/'
    files = sorted(glob.glob(twopho_path+'*.p'))
    nfile = len(files)
    for ifile in range(nfile):
        print('d. reading ', files[ifile])
        twophorate_dict = pickle.load(open(files[ifile], "rb"))
        ni = twophorate_dict['abundance']*np.array(twophorate_dict['eqfraction'])   
        twophorate += twophorate_dict['rate']*ni

    # e. Sum = ff + fb + bb + twopho
    rate_table = ffrate + fbrate + bbrate + twophorate    
    return rate_table

# ==============================================================================
# Start Main
# ==============================================================================
#
# Define element list
#
minabund = 2.0e-5
if (minabund >= 1.0e-4):
    elemlist = ['h', 'he', 'c', 'o', 'ne']
    izlist = [1, 2, 6, 8, 10]
elif (minabund >= 2.0e-5):
    elemlist = ['h', 'he', 'c', 'n', 'o', 'ne', 'mg', 'si', 's', 'fe']
    izlist = [1, 2, 6, 7, 8, 10, 12, 14, 16, 26]
else:
    elemlist = (['h', 'he', 'c', 'n', 'o', 'ne', 'na','mg', 'al', 'si', 's', 'ar', 'ca', 'fe', 'ni'])
    izlist = [1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28]
print(minabund)
print(elemlist)
print(izlist)

# 
# Read tables
# 
root_path = '/Users/chengcai/Works/Project/ionization/ChiantiPy_subs/Radiativeloss/'

#
# Compare it with Radiative loss rate (in ChiantiPy Radloss)
#
eDensity = 1.0e9
rl_file = "chiantipy_totalloss.p"
if (os.path.exists(rl_file)):
    rl_dict = pickle.load(open(rl_file, "rb"))
else:
    temp = np.logspace(4.0, 8.0, 41)
    rl = ch.radLoss(temp, eDensity, elementList=elemlist, \
        abundance='sun_coronal_1992_feldman_ext', verbose=1)
    rate_ch = rl.RadLoss['rate']
    rl_dict = {"rate": rate_ch, "temperature":temp}
    pickle.dump(rl_dict, open(rl_file, "wb"))

#
# Read all save 2d tables
#
nelem = len(elemlist)
for ielem in range(nelem):
    elem = elemlist[ielem]
    file = "{0:s}_radloss.p".format(elem)
    print(file)
    dict = pickle.load(open(file, "rb"))
    if (ielem == 0):
        nte = len(dict["temperature"])
        rate_table_list = np.zeros((nelem, nte))
        rate_ratio = np.zeros((nelem, nte))
        rate_table_sum = np.zeros(nte)
    rate_c = np.sum(dict["rate"]*dict["eifraction"]*dict["abundance"], axis=0)
    rate_table_sum += rate_c
    rate_table_list[ielem, 0:nte] = rate_c[0:nte]

#  
# Get the relative contribution of elements
#
for ielem in range(nelem):
    rate_ratio[ielem, 0:nte] = rate_table_list[ielem, 0:nte]\
        /rate_table_sum[0:nte]

#
# Plot figures
#     
fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(1,2,1)

# show each element
for ielem in range(nelem):
    elem_str = elemlist[ielem].capitalize()
    ax1.loglog(dict["temperature"], rate_table_list[ielem, 0:nte], \
        label=elem_str)
ax1.loglog(rl_dict["temperature"], rl_dict["rate"], \
    'bo', label='Chianti')
ax1.loglog(dict["temperature"], rate_table_sum, \
    c='black', label='Sum')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Radiative Loss Rate (ergs cm$^{3}$ s$^{-1}$)')
ax1.set_ylim(1.0e-25, 1.0e-21)
ax1.legend()

ax2 = fig.add_subplot(1,2,2)
# show each element
for ielem in range(nelem):
    elem_str = elemlist[ielem].capitalize()
    ax2.plot(dict["temperature"], rate_ratio[ielem, 0:nte], \
        label=elem_str)
# show [he+o+mg+fe]
rate_ratio_o5 = rate_ratio[1, 0:nte] + rate_ratio[4, 0:nte]\
							+ rate_ratio[6, 0:nte] + rate_ratio[9, 0:nte]
#ax2.plot(dict["temperature"], rate_ratio_o5, \
#        label='He+O+Mg+Fe', c='black')

ax2.set_xscale('log')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Normalized Ratio')
ax2.set_ylim(1.0e-3, 1.0)
ax2.legend()
plt.savefig("fig_all_elements.png", dpi=300)
plt.close('all')
