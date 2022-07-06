#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os.path import join
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
# setting-up the paths
root_dir = os.getcwd()
# insert path where the local libraries are located
sys.path.insert(0,'/net/home/h01/nvalient/nvalient-python/')
# import local library
import pyww3.plot.ver_plot_cmems as vpc
import pyww3.obs_funs.obs_reader as rdobs

"""
Created on Mon Nov 22 12:30:24 2021
EXAMPLE SCRIPT TO EVALUATE OPERATIONAL MODEL PERFORMANCE:
     RETURNS SCATTER AND QQ PLOTS WITH METRICS FOR A LOCATION AND/OR AREA 
Wave model accuracy across the English Channel. From March 2021 to 22nd November 2021
- INPUT: operational matchup files.  
- OUTPUT: .png scatter plot with confidence intervals, R, bias and RMSD.

@author: nvalient
"""


DIR_IN = '/project/ofrd/waves/verification/matchup/uk/'
# MM_ini = int(202112)
# MM_end = int(202203)

MM_ini = int(202012)
MM_end = int(202103)

# LOC_ID = ['6201053','6201025','6201005','6201006','6201011']
# LOC_NAME = ['Scilly_Island_2','Looe_Bay','Chesil_Beach_1','Chesil_Beach_2','Hayling_island']
# LOC_TITLE = ['CCO Isles of Scilly','CCO Looe Bay','CCO West Bay','CCO Chesil',
#              'CCO Hayling Island']

LOC_ID = ['6201054','6201053','6201025','6201005','6201006','6201011']
LOC_NAME = ['Scilly_Island_1','Scilly_Island_2','Looe_Bay','Chesil_Beach_1','Chesil_Beach_2','Hayling_island']
LOC_TITLE = ['Cefas SW Isles of Scilly','CCO Isles of Scilly','CCO Looe Bay','CCO West Bay','CCO Chesil',
              'CCO Hayling Island']
outdir = os.path.join(root_dir,'img/')

# check if dirout exists
if not os.path.exists(os.path.join(root_dir,'img/')):
    os.makedirs(os.path.join(root_dir,'img/'))


# -------------------------------------------------------------------------
    # START CODE TO GET THE DATA AND PLOT

for ii,ids in enumerate(LOC_ID):
    print('Read files for ID: '+ids)
    var_obs = []
    var_mod = []
    MM_ini = int(202012) #int(202012)
    while MM_ini <= MM_end:
        # List all the files in the directory
        dir_files = listdir(join(DIR_IN,str(MM_ini)))
        dir_files.sort()
        for file in dir_files:
            print('[INFO] Reading data from '+file)
    
            latsarr, lonsarr, vtarr, myvarobs, myvarmod = rdobs.read_collocation_files(join(DIR_IN,str(MM_ini),file),ids)
            var_obs = np.append(var_obs,myvarobs)
            var_mod = np.append(var_mod,myvarmod)
            
        if '12' in str(MM_ini):
            MM_ini = int(str(int(str(MM_ini)[0:4])+1)+'01')
        else:
            MM_ini = MM_ini+1
            
    # eliminate bad obs where differential to model is too large
    if not(var_obs is None):
        var_obs = np.array(var_obs)
        var_mod = np.array(var_mod)
        modstd = np.nanstd(var_mod)
        goodobs = np.where(np.abs(var_obs-var_mod)<5.0*modstd)
        rejectcount = len(var_obs)-len(goodobs[0])
        print('[INFO] Rejected %d' %rejectcount + ' bad observations')
        obsforver = var_obs[goodobs]
        modforver = var_mod[goodobs]
        #print('var_obs are '+str(var_obs))
    
    xdata = modforver 
    ydata = obsforver
    
    # TEST PLOT
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatter( xdata, ydata, axisres=1.0, hexbin=None, linfit=True, grid=False, showplt=False ) 
    #sns.regplot(xdata, ydata)
    ax.set_ylabel('$H_s$ observations [m]')
    ax.set_xlabel('$H_s$ model [m]')
    plt.title('UK model verification - '+LOC_TITLE[ii])
    rect = vpc.TextStats(xdata, ydata, units=None, linfit=True, errorstats=True, basicstats=True, forplot=True, ptloc=[0.1,max(ydata)], font='x-small',
                          tolerance=True, dirn=False)
    
    ax.add_patch(rect)
                 
    # save the plot
    out_name = os.path.join(outdir,LOC_NAME[ii]+'_scatter_W2020_2021.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
             
    
    # TEST PLOT
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatQQ( xdata, ydata, axisres=1.0, hexbin=True, linfit=True, grid=True, showplt=False ) 
    ax.set_ylabel('$H_s$ observations [m]')
    ax.set_xlabel('$H_s$ model [m]')
    plt.title('UK model verification - '+LOC_TITLE[ii])
    #vpc.TextStats(xdata, ydata, units=None, linfit=True, errorstats=True, forplot=True, dirn=False)
    # save the plot
    #out_name = os.path.join(outdir,'scatterQQ-all.png')
    out_name = os.path.join(outdir,LOC_NAME[ii]+'_QQ_W2020_2021.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
    
    print('-----------------------------------------------')
    print('Saved figures for ID: '+ids)

    
    
    
    
    
