#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import os

#import VerCalls as vc
#import VerLoad as vl

# setting-up the paths
root_dir = os.getcwd()
# insert path where the local libraries are located
sys.path.insert(0,'/net/home/h01/nvalient/nvalient-python/')
# import local library
import pyww3.plot.ver_plot_cmems as vpc
import pyww3.obs_funs.read_obs_insitu as rdobs

"""
Created on Thu Feb 11 11:41:30 2021

 SCRIPT EXAMPLE TO PLOT SCATTER QQ PLOT model vs obs

    Inputs: - xdata: model timeseries
            - ydata: observation timeseries

@author: nvalient
"""

outdir = os.path.join(root_dir,'img/')

# Obs data
IDs = ['62023'] 
obs_name = 'WFVS'
year=2013
month=11
day=1
fcday=0
cycle=datetime.datetime(year,month,day,0)
days_lead = 4 # months
varname = 'hs'
obsforver0 = rdobs.create_timeseries(IDs[0], obs_name, varname, cycle, days_lead)
obsforver  = obsforver0.data[:]

# Load model data - e.g., wiht timeseries previously extracted
modFile = '/data/users/nvalient/T2013_14/TIMESERIES/UKC4aow-ST4_timeseries_WFVS_pnts.nc'
# e.g., for Hs
d = nc.Dataset(modFile)
iID = d.variables['point'][:]=="62023"
modforver = d.variables['hs'][:,iID]
modforver = modforver[:,0]
t         = d.variables['time']
mod_date  = nc.num2date(t[:],t.units,calendar='gregorian',only_use_cftime_datetimes=False,only_use_python_datetimes=True)
# Load model data - e.g., using nc2pnt.py 

# Add matchup for time!
date_mask = np.in1d(mod_date[:],obsforver0.times[:])
modforver = modforver[date_mask]

# # ensure observations and model have the same data mask
# modforver = np.ma.masked_where(obsforver.mask,modforver)
# obsforver = np.ma.masked_where(modforver.mask,obsforver)
# obsforver = np.ma.compressed(obsforver)
# modforver = np.ma.compressed(modforver)

# eliminate bad obs where differential to model is too large
if not(obsforver is None):
    modstd = np.std(modforver)
    goodobs = np.where(np.abs(obsforver-modforver)<5.0*modstd)
    rejectcount = len(obsforver)-len(goodobs[0])
    print('[INFO] Rejected %d' %rejectcount + ' bad observations')
    obsforver = obsforver[goodobs]
    modforver = modforver[goodobs]

xdata = modforver 
ydata = obsforver   

# Plot - Hs
fig = plt.figure()
ax  = fig.add_subplot(111)
pl = vpc.PlotScatQQ( xdata, ydata, axisres=1.0, qqplt=1, hexbin=True, linfit=True, grid=True, showplt=True )     
ax.set_ylabel('$H_s$ observations [m]')
ax.set_xlabel('$H_s$ model [m]')

# save the plot
out_name = os.path.join(outdir,'scatterQQ-test.png')
plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
         