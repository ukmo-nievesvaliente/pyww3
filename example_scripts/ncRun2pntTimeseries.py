#!/usr/bin/env python

#import subprocess
import matplotlib as mpl
mpl.use('Agg')
import sys
import os
import glob
import csv
from os.path import join, abspath
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import configparser as cfg
#from optparse import OptionParser
from datetime import datetime
#from scipy.interpolate import interp2d, RectBivariateSpline, griddata
import ver.collocation as col
import point_extract as pnt

def get_model_indices(model_ncfile, obs_lon, obs_lat):
    """ Return a set of model indices (either into a regular or
        SMC grid) for the nearest grid cell to the supplied
        observation lat/lons.
    """

    with nc.Dataset(model_ncfile, 'r') as d:

        # SMC or regular grid?
        is_smc = "seapoint" in d.dimensions

        if is_smc:
            # Uses a KdTree search algorithm:
            ind,dist = col.smc_nn_indices(d, obs_lon, obs_lat, max_dist=10000)
            ind = (ind,)
        else:
            # Uses a simple regular grid lookup algorithm:
            ind = col.reg_nn_indices(d, obs_lon, obs_lat)
            ind = ind[::-1]    # reverse order so is (iy, ix)

    return ind

# SCRIPT FOR EXTRACTING POINT TIMESERIES FROM MODEL NETCDF FILES USING Chris' routine in
#                 the python library "~frwave/FCM/python/"
# --------------------------------------------------------------------------

# 1. Define site locations and names
# Extracting Ids using matchup Dataset
obsDir_new = '/data/users/nvalient/trials2013_14/matchup/UKC4aow'
obstype    = 'WFVS'#'SHPSYN'# WFVS; WAVENET
Jfile      = join(obsDir_new,'match_'+obstype+'_UKC4aow_20131215.nc')
obs        = nc.Dataset(Jfile,"r")
robs_lat   = np.array(obs.variables["latitude_obs"])
robs_lon   = np.array(obs.variables["longitude_obs"])
robs_id    = np.array(obs.variables["station_id"])
rloc_N     = np.array(obs.variables['station_name'])
print('Extracting IDs')
obs.close()
# Map unique values in obs_id
mapping = {}
for i in np.unique(robs_id):
    mapping[i] = np.where(robs_id == i)[0]
# Extract indexes; get the first value in the dictionary
myList  = [mapping [i][0] for i in sorted(mapping.keys()) ]

obs_id  = np.unique(robs_id)
obs_lon = robs_lon[myList]
obs_lat = robs_lat[myList]

# 2. Define your model files and variables you want to extract:
RunFolder = '/data/users/nvalient/T2013_14/RUN_RESULTS'
# Runs  = ['UKC4ow']# folder run name
# Runs1 = ['UKC4owg']# prefix file MODEL RUNS
# Runs_name = ['UKC4owg']
Runs  = ['UKC4aow_st6JE']# folder run name
Runs1 = ['UKC4aow']# prefix file MODEL RUNS
Runs_name = ['UKC4aow-ST6mod']#,'UKW4h-ST4','UKW4h-ST6'] # name for output

out_dir = '/data/users/nvalient/T2013_14/TIMESERIES/'
varnames = ['hs','t01','fp','wspd','cha','utaw','vtaw','uuss','vuss','utwo','vtwo','utoc','vtoc','foc','wlv','sdc','wnm','twf','uust','vust']

# ---------------------------------------------------------------------------
# 3. Get the indices into the model grid for the observation locations:
# NB: ind will be a tuple of size 1 for SMC grids and size 2 for regular grids.
ncfile = join(RunFolder,Runs[0],Runs1[0]+'_2013120500.nc') #
ind = get_model_indices(ncfile[0] if isinstance(ncfile, (list,tuple)) else ncfile,
        obs_lon, obs_lat)

# Remove locations that are outside the grid, or on land (ind = -1)
m = ind[0] != -1
if not m.all():
    print("Removing {} of {} observations (outside grid)".format(
        np.count_nonzero(~m), len(m)))

    ind = tuple(i[m] for i in ind)
    print(ind)
    obs_lat = obs_lat[m]
    obs_lon = obs_lon[m]
    obs_id = obs_id[m]
##--
# --------------------------------------------------------------------------
# 4. Call the extraction routine to extract the points timeseries:
for i in range(len(Runs)):
    rncfile =os.listdir(join(RunFolder,Runs[i])) #
    rncfile.sort()
    ncfile = list()
    # Get the full path
    for s in range(len(rncfile)):
        ncfile.append(join(RunFolder,Runs[i],rncfile[s]))
    ncout = join(out_dir,Runs_name[i]+'_timeseries_'+obstype+'_pnts.nc')
    pnt.extract_point_timeseries_nc(ncfile, ind, obs_id,
            ncout, obs_lat=obs_lat, obs_lon=obs_lon, variables=varnames)
