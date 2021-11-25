import netCDF4 as nc
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
from os.path import join, abspath
#from collections import OrderedDict
import sys
sys.path.insert(0,'/home/h01/nvalient/fcm_python/ver_trials_scripts/matchup_for_trials/')
import waves_matchup_rw as wamrw

"""
Created on Mon Feb  8 16:07:19 2021
Module with functions to read and concatenate the different observations 

@author: nvalient
"""

# Common folder for observations
dirIn = '/project/ofrd/waves/ObsForTrials/'

##- Define individual functions for each type of observations
 
def read_SHPSYN(cycle, fcday=0):
    dirInF = join(dirIn,'MetDB_SHPSYN')
    cycstr = cycle.strftime('%Y%m%d')
    fname = 'SHPSYN_'+ cycstr + '.nc'
    ncfile = dirInF + '/' + fname
    
    return ncfile

def read_WFVS(cycle, fcday=0):
    dirInF = join(dirIn,'JCOMM_WFVS')
    cycstr = cycle.strftime('%m%Y')
    fname = 'waves_'+ cycstr[4:6] + cycstr[0:2] + '_t000'
    ncfile = dirInF + '/' + fname
    
    return ncfile    
    
def read_WAVENET(cycle, fcday=0):
    """Generates a Met Office model output file name"""      
    dirInF = join(dirIn,'MetDB_WAVENET')
    cycstr = cycle.strftime('%Y%m%d')
    fname = 'WAVENET_'+ cycstr + '.nc'
    ncfile = dirInF + '/' + fname

    return ncfile

# Function specific to JCOMM-WFVS
def var_WFVS(varname,hsobs,wsobs,wdobs,tpobs,tzobs,teobs):
    if varname == 'hs':
        datall = hsobs
    elif varname == 'ws':
        datall = wsobs
    elif varname == 'wdir':
        datall = wdobs
    elif varname == 'tp':
        datall = tpobs
    elif varname == 't02':
        datall = tzobs
    elif varname == 'te':
        datall = teobs
    
    #wsobs, wdobs, hsobs, tpobs, tzobs, teobs
    return datall

##- classes and methods to read data

# variable grid class
class varTimeseries:

    def __init__(self):
        self.varname     = None  # variable name
        self.longname    = None  # variable long name
        self.lat         = None  # latitude of model grid cell
        self.lon         = None  # longitude of model grid cell
        self.times       = None  # datetime object for forecast validity times
        self.data        = None  # variable data
        self.station_id  = None  # observation station id
        self.sites       = None  # observation station name

    def loadTimeseries(self, ncfile, varname, IDs, obs_name):
        """Load data from observation files"""
        print('[INFO] Loading in data from %s; %s' %(varname,ncfile))
        self.varname        = varname
        if obs_name == 'WFVS':
            buoylist    = '/data/cr1/frxs/waves_python/r1272_79/ver_trials_scripts/platform_lists/proposed_buoy_list_intercomparison_October2013'
            sitesobs, idsobs, latsobs, lonsobs, vtobs, wsobs, wdobs, hsobs, tpobs, tzobs, teobs = wamrw.read_wfvs(ncfile, buoylist)
            imask = idsobs == IDs
            ia = np.where(imask==1)[0][0]
            self.lat        = latsobs[imask]
            self.lon        = lonsobs[imask]
            #t = d.variables['time']
            self.times      = vtobs 
            datall          = var_WFVS(varname,hsobs,wsobs,wdobs,tpobs,tzobs,teobs)
            self.data       = np.array(datall[ia,:,0]) # Add the mask for the IDs
            self.data[self.data == -99.99]    = np.nan
            #self.units      = d.variables[var].units
            self.station_id = idsobs[imask] # to confirm we are extracting the right one
            self.sites      = sitesobs[imask]   
            
        else:
            d = nc.Dataset(ncfile)
            if obs_name == 'WAVENET':
                obs_id0 = d.variables['station_id'][:]
                obs_id  = np.char.mod('%d',obs_id0) # convert to string
            else:
                obs_id  = d.variables['station_id'][:]
            imask = obs_id == IDs
            self.longname   = d.variables[varname].long_name
            self.lat        = np.array(d.variables['latitude'][imask])
            self.lon        = np.array(d.variables['longitude'][imask])
            t = d.variables['time']
            self.times      = nc.num2date(t[:], t.units,calendar='gregorian',only_use_cftime_datetimes=False,only_use_python_datetimes=True)
            self.data       = np.array(d.variables[varname][imask,:])[0][:] # Add the mask for the IDs
            self.data[self.data == -3.2768000e+04]    = np.nan
            self.units      = d.variables[varname].units
            self.station_id = obs_id[imask] # to confirm we are extracting the right one
            self.sites      = d.variables['sites'][imask]     

            d.close()

    def concatenate_time(self, var2):
        """Concatenate two variables along the time axis"""
        print('[INFO] Concatenating time axis data')
        self.times      = np.concatenate((self.times,var2.times), axis=0)
        self.data = np.ma.concatenate((self.data,var2.data), axis=0)

    
def create_timeseries(IDs, obs_name, varname, cycle, days_lead):
    # loop for days/ months (for JCOMM_WFVS) to generate timeseries
    for i in range(days_lead):
        
        if obs_name == 'WFVS':
            delta= relativedelta(months = +1)
            print('[INFO] JCOMM WFVS are monthly')
            ncfile = read_WFVS(cycle, fcday=0)            
        elif obs_name == 'WAVENET':
            delta= dt.timedelta(days = 1)
            ncfile = read_WAVENET(cycle, fcday=0)
        elif obs_name == 'SHPSYN':
            delta= dt.timedelta(days = 1)
            ncfile = read_SHPSYN(cycle, fcday=0)

        if i == 0:
            try:
                if obs_name == 'WFVS':
                    print('[INFO] Entering varTimeseries fun - first month from %s' %ncfile)
                else:
                    print('[INFO] Entering varTimeseries fun - first day from %s' %ncfile)
                var = varTimeseries()
                var.loadTimeseries(ncfile, varname, IDs, obs_name)
            except:
                print('[ERROR] Data not available for these times from %s' %ncfile)
        else:
            try:
                if obs_name == 'WFVS':
                    print('Entering varTimeseries fun - concatenate month')
                else:
                    print('Entering varTimeseries fun - concatenate day')
                tmp = varTimeseries()
                tmp.loadTimeseries(ncfile, varname, IDs, obs_name)                
                var.concatenate_time(tmp)
            except:
                print('[ERROR] Data not available for these times from %s' %ncfile)
        # Loop over the number of days you want to concatenate
        cycle += delta
        
    return  var