#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:41:56 2021

@author: nvalient

------- SET OF FUNCTIONS TO OBTAIN THE EXTREMES VERIFICATION -----------
* The fucntions deal with several variables depending on the observation type.
Observation type does not include Merge Altimeter data.

INPUT:  -  daily matchup files
OUTPUT: - .csv files per in-situ location

"""

#import subprocess
import os.path
import matplotlib as mpl
mpl.use('Agg')
import csv
#from os.path import join
import numpy as np
import netCDF4 as nc4
import datetime 
from scipy.interpolate import griddata
from scipy.stats.stats import pearsonr
     
# ------------------------------------------------

def var_inObs(obstype,all_var=False):
    if all_var == True:
        if obstype == 'WFVS' or obstype == 'SHPSYN':
            #VAR = ['hs','ws','wdir','tp','t02']  #sometimes we have t0m1 instead of Tp
            #VAR = ['hs','ws','wdir','t02'] 
            VAR = ['hs','t02']
            #VAR = ['hs','ws','wdir','tp']
        elif obstype == 'WAVENET':
            VAR = ['hs','tp','t02','dir','spr']  #sometimes we have t0m1 instead of Tp
            #VAR = ['hs','t02','dir','spr']
            #VAR = ['t02']
            #VAR = ['hs','tp','dir']
    else:
        if obstype == 'WFVS' or obstype == 'SHPSYN':
            VAR = ['hs','ws','wdir','t02']
        elif obstype == 'WAVENET':
            VAR = ['hs','t02']
    
    return VAR

def get_NCfile(obsDir,obstype,t,run,run_folder=None):
    if run_folder is not None:
        run_dir = run_folder
    else:
        run_dir = run
        
    FILEN = os.path.join(obsDir,run_dir,'match_'+obstype+'_'+run+'_'+t+'.nc')
    return FILEN

def get_NCvar(FNAME,u_ID,VAR):  
    
    #variables to extract depend on type of obs: time, station_id, hs_obs, hs_Hx, tp_obs, tp_Hx
    nc    = nc4.Dataset(FNAME)
    nvar  = nc.variables
    units = nvar['time'].units    
    idN = np.array(nvar['station_name'][:])
    m1 =  idN == u_ID
    #m = m1.T
    m = m1

    sidJ = np.array(nvar['station_id'][:])[m]
    
    #idN = np.array(nvar['station_name'][:])[m]             
    tJ = np.array(nvar['time'][:])[m]
    COORD_ID = np.column_stack((nvar["longitude_obs"][:],nvar["latitude_obs"][:]))
    coordinate = np.array(COORD_ID); coordinate = coordinate[m,:]; 
    # get dimension of the variable to extract 
    var_dim = len(tJ)
    if var_dim == 0:
        print('[WARNING] No observations found!')        
        nc.close()
        var_obs = []
        var_mod = []
    elif var_dim != 0:
        #print('[Debug] Dimension of the variable is = '+str(var_dim))
        var_o   = np.empty((len(VAR),var_dim)); var_obs = np.full_like(var_o,np.nan,dtype=np.double)
        var_m = np.empty((len(VAR),var_dim)); var_mod = np.full_like(var_m,np.nan,dtype=np.double)
        for ii, var in enumerate(VAR):
            #print('[Debug] Which variable? '+var)
            ao = np.array(nvar[var+'_obs'][:])[m] 
            is_empty = ao.size == 0 
            if is_empty == True:
                print('[WARNING] No observations found for '+ var)
            elif is_empty == False:
                print('[INFO] Observations found for '+ var)
                var_obs[ii] = np.array(nvar[var+'_obs'][:])[m] # OBS
                var_mod[ii] = np.array(nvar[var+'_Hx'][:])[m]  # MODEL OUTPUT
        nc.close()
    
    return tJ,var_obs,var_mod,sidJ,idN,units,var_dim,coordinate


def get_depth_offshelf(obstype,ndir,run,tini,run_folder=None):
    
    """Function to obtain the index (I_IN AND I_OFF) of the observations that are located on-shelf and off-shelf"""
    
    # Reading bathymetry
    bathy_file = '/data/users/nvalient/bathymetry/NWshelf_amm15_bathy.nc'
    print("Reading "+bathy_file)
    fbathy     = nc4.Dataset(bathy_file,"r")
    bat        = np.array(fbathy.variables["Bathymetry"])
    lat        = np.array(fbathy.variables["lat"])
    lon        = np.array(fbathy.variables["lon"])
    fbathy.close()

    # ------------------------------------------------
    # Try extracting Ids using a matchup Dataset
    Jfile      = get_NCfile(ndir,obstype,tini,run,run_folder)
    obs        = nc4.Dataset(Jfile,"r")
    lat_new    = np.array(obs.variables["latitude_obs"])
    lon_new    = np.array(obs.variables["longitude_obs"])
    IDs        = np.array(obs.variables["station_id"])
    LOC_N      = np.array(obs.variables['station_name'])
    
    # ------------------------------------------------
    # Estimate on/off shelf
    points = np.array( (lon.flatten(), lat.flatten()) ).T
    values = bat.flatten()

    print('Interpolating obs location to nearest bathy points')
    depth_ID = griddata(points, values, (lon_new, lat_new), method='nearest')

    I_IN  = depth_ID <= 200
    I_OFF = depth_ID > 200
    # Sanity check; number of observations during the storms
    n_in  = sum(I_IN) # 396
    n_off = sum(I_OFF) # 64
    COORD_ID = np.column_stack((obs.variables["longitude_obs"][:],obs.variables["latitude_obs"][:]))
    obs.close()
    print('Number of observations on-shelf = '+str(n_in))
    print('Number of observations off-shelf = '+str(n_off))
    
    return I_IN, I_OFF, IDs, LOC_N, COORD_ID

# ---------------------------------------------------------------------------
def get_obs_ID(obstype,ndir,run,tini,run_folder=None):
    
    """Function to obtain Ids, coordinates and names of the observations """
    
    # ------------------------------------------------
    # Try extracting Ids using a matchup Dataset
    Jfile      = get_NCfile(ndir,obstype,tini,run,run_folder)
    obs        = nc4.Dataset(Jfile,"r")
    IDs        = np.array(obs.variables["station_id"])
    LOC_N      = np.array(obs.variables['station_name'])
    
    # ------------------------------------------------
    COORD_ID = np.column_stack((obs.variables["longitude_obs"][:],obs.variables["latitude_obs"][:]))
    obs.close()
    
    return IDs, LOC_N, COORD_ID

# ---------------------------------------------------------------------------
# Define function to compute RMSE
def rmse(predictions, targets):
    a = (predictions - targets) ** 2
    return np.sqrt(np.nanmean(a))

# --------------------------------------------------------------------------- 
    
def get_timeseries(ndir,obstype,run,TINI,TEND,u_ID,VAR,run_folder=None):
   
    """Function to build the timeseries per location (IDs)"""
       
      
    # Iterate through the days you are interested in              
    start_date = datetime.date(int(TINI[0:4]),int(TINI[4:6]),int(TINI[6:8]))
    end_date   = datetime.date(int(TEND[0:4]),int(TEND[4:6]),int(TEND[6:8]))
    delta = datetime.timedelta(days=1)
            
    while start_date <= end_date:
        filedate = start_date.strftime("%Y%m%d")
        print(['[INFO] Extracting date = '+filedate])
        # build file name, e.g., match_SHPSYN_UKC4aow-st6_20140208.nc
        fileName = get_NCfile(ndir,obstype,filedate,run,run_folder)
            
        # variables to extract: time, station_id, hs_obs, hs_Hx, tp_obs, tp_Hx
        if filedate == TINI:
            tJ,var_obs,var_mod,sidJ,idN,units,var_dim,coordinates = get_NCvar(fileName,u_ID,VAR)
        else:
            tJi,var_obsi,var_modi,sidJi,idNi,units,var_dim,coordinatesi = get_NCvar(fileName,u_ID,VAR)
            if var_dim != 0:
                tJ = np.hstack((tJ,tJi))
                var_obs = np.concatenate((var_obs,var_obsi),axis=1)
                var_mod = np.concatenate((var_mod,var_modi),axis=1)

        #print ([start_date.strftime("%Y%m%d") +' - '+ np.str(len(hs_obsJ))])
        start_date += delta
    
    datesJ = nc4.num2date(tJ,units,'gregorian')
    
    # Include NaNs in fill values if any
    Nm = var_obs == -32768
    var_obs[Nm] = np.NaN
    station_ID = sidJ[0]
    name_ID = idN[0]
    #print('[Debug] coordinates are ='+str(coordinates))
    COORD = coordinates[0]
    
    return var_mod, var_obs, datesJ, station_ID, name_ID, COORD


def get_extremes_CSV(ndir,out_dir,obstype,COORD_ID,LOC_N,RUN,TINI,TEND,Q1,Q2,VAR,run_folder=None,opt=None):
    """Function to call the timeseries function per location (IDs) and get the extremes"""
    # BUILD TIMESERIES AND GET EXTREMES
    
    var = VAR
    q1 = str(Q1)
    q2 = str(Q2)
    
    # Get unique locations and their indices
    #print('[Debug] Location name is : '+LOC_N)
    u, indices = np.unique(LOC_N,return_index=True)
    #print('[Debug] Location name UNIQUE is : '+u)
    #u, indices = np.unique(LOC_N,axis=0,return_index=True)
    
    # Create empty values for stats: BIAS, PIERSON CORR COEF, RMSE 
    bs   = np.empty([len(var),len(indices)]); bstorm = np.full_like(bs,np.nan,dtype=np.double)   
    be   = np.empty([len(var),len(indices)]); bext = np.full_like(be,np.nan,dtype=np.double)
    rs   = np.empty([len(var),len(indices)]); rstorm = np.full_like(rs,np.nan,dtype=np.double)   
    re   = np.empty([len(var), len(indices)]); rext = np.full_like(re,np.nan,dtype=np.double)
    Es   = np.empty([len(var), len(indices)]); Estorm = np.full_like(Es,np.nan,dtype=np.double)   
    Ee   = np.empty([len(var), len(indices)]); Eext = np.full_like(Ee,np.nan,dtype=np.double)
    
    #print('[Debug] the shape of bstorm is ' + str(bstorm.shape))
    
    var_limS = [] 
    var_limE = []
    loc_n  = []
    loc_id = []
    Y      = []
    X      = []   
    for kk in range(len(indices)):
        # coordinate of each in-situ location
        u_ID = u[kk]   
        print('[INFO] Computation of '+u_ID+' for '+RUN+' - COORD = '+str(COORD_ID[indices[kk]])) 
        
        var_mod, var_obs, datesJ, sidJ, idN, COORD = get_timeseries(ndir,obstype,RUN,TINI,TEND,u_ID,VAR,run_folder)
        # TO DO: Possibility to add a condition for when we want to plot the timeseries?
        
        # ----------------------------------------------------------------------
        # USE the HS in order to get the extremes ([0] in var_mod and var_obs)
        # for OBS
        var_storm = np.nanpercentile(var_obs[0], Q1) # quantile for Hs; 
        var_ext = np.nanpercentile(var_obs[0], Q2) 
        # for MODEL
        var_stormM = np.nanpercentile(var_mod[0], Q1)
        var_extM = np.nanpercentile(var_mod[0], Q2)

        # Loop to get stats per location  
        if np.isnan(var_obs).all() == True:
            print('No obs for '+u_ID)
        #if 
        elif np.isnan(var_storm) == False: # when observations          
            if opt == 1 or opt == None:
                # obtain mask; obs + model condition > Q
                print('[INFO] the threshold for Hs,q75 is '+str(var_storm)+' [m]') 
                mstorm     = np.logical_and((var_obs[0] >= var_storm),(var_mod[0] >= var_stormM)) # Q75
                mext       = np.logical_and((var_obs[0] >= var_ext),(var_mod[0] >= var_extM)) # Q90
            elif opt == 2:
                # obtain mask; smaller threshold to use as obs and model condition > Q
                print('[INFO] the threshold for Hs,q75 is '+str(var_storm)+' [m]') 
                min_storm = min([var_storm,var_stormM])
                min_ext   = min([var_ext,var_extM]) 
                mstorm     = np.logical_and((var_obs[0] >= min_storm),(var_mod[0] >= min_storm)) # Q75
                mext       = np.logical_and((var_obs[0] >= min_ext),(var_mod[0] >= min_ext)) # Q90
                  
            # Do some QC for the extremes; is the sample statistically significant? only if p-value < 0.1
            print('[Debug] var_mod is :'+str(var_mod))
            print('[Debug] var_mod is :'+str(var_obs))
            QC = pearsonr(var_mod[0,mstorm],var_obs[0,mstorm])
            print('[INFO] The p-value for Hs,q75 is ' +str(QC[1]))
            if QC[1] <= 0.1:                
                for var_i in range(len(var)):
                    bstorm[var_i,kk] = np.nanmean(var_mod[var_i,mstorm]-var_obs[var_i,mstorm]) # bias Q75
                    bext[var_i,kk]   = np.nanmean(var_mod[var_i,mext]-var_obs[var_i,mext]) # bias Q90
                    rext[var_i,kk]   = np.ma.corrcoef(np.ma.masked_invalid(var_mod[var_i,mext]),\
                                                      np.ma.masked_invalid(var_obs[var_i,mext]))[0,1] # r, pierson coef.
                    rstorm[var_i,kk] = np.ma.corrcoef(np.ma.masked_invalid(var_mod[var_i,mstorm]),\
                                                      np.ma.masked_invalid(var_obs[var_i,mstorm]))[0,1] # r, pierson coef.
                    Estorm[var_i,kk] = rmse(var_mod[var_i,mstorm],var_obs[var_i,mstorm])
                    Eext[var_i,kk]   = rmse(var_mod[var_i,mext],var_obs[var_i,mext])                     
                
        # save values per location bias, r, RMSE
        var_limS.append(var_storm)
        var_limE.append(var_ext)
        loc_n.append(sidJ)
        loc_id.append(u_ID)
        Y.append(COORD[1])
        X.append(COORD[0])   
        print('------------------------------------------------------------')
    
    # Create list and write .csv file per variable
    for iv, var in enumerate(VAR):
        var_nq1 = var+'_'+q1
        var_nq2 = var+'_'+q2
        header = ['ID','Name','Latitude','Longitude',var_nq1,'bias_'+q1,'R_'+q1,'RMSD_'+q1,var_nq2,'bias_'+q2,'R_'+q2,'RMSD_'+q2] 
              
        with open(os.path.join(out_dir,obstype,var+'_'+RUN+'_quantileStats.csv'), mode='w+') as f_75:
            f_75w = csv.writer(f_75, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            f_75w.writerow(header)
            f_75w.writerows(np.transpose([loc_id,loc_n,Y,X,var_limS,bstorm[iv,:],rstorm[iv,:],Estorm[iv,:],\
                                          var_limE,bext[iv,:],rext[iv,:],Eext[iv,:]]))
    
        print('.CSV File saved for '+RUN+' - Variable '+var)  
    
    return

def get_mean_extremes(stat):
    return np.nanmean(np.array(stat))

def read_extremes(out_dir,obstype,RUN,Q1,Q2,var):
    """Function to read the .CSV files with the stats for the extremes;
    Stats are per location (IDs)
    
    Average stats are printed"""
    
    path_to_file = os.path.join(out_dir,obstype,var+'_'+RUN+'_quantileStats.csv')
    
    if os.path.exists(path_to_file) == True:     
        lat = []
        lon = []
        b1  = []
        b2  = []
        rmse1 = []
        rmse2 = []
        r1  = []
        r2  = []
        
        q1 = str(Q1)
        q2 = str(Q2)
        
        # Read EXTREMES.CSV file per variable
      
        var_nq1 = var+'_'+q1
        var_nq2 = var+'_'+q2
        header = ['ID','Name','Latitude','Longitude',var_nq1,'bias_'+q1,'R_'+q1,'RMSD_'+q1,var_nq2,'bias_'+q2,'R_'+q2,'RMSD_'+q2] 
                  
        with open(os.path.join(out_dir,obstype,var+'_'+RUN+'_quantileStats.csv'), mode='r') as csvfile:
            print('[Info] Reading EXTREMES.CSV File for '+RUN+' - Variable '+var)
            csvreader = csv.DictReader(csvfile, delimiter=',')

            for lines in csvreader:
                b1.append(lines['bias_'+q1])
                b2.append(lines['bias_'+q2])
                rmse1.append(lines['RMSD_'+q1])
                rmse2.append(lines['RMSD_'+q2])
                lat.append(lines['Latitude'])
                lon.append(lines['Longitude'])
                r1.append(lines['R_'+q1])
                r2.append(lines['R_'+q2])
        
        for i in range(len(b1)):
            if b1[i] == 'nan':
                b1[i] = np.nan
                b2[i] = np.nan
                r1[i] = np.nan
                r2[i] = np.nan
                rmse1[i] = np.nan
                rmse2[i] = np.nan
            else:
                b1[i] = float(b1[i])
                b2[i] = float(b2[i])
                r1[i] = float(r1[i])
                r2[i] = float(r2[i])
                rmse1[i] = float(rmse1[i])
                rmse2[i] = float(rmse2[i])
        lat = list(map(float,lat))
        lon = list(map(float,lon))
        
        #print('[Debug] bias = ' + str(b1))            
        Mstats = [get_mean_extremes(b1),get_mean_extremes(r1),get_mean_extremes(rmse1),get_mean_extremes(b2),get_mean_extremes(r2),\
                  get_mean_extremes(rmse2)] 
        print(obstype+' mean stats for '+var+' - [b75,r75,rmse75,b90,r90,rmse90] = ')
        print(str(Mstats))
        
    elif os.path.exists(path_to_file) == False:
        print('.csv for '+var+' - '+obstype+' does not exist!')
        
    return lat, lon, b1, b2, rmse1, rmse2, r1, r2 
