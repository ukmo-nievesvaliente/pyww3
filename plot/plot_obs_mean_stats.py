#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:14:05 2021

    SET OF FUNCTIONS TO PLOT MEAN STATISTICS OF MODEL RUNS AGAINST  IN-SITU OBSERVATIONS
    ADDITIONAL FUNCTION TO PLOT FCST VERIFICATION USING .csv FILES 
    TO dO: add the possibility to compute basic stats from MA data (including the gridding option)

@author: nvalient
"""

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from pathlib import Path
import obs_funs.obs_reader as obsread

transform = ccrs.PlateCarree()

def fill_land_cfeature(color='silver'):
    """Fill land when using Cartopy
       Inputs:
           color: chosen color to fill the land; silver (default)
       Output:
           land_50; to be usesd as axes.add_feature(land_50)"""
           
    if color.lower() != 'silver':
        land_50 = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='none',
                                        facecolor=color) #facecolor=cfeature.COLORS['land'])
    else:
        
        land_50 = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='none',
                                        facecolor='silver') #facecolor=cfeature.COLORS['land'])
    return land_50

def get_contours_NWS_bat():
        
    # Reading bathymetry
    bathy_file = '/data/users/nvalient/bathymetry/NWshelf_amm15_bathy.nc'
    print("Reading "+bathy_file)
    fbathy     = nc.Dataset(bathy_file,"r")
    bat        = np.array(fbathy.variables["Bathymetry"])
    lat        = np.array(fbathy.variables["lat"])
    lon        = np.array(fbathy.variables["lon"])
    fbathy.close()
    
    return lon, lat, bat

def get_contours_GBL_bat():
        
    # Reading bathymetry
    bathy_file = '/data/users/nvalient/GMD_paper/data/gbl_dpt_4xres.nc'
    print("Reading "+bathy_file)
    fbathy     = nc.Dataset(bathy_file,"r")
    bathy      = np.array(fbathy.variables["dpt"][0,:,:])
    lat        = np.array(fbathy.variables["standard_latitude"])
    lon        = np.array(fbathy.variables["standard_longitude"])
    fbathy.close()
    bat = np.ma.masked_where(bathy < 200., bathy)
    
    return lon, lat, bat

def var_inObs(all_var=False):
    if all_var == True:
        VAR_PLT = ['Hs','Ws','Wdir','Tp','T02','Dir','Spr']            
    else:
        #VAR_PLT = ['Hs','T02']
        VAR_PLT = ['Hs','Ws','T02']
    
    return VAR_PLT


def get_dict():
    
    """
    Created 22 June 15:17:23 2022
    
    Function to create dictionary with some of the mean variables that can be read from .csv file
    Values can be extracted using previously the  read_summary_csv function
    E.g.: import obs_funs.obs_reader as obsread 
          obsread.read_summary_csv(fileIn) 
    
        
    @author: nvalient
    """
    # ------------------ MEAN -----------------------
    VAR = {}
    VAR['mean_ws']={}
    VAR['mean_ws']['short_n'] = 'U10'
    VAR['mean_ws']['colorbar'] = 'gist_ncar'
    VAR['mean_ws']['limits'] = [0,16]
    VAR['mean_ws']['limits_diff'] = [-4,4]
    
    VAR['mean_wdir']={}
    VAR['mean_wdir']['short_n'] = 'U10 dir'
    VAR['mean_wdir']['colorbar'] = 'PiYG_r'
    VAR['mean_wdir']['limits'] = [0,360]
    VAR['mean_wdir']['limits_diff'] = [-4,4]
    
    VAR['mean_hs']={}
    VAR['mean_hs']['short_n'] = 'Hs'
    VAR['mean_hs']['colorbar'] = 'nipy_spectral'
    VAR['mean_hs']['limits'] = [0,7]
    VAR['mean_hs']['limits_diff'] = [-2,2]
    
    VAR['mean_tp']={}
    VAR['mean_tp']['short_n'] = 'Tp'
    VAR['mean_tp']['colorbar'] = 'hsv'
    VAR['mean_tp']['limits'] = [0,12]
    VAR['mean_tp']['limits_diff'] = [-2,2]
    
    VAR['mean_t02']={}
    VAR['mean_t02']['short_n'] = 'T02'
    VAR['mean_t02']['colorbar'] = 'hsv'
    VAR['mean_t02']['limits'] = [0,12]
    VAR['mean_t02']['limits_diff'] = [-2,2]  
    # ------------------ BIAS -----------------------
    
    VAR['bias_ws']={}
    VAR['bias_ws']['short_n'] = 'U10 bias'
    VAR['bias_ws']['colorbar'] = 'PiYG_r'
    VAR['bias_ws']['limits'] = [-2,2]
    VAR['bias_ws']['limits_diff'] = [-4,4]

    VAR['bias_wdir']={}
    VAR['bias_wdir']['short_n'] = 'U10 dir bias'
    VAR['bias_wdir']['colorbar'] = 'PiYG_r'
    VAR['bias_wdir']['limits'] = [-15,15]
    VAR['bias_wdir']['limits_diff'] = [-4,4]
    
    VAR['bias_hs']={}
    VAR['bias_hs']['short_n'] = 'Hs bias'
    VAR['bias_hs']['colorbar'] = 'PiYG_r'
    VAR['bias_hs']['limits'] = [-0.4,0.4]
    VAR['bias_hs']['limits_diff'] = [-2,2]
    
    VAR['bias_tp']={}
    VAR['bias_tp']['short_n'] = 'Tp bias'
    VAR['bias_tp']['colorbar'] = 'PiYG_r'
    VAR['bias_tp']['limits'] = [-3,3]
    VAR['bias_tp']['limits_diff'] = [-2,2]
    
    VAR['bias_t02']={}
    VAR['bias_t02']['short_n'] = 'T02 bias'
    VAR['bias_t02']['colorbar'] = 'PiYG_r'
    VAR['bias_t02']['limits'] = [-2,2]
    VAR['bias_t02']['limits_diff'] = [-2,2]    
    
    # ------------------ RMSD -----------------------

    VAR['rmse_ws']={}
    VAR['rmse_ws']['short_n'] = 'U10 RMSD'
    VAR['rmse_ws']['colorbar'] = 'gist_stern_r'
    VAR['rmse_ws']['limits'] = [0,3]
    VAR['rmse_ws']['limits_diff'] = [-1,1]
    
    VAR['rmse_wdir']={}
    VAR['rmse_wdir']['short_n'] = 'U10 dir RMSD'
    VAR['rmse_wdir']['colorbar'] = 'gist_stern_r'
    VAR['rmse_wdir']['limits'] = [0,60]
    VAR['rmse_wdir']['limits_diff'] = [-1,1]
  
    VAR['rmse_hs']={}
    VAR['rmse_hs']['short_n'] = 'Hs RMSD'
    VAR['rmse_hs']['colorbar'] = 'gist_stern_r'
    VAR['rmse_hs']['limits'] = [0,0.8]
    VAR['rmse_hs']['limits_diff'] = [-1,1]

    VAR['rmse_tp']={}
    VAR['rmse_tp']['short_n'] = 'Tp RMSD'
    VAR['rmse_tp']['colorbar'] = 'gist_stern_r'
    VAR['rmse_tp']['limits'] = [0,2.5]
    VAR['rmse_tp']['limits_diff'] = [-1,1]
 
    VAR['rmse_t02']={}
    VAR['rmse_t02']['short_n'] = 'T02 RMSD'
    VAR['rmse_t02']['colorbar'] = 'gist_stern_r'
    VAR['rmse_t02']['limits'] = [0,2]
    VAR['rmse_t02']['limits_diff'] = [-1,1]   
    
    # ------------------ Error STD -----------------------
    VAR['error_std_t02']={}
    VAR['error_std_t02']['short_n'] = 'T02 StdE'
    VAR['error_std_t02']['colorbar'] = 'brg'
    VAR['error_std_t02']['limits'] = [0,2]
    VAR['error_std_t02']['limits_diff'] = [-1,1]
    
    VAR['error_std_hs']={}
    VAR['error_std_hs']['short_n'] = 'Hs StdE'
    VAR['error_std_hs']['colorbar'] = 'brg'
    VAR['error_std_hs']['limits'] = [0,0.8]
    VAR['error_std_hs']['limits_diff'] = [-1,1] 
    
    VAR['error_std_ws']={}
    VAR['error_std_ws']['short_n'] = 'U10 StdE'
    VAR['error_std_ws']['colorbar'] = 'brg'
    VAR['error_std_ws']['limits'] = [0,3]
    VAR['error_std_ws']['limits_diff'] = [-1,1]
    
    # ------------------ PIERSON corr coef -----------------------
    
    VAR['pierson_t02']={}
    VAR['pierson_t02']['short_n'] = 'T02 r'
    VAR['pierson_t02']['colorbar'] = 'terrain_r'
    VAR['pierson_t02']['limits'] = [0.75,1.]
    VAR['pierson_t02']['limits_diff'] = [-1,1]
    
    VAR['pierson_hs']={}
    VAR['pierson_hs']['short_n'] = 'Hs r'
    VAR['pierson_hs']['colorbar'] = 'terrain_r'
    VAR['pierson_hs']['limits'] = [0.9,1.]
    VAR['pierson_hs']['limits_diff'] = [-1,1] 
    
    VAR['pierson_ws']={}
    VAR['pierson_ws']['short_n'] = 'U10 r'
    VAR['pierson_ws']['colorbar'] = 'brg'
    VAR['pierson_ws']['limits'] = [0.85,1.]
    VAR['pierson_ws']['limits_diff'] = [-1,1] 
        
    return VAR

def var_inObs4plot(var):
    if var == 'Hs':
        VAR_4PLOT = ['mean_hs','bias_hs','rmse_hs','error_std_hs','pierson_hs']
    elif var == 'Ws':
        VAR_4PLOT = ['mean_ws','bias_ws','rmse_ws','error_std_ws','pierson_ws']
    elif var == 'Tp':
        VAR_4PLOT = ['bias_tp','rmse_tp']
    elif var == 'T02':
        VAR_4PLOT = ['mean_t02','bias_t02','rmse_t02','error_std_t02','pierson_t02']
    elif var == 'Wdir':
        VAR_4PLOT = ['bias_wdir','rmse_wdir']
    
    return VAR_4PLOT

def plot_obs_stats(out_dir,lon_stat,lat_stat,var,val_stat,run,nwshelf=True):
    
    land_50 = fill_land_cfeature()
    VAR_4PLOT = var_inObs4plot(var)
    if nwshelf is True:
        lon, lat, bat = get_contours_NWS_bat()
        levels = [40,80,120,240]
        print('[INFO] Domain is NWshelf')
    if nwshelf is not True:
        lon, lat, bat = get_contours_GBL_bat()
        levels = [200,500,1000,2000,3000,4000,5000]
        print('[INFO] Domain is GBL')
        
    VARs = get_dict()
    
    #print('[Debug] vars 4 plot are '+str(VAR_4PLOT))
    for ii,stat in enumerate(VAR_4PLOT):
        # Iterate over the different stats that can be plotted
        rr = np.array(val_stat[ii])
        
        fig2 = plt.figure()
        #axes = fig2.add_subplot(111,projection=ccrs.PlateCarree())
        axes = fig2.add_subplot(111,projection=ccrs.Mercator(central_longitude=0))
        a = axes.contour(lon,lat,bat, levels, colors='grey',linewidths=.25,transform = transform)
        axes.coastlines(resolution='50m', color='black', linewidth=1)
        axes.add_feature(land_50,zorder=1)
        if nwshelf is True:
            axes.clabel(a, inline=True, fmt = '%3.0f', fontsize=6)
            e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*16,cmap=VARs[stat]['colorbar'],\
                             vmin=VARs[stat]['limits'][0],vmax=VARs[stat]['limits'][1],transform = transform, zorder=2)
            cbar = fig2.colorbar(e,extend='both')
            axes.scatter(lon_stat,lat_stat,c='k',marker="x",s=0.6,transform = transform,zorder=3)
        if nwshelf is not True:
            e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*2,cmap=VARs[stat]['colorbar'],\
                             vmin=VARs[stat]['limits'][0],vmax=VARs[stat]['limits'][1],transform = transform,zorder=2)
            cbar = fig2.colorbar(e,shrink=0.7,extend='both')
            #axes.scatter(lon_stat,lat_stat,c='k',marker="x",s=0.4,zorder=3)
        #cbar.set_label(VAR[stat]['short_n'])
        #axes.gridlines()
        axes.set_axisbelow(True)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axes.xaxis.set_major_formatter(lon_formatter)
        axes.yaxis.set_major_formatter(lat_formatter)
        if nwshelf is True:
            axes.set_xticks([-12,-6,0,6],crs=ccrs.PlateCarree())
            axes.set_yticks([48,54,60],crs=ccrs.PlateCarree())
            axes.set_extent([-14,11,45,63],crs=ccrs.PlateCarree())
            # axes.set_xlim([-13,9])
            # axes.set_ylim([45,63])
        plt.title(VARs[stat]['short_n']+' for '+run,fontsize=8)
        out_name2 = join(out_dir,run+'_'+stat+'.png')
        print("Saving figure " +out_name2)
        plt.savefig(out_name2,bbox_inches="tight", pad_inches=0.1,dpi=300)
        plt.close()
    
    return

def plot_obs_stats_r(out_dir,lon_stat,lat_stat,var,val_stat,run,nwshelf=True):
    
    land_50 = fill_land_cfeature()
    VAR_4PLOT = var_inObs4plot(var) # includes the mean values; not relevant to plot the relative change as it should be =0
    
    if nwshelf is True:
        lon, lat, bat = get_contours_NWS_bat()
        levels = [40,80,120,240]
        print('[INFO] Domain is NWshelf')
    if nwshelf is not True:
        lon, lat, bat = get_contours_GBL_bat()
        levels = [200,500,1000,2000,3000,4000,5000]
        print('[INFO] Domain is GBL')
    
    VAR = get_dict()
    
    for ii,stat in enumerate(VAR_4PLOT):
        
        if len(VAR_4PLOT) == 5 and ii == 0: # jump plotting relative mean values
            continue
        else:           
            
            # Iterate over the different stats that can be plotted
            rr = val_stat[ii]
            
            fig2 = plt.figure()
            #axes = fig2.add_subplot(111,projection=ccrs.PlateCarree())
            axes = fig2.add_subplot(111,projection=ccrs.Mercator())
            a = axes.contour(lon,lat,bat, levels, colors='grey',linewidths=.25,transform = transform)
            axes.coastlines(resolution='50m', color='black', linewidths=1)
            axes.add_feature(land_50, zorder=1)
            if nwshelf is True:
                axes.clabel(a, inline=True, fmt = '%3.0f', fontsize=6)
                e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*12,cmap='seismic',\
                                 vmin=VAR[stat]['limits_diff'][0],vmax=VAR[stat]['limits_diff'][1], transform = transform,zorder=2)
                cbar = fig2.colorbar(e,extend='both')
                axes.scatter(lon_stat,lat_stat,c='k',marker="x",s=0.6,transform = transform,zorder=3)
            if nwshelf is not True:
                e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*2,cmap='seismic',\
                                 vmin=VAR[stat]['limits_diff'][0],vmax=VAR[stat]['limits_diff'][1], transform = transform,zorder=2)     
                cbar = fig2.colorbar(e,shrink=0.7,extend='both')
                #axes.scatter(lon_stat,lat_stat,c='k',marker="x",s=0.4,zorder=3)
            #cbar.set_label(VAR[stat]['short_n'])
            #axes.gridlines()
            axes.set_axisbelow(True)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            axes.xaxis.set_major_formatter(lon_formatter)
            axes.yaxis.set_major_formatter(lat_formatter)
            if nwshelf is True:
                axes.set_xticks([-12,-6,0,6],crs=ccrs.PlateCarree())
                axes.set_yticks([48,54,60],crs=ccrs.PlateCarree())
                # axes.set_xlim([-13,9])
                # axes.set_ylim([45,63])
                axes.set_extent([-14,11,45,63],crs=ccrs.PlateCarree())
            plt.title(VAR[stat]['short_n']+' for '+run+' relative to ctr',fontsize=8)
            out_name2 = join(out_dir,run+'_relative2ctrl_'+stat+'.png')
            print("Saving figure " +out_name2)
            plt.savefig(out_name2,bbox_inches="tight", pad_inches=0.1,dpi=300)
            plt.close()
    
    return

def get_var_from_obstype(obs_type):
    
    if obs_type == 'WFVS':
        var_obstype = ['Hs','T02','Tp','Wdir','Ws']
        var_title = ['$H_s$ [$m$]','$T_{02}$ [$s$]','$T_p$ [$s$]' ,'$U_{10}\ dir.$ [$^\circ$]','$U_{10}$ [$ms^{-1}$] ']
    elif obs_type == 'WAVENET':
        var_obstype = ['Hs','T02','Tp','Spr','Dirn']
        var_title = ['$H_s$ [$m$]','$T_{02}$ [$s$]','$T_p$ [$s$]','Spr','$Dir.$ [$^\circ$]']
    elif obs_type == 'SHPSYN':
        var_obstype = ['Hs','T02','Tp','Wdir','Ws']
        var_title = ['$H_s$ [$m$]','$T_{02}$ [$s$]','$T_p$ [$s$]' ,'$U_{10}\ dir.$ [$^\circ$]','$U_{10}$ [$ms^{-1}$] ']
    elif obs_type == 'MA_SUP03':
        var_obstype = ['Hs','Ws']
        var_title = ['$H_s$ [$m$]','$U_{10}$ [$ms^{-1}$] ']
    
    return var_obstype,var_title

def plot_FCST(inFolder,inFolderCTRL,date_str,FCST_LEN,obs_type,OUT_DIR):
    
    """"
    Function to plot the mean bias and RMSD of model run versus control run
        Input parameters:
            - inFolder = complete path to verification including the running folder
            - inFolderCTRL = complete path to verification including the running folder CTRL
            - date_str = string with STARTDATE_ENDDATE e.g., '20191204_20200125'
            - obs_type = string with the observations name; i.e., WAVENET, SHPSYN, WFVS, MA_SUP03
            - OUT_DIR = complete path where the plot should be stored
            - FCST_LEN = array with strings with the name of the lead time as per folder in verification  e.g., ['T+24','T+48']
            
        Output: .png with two subplots including stats (bias and RMSD) over FCST lead time 
    """
    
    var_obstype,var_title = get_var_from_obstype(obs_type)
    print('Observations are '+obs_type+' with variables '+str(var_obstype))
    
    
    for indvar, ivar in enumerate(var_obstype):
        biasActr = []
        biasA    = []
        RMSDActr = []
        RMSDA    = []
        
        for ia, ifcst in enumerate(FCST_LEN):
            # CONTROL (CTRL)
            filenamectr = join(inFolderCTRL,'plots',ifcst,'SummaryStats_'+ifcst+'_'+obs_type+'_'+date_str+'_'+ivar+'.csv')
            areactr, BIASctr, RMSDctr, RVALUEctr, STDERRORctr = obsread.read_summary_csv(filenamectr)
            print('Reading data from '+filenamectr)
            # NEW RUN/PS 
            filename = join(inFolder,'plots',ifcst,'SummaryStats_'+ifcst+'_'+obs_type+'_'+date_str+'_'+ivar+'.csv')
            area, BIAS, RMSD, RVALUE, STDERROR = obsread.read_summary_csv(filename)
            print('Reading data from '+filename)
            
            # concatenate the values per variable as areas are expected to match between trials matchup 
            # At the moment only for BIAS and RMSD        
            
            area_KEY = area
            area_KEYctr = areactr
            # Check if areas are the same in both runs
    
            if len(area_KEY) != len(area_KEYctr):
                 if len(area_KEY) > len(area_KEYctr):
                    imatch = [area_KEY.index(iw) for iw in area_KEYctr]
                    biasA0=np.array(BIAS,dtype=float)[imatch]
                    biasActr0 = BIASctr
                   
                    RMSDA0=np.array(RMSD,dtype=float)[imatch]
                    RMSDActr0 = RMSDctr
                    area_KEY0 = area_KEYctr
                 elif len(area_KEYctr) > len(area_KEY):
                    imatch = [area_KEYctr.index(iw) for iw in area_KEY]
                    #print('imatch is'+str(imatch))
                    biasActr0=np.array(BIASctr,dtype=float)[imatch]
                    #print('[DEBUG] '+str(biasActr0.size))  
                    biasA0 = BIAS
                        
                    RMSDActr0 = np.empty(len(area_KEY))
                    RMSDActr0[:] = np.NaN
                    RMSDActr0=np.array(RMSDctr,dtype=float)[imatch]
                    RMSDA0 = RMSD 
                    area_KEY0 = area_KEY
            else:
                # Same areas in both runs
                biasActr0 = BIASctr
                biasA0    = BIAS
                RMSDActr0 = RMSDctr
                RMSDA0    = RMSD                
                area_KEY0 = area
                
            if ia ==0:
                biasActr = biasActr0
                biasA    = biasA0
                RMSDActr = RMSDActr0
                RMSDA    = RMSDA0 
                area_KEY_final = area_KEY0
            if ia != 0:
                # Need to make sure that areas are the same for the =! lead times:
                imatch=[]
                # Search for area
                for iw in area_KEY_final:
                    try:
                        imatch.append(area_KEY0.index(iw))
                    except ValueError:
                        print('[WARNING] '+iw+' is not present')
                                  
                # imatch = [area_KEY0.index(iw) for iw in area_KEY] # index = mylist.index(element)
                #-----------------------------------------------------------------
                # Create nan array of length areas key
                arrNan = np.empty(len(area_KEY_final))
                arrNan[:] = np.NaN
                arrNanctr = np.copy(arrNan)
                arrNan[imatch]=np.array(biasA0,dtype=float)        
                arrNanctr[imatch]=np.array(biasActr0,dtype=float)
                           
                biasActr = np.vstack((biasActr,arrNanctr))
                biasA    = np.vstack((biasA,arrNan))
                
                arrNan[imatch]=np.array(RMSDA0,dtype=float)  
                arrNanctr[imatch]=np.array(RMSDActr0,dtype=float)
                RMSDActr = np.vstack((RMSDActr,arrNanctr))
                RMSDA    = np.vstack((RMSDA,arrNan))
                
        # transpose the array
        biasActr = np.transpose(biasActr)
        biasA    = np.transpose(biasA)
        RMSDActr = np.transpose(RMSDActr)
        RMSDA    = np.transpose(RMSDA)
        
        for i,area_name in enumerate(area_KEY_final):
            print('Plotted area is '+str(area_name))
            # setup figure
            f = plt.figure(figsize=(7,5))
            
            axs = f.add_subplot(211)
            axs.plot(range(len(FCST_LEN)),biasActr[i], label='CTRL',
                           c ='k', marker='o', ls='-', lw=2, ms=6)
            axs.plot(range(len(FCST_LEN)),biasA[i],label='PS45',
                           c ='g', marker='o', ls='-', lw=2, ms=6)
            axs.set_xticks(range(len(FCST_LEN))) 
            axs.set_xticklabels(FCST_LEN)
            plt.setp(axs.get_yticklabels(), fontsize=12)
            plt.setp(axs.get_xticklabels(), fontsize=12)
            axs.set_ylabel('Bias',fontsize=12)
            axs.legend()
            axs.set_title(var_title[indvar], fontsize =12)
            axs.grid(True, lw=0.5, ls=':', c='gray')
                      
            axs = f.add_subplot(212)
            axs.plot(range(len(FCST_LEN)),RMSDActr[i],label='CTRL',
                           c ='k', marker='o', ls='-', lw=2, ms=6)
            axs.plot(range(len(FCST_LEN)),RMSDA[i],label='PS45',
                           c ='g', marker='o', ls='-', lw=2, ms=6)
            axs.set_xticks(range(len(FCST_LEN))) 
            axs.set_xticklabels(FCST_LEN)
            plt.setp(axs.get_yticklabels(), fontsize=12)
            plt.setp(axs.get_xticklabels(), fontsize=12)
            axs.set_ylabel('RMSD',fontsize=12) 
            axs.grid(True, lw=0.5, ls=':', c='gray')
            
            out_name = join(OUT_DIR,ivar+'_'+area_name+'_FCST_verification.png')
            plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=300)
            plt.close()   
            print('Figure FCST evaluation saved')
            
    return