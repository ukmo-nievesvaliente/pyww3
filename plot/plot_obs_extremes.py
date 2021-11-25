#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:14:05 2021

    SET OF FUNCTIONS TO PLOT OBSERVATIONS

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

def var_inObs(all_var=False):
    if all_var == True:
        VAR_PLT = ['hs','ws','wdir','tp','t02','dir','spr']            
    else:
        VAR_PLT = ['hs','ws','wdir','t02']
    
    return VAR_PLT

def get_dict():
    
    """
    Created on 30 September 15:17:23 2021
    
    Function to create dictionary with some of the variables and stats that can be read from the extremes.csv file   
        
    @author: nvalient
    """
    # ------------------ BIAS -----------------------
    VAR = {}
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
    VAR['bias_t02']['limits'] = [-3,3]
    VAR['bias_t02']['limits_diff'] = [-2,2]    
    
    # ------------------ RMSD -----------------------

    VAR['rmse_ws']={}
    VAR['rmse_ws']['short_n'] = 'U10 RMSD'
    VAR['rmse_ws']['colorbar'] = 'gist_stern_r'
    VAR['rmse_ws']['limits'] = [0,4]
    VAR['rmse_ws']['limits_diff'] = [-1,1]
    
    VAR['rmse_wdir']={}
    VAR['rmse_wdir']['short_n'] = 'U10 dir RMSD'
    VAR['rmse_wdir']['colorbar'] = 'gist_stern_r'
    VAR['rmse_wdir']['limits'] = [0,60]
    VAR['rmse_wdir']['limits_diff'] = [-1,1]
  
    VAR['rmse_hs']={}
    VAR['rmse_hs']['short_n'] = 'Hs RMSD'
    VAR['rmse_hs']['colorbar'] = 'gist_stern_r'
    VAR['rmse_hs']['limits'] = [0,1.5]
    VAR['rmse_hs']['limits_diff'] = [-1,1]

    VAR['rmse_tp']={}
    VAR['rmse_tp']['short_n'] = 'Tp RMSD'
    VAR['rmse_tp']['colorbar'] = 'gist_stern_r'
    VAR['rmse_tp']['limits'] = [0,2.5]
    VAR['rmse_tp']['limits_diff'] = [-1,1]
 
    VAR['rmse_t02']={}
    VAR['rmse_t02']['short_n'] = 'T02 RMSD'
    VAR['rmse_t02']['colorbar'] = 'gist_stern_r'
    VAR['rmse_t02']['limits'] = [0,2.5]
    VAR['rmse_t02']['limits_diff'] = [-1,1]   
    
    return VAR

def var_inObs4plot(var):
    if var == 'hs':
        VAR_4PLOT = ['bias_hs','rmse_hs']
    elif var == 'ws':
        VAR_4PLOT = ['bias_ws','rmse_ws']
    elif var == 'tp':
        VAR_4PLOT = ['bias_tp','rmse_tp']
    elif var == 't02':
        VAR_4PLOT = ['bias_t02','rmse_t02']
    elif var == 'wdir':
        VAR_4PLOT = ['bias_wdir','rmse_wdir']
    
    return VAR_4PLOT

def plot_obs_stats(out_dir,lon_stat,lat_stat,var,val_stat,run,Q):
    
    land_50 = fill_land_cfeature()
    VAR_4PLOT = var_inObs4plot(var)
    
    lon, lat, bat = get_contours_NWS_bat()
    VARs = get_dict()
    levels = [40,80,120,240]
    
    #print('[Debug] vars 4 plot are '+str(VAR_4PLOT))
    for ii,stat in enumerate(VAR_4PLOT):
        # Iterate over the different stats that can be plotted
        rr = np.array(val_stat[ii])
        
        fig2 = plt.figure()
        axes = fig2.add_subplot(111,projection=ccrs.PlateCarree())
        a = axes.contour(lon,lat,bat, levels, colors=(.4, .4, .4),linewidths=.25)
        axes.coastlines(resolution='50m', color='black', linewidth=1)
        axes.add_feature(land_50,zorder=1)
        axes.clabel(a, inline=True, fmt = '%3.0f', fontsize=6)
        e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*10,cmap=VARs[stat]['colorbar'],\
                         vmin=VARs[stat]['limits'][0],vmax=VARs[stat]['limits'][1],zorder=2)
        cbar = fig2.colorbar(e,extend='both')
        #cbar.set_label(VAR[stat]['short_n'])
        axes.set_xticks([-12,-6,0,6],crs=ccrs.PlateCarree())
        axes.set_yticks([48,54,60],crs=ccrs.PlateCarree())
        #axes.gridlines()
        axes.set_axisbelow(True)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axes.xaxis.set_major_formatter(lon_formatter)
        axes.yaxis.set_major_formatter(lat_formatter)
        axes.set_xlim([-13,9])
        axes.set_ylim([45,63])
        plt.title(VARs[stat]['short_n']+' for '+run+' - '+Q,fontsize=8)
        out_name2 = join(out_dir,run+'_'+stat+'_'+Q+'.png')
        print("Saving figure " +out_name2)
        plt.savefig(out_name2,bbox_inches="tight", pad_inches=0.1,dpi=300)
        plt.close()
    
    return

def plot_obs_stats_r(out_dir,lon_stat,lat_stat,var,val_stat,run,Q):
    
    land_50 = fill_land_cfeature()
    VAR_4PLOT = var_inObs4plot(var)
    
    lon, lat, bat = get_contours_NWS_bat()
    VAR = get_dict()
    levels = [40,80,120,240]
    
    for ii,stat in enumerate(VAR_4PLOT):
        # Iterate over the different stats that can be plotted
        rr = val_stat[ii]
        
        fig2 = plt.figure()
        axes = fig2.add_subplot(111,projection=ccrs.PlateCarree())
        a = axes.contour(lon,lat,bat, levels, colors=(.4, .4, .4),linewidths=.25)
        axes.coastlines(resolution='50m', color='black', linewidths=1)
        axes.add_feature(land_50, zorder=1)
        axes.clabel(a, inline=True, fmt = '%3.0f', fontsize=6)
        e = axes.scatter(lon_stat,lat_stat,c=rr,s=(np.ones((1,len(rr))))*10,cmap='seismic',\
                         vmin=VAR[stat]['limits_diff'][0],vmax=VAR[stat]['limits_diff'][1], zorder=2)
        cbar = fig2.colorbar(e,extend='both')
        #cbar.set_label(VAR[stat]['short_n'])
        axes.set_xticks([-12,-6,0,6],crs=ccrs.PlateCarree())
        axes.set_yticks([48,54,60],crs=ccrs.PlateCarree())
        axes.gridlines()
        axes.set_axisbelow(True)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axes.xaxis.set_major_formatter(lon_formatter)
        axes.yaxis.set_major_formatter(lat_formatter)
        axes.set_xlim([-13,9])
        axes.set_ylim([45,63])
        plt.title(VAR[stat]['short_n']+' for '+run+' relative to ctr'+' - '+Q,fontsize=8)
        out_name2 = join(out_dir,run+'_relative2ctrl_'+stat+'_'+Q+'.png')
        print("Saving figure " +out_name2)
        plt.savefig(out_name2,bbox_inches="tight", pad_inches=0.1,dpi=300)
        plt.close()
    
    return