#!/usr/bin/env python3

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import general_funs.wavemaths_fun as wmf
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

def get_dict():
    
    """
    Created on 30 September 15:17:23 2021
    
    Function to create dictionary with some of the variables that can be read from a model output .nc
        The function deals with 1D and 2D variables    
        
    @author: nvalient
    """
    
    VAR = {}
    VAR['hs']={}
    VAR['hs']['short_n'] = 'hs'
    VAR['hs']['cbarn'] = '$H_s\ [m]$'
    VAR['hs']['colorbar'] = 'nipy_spectral'
    VAR['hs']['limits'] = [0,12]
    VAR['hs']['limits_diff'] = [-1,1]
    
    VAR['hmaxe']={}
    VAR['hmaxe']['short_n'] = 'hmaxe'
    VAR['hmaxe']['cbarn'] = '$H_{max}\ [m]$'
    VAR['hmaxe']['colorbar'] = 'magma'
    VAR['hmaxe']['limits'] = [0,16]
    VAR['hmaxe']['limits_diff'] = [-2,2]
    
    VAR['sdc']={}
    VAR['sdc']['short_n'] = 'sdc'
    VAR['sdc']['cbarn'] = '$C_d$'
    VAR['sdc']['colorbar'] = 'pink_r'
    VAR['sdc']['limits'] = [0,0.0032]
    VAR['sdc']['limits_diff'] = [-0.001,0.001]
    
    VAR['twf']={}
    VAR['twf']['short_n'] = 'twf'
    VAR['twf']['cbarn'] =  '$tauoc\ [W\ m^{-2}]$'
    VAR['twf']['colorbar'] = 'jet'
    VAR['twf']['limits'] = [0,0.]
    VAR['twf']['limits_diff'] = [-0.2,0.2]
    
    VAR['cha']={}
    VAR['cha']['short_n'] = 'cha'
    VAR['cha']['cbarn'] =  '$\\alpha$'
    VAR['cha']['colorbar'] = 'hot'
    VAR['cha']['limits'] = [0,0.035]
    VAR['cha']['limits_diff'] = [-0.011,0.011]
    
    VAR['fp']={}
    VAR['fp']['short_n'] = 'fp'
    VAR['fp']['cbarn'] =  '$f_p\ [s^{-1}]$'
    VAR['fp']['colorbar'] = 'gist_stern_r'
    VAR['fp']['limits'] = [0,0.05]
    VAR['fp']['limits_diff'] = [-0.01,0.01]
    
    VAR['t01']={}
    VAR['t01']['short_n'] = 't01'
    VAR['t01']['cbarn'] =  '$T_{01}\ [s]$'
    VAR['t01']['colorbar'] = 'gist_stern_r'
    VAR['t01']['limits'] = [0,18]
    VAR['t01']['limits_diff'] = [-2,2]
    
    VAR['t02']={}
    VAR['t02']['short_n'] = 't01'
    VAR['t02']['cbarn'] =  '$T_{01}\ [s]$'
    VAR['t02']['colorbar'] = 'gist_stern_r'
    VAR['t02']['limits'] = [0,18]
    VAR['t02']['limits_diff'] = [-2,2]
    
    VAR['u10']={}
    VAR['u10']['short_n'] = ['uwnd','vwnd']
    VAR['u10']['colorbar'] = plt.cm.magma
    VAR['u10']['cbarn'] =  '$U_{10}\ [ms^{-1}]$'
    VAR['u10']['limits'] = [0,28]
    VAR['u10']['limits_diff'] = [-4,4]
    
    VAR['taw']={}
    VAR['taw']['short_n'] = ['utaw','vtaw']
    VAR['taw']['colorbar'] = 'gist_heat'
    VAR['taw']['cbarn'] =  '$\\tau_{aw}\ [Nm^{-2}]$'
    VAR['taw']['limits'] = [0,3]
    VAR['taw']['limits_diff'] = [-0.4,0.4]
    VAR['taw']['density'] = [1000.]
    
    VAR['ust']={}
    VAR['ust']['short_n'] = ['uust','vust']
    VAR['ust']['colorbar'] = 'gist_rainbow'
    VAR['ust']['cbarn'] =  '$u_{*}\ [ms^{-1}]$'
    VAR['ust']['limits'] = [0,2]
    VAR['ust']['limits_diff'] = [-0.2,0.2]
    
    VAR['tauaww3']={}  # This is the variable with three dimensions instead of 2
    VAR['tauaww3']['short_n'] = ['uust','vust','rhoa']
    VAR['tauaww3']['colorbar'] = 'jet'
    VAR['tauaww3']['cbarn'] =  '$\\tau_{a}\ [Nm^{-2}]$'
    VAR['tauaww3']['limits'] = [0,3]
    VAR['tauaww3']['limits_diff'] = [-0.2,0.2]
    
    VAR['taua']={}
    VAR['taua']['short_n'] = ['utaua','vtaua']
    VAR['taua']['colorbar'] = 'gist_heat'#'jet'
    VAR['taua']['cbarn'] =  '$\\tau_{a}\ [Nm^{-2}]$'
    VAR['taua']['limits'] = [0,3]
    VAR['taua']['limits_diff'] = [-0.2,0.2]
    
    VAR['uss']={}
    VAR['uss']['short_n'] = ['uuss','vuss']
    VAR['uss']['colorbar'] = 'inferno' # or gnuplot
    VAR['uss']['cbarn'] =  '$U_{ss}\ [ms^{-1}]$'
    VAR['uss']['limits'] = [0,2]
    VAR['uss']['limits_diff'] = [-0.2,0.2]
    
    VAR['two']={}
    VAR['two']['short_n'] = ['utwo','vtwo']
    VAR['two']['colorbar'] = 'PiYG'
    VAR['two']['cbarn'] =  '$\\tau_{wo}\ [Nm^{-2}]$'
    VAR['two']['limits'] = [0,4]
    VAR['two']['limits_diff'] = [-0.4,0.4]
    VAR['two']['density'] = [1000.]
    
    VAR['Uc']={}
    VAR['Uc']['short_n'] = ['ucur','vcur']
    VAR['Uc']['colorbar'] = 'hot'
    VAR['Uc']['cbarn'] =  '$U_{cur}\ [ms^{-1}]$'
    VAR['Uc']['limits'] = [0,1.2]
    VAR['Uc']['limits_diff'] = [-0.4,0.4]
    
    VAR['rhoa']={}
    VAR['rhoa']['short_n'] = 'rhoa'
    VAR['rhoa']['colorbar'] = 'viridis'
    VAR['rhoa']['cbarn'] =  '$\\rho_a [kgm^{3}]$'
    VAR['rhoa']['limits'] = [1.1,1.27]
    VAR['rhoa']['limits_diff'] = [-0.1,0.1]
    
    VAR['wave_age']={}
    VAR['wave_age']['short_n'] =  ['uust','vust']
    VAR['wave_age']['colorbar'] = 'Set1'
    VAR['wave_age']['cbarn'] =  '$u_*\ /\ c_p$'
    VAR['wave_age']['limits'] = [0,0.1] # For inverse wave age
    VAR['wave_age']['limits_diff'] = [-0.04,0.04] # For inverse wave age
    
    return VAR    

def get_limits_plots(matrix):
    """ Get min and max values of a variable in order to get a neat plot
    Input: matrix with variable values 
    """
    dup = []
    for k in matrix:
        for i in k:
            dup.append(i)
    
    return max(dup), min(dup)

def get_snapshots(filer, ST, var,  dimension, out_name, title_ini):
    
    """
    Created on 30 September 15:17:23 2021
    
    QUICK PLOT/ snapshots of model output over a day
        The function deals with 1D and 2D variables    
        
    @author: nvalient
    """
    
    print('Plotting snapshots')
    print('Using file '+filer)
    
    land_50 = fill_land_cfeature()
    
    VAR   = get_dict()
    d     = nc.Dataset(filer)    
    
    if dimension == '1D':        
        hs   = d.variables[var]
        
    elif dimension == '2D':
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        hs    = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        if var == 'taw':
            #print('Using rho water to compute Taw: '+str(VAR[var]['density'][:]))
            hs = hs*VAR[var]['density'][:]
        elif var == 'wave_age':
            dp1 = d.variables['dpt']
            fp1 = d.variables['fp']
            tp1 = 1./fp1[:,:,:]
            hs,wage1 = wmf.calc_wave_age(dp1[:,:,:],tp1[:,:,:],u11[:,:,:],u12[:,:,:])
    
    elif dimension == '3D': # valid for the computation of the atmospheric stress from tau=rho*Ufric**2
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        Ufric = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        if 'rhoa' in d.variables:
            rho   = d.variables[VAR[var]['short_n'][2]][:]
        else:
            rho = 1.225 # cte atmospheric density as per WW3 code
        hs    = rho*Ufric**2
        

    units_storm = d.variables['time'].units
    date_storm = np.array(d.variables['time'])
    t_storm    = nc.num2date(date_storm[:], units_storm,'gregorian')
    lat        = d.variables['latitude'][:]
    lon        = d.variables['longitude'][:]
 
    for i in range(len(t_storm)):#range(0,4):
                   
        # PLOT 
        fig2 = plt.figure(figsize=(10, 5))
        axes = fig2.add_subplot(111,projection=ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5))
        # remove a margin around the data
        axes.set_xmargin(0)
        axes.set_ymargin(0)
        cmap = VAR[var]['colorbar']
        pc = axes.pcolormesh(lon, lat, hs[i,:,:], cmap = cmap, vmin = VAR[var]['limits'][0], vmax = VAR[var]['limits'][1])
        cbar = fig2.colorbar(pc)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(VAR[var]['cbarn'],size=12)
        axes.coastlines(resolution='50m', color='black', linewidth=1)
        axes.add_feature(land_50)
        if VAR[var]['short_n'][0] == 'uwnd':
            U10         = u11[i,:,:]
            V10         = u12[i,:,:]
            s = 35
            QV = axes.quiver(lon[::s], lat[::s], U10[::s,::s], V10[::s,::s],
                         units='xy',
                         angles='xy',
                         scale=190,
                         scale_units='inches',
                         color='black') 
            plt.quiverkey(QV, 0.66, 0.2, 20, "20 $ms^{-1}$", labelpos = "S", coordinates='figure')
        
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axes.xaxis.set_major_formatter(lon_formatter)
        axes.yaxis.set_major_formatter(lat_formatter)
        
        date_time = t_storm[i]
        date_title = date_time.strftime("%d/%m/%Y, %H:%M:%S")
        plt.title(title_ini+' '+date_title)        

        out_name_end = join(out_name+'_T'+str(i)+'.png')
        plt.savefig(out_name_end,bbox_inches="tight", pad_inches=0.1, dpi=150)
        print ('Saving '+var+' png ' +'snapshot' +' time='+str(i))
        plt.close("all")
    
    return

def get_snapshots_diff(filer1, filer2, ST, var,  dimension, out_name, title_ini):
    
    """
    Created on 30 September 15:17:23 2021
    
    QUICK PLOT/ snapshots of the difference between two models over a day
        The function deals with 1D and 2D variables    
        
    @author: nvalient
    """
    
    land_50 = fill_land_cfeature()
    
    VAR   = get_dict()
    d     = nc.Dataset(filer1)    
    d2    = nc.Dataset(filer2)
    
    print('Plotting snapshots')
    print('Difference between file '+filer1)
    print('and')
    print('file '+filer2)
    
    if dimension == '1D':    
        if var in d2.variables:
            hs1   = d.variables[var][:,:,:]
            hs2   = d2.variables[var][:,:,:]
            
        else:
            print('[WARNING] Does variable exist?')
            hs1   = d.variables[var][:,:,:]
            if var == 'rhoa':
                print('[WARNING] Using cte for density')
                hs2   = 1.225
            else:
                print('[ERROR] variable not found!')
                exit()
        hs = hs1-hs2
        
    elif dimension == '2D':
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        U1 = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        u21   = d2.variables[VAR[var]['short_n'][0]]
        u22   = d2.variables[VAR[var]['short_n'][1]]
        U2 = np.sqrt(u21[:,:,:]**2+u22[:,:,:]**2)
        
        if var == 'taw':
            #print('Using rho water to compute Taw: '+str(VAR[var]['density'][:]))
            hs = hs*VAR[var]['density'][:]
        elif var == 'wave_age':
            dp1 = d.variables['dpt']
            fp1 = d.variables['fp']
            tp1 = 1./fp1[:,:,:]
            U1,wage1 = wmf.calc_wave_age(dp1[:,:,:],tp1[:,:,:],u11[:,:,:],u12[:,:,:])
            dp2 = d2.variables['dpt']
            fp2 = d2.variables['fp']
            tp2 = 1./fp2[:,:,:]
            U2,wage2 = wmf.calc_wave_age(dp2[:,:,:],tp2[:,:,:],u21[:,:,:],u22[:,:,:])
        
        hs = U1-U2
        
    elif dimension == '3D': # valid for the computation of the atmospheric stress from tau=rho*Ufric**2
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        Ufric = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        u21   = d2.variables[VAR[var]['short_n'][0]]
        u22   = d2.variables[VAR[var]['short_n'][1]]
        Ufric2 = np.sqrt(u21[:,:,:]**2+u22[:,:,:]**2)
        if 'rhoa' in d2.variables:
            rho   = d.variables[VAR[var]['short_n'][2]][:]
            rho2  = d2.variables[VAR[var]['short_n'][2]][:]
        else:
            rho   = d.variables[VAR[var]['short_n'][2]][:]
            # cte atmospheric density as per WW3 code
            rho2= 1.225
        hs    = rho*Ufric**2-rho2*Ufric2**2

    units_storm = d.variables['time'].units
    date_storm = np.array(d.variables['time'])
    t_storm    = nc.num2date(date_storm[:], units_storm,'gregorian')
    lat        = d.variables['latitude'][:]
    lon        = d.variables['longitude'][:]
 
    for i in range(len(t_storm)):
                   
        # PLOT 
        fig2 = plt.figure(figsize=(10, 5))
        axes = fig2.add_subplot(111,projection=ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5))
        # remove a margin around the data
        axes.set_xmargin(0)
        axes.set_ymargin(0)
        cmap = 'seismic'
        pc = axes.pcolormesh(lon, lat, hs[i,:,:], cmap = cmap, vmin = VAR[var]['limits_diff'][0], vmax = VAR[var]['limits_diff'][1])
        cbar = fig2.colorbar(pc)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(VAR[var]['cbarn'],size=12)
        axes.coastlines(resolution='50m', color='black', linewidth=1)
        axes.add_feature(land_50)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axes.xaxis.set_major_formatter(lon_formatter)
        axes.yaxis.set_major_formatter(lat_formatter)
        
        date_time = t_storm[i]
        date_title = date_time.strftime("%d/%m/%Y, %H:%M:%S")
        plt.title(title_ini+' '+date_title)        

        out_name_end = join(out_name+'_'+str(i)+'.png')
        plt.savefig(out_name_end,bbox_inches="tight", pad_inches=0.1, dpi=150)
        print ('Saving '+var+' png ' +'snapshot_diff time='+str(i))
        plt.close("all")
    
    return

def get_mean_diff(filer1, filer2, ST, var,  dimension, out_name, title_ini):
    
    """
    Created on 11 October 15:17:23 2021
    
    QUICK PLOT/ mean difference between two models over a day
        The function deals with 1D and 2D variables    
        
    @author: nvalient
    """
    
    print('Plotting mean difference')
    print('Difference between file '+filer1)
    print('and')
    print('file '+filer2)
    
    land_50 = fill_land_cfeature() 
    
    VAR   = get_dict()
    d     = nc.Dataset(filer1)    
    d2    = nc.Dataset(filer2)
    
    if dimension == '1D':     
        if var in d2.variables:
            hs1   = d.variables[var]
            hs2   = d2.variables[var]
            hs = np.nanmean(hs1[:,:,:]-hs2[:,:,:],axis=0)
        else:
            print('[WARNING] Does variable exist?')
            hs1   = np.nanmean(d.variables[var][:,:,:],axis=0)
            if var == 'rhoa':
                print('[WARNING] Using cte for density')
                hs2   = 1.225
                hs = hs1-hs2
            else:
                print('[ERROR] variable not found!')
                exit()
        
    elif dimension == '2D':
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        U1 = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        u21   = d2.variables[VAR[var]['short_n'][0]]
        u22   = d2.variables[VAR[var]['short_n'][1]]
        U2 = np.sqrt(u21[:,:,:]**2+u22[:,:,:]**2)
        hs = np.nanmean(U1[:,:,:]-U2[:,:,:],axis=0)
        if var == 'taw':
            #print('Using rho water to compute Taw: '+str(VAR[var]['density'][:]))
            hs = hs*VAR[var]['density'][:]
        
    elif dimension == '3D': # valid for the computation of the atmospheric stress from tau=rho*Ufric**2
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        Ufric = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
        u21   = d2.variables[VAR[var]['short_n'][0]]
        u22   = d2.variables[VAR[var]['short_n'][1]]
        Ufric2 = np.sqrt(u21[:,:,:]**2+u22[:,:,:]**2)
        if 'rhoa' in d2.variables:
            rho   = d.variables[VAR[var]['short_n'][2]][:]
            rho2  = d2.variables[VAR[var]['short_n'][2]][:]
        else:
            rho   = d.variables[VAR[var]['short_n'][2]][:]
            # cte atmospheric density as per WW3 code
            rho2= 1.225
        U1 = np.nanmean(rho*Ufric**2,axis=0)
        U2 = np.nanmean(rho2*Ufric2**2,axis=0)
        hs = U1-U2
        
    #units_storm = d.variables['time'].units
    #date_storm = np.array(d.variables['time'])
    #t_storm    = nc.num2date(date_storm[:], units_storm,'gregorian')
    lat        = d.variables['latitude'][:]
    lon        = d.variables['longitude'][:]
 
                     
    # PLOT 
    fig2 = plt.figure(figsize=(10, 5))
    axes = fig2.add_subplot(111,projection=ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5))
    # remove a margin around the data
    axes.set_xmargin(0)
    axes.set_ymargin(0)
    cmap = 'seismic'
    pc = axes.pcolormesh(lon, lat, hs, cmap = cmap, vmin = VAR[var]['limits_diff'][0], vmax = VAR[var]['limits_diff'][1])
    cbar = fig2.colorbar(pc)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(VAR[var]['cbarn'],size=12)
    axes.coastlines(resolution='50m', color='black', linewidth=1)
    axes.add_feature(land_50)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    axes.xaxis.set_major_formatter(lon_formatter)
    axes.yaxis.set_major_formatter(lat_formatter)
        
    plt.title(title_ini)        

    out_name_end = join(out_name+'mean_difference.png')
    plt.savefig(out_name_end,bbox_inches="tight", pad_inches=0.1, dpi=150)
    print ('Saving '+var+' png ' +'mean_diff')
    plt.close("all")
    
    return

def get_mean(filer, ST, var,  dimension, out_name, title_ini):
    
    """
    Created on 30 September 15:17:23 2021
    
    QUICK PLOT/ average values of model output over a day
        The function deals with 1D and 2D variables    
        
    @author: nvalient
    """
    
    print('Plotting mean')
    print('Using file '+filer)
    
    land_50 = fill_land_cfeature()
    
    VAR   = get_dict()
    d     = nc.Dataset(filer)    
    
    if dimension == '1D':        
        hs   = d.variables[var]
        
    elif dimension == '2D':
        u11   = d.variables[VAR[var]['short_n'][0]]
        u12   = d.variables[VAR[var]['short_n'][1]]
        hs = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)

    lat        = d.variables['latitude'][:]
    lon        = d.variables['longitude'][:]
                   
    # PLOT 
    fig2 = plt.figure(figsize=(10, 5))
    axes = fig2.add_subplot(111,projection=ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5))
    # remove a margin around the data
    axes.set_xmargin(0)
    axes.set_ymargin(0)
    cmap = VAR[var]['colorbar']
    pc = axes.pcolormesh(lon, lat, np.nanmean(hs[:,:,:],axis=0), cmap = cmap, vmin = VAR[var]['limits'][0], vmax = VAR[var]['limits'][1])
    cbar = fig2.colorbar(pc)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(VAR[var]['cbarn'],size=12)
    axes.coastlines(resolution='50m', color='black', linewidth=1)
    axes.add_feature(land_50)
    if VAR[var]['short_n'][0] == 'uwnd':
        U10         = np.nanmean(u11[:,:,:],axis=0)
        V10         = np.nanmean(u12[:,:,:],axis=0)
        s = 35
        QV = axes.quiver(lon[::s], lat[::s], U10[::s,::s], V10[::s,::s],
                     units='xy',
                     angles='xy',
                     scale=190,
                     scale_units='inches',
                     color='black') 
        plt.quiverkey(QV, 0.66, 0.2, 20, "20 $ms^{-1}$", labelpos = "S", coordinates='figure')
        
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    axes.xaxis.set_major_formatter(lon_formatter)
    axes.yaxis.set_major_formatter(lat_formatter)
  

    out_name_end = join(out_name+'_mean.png')
    plt.savefig(out_name_end,bbox_inches="tight", pad_inches=0.1, dpi=150)
    print ('Saving '+var+' png ' +'snapshot mean')
    plt.close("all")
    
    return


# def get_snapshots_diff_2subplots(filer1, filer2, ST, subplot, VARS, dimensions, out_name, title_ini):
    
#     """
#     Created on 30 September 15:17:23 2021
    
#     QUICK PLOT/ snapshots of the difference between two models over a day
#         This function can deal with subplots. For this, the different variables and the respective dimensions should be included 
#         as a list. 
#         The function deals with 1D and 2D variables    
        
#     @author: nvalient
#     """
    
#     land_50 = fill_land_cfeature()
    
#     VAR   = get_dict()
#     d     = nc.Dataset(filer1)    
#     d2    = nc.Dataset(filer2)
    
#     units_storm = d.variables['time'].units
#     date_storm = np.array(d.variables['time'])
#     t_storm    = nc.num2date(date_storm[:], units_storm,'gregorian')
#     lat        = d.variables['latitude'][:]
#     lon        = d.variables['longitude'][:]
        
#     for dd in range(len(subplot)):
#         dimension = dimensions[dd]
#         var = VARS[dd]
        
#         if dimension == '1D':
#             hs1   = d.variables[var]
#             hs2   = d2.variables[var]
#             hs = hs1[:,:,:]-hs2[:,:,:]
        
#         elif dimension == '2D':
#             u11   = d.variables[VAR[var]['short_n'][0]]
#             u12   = d.variables[VAR[var]['short_n'][1]]
#             U1 = np.sqrt(u11[:,:,:]**2+u12[:,:,:]**2)
#             u21   = d2.variables[VAR[var]['short_n'][0]]
#             u22   = d2.variables[VAR[var]['short_n'][1]]
#             U2 = np.sqrt(u21[:,:,:]**2+u22[:,:,:]**2)
#             hs = U1-U2
            
#         if dd == 1:
#             hsa = hs
#         else:
#             hsa = np.vstack((hsa,hs)

#     # Start plotting routine
#     for i in range(len(t_storm)):
           
#         # PLOT 
#         if subplot == 0:
#             fig, axes = plt.subplots()
            
#         elif subplot != 0:
#             if subplot == 2:
#                 fig, axes = plt.subplots(1,2)                
#             elif subplot == 3:
#                 fig, axes = plt.subplots(1,3)                
#             elif subplot == 4:
#                 fig, axes = plt.subplots(2,2)
#             elif subplot == 6:
#                 fig, axes = plt.subplots(2,3)                
#             elif subplot == 9:
#                 fig, axes = plt.subplots(3,3)

#         axes = fig2.add_subplot(subplot,1dd,projection=ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5))
#         # remove a margin around the data
#         axes.set_xmargin(0)
#         axes.set_ymargin(0)
#             cmap = 'seismic'
#             pc = axes.pcolormesh(lon, lat, hs[i,:,:], cmap = cmap, vmin = VAR[var]['limits_diff'][0], vmax = VAR[var]['limits_diff'][1])
#             cbar = fig2.colorbar(pc)
#             cbar.ax.tick_params(labelsize=12)
#             cbar.set_label(VAR[var]['cbarn'],size=12)
#             axes.coastlines(resolution='50m', color='black', linewidth=1)
#             axes.add_feature(land_50)
#             lon_formatter = LongitudeFormatter(zero_direction_label=True)
#             lat_formatter = LatitudeFormatter()
#             axes.xaxis.set_major_formatter(lon_formatter)
#             axes.yaxis.set_major_formatter(lat_formatter)
            
#             date_time = t_storm[i]
#             date_title = date_time.strftime("%d/%m/%Y, %H:%M:%S")
#             plt.title(title_ini+' '+date_title)        
    
#             out_name_end = join(out_name+'_'+str(i)+'.png')
#             plt.savefig(out_name_end,bbox_inches="tight", pad_inches=0.1, dpi=150)
#             print ('Saving '+var+' png ' +'snapshot_diff time='+str(i))
#             plt.close("all")
    
#     return
