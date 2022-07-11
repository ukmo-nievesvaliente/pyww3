#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 09:14:05 2022

    EXAMPLE SCRIPT - TESTING HMAX VERIFICATION USING FORRISTAL/RAYLEIGH DISTRIBUTION
    From 
        # Forristall constants
            alpha = 2.126
            beta = 8.42
        # Rayleigh constants
            alpha = 2.0
            beta  = 8.0
        p = 0.99
        
        From Forristall, G.Z. (1978) 
                h = ((np.log(1.0-p) * beta * -1.0) ** (1.0 / alpha) ) 
                h = h * hs / 4.0

@author: nvalient
"""

# import standard libraries
import datetime
import numpy as np
import matplotlib.pyplot
import sys
import os
import VerLoad as vl
import matplotlib.pyplot as plt
sys.path.insert(0,'/net/home/h01/nvalient/nvalient-python/')
import pyww3.plot.ver_plot_cmems as vpc

sys.path.insert(0,'/net/home/h01/nvalient/nvalient-python/MOLevel1-wavetools/')
import wavetools.seastates.seastates as sst

# get these fields from matchup files for the testing
# ---------------------------------------------------
#DIR_IN = '/data/users/nvalient/verification/GMD_paper/matchup/T1uk'
DIR_IN = '/data/users/nvalient/verification/test_hmax/matchup/ver_hmax-uk-T02'
DIR_OUT = '/scratch/nvalient/img_testing'
MOD = 'ver_hmax-uk'
MM_start = '20191201'#int(202001)
MM_end   = '20191231'#int(202112)
OBSTYPE = 'WAVENET'
VAR = ['Hs','Tp','T02','Dir','Spr','Hmax']
LOC_ID = ['6201002']
LOC_NAME = ['Start_bay']
LOC_TITLE = ['Start bay (English Channel) - 6201002']

ids, locns, paras, obsdata, modeldata, dates = vl.LoadMatchUpInSituWithTime(DIR_IN, MOD, OBSTYPE, 
                                                                     MM_start, MM_end, hsthresh=False)

print('[INFO] Parameters included are :'+str(paras))
# Loop to obtain obs and model in a particular location
for ll,locs in enumerate(LOC_ID):        
    #only values for idtype == wavenet and buoyB
    # Obs for location id 6200288
    mask = ids == locs            
    
    #vtarr   = vtarr1[mask,:]    
    latsarr = locns[mask,0]
    lonsarr = locns[mask,1]
    obsall = obsdata[mask,:] # [N location,N variables] N variables depend on the OBSTYPE
    modall = modeldata[mask,:] 
    tdate = dates[mask]
    
    hs_obs = obsall[:,0]
    hs_mod = modall[:,0]
    tp_obs = obsall[:,1]
    tp_mod = modall[:,1]
    t02_obs = obsall[:,2]
    t02_mod = modall[:,2]
    
    try:
        hmax_mod = modall[:,5]
    except:
        hmax_mod = None
    
    if not(hs_obs is None):
        modstd = np.nanstd(hs_mod)
        goodobs = np.where(np.abs(hs_obs-hs_mod)<5.0*modstd)
        rejectcount = len(hs_obs)-len(goodobs[0])
        print('[INFO] Rejected %d' %rejectcount + ' bad observations')
        hsobsforver = hs_obs[goodobs]
        hsmodforver = hs_mod[goodobs]
        tpobsforver = tp_obs[goodobs]
        tpmodforver = tp_mod[goodobs]
        t02obsforver = t02_obs[goodobs]
        t02modforver = t02_mod[goodobs]
        tvector = tdate[goodobs]
        if hmax_mod is not None:
            hmaxmodforver = hmax_mod[goodobs]
        
    # -------------------------------------------------  
    # Compute HMAX from obs
    p, hmaxFobst02 = sst.calcHmax(hsobsforver, ptype='dist', p=0.99, tp=t02obsforver, window=3600., expected=False, rayleigh=False)
    # if rayleigh=True; Rayleigh distribution is used
    p, hmaxRobst02 = sst.calcHmax(hsobsforver, ptype='dist', p=0.99, tp=t02obsforver, window=3600., expected=False, rayleigh=True)
    
    # Compute HMAX from obs
    p, hmaxFobs = sst.calcHmax(hsobsforver, ptype='dist', p=0.99, tp=tpobsforver, window=3600., expected=False, rayleigh=False)
    # if rayleigh=True; Rayleigh distribution is used
    p, hmaxRobs = sst.calcHmax(hsobsforver, ptype='dist', p=0.99, tp=tpobsforver, window=3600., expected=False, rayleigh=True)  
    
    if hmax_mod is not None:
        hmaxFmod = hmaxmodforver
        hmaxRmod = hmaxmodforver
        print('[INFO] Hmax in matchup file')
    else:
        # Compute HMAX from model bulk stats
        print('[INFO] Hmax computed using bulk stats')
        # if rayleigh=False; Forristal distribution is used
        p, hmaxFmod = sst.calcHmax(hsmodforver, ptype='dist', p=0.99, tp=tpmodforver, window=3600., expected=False, rayleigh=False)
        # if rayleigh=True; Rayleigh distribution is used
        p, hmaxRmod = sst.calcHmax(hsmodforver, ptype='dist', p=0.99, tp=tpmodforver, window=3600., expected=False, rayleigh=True)
        #--------------- t02
        p, hmaxFmodt02 = sst.calcHmax(hsmodforver, ptype='dist', p=0.99, tp=t02modforver, window=3600., expected=False, rayleigh=False)
        # if rayleigh=True; Rayleigh distribution is used
        p, hmaxRmodt02 = sst.calcHmax(hsmodforver, ptype='dist', p=0.99, tp=t02modforver, window=3600., expected=False, rayleigh=True)
        

    ptloc=[0.1,max(hmaxFobs)]    
    
    # Get model accuracy
    # -------------------------------------------------
    # Plotting routine
    
    # 1) HMAX from Model bulk stats - using T02
    # 1A Forristall
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatQQ( hmaxFmod, hmaxFobst02, axisres=1.0, hexbin=True, linfit=True, grid=True, showplt=False ) 
    ax.set_ylabel('HMAX observations [m]')
    ax.set_xlabel('HMAX model [m]')
    plt.title('Forristall distribution - T02')
                
    rect = vpc.TextStats(hmaxFmod, hmaxFobst02, units=None, linfit=True, errorstats=True, basicstats=False, forplot=True, ptloc=ptloc, font='small',
                                      tolerance=True, dirn=False)    
    ax.add_patch(rect)
    # save the plot
    out_name = os.path.join(DIR_OUT,MOD+'_withT02_'+LOC_ID[ll]+'_QQ_HMAX_Forristall_fromMatchup.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
    
    # 1B Rayleigh
    ptloc=[0.1,max(hmaxRobs)]
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatQQ( hmaxRmod, hmaxRobst02, axisres=1.0, hexbin=True, linfit=True, grid=True, showplt=False ) 
    ax.set_ylabel('HMAX observations [m]')
    ax.set_xlabel('HMAX model [m]')
    plt.title('Rayleigh distribution - T02')
                
    rect = vpc.TextStats(hmaxRmod, hmaxRobst02, units=None, linfit=True, errorstats=True, basicstats=False, forplot=True, ptloc=ptloc, font='small',
                                      tolerance=True, dirn=False)    
    ax.add_patch(rect)
    # save the plot
    out_name = os.path.join(DIR_OUT,MOD+'_withT02_'+LOC_ID[ll]+'_QQ_HMAX_Rayleigh_fromMatchup.png')
    #out_name = os.path.join(DIR_OUT,'T1UK_'+LOC_ID[ll]+'_QQ_HMAX_Rayleigh_fromMatchup.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
            
    print('-----------------------------------------------')
    print('Saved figures for ID: '+ids)
    
    
    # 2) HMAX from Model bulk stats - using Tp
    # 2A Forristall
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatQQ( hmaxFmod, hmaxFobs, axisres=1.0, hexbin=True, linfit=True, grid=True, showplt=False ) 
    ax.set_ylabel('HMAX observations [m]')
    ax.set_xlabel('HMAX model [m]')
    plt.title('Forristall distribution - Tp')
                
    rect = vpc.TextStats(hmaxFmod, hmaxFobs, units=None, linfit=True, errorstats=True, basicstats=False, forplot=True, ptloc=ptloc, font='small',
                                      tolerance=True, dirn=False)    
    ax.add_patch(rect)
    # save the plot
    out_name = os.path.join(DIR_OUT,MOD+'_withTp_'+LOC_ID[ll]+'_QQ_HMAX_Forristall_fromMatchup.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
    
    # 2B Rayleigh
    ptloc=[0.1,max(hmaxRobs)]
    fig = plt.figure()
    ax  = fig.add_subplot(111)  
    pl = vpc.PlotScatQQ( hmaxRmod, hmaxRobs, axisres=1.0, hexbin=True, linfit=True, grid=True, showplt=False ) 
    ax.set_ylabel('HMAX observations [m]')
    ax.set_xlabel('HMAX model [m]')
    plt.title('Rayleigh distribution - Tp')
                
    rect = vpc.TextStats(hmaxRmod, hmaxRobs, units=None, linfit=True, errorstats=True, basicstats=False, forplot=True, ptloc=ptloc, font='small',
                                      tolerance=True, dirn=False)    
    ax.add_patch(rect)
    # save the plot
    out_name = os.path.join(DIR_OUT,MOD+'_withTp_'+LOC_ID[ll]+'_QQ_HMAX_Rayleigh_fromMatchup.png')
    #out_name = os.path.join(DIR_OUT,'T1UK_'+LOC_ID[ll]+'_QQ_HMAX_Rayleigh_fromMatchup.png')
    plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1, dpi=150)
            
    print('-----------------------------------------------')
    print('Saved figures for ID: '+ids) 