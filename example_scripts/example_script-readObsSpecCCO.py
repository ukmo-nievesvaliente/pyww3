#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
sys.path.insert(0,'/home/h01/nvalient/nvalient-python/pyww3')
import  obs_funs.read_obs_spec as obspec

"""
Created on Thu Apr 13 15:12:54 2023

EXAMPLE SCRIPT TO READ AND PLOT SPECTRAL DATA FROM OBS - CCO (.SPT format)

@author: nvalient

"""
 
# read obs data
dirIn = '/data/users/nvalient/GMD_paper/data/spec_obs_Perranporth/'
filename = join(dirIn,'Perranporth}2019-08-10T11h45Z.spt')#'Perranporth}2019-08-10T05h45Z.spt')

# Call to local library to extract and get observed spectral data from CCO wave buoys (.SPT format)
xi, yi, overlapping_z, orig_frequencies, orig_norm_psds = obspec.get_spec_obs(filename)


# PLOT - 1D
plt.figure(figsize=(5,4))
ax=plt.subplot(1,1,1)
ax.plot(1/orig_frequencies,orig_norm_psds,'.-k',label='obs')
ax.set_ylabel('Energy Density (m$^2$/s)')
ax.set_xlabel('Period (s)')
plt.savefig('PPT_1Dspec_example.png')

# PLOT - 2D
# Directional spectrum
fig=plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location("N")
labels=[ 'N','','','W','','','S','','','E','' ,'']
ax.set_xticks(np.pi/180. * np.linspace(0,  360, 12, endpoint=False),labels=labels)
ax.set_ylabel('Frequency (Hz)')
# ax.set_xlabel('Direction (degrees)')
ax.yaxis.set_label_coords(0.35, 0.75)

# contour color levels:
levels=[0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0]

# ax.set_rlim([0,0.65])
# ax.set_rscale('symlog')
# ax.set_rticks([0.04,0.1,0.3,0.6])

min_index = 180
max_index = 360

# ax_im = ax.contourf(xi[min_index:max_index], yi, overlapping_z, levels=levels,cmap = 'bone_r', extend='min')#, vmin=0.1,vmax=200)

ax_im = ax.contourf(xi[min_index:max_index], np.log10(yi/0.02), overlapping_z, levels=levels, locator=ticker.LogLocator(),cmap = 'viridis')

# Get labels in log scale
rr=ax.yaxis.get_ticklabels()
# set default freq values, if none specified:
freq_labels = []
freq_pos = []
for label in ax.yaxis.get_ticklabels():
    next = 10**float(label.get_text())*0.02
    freq_pos.append(float(label.get_text()))
    freq_labels.append(str(round(next,3)))
    
# Replace scale with real freq labels
ax.set_yticks(freq_pos,labels=freq_labels)
ax.set_ylabel('Frequency (Hz)')
ax.yaxis.set_label_coords(0.35, 0.75)
cb=plt.colorbar(ax_im, orientation='horizontal')
cb.ax.set_xlabel('m$^2$/Hz/Rad')
cb.ax.xaxis.set_label_position('top')

plt.savefig('PPT_2Dspec_example.png')
