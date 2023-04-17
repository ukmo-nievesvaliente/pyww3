#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import specwave.specio as spcio
#import specwave.specplot as splt
#import matplotlib.pyplot as plt
#from os.path import isfile, join
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from scipy import pi
from scipy.special import erf

"""
Created on Thu Apr 13 12:14:14 2023

WORK IN PROGRESS

Call as: 
    orig_norm_psds, orig_frequencies,  orig_direction, orig_spread, orig_skew = read_spt_files(filename)
    # filename must include the path 
    # E.g. filename = join(dirIn,'Perranporth}2019-08-10T05h45Z.spt')
@author: nvalient
"""

total_degrees = 540
min_degrees = -360
max_degrees = 720

def pdf(norm_psd, x):
    return (100*norm_psd)/np.sqrt(2*pi) * np.exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/np.sqrt(2))) / 2

def skew(x, e=0, w=1, a=0, norm_psd = 1):
    t = (x-e) / w
    return 2 / w * pdf(norm_psd,t) * cdf(a*t)
    # You can of course use the scipy.stats.norm versions
    # return 2 * norm.pdf(t) * norm.cdf(a*t)

def distribution(direction, spread, norm_psd, skewness):

    x = np.linspace(min_degrees, max_degrees, total_degrees) 

    p = skew(x, direction, spread, skewness, norm_psd)
    return p


def read_spt_files(spectra_file):
    
    # Read files from CCO
    # file format is:
    #         Spectrum file (64 x 6 array)
    #         f, S(f) / Smax, Dir(f), Spr(f), Skew(f), Kurt(f) 
    #         f: wave frequency [Hz]
    #         S(f) / Smax: relative psd (power spectral density) [-]
    #         Dir(f): wave direction [°]
    #         Spr(f): directional spread [°]
    #         Skew(f): skewness of the directional distribution [-]
    #         Kurt(f): kurtosis of the directional distribution [-] 
    
    max_psd = float(open(spectra_file).readlines()[3])
    print('max psd is '+str(max_psd))
    spt = np.genfromtxt(spectra_file, delimiter=',', skip_header=12)
    orig_frequencies = spt[:,0]
    orig_norm_psds = spt[:,1] * max_psd
    orig_direction = spt[:,2]
    orig_spread = spt[:,3]
    orig_skew = spt[:,4]
    return orig_norm_psds, orig_frequencies, orig_direction, orig_spread, orig_skew


def spec_reader(filename):
    # Read spectral file from obs
    df=pd.read_csv(filename,
           delimiter=',',header=12,engine='python')
    return df

def get_spec_obs(filename):
    # filename must include the path 
    # E.g. filename = join(dirIn,'Perranporth}2019-08-10T05h45Z.spt')
    orig_norm_psds, orig_frequencies,  orig_direction, orig_spread, orig_skew = read_spt_files(filename)
    
    #y is frequency 
    yi = np.linspace(min(orig_frequencies), max(orig_frequencies),123)
    #x is direction
    xi = np.array(np.radians(np.linspace(min_degrees, max_degrees, 540)))
    X,Y = np.meshgrid(xi,yi)
    all_spreads = np.array([[],[],[]])
    for index, frequency in enumerate(orig_frequencies):
        norm_psd = orig_norm_psds[index]
        direction = orig_direction[index]
        spread = orig_spread[index]
        skewness = orig_skew[index]
        spread_distribution = np.array(distribution(direction, spread, norm_psd, skewness))
        spread_plus_dir = np.array([spread_distribution, xi, np.linspace(frequency,frequency,total_degrees)])
        all_spreads = np.concatenate((all_spreads, spread_plus_dir),axis=1)    
    # Grid spectrum
    Z = griddata((all_spreads[1], all_spreads[2]), all_spreads[0], (X, Y), method='linear') # griddata((x,y),z,(xi,yi))
    min_index = 180
    max_index = 360
    subset_Z = Z[:,np.arange(min_index,max_index)]
    Z_less_than_0 = Z[:,np.arange(0,180)]
    z_over_360 = Z[:,np.arange(360,540)]
    overlapping_z = subset_Z + Z_less_than_0 + z_over_360
    overlapping_z= np.ma.masked_where( overlapping_z <= 0.1, overlapping_z)
    
    return xi, yi, overlapping_z, orig_frequencies, orig_norm_psds

