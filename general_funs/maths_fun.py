#!/usr/bin/env python3

import numpy as np
from scipy.signal import butter,filtfilt
import scipy.spatial as sp
import scipy.interpolate as interpolate
from collections import Counter
from scipy.stats import binned_statistic_2d

# Function to bin/counts data for pcolor
def get_counts_binnned_data(VAR1,VAR2,BINS):
    """
    FUNCTION TO COMPUTE DENSITY/COUNTS FROM TWO VARIABLES
    
    Parameters
    ----------
    VAR1: np array that can be masked; can be 2D
    VAR2: np array that can be masked; can be 2D
    BINS: number of bins that you want to use to bin the data
    """
    
    if len(VAR1.shape) >= 2:
        VAR1 = VAR1.ravel(order='C')
        VAR2 = VAR2.ravel(order='C')
        
    if  np.ma.is_masked(VAR1) == 'True':
        VAR2 = VAR2[~VAR2.mask]
        VAR1 = VAR1[~VAR1.mask]
    
    # Obtain counts
    stats_counts = binned_statistic_2d(x=VAR2, y=VAR1, values=VAR1, statistic='count', bins=BINS)
    
    z_sigma, x_edges, y_edges = stats_counts[0], stats_counts[1], stats_counts[2]
    z_rot = np.rot90(z_sigma)  # rotate and flip to properly display...
    z_rot_flip = np.flipud(z_rot)
    # -----------------------------------------------------------------------
    # Estimate density function and probabiblity 
    VAR_density = VAR2[~VAR2.mask].data
    c = Counter(VAR_density)
    total = sum(c.values())
    probability_VAR = {k:v/total for k,v in c.items()}
    VAR_densityList = np.ndarray.tolist(VAR_density)
    
    probabilityCounter = []
    for i in range(len(VAR_densityList)):
        probabilityCounter.append(probability_VAR.get(VAR_densityList[i],0))

    probabilityCounter = np.array(probabilityCounter)*10000.
    
    return x_edges, y_edges, z_rot_flip

# moving average filter
def mov_av(N,mylist):
    """
    FUNCTION TO COMPUTE MOVING AVERAGE 
    
    Parameters
    ----------
    N : number of elements for running average.
    mylist : list of elements to apply the filter to.

    Returns
    -------
    moving_aves : list of elements after the application of
            the running average.

    """
    cumsum, moving_aves = [0], []

    for i, x in enumerate(mylist, 1):
        cumsum.append(cumsum[i-1] + x)
        if i>=N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            #can do stuff with moving_ave here
            moving_aves.append(moving_ave)
        
       # out   = np.ones(len(mylist))*np.nan
    return moving_aves

def butter_lowpass_filter(data, cutoff, T, fs, order):
    
    """
    # Filter requirements.
    T         # Sample Period
    fs        # sample rate, Hz
    cutoff    # desired cutoff frequency of the filter, Hz ,       
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
    """
   
    nyq = 0.5 * fs  # Nyquist Frequency

        
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    
    return y

def eoshift(array,shift,boundary,dim):

    a = np.roll(array, shift, axis=dim)
    if dim == 1:
       b = np.reshape(boundary,(boundary.shape[0],1))
    else:
       b = boundary

    if shift > 0:
       if dim == 0:
          a[0:shift,:] = b

       if dim == 1:
          a[:,0:shift] = b

    elif  shift < 0:
       if dim == 0:
          a[shift:,:] = b

       if dim == 1:
          a[:,shift:] = b

    return a

def seaoverland(c):
    """
    c = 2D field with np.nan land values  
    """
    #print '! SEAOVERLAND !'
    a = np.copy(c)
    b = np.copy(c)

    [nj, ni] = a.shape

    mat8 = eoshift(a   ,  1, boundary = a[:,0]      , dim = 1)
    mat1 = eoshift(mat8,  1, boundary = mat8[0,:]   , dim = 0)
    mat2 = eoshift(a   ,  1, boundary = a[0,:]      , dim = 0)
    mat4 = eoshift(a   , -1, boundary = a[:,ni-1]   , dim = 1)
    mat3 = eoshift(mat4,  1, boundary = mat4[0,:]   , dim = 0)
    mat5 = eoshift(mat4, -1, boundary = mat4[nj-1,:], dim = 0)
    mat6 = eoshift(a   , -1, boundary = a[nj-1,:]   , dim = 0)
    mat7 = eoshift(mat8, -1, boundary = mat8[nj-1,:], dim = 0)

    S = np.array([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])

    mask_logic = np.isnan(a) # TRUE for land point

    try:
       SS = np.nanmean(S, axis=0)
    except AttributeError: # added this fix to make it compatible with numpy < 1.8
       from scipy.stats import nanmean
       SS = nanmean(S, axis=0)

    a[mask_logic] = SS[mask_logic]
    a[~mask_logic] = b[~mask_logic]

    return a

