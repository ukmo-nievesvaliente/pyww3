#!/usr/bin/env python3

import numpy as np
from scipy.signal import butter,filtfilt
import scipy.spatial as sp
import scipy.interpolate as interpolate


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

