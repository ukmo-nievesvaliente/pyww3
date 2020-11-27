#!/usr/bin/env python3
from scipy.signal import butter,filtfilt

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