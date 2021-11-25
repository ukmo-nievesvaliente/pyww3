#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.stats import t, norm, nct
#import iris
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

"""
Created on Thu Nov 25 11:11:17 2021

SET OF FUNCTIONS TO COMPUTE STANDARD MODEL PERFORMANCE STATS  

@author: nvalient
"""
def MyValStr( value, units=None ):
    """create a value string"""

    myvalstr = '%8.3f' %value
    if units is not None:
        myvalstr = myvalstr + units
    myvalstr = myvalstr.strip(' ')

    return myvalstr

def get_tolerance(xdata,ydata,SI,rmse):
    """function to obtain the confidence interval and the tolerance interval for scatter index and RMSE between
    observations and model data"""
    
    # Two-sided inverse Students t-distribution
    # p - probability, df - degrees of freedom            
    tinv = lambda p, df: abs(t.ppf(p/2, df))
    ts = tinv(0.05, len(xdata)-2) # 95% confidence
            
    # confidence intervals; Average t*Stdev*(1/sqrt(n))
    #plt.plot(axlims,[fitvalmin+ts*fitvalmin,fitvalmax+ts*fitvalmax],'k..', linewidth=1, label='95% confidence')
    print('[Information] The 95% confidence interval (scatter index) is = '+str(SI*ts))
    print('[Information] The 95% confidence interval (RMSD) is = '+str(ts*rmse))
    # prediction intervals; t*StDev*(sqrt(1+(1/n)))
           
    # Tolerance intervals; Average k*StDev
    # sample size
    n=len(xdata)
    # Percentile for the TI to estimate
    p=0.99
    # confidence level
    g = 0.95
    # (100*p)th percentile of the standard normal distribution
    zp = norm.ppf(p)
    # gth quantile of a non-central t distribution
    # with n-1 degrees of freedom and non-centrality parameter np.sqrt(n)*zp
    tt = nct.ppf(g, df=n-1., nc=np.sqrt(n)*zp)
    # k factor from Young et al paper
    k = tt / np.sqrt(n)

    print('[Information] The 99% tolerance interval for a 95% confidence (scatter index) is = '+str(k*SI))
    print('[Information] The 99% tolerance interval for a 95% confidence (RMSD) is = '+str(k*rmse))
    
    return SI*ts, ts*rmse, k*SI, k*rmse

def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b) 

def get_comp_CIandPIbands(xdata,ydata):
    """"Return the confidence intervals of a 1D polynomial."""
    
    p, cov = np.polyfit(xdata, ydata, 1, cov=True)             # parameters and covariance from of the fit of 1-D polynom.
    y_model = equation(p, xdata)                               # model using the fit parameters; NOTE: parameters here are coefficients
    
    # Statistics
    n = xdata.size                                             # number of observations
    m = p.size                                                 # number of parameters
    dof = n - m                                                # degrees of freedom
    t = stats.t.ppf(0.975, n - m)                              # used for CI and PI bands
    
    # Estimates of Error in Data/Model
    resid = ydata - y_model                           
    chi2 = np.sum((resid / y_model)**2)                        # chi-squared; estimates error in data
    chi2_red = chi2 / dof                                      # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid**2) / dof)                    # standard deviation of the error
    
    x2 = np.linspace(np.min(xdata), np.max(xdata), 100)
    y2 = equation(p, x2)
    
    # Prediction Interval
    pi = t * s_err * np.sqrt(1 + 1/n + (x2 - np.mean(xdata))**2 / np.sum((xdata - np.mean(xdata))**2))

    return t, resid, s_err, y_model, x2, y2, pi


def TextStats(xdata, ydata, dirn=False, Wprint=True):
    """Generates text strings for standard metrics.
       Output strings are either to be printed in the linux console, or
       to be passed back for write to text file 
       """

    # calculate the stats (might need to change the mean for medians for dirn data??)
    m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
    xmean = np.mean(xdata)
    ymean = np.mean(ydata)
    if not dirn:
        xstd = np.std(xdata)
        ystd = np.std(ydata)
    else:
        xtmp = xdata - xmean
        xtmp[xtmp > 180.0] = xtmp[xtmp > 180.0] - 360.0
        xtmp[xtmp < -180.0] = xtmp[xtmp < -180.0] + 360.0
        xstd = np.std(xtmp)
        ytmp = ydata - ymean
        ytmp[ytmp > 180.0] = ytmp[ytmp > 180.0] - 360.0
        ytmp[ytmp < -180.0] = ytmp[ytmp < -180.0] + 360.0
        ystd = np.std(ytmp)

    errors = xdata - ydata
    if dirn:
        errors[errors > 180.0] = errors[errors > 180.0] - 360.0
        errors[errors < -180.0] = errors[errors < -180.0] + 360.0
    bias = np.mean(errors)
    rmse = np.sqrt( np.mean( (errors)**2.0 ) )
    estd = np.std( (errors) )
    # note scatter index defn relative to observed std,
    # i.e. 1.0 when si same as for observed mean as predictor
    sind = estd / ystd          # assume observed data will be y value
    syms = (xstd / ystd)**2.0   # assume observed data will be y value
    
    tsind, trmse, ksind, krmse = get_tolerance(xdata,ydata,sind,rmse)
    
    outstr = '%d' %len(xdata)
    outstr = outstr + ',' + MyValStr( xmean )
    outstr = outstr + ',' + MyValStr( xstd )
    outstr = outstr + ',' + MyValStr( ymean )
    outstr = outstr + ',' + MyValStr( ystd )

    outstr = outstr + ',' + MyValStr( bias )
    outstr = outstr + ',' + MyValStr( rmse )
    outstr = outstr + ',' + MyValStr( estd )
    outstr = outstr + ',' + MyValStr( sind )
    outstr = outstr + ',' + MyValStr( syms )

    outstr = outstr + ',' + MyValStr( r_value )
    outstr = outstr + ',' + MyValStr( m )
    outstr = outstr + ',' + MyValStr( c )
    
    outstr = outstr + ',' + MyValStr( tsind )
    outstr = outstr + ',' + MyValStr( trmse)
    outstr = outstr + ',' + MyValStr( ksind )
    outstr = outstr + ',' + MyValStr( krmse)
    
    if Wprint:
        print('---  X-Y Statistics ---')
        print('No. data = ' + str(len(xdata)))
        print('Xmean = '+str(xmean) )
        print('Xstdev = '+str(xstd))
        print('Ymean = '+str(ymean) )
        print('Ystdev = '+str(ystd))
        
        print('---  X-Y Errors ---')
        print('Bias = '+str(bias) )
        print('RMSD = '+str(rmse))
        print('stdE = '+str(estd) )
        print('SI = '+str(sind))
        print('Sym. Slope= '+str(syms) )
        
        print('---  X-Y Linear Fit ---')
        print('R = '+str(r_value) )
        print('Slope = '+str(m))
        print('Offset = '+str(c) )
        
        print('---  Confidence Intervals ---')
        print('CI-SI = '+str(tsind))
        print('TI-SI = '+str(ksind))
        print('CI-RMSD = '+str(trmse))
        print('TI-RMSD = '+str(krmse))
        
    return outstr
            
