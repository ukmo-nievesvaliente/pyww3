#! /usr/bin/env python

# plotting routines to be used for CMEMS wave verification
# functions inherited from MOLevel1-wavetools by A. Saulter See https://github.com/ukmo-ansaulter/MOLevel1-wavetools.git
# Additions and modifications for more functionality by N.G. Valiente   

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.stats import t, norm, nct
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import pyww3.general_funs.stats_fun as sf


def MyValStr( value, units=None ):
    """create a value string"""

    myvalstr = '%8.2f' %value
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


def RoundAxLims(data,rndval=False,pad=False):
    """set axis limits based on rounding value and input data"""

    datamax = np.max(data)
    datamin = np.min(data)

    if pad:
        datamax = datamax + pad
        datamin = datamin - pad

    # automatically estimate rounding value
    if not rndval:
        if datamax - datamin > 100.0:
            rndval = 100.0
        elif datamax - datamin > 50.0:
            rndval = 50.0
        elif datamax - datamin > 20.0:
            rndval = 20.0
        elif datamax - datamin > 10.0:
            rndval = 10.0
        elif datamax - datamin > 5.0:
            rndval = 5.0
        elif datamax - datamin > 2.0:
            rndval = 2.0
        elif datamax - datamin > 1.0:
            rndval = 1.0
        elif datamax - datamin > 0.5:
            rndval = 0.5
        elif datamax - datamin > 0.2:
            rndval = 0.2
        else: rndval = 0.1

    axmax = np.ceil(datamax / rndval) * rndval
    axmin = np.floor(datamin / rndval) * rndval

    lims = [axmin, axmax]

    return lims


def PlotScatQQ( xdata, ydata, axisres=1.0, qqplt=1, hexbin=None, linfit=False, grid=True, showplt=False ):
    """plot scatter data and qq data"""

    # set axis limits
    axlims = RoundAxLims( np.array([np.min(xdata),np.min(ydata),np.max(xdata),np.max(ydata)]), axisres)

    # calculate a least squares linear fit
    if linfit:
        m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
        fitvalmin = axlims[0] * m + c
        fitvalmax = axlims[1] * m + c

    # plot the data

    # 1:1 line
    plt.plot( axlims, axlims, 'c-' )

    # scatter or hexbin
    if hexbin is None:
        plt.scatter( xdata, ydata, s=2, color='grey' )
    else:
        plt.hexbin( xdata, ydata, gridsize=np.int(axlims[1]/hexbin), mincnt=1, cmap=plt.cm.Set2_r, edgecolor='white')
        cb = plt.colorbar()
        cb.set_label('Data Frequency')

    # qq data
    # calculate and plot qq data - use percentile factors of 10 associated with different markers
    if qqplt is not None:
        markers = ['o', 's', '^', '+', '*', '-']

        # calculate qq percentiles
        qqlog = np.floor( np.log10( len(xdata) ) ) - 1
        pcstart = 0.0
        for i in range( np.int(qqlog) ):
            xdata_dist = []
            ydata_dist = []
            pcstep = 1.0 / 10.0**i
            if qqplt == 2:
                for j in np.arange(0.0, 1/10.0**(i-1), pcstep):
                    xdata_dist.append( stats.scoreatpercentile(xdata,j+pcstep) )
                    ydata_dist.append( stats.scoreatpercentile(ydata,j+pcstep) )
            for j in np.arange(pcstart, 100.0-2.0*pcstep, pcstep):
                xdata_dist.append( stats.scoreatpercentile(xdata,j+pcstep) )
                ydata_dist.append( stats.scoreatpercentile(ydata,j+pcstep) )
            pcstart = j+pcstep

            # plot the data
            if i == 0:
                plt.scatter( xdata_dist, ydata_dist, color='black', marker=markers[i] , s=25, label='QQ data')
            else:
                plt.scatter( xdata_dist, ydata_dist, color='black', marker=markers[i] , s=25)


    # linear fit line
    if linfit:
        plt.plot(axlims,[fitvalmin,fitvalmax],'r--', linewidth=2, label='Linear Fit')

    # apply axis limits
    plt.xlim( axlims )
    plt.ylim( axlims )

    # legend
    plt.legend(loc='lower right', fontsize='small')

    if grid:
        plt.grid()

    if showplt:
        plt.show()

    return


def PlotScatter( xdata, ydata, axisres=1.0, hexbin=None, linfit=False, grid=True, showplt=False ):
    """plot scatter data"""

    # set axis limits
    axlims = RoundAxLims( np.array([np.min(xdata),np.min(ydata),np.max(xdata),np.max(ydata)]), axisres)

    # calculate a least squares linear fit
    if linfit:
        m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
        fitvalmin = axlims[0] * m + c
        fitvalmax = axlims[1] * m + c

    # plot the data
    # 1:1 line
    plt.plot( axlims, axlims, 'c-' )

    # scatter or hexbin
    if hexbin is None:
        plt.scatter( xdata, ydata, s=2, color='grey' )
    else:
        plt.hexbin( xdata, ydata, gridsize=hexbin, mincnt=1, cmap=plt.cm.Set2_r, edgecolor='white')
        cb = plt.colorbar()
        cb.set_label('Data Frequency')

    # linear fit line
    if linfit:
        plt.plot(axlims,[fitvalmin,fitvalmax],'r--', linewidth=2, label='Linear Fit')
    
    # # apply axis limits
    # plt.xlim( RoundAxLims( np.array([np.min(xdata),np.max(xdata)]), axisres) )
    # plt.ylim( RoundAxLims( np.array([np.min(ydata),np.max(ydata)]), axisres) )
    
    # apply axis limits
    plt.xlim( axlims )
    plt.ylim( axlims )

    # legend
    plt.legend(loc='lower right', fontsize='small')

    if grid:
        plt.grid(linewidth=0.5,zorder=1)

    if showplt:
        plt.show()

    return


def TextStats(xdata, ydata, units=None, linfit=True, errorstats=True, basicstats=False, forplot=False, ptloc=None, font=None,
              tolerance=False, dirn=False):
    """Generates text strings for standard metrics.
       Output strings are either to be used with scatter plot layout (forplot=True), or
       to be passed back for write to text file (forplot=False)
       Note that if forplot=True, the location of the text and the fontsize can be specified"""

    # calculate the stats (might need to change the mean for medians for dirn data??)
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

    if errorstats:
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
        if tolerance:
            tsind, trmse, ksind, krmse = get_tolerance(xdata,ydata,sind,rmse)

    if linfit:
        m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
            
    # plot the text, assumed in axis coordinates
    if forplot:
        if not ptloc:
            xpt = 0.1
            ypt = 1.00
            dy  = 0.05
            ymax= ypt+0.05
        else:
            xpt = ptloc[0]
            ypt = ptloc[1]
            if ypt < 5:
                ypt = 4.7
            if max(ydata)>=9 and max(ydata)<=16:
                dy = 0.36
                dy  = 0.4 #STB
            elif max(ydata)>=16 and max(ydata)<40:
                dy  = 0.58
            elif max(ydata)>40 and max(ydata)<250:
                dy  = 2.
            elif max(ydata)>250:
                dy  = 15.
            elif max(ydata) >= 6 and max(ydata) < 9:
                dy  = 0.22
            else:
                dy  = 0.18
            ymax= ypt+0.05
        if not font:
            fonts = 'medium'
        else:
            fonts = font
        
        if basicstats:    
            ypt = ypt - dy
            plt.text(xpt,ypt,'X-Y Statistics',fontsize=fonts)
            ypt = ypt - dy
            myvalstr = '%d' %len(xdata)
            plt.text(xpt,ypt,'No. data = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( xmean, units=units )
            plt.text(xpt,ypt,'Xmean = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( xstd, units=units )
            plt.text(xpt,ypt,'Xstdev = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( ymean, units=units )
            plt.text(xpt,ypt,'Ymean = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( ystd, units=units )
            plt.text(xpt,ypt,'Ystdev = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            plt.text(xpt,ypt,'----' )
        if errorstats:            
            ypt = ypt - dy
            plt.text(xpt,ypt,'X-Y Errors',fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( bias, units=units )
            plt.text(xpt,ypt,'Bias = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( rmse, units=units )
            plt.text(xpt,ypt,'RMSD = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( estd, units=units )
            plt.text(xpt,ypt,'StdE = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( sind )
            plt.text(xpt,ypt,'SI = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( syms )
            plt.text(xpt,ypt,'Sym. Slope = ' + myvalstr, fontsize=fonts )
        if linfit:
            ypt = ypt - dy
            plt.text(xpt,ypt,'----' )
            ypt = ypt - dy
            plt.text(xpt,ypt,'X-Y Linear Fit', fontsize=fonts)
            ypt = ypt - dy
            myvalstr = MyValStr( r_value )
            plt.text(xpt,ypt,'r = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( m )
            plt.text(xpt,ypt,'Slope = ' + myvalstr, fontsize=fonts )
            ypt = ypt - dy
            myvalstr = MyValStr( c )
            plt.text(xpt,ypt,'Offset = ' + myvalstr, fontsize=fonts )
        # Put stats in a box
        if max(xdata) <=5:
            xmax = max(xdata)/4 + 0.5
            xmax = max(xdata)/3 + 0.7 #STB
        else:
            xmax = max(xdata)/4
            xmax = max(xdata)/3 #STB
        rect = patches.Rectangle((xpt-0.05,ypt-0.05),xmax-0.05,ymax-ypt, edgecolor='black', alpha=0.5, facecolor='white',zorder=3 )
        
        return rect
    # otherwise pass data out as a string
    else:
        outstr = '%d' %len(xdata)
        outstr = outstr + ',' + MyValStr( xmean )
        outstr = outstr + ',' + MyValStr( xstd )
        outstr = outstr + ',' + MyValStr( ymean )
        outstr = outstr + ',' + MyValStr( ystd )
        if errorstats:
            outstr = outstr + ',' + MyValStr( bias )
            outstr = outstr + ',' + MyValStr( rmse )
            outstr = outstr + ',' + MyValStr( estd )
            outstr = outstr + ',' + MyValStr( sind )
            outstr = outstr + ',' + MyValStr( syms )
        if linfit:
            outstr = outstr + ',' + MyValStr( r_value )
            outstr = outstr + ',' + MyValStr( m )
            outstr = outstr + ',' + MyValStr( c )

        return outstr

    return

def ParaStrings():
    """Dictionary to define strings for the various variables we could be verifying.
       Defined strings are:
       - long name
       - acronym
       - CMEMS parameter name (guessed at where not actually in catalogue!)
       - units"""

    paralist = { 'hs':['Significant wave height','Hs','VHM0','m'],
                 'tp':['Wave peak period','Tp','VTPK','s'],
                 't02':['Mean zero-upcrossing period','T02','VTM02','s'],
                 't01':['Wave mean period','T01','VTM01','s'],
                 't0m1':['Wave energy period','Tm10','VTM10','s'],
                 'dir':['Wave direction','Dirn','VMDR','deg'],
                 'spr':['Wave directional spread','Spr','VSPR','deg'],
                 'ws':['Wind Speed at 10m asl','Ws','VWS10','m/s'],
                 'wdir':['Wind Direction (from) at 10m asl','Wdir','VWD10','deg'] }

    return paralist


def ParaRanges(parastr):
    """Returns default ranges for the various variables we could be verifying
       Ranges are given for mean, std, bias, si"""

    paralist = { 'hs':[[0.0,4.0],[0.0,3.0],[-0.3,0.3],[0.05,0.75]],
                 'tp':[[3.0,15.0],[0.0,5.0],[-2.0,2.0],[0.05,0.75]],
                 't02':[[2.0,12.0],[0.0,5.0],[-2.0,2.0],[0.05,0.75]],
                 't01':[[2.0,12.0],[0.0,5.0],[-2.0,2.0],[0.05,0.75]],
                 't0m1':[[2.0,12.0],[0.0,5.0],[-2.0,2.0],[0.05,0.75]],
                 'dir':[[0.0,360.0],[0.0,120.0],[-30.0,30.0],[0.05,0.75]],
                 'spr':[[5.0,45.0],[0.0,45.0],[-15.0,15.0],[0.05,0.75]],
                 'ws':[[2.0,12.0],[0.0,6.0],[-2.0,2.0],[0.05,0.75]],
                 'wdir':[[0.0,360.0],[0.0,120.0],[-15.0,15.0],[0.05,0.75]] }

    meanrng = paralist[parastr][0]
    stdrng  = paralist[parastr][1]
    biasrng = paralist[parastr][2]
    sirng   = paralist[parastr][3]

    return meanrng, stdrng, biasrng, sirng


def ConvertDirn(dirseries):
    """Converts directions from 0-360 convention to -180 to +180 relative to central value"""

    # use median as central direction value
    dircentre = np.median(dirseries)

    # calculate differentials and correct to within +/-180
    dirchk = dirseries - dircentre
    dirseries[dirchk > 180.0] = dircentre + dirchk[dirchk > 180.0] - 360.0
    dirseries[dirchk < -180.0] = dircentre + dirchk[dirchk < -180.0] + 360.0

    return dirseries


def SetVerMap(extent, meridians, parallels, resolution='50m', projection=ccrs.PlateCarree()):

    map = plt.axes(projection=projection)
    map.coastlines(resolution=resolution)
    map.set_xticks(meridians, crs=projection)
    map.set_yticks(parallels, crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    map.xaxis.set_major_formatter(lon_formatter)
    map.yaxis.set_major_formatter(lat_formatter)
    map.gridlines(projection, draw_labels=False, xlocs=meridians, ylocs=parallels)
    map.set_extent(extent, projection) 

    return map

def plot_ci_manual(t, s_err, n, x, x2, y2, ax=None):
    """Return an axes of confidence bands using a simple approach.
    
    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}
    
    References
    ----------
    .. [1] M. Duarte.  "Curve fitting," Jupyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb
    
    """
    if ax is None:
        ax = plt.gca()
    
    ci = t * s_err * np.sqrt(1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    ax.fill_between(x2, y2 + ci, y2 - ci, color="#b9cfe7", edgecolor="")

    return ax

def plot_ci_bootstrap(xs, ys, resid, nboot=500, ax=None):
    """Return an axes of confidence bands using a bootstrap approach.

    Notes
    -----
    The bootstrap approach iteratively resampling residuals.
    It plots `nboot` number of straight lines and outlines the shape of a band.
    The density of overlapping lines indicates improved confidence.

    Returns
    -------
    ax : axes
        - Cluster of lines
        - Upper and Lower bounds (high and low) (optional)  Note: sensitive to outliers

    References
    ----------
    .. [1] J. Stults. "Visualizing Confidence Intervals", Various Consequences.
       http://www.variousconsequences.com/2010/02/visualizing-confidence-intervals.html

    """ 
    if ax is None:
        ax = plt.gca()

    bootindex = sp.random.randint

    for _ in range(nboot):
        resamp_resid = resid[bootindex(0, len(resid) - 1, len(resid))]
        # Make coeffs of for polys
        pc = sp.polyfit(xs, ys + resamp_resid, 1)                   
        # Plot bootstrap cluster
        ax.plot(xs, sp.polyval(pc, xs), "b-", linewidth=2, alpha=3.0 / float(nboot))

    return ax

def PlotScatterCI( xdata, ydata, CI=None):
    """plot scatter data"""
     
    axisres=1
    # set axis limits
    axlims = RoundAxLims( np.array([np.min(xdata),np.min(ydata),np.max(xdata),np.max(ydata)]), axisres)

    # calculate a least squares linear fit
    m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
    fitvalmin = axlims[0] * m + c
    fitvalmax = axlims[1] * m + c

    # plot the data
    fig, ax = plt.subplots()
    # 1:1 line
    ax.plot( axlims, axlims, 'c-' )
    
    ax.scatter( xdata, ydata, s=2, color='grey' )
    
    # Call function to obtain CI and PI bands
    t, resid, s_err, y_model, x2, y2, pi = sf.get_comp_CIandPIbands(xdata,ydata)
    # Fit
    ax.plot(xdata, y_model, "-", color="0.1", linewidth=1.5, alpha=0.5, label="Fit")  
    
    # Plot confidence Interval (select one)
    if CI:
        plot_ci_bootstrap(xdata, ydata, resid, ax=ax)
    else:    
        plot_ci_manual(t, s_err, len(xdata), xdata, x2, y2, ax=ax)

    # Plot prediction interval
    ax.fill_between(x2, y2 + pi, y2 - pi, color="None", linestyle="--")
    ax.plot(x2, y2 - pi, "--", color="0.5", label="95% Prediction Limits")
    ax.plot(x2, y2 + pi, "--", color="0.5")

    return