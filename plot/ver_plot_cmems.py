#! /usr/bin/env python

# plotting routines to be used for CMEMS wave verification

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
#import iris
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def MyValStr( value, units=None ):
    """create a value string"""

    myvalstr = '%8.3f' %value
    if units is not None:
        myvalstr = myvalstr + units
    myvalstr = myvalstr.strip(' ')

    return myvalstr


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

    # apply axis limits
    plt.xlim( RoundAxLims( np.array([np.min(xdata),np.max(xdata)]), axisres) )
    plt.ylim( RoundAxLims( np.array([np.min(ydata),np.max(ydata)]), axisres) )

    # legend
    plt.legend(loc='lower right', fontsize='small')

    if grid:
        plt.grid()

    if showplt:
        plt.show()

    return


def TextStats(xdata, ydata, units=None, linfit=True, errorstats=True, forplot=True, dirn=False):
    """Generates text strings for standard metrics.
       Output strings are either to be used with scatter plot layout (forplot=True), or
       to be passed back for write to text file (forplot=False)"""

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

    if linfit:
        m, c, r_value, p_value, stderr = stats.linregress(xdata, ydata)
    
    # plot the text, assumed in axis coordinates
    if forplot:
        ypt = 1.00
        dy  = 0.05
        ypt = ypt - dy
        plt.text(0.1,ypt,'X-Y Statistics')
        ypt = ypt - dy
        myvalstr = '%d' %len(xdata)
        plt.text(0.1,ypt,'No. data = ' + myvalstr, fontsize='medium' )
        ypt = ypt - dy
        myvalstr = MyValStr( xmean, units=units )
        plt.text(0.1,ypt,'Xmean = ' + myvalstr, fontsize='medium' )
        ypt = ypt - dy
        myvalstr = MyValStr( xstd, units=units )
        plt.text(0.1,ypt,'Xstdev = ' + myvalstr, fontsize='medium' )
        ypt = ypt - dy
        myvalstr = MyValStr( ymean, units=units )
        plt.text(0.1,ypt,'Ymean = ' + myvalstr, fontsize='medium' )
        ypt = ypt - dy
        myvalstr = MyValStr( ystd, units=units )
        plt.text(0.1,ypt,'Ystdev = ' + myvalstr, fontsize='medium' )
        if errorstats:
            ypt = ypt - dy
            plt.text(0.1,ypt,'----' )
            ypt = ypt - dy
            plt.text(0.1,ypt,'X-Y Errors' )
            ypt = ypt - dy
            myvalstr = MyValStr( bias, units=units )
            plt.text(0.1,ypt,'Bias = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( rmse, units=units )
            plt.text(0.1,ypt,'RMSD = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( estd, units=units )
            plt.text(0.1,ypt,'StdE = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( sind )
            plt.text(0.1,ypt,'SI = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( syms )
            plt.text(0.1,ypt,'Sym. Slope = ' + myvalstr, fontsize='medium' )
        if linfit:
            ypt = ypt - dy
            plt.text(0.1,ypt,'----' )
            ypt = ypt - dy
            plt.text(0.1,ypt,'X-Y Linear Fit' )
            ypt = ypt - dy
            myvalstr = MyValStr( r_value )
            plt.text(0.1,ypt,'R = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( m )
            plt.text(0.1,ypt,'Slope = ' + myvalstr, fontsize='medium' )
            ypt = ypt - dy
            myvalstr = MyValStr( c )
            plt.text(0.1,ypt,'Offset = ' + myvalstr, fontsize='medium' )

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
