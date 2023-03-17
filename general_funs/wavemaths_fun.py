#!/usr/bin/env python3

import numpy as np
import math
from scipy import optimize
import netCDF4 as nc4
from scipy.interpolate import griddata
from os.path import join
import pandas as pd

# ---------------------------------------------------------------------------
# Funtions for wave direction correction and return of x,y components 
# ---------------------------------------------------------------------------
def MdirCopernicus2zonal(Mdir):  
    """ Copernicus wave direction is Mean wave direction FROM (Mdir) wrt the N (we call it phi)
    
    In order to get the x,y components (for quiver) we need theta (wrt the zonal direction; i.e., E = 0 deg angle).
    Direction angles gamma are wrt True North: the angle wrt the zonal direction theta is 
                                               theta = 90.-gamma
    Direction angles theta gives where waves are going: the angles where waves are coming from phi will be gamma+180; therefore 
    gamma = phi - 180;
    Combining the two we have that the angle theta we want is
    theta = 90 - (phi - 180) """
    theta = 270. - Mdir
    
    return theta

def getXY_MdirCopernicus(theta):
    """ To use in quiver - Get meridional and zonal components from Mdir (WAV Copernicus products)
        Inputs: Theta, Mdir corrected wrt zonal direction + direction To (not from)
        Outputs: x_hs; zonal component for Hs direction
                 y_hs; meridional component for Hs direction """
    x_hs = 1. * np.cos((theta)*math.pi/180.)
    y_hs = 1. * np.sin((theta)*math.pi/180.)
    
    return x_hs, y_hs

# ---------------------------------------------------------------------------
# Function to get depth at in-situ locations 
# ---------------------------------------------------------------------------
def get_dpt(lon_i,lat_i):
    """ Interpolate depth at lon_i lat_i location; valid for UK regional domain"""
    # -------------------------------------------------------------------------
# Reading bathymetry
    bathy_file = '/data/users/nvalient/bathymetry/NWshelf_amm15_bathy.nc'
    print("Reading "+bathy_file)
    fbathy     = nc4.Dataset(bathy_file,"r")
    bat        = np.array(fbathy.variables["Bathymetry"])
    lat        = np.array(fbathy.variables["lat"])
    lon        = np.array(fbathy.variables["lon"])
    fbathy.close()

# ---------------------------------------------------------------------------
# Estimate depth
    points = np.array( (lon.flatten(), lat.flatten()) ).T
    values = bat.flatten()

    print('Interpolating obs location to nearest bathy points')
    depth_ID = griddata(points, values, (lon_i, lat_i), method='nearest')
    
    return depth_ID

# ---------------------------------------------------------------------------
# Wave basic funtions 
# ---------------------------------------------------------------------------
def WVLFenton( t, depth, return_deptype=False):
     """Calculate wavelength for intermediate wave depths
        according to Fenton & McKee 1990, 
        'On calculating the lengths of water waves'
 
        Expected inputs:
        t     - period; array
        depth - depth; scalar or array with appropriate dimensions
 
        Outputs:
        wvint - wavelength value
        deptype - water depth type, 0:shallow, 1:intermediate, 2:deep"""

     # set up the output arrays
     wvint   = np.zeros( np.shape(t) )
     deptype = np.ones( np.shape(t), dtype=int )

     # use missing data array for calcs to ensure the data is masked when dividing by small numbers
     tma = np.ma.masked_less( t, 0.5 )

     # calculate the shallow, deep and intermediate wavelength options
     wlshallow = np.sqrt(9.81 * depth) * tma
     wrshallow = depth / wlshallow
     wldeep    = 9.81 * tma**2.0 / (2.0 * np.pi)
     wrdeep    = depth / wldeep
     wlint     = wldeep * np.tanh( (2.0*np.pi*depth/wldeep)**0.75 )**(2.0/3.0)

     # put the correct option into the wvint and deptype arrays
     wvint[ np.ma.MaskedArray.nonzero(wlint) ] = wlint[ np.ma.MaskedArray.nonzero(wlint) ]
     wvint[ np.where(wrshallow < 0.05) ] = wlshallow[ np.where(wrshallow < 0.05) ]
     deptype[ np.where(wrshallow < 0.05) ] = 0
     wvint[ np.where(wrdeep > 0.5) ] = wldeep[ np.where(wrdeep > 0.5) ]
     deptype[ np.where(wrdeep > 0.5) ] = 2

     if return_deptype:
         return wvint, deptype
     else:
         return wvint

def calc_wave_age(depth,tp,ust,vst):
    '''
    Function to compute wave age or its inverse
    
    Parameters
    ----------
    depth : depth.
    tp  : peak period (most desirable would be peak period from wind sea partition)
    ust : u-component friction velocity
    vst : v-component friction velocity

    Returns
    -------
    INVWAGE : Inverse of wave age; U_*/Cp.
    WAGE : wave age; Cp/U_*
    '''
    UST = np.sqrt(ust**2+vst**2)
    LP  = WVLFenton(tp, depth)
    TP  = np.maximum(tp,0.001)
    CP  = LP / TP
    WAGE  = CP / UST 
    INVWAGE = UST / CP
    
    return INVWAGE, WAGE

def dimensionless_fetch(U10,tvals):
    '''
    Compute dimensionless fetch value

    Parameters
    ----------
    U10 : module for wind speed.
    tvals : time vector; e.g. resn = 5
            tvals = np.arange(resn,201*resn-0.1,resn).

    Returns
    -------
    tdim : dimensionless fetch value.
    '''
    tdim = 9.81 * 1000.0 * tvals / np.float(U10)**2.0
    
    return tdim

def dimensionless_growth(U10,HS):
    '''
    Compute dimensionless wave growth

    Parameters
    ----------
    U10 : module for wind speed.
    HS  : significant wave height from wave model.

    Returns
    -------
    wgrowth : dimensionless wave growth (hs/u10)
    '''
    wgrowth = 9.81 * HS / U10**2.0
    
    return wgrowth
 
def YoungVerhagen(ws,tvals,depth=None,dimensionless=True):
    '''
    Function to compute dimensionless growth (hs/ws) values using the theoretical curve of
    Young & Verhagen (1996)  
    '''

    if depth is not None:
        dcorr = np.tanh(0.343 * (9.81 * depth / np.float(ws)**2.0)**1.14)
        hs_yv = (np.float(ws)**2.0 * 0.24 / 9.81) * (dcorr * np.tanh(0.000414 * (9.81 * tvals * 1000. / np.float(ws)**2.0)**0.79 / dcorr))** 0.572
    else:
        hs_yv = (np.float(ws)**2.0 * 0.24 / 9.81) * np.tanh(0.000414 * (9.81 * tvals * 1000. / np.float(ws)**2.0)**0.79) ** 0.572

    if dimensionless:
        hs_yv = 9.81 * hs_yv / np.float(ws)**2.0

    return hs_yv
       
def growthDecay_coef(depth,tp,ust,vst):
    
    """
    Function to compute the growth/ decay coefficient function of wave age
    """
    
    INVWAGE, WAGE = calc_wave_age(depth,tp,ust,vst)
    
    Cb = np.ones(WAGE.shape)
    Cb[WAGE <= 20.] = 32.
    Cb[WAGE > 20.] = -30.
    #np.where(WAGE <= 20.,Cb,Cb*32.)
    #np.where(WAGE > 20.,Cb,Cb*(-30.))
    
    return Cb


def growthDecay_rate(depth,tp,k,ust,vst):
    """
    Function to compute growth/ decay rate = Cb*k*rhoa/rhow*U_*^2/CP
    (Belcher and Hunt 1993)
    Inputs:
        depth : depth
        tp  : peak period (most desirable would be peak period from wind sea partition)
        Cb = wave growth/decay coef., fc. of wave age obtained using growthDecay_coef
        k = wave number
        ust : u-component friction velocity
        vst : v-component friction velocity

    Returns
    -------
    GDrate = growth/ decay rate.

    """
    rhoa      = 1.225       # [kg/m3] air density
    rhow      = 1023        # [kg/m3] water density
    
    LP     = WVLFenton(tp, depth)
    CP     = LP / tp
    #Cb     = growthDecay_coef(depth,tp,ust,vst)
    Cb     = 30
    Ust    = np.sqrt(ust**2+vst**2)
    GDrate = Cb*(2*np.pi/LP)*(rhoa/rhow)*Ust**2/CP
    fp = 1/tp
    GDrate_nonD = GDrate/fp
    
    return GDrate, GDrate_nonD

def waveSteepness(hs,depth,tp):
    """
    Function to compute wave steepness = Hs/Lp
    Inputs:
        depth : depth
        tp  : peak period (most desirable would be peak period from wind sea partition)

    Returns
    -------
    WS = wave steepness.

    """
    
    LP     = WVLFenton(tp, depth)
    WS     = hs / LP
    
    return WS

# ----------------------------------------------------------------------
# Functions for Relation between Cd and U10
# ----------------------------------------------------------------------
def func_CdHwang(U10):
    '''
    Function to compte Cd following Hwang (2011); as per ST6 WW3 physics.  
    http://doi.org/10.1175/2010JTECHO812.1 and
    http://doi.org/10.1175/JTECH-D-11-00092.1
    '''
    cd_st6 =  (8.058+0.967*U10-0.016*U10**2)/1e4
    
    return cd_st6

def func_cd(xx):
    '''
    Compute the drag coefficient Cd as a function of 
    the sea surface roughness
    xx : sea surface roughness
    
    How to use: 
            z0_0 = first_guess_z0(U10,charn)
            z0   = z_0(func_z0, z0_0, fder_z0,charn,U10)
            cd.append(func_cd(z0))
    '''
    kk        = 0.4          # von Karman constant
    
    return kk**2. / (np.log(1. + 10./xx))**2.

def first_guess_z0(U10,charn):
    '''
    Function to determing the first guess for z0.
    From https://www.ecmwf.int/sites/default/files/elibrary/2010/9875-sea-surface-roughness-and-drag-coefficient-function-neutral-wind-speed.pdf
    '''
    g         = 9.81         # [m2/s] gravitational acceleration
    cva       = 1.54e-6      # kinematic viscosity of air times
    # guess for ustar
    ustar = U10 / 25.
    
    # guess for z0
    z0 = (charn*ustar**2 / g) + (cva / max(ustar,0.0001))

    return z0

def func_z0(x,charn,U10):
    '''
    Callable function for the sea surface roughness.

    In takes the following inputs:
       x   : the unknown sea surface roughness
       charn: the Charnock parameter
       u10  : the module of the wind speed at 10m

    The function results from combining:

    1) z0 = (charn*ustar**2 / g) + (cva / ustar)   
    
    2) ustar = (kk * u10) / ln(1+10/z0)
    
    How to use: 
            z0_0 = first_guess_z0(U10,charn)
            z0   = z_0(func_z0, z0_0, fder_z0,charn,U10)
            cd.append(func_cd(z0))
    '''
# PHYSICAL variables ---------------------------------------------------

    g         = 9.81         # [m2/s] gravitational acceleration
    kk        = 0.4          # von Karman constant
    cva       = 1.54e-6      # kinematic viscosity of air times
                         # gustiness constant for z0 (=0.11)
    A  = g * kk * U10
    B  = charn * kk**3 * U10**3
    D  = cva * g
    ln = np.log(1. + 10./x)
    #print x
    
    return A*x*ln**2 - D*ln**3 - B

def fder_z0(x,charn,U10):
    '''
    Callable function for the first derivative of the 
    sea surface roughness.

    In takes the following inputs:
       x or z0   : the unknown sea surface roughness
       charn: the Charnock parameter
       u10  : the module of the wind speed at 10m
       
    How to use: 
            z0_0 = first_guess_z0(U10,charn)
            z0   = z_0(func_z0, z0_0, fder_z0,charn,U10)
            cd.append(func_cd(z0))
    '''
    g         = 9.81         # [m2/s] gravitational acceleration
    kk        = 0.4          # von Karman constant
    cva       = 1.54e-6      # kinematic viscosity of air times
    
    A   = g * kk * U10
    D   = cva * g
    ln  = np.log(1. + 10./x)
    ra1 = x*(1.+ 10./x)
    ra2 = x*x*(1.+ 10./x)

    return A*ln**2 - (20./ra1)*A*ln - (30./ra2)*D*ln**2

def z_0(f_z0,z0_0,f1_z0,charn,U10):
    '''
    This function return the sea surface roughness
    approximated with the Newton-Raphson iterative 
    method. 

    It takes the following inputs:
        f_z0   : from the function
        z0_0   : from the funtion
        f1_z0  : method
        charn  : the Charnock parameter
        u10    : the module of the wind speed at 10m
       
    How to use: 
            z0_0 = first_guess_z0(U10,charn)
            z0   = z_0(func_z0, z0_0, fder_z0,charn,U10)
            cd.append(func_cd(z0))    
    '''

    if f1_z0: # Newton-Raphson method
       root = optimize.newton(f_z0, z0_0, fprime=f1_z0, maxiter=100, args=(charn,U10))
    else:     # Secant method
       root = optimize.newton(f_z0, z0_0, fprime=None , maxiter=100, args=(charn,U10))

    return root

def func_Donelan():
    '''
    Function to plot the values in fig. 6 of Donelan 2018 "On the Decrease of the Oceanic Drag
    Coefficient in High Winds" see https://doi.org/10.1002/2017JC013394
    '''
    
    dataIn='/data/users/nvalient/tabu_stress/Donelan2018'
    datasets=['Donelan_2018_modelled_50km_fetch.txt','Donelan_2018_modelled_230km_fetch.txt','obs_edson_et_al_2013_Jarosz_el_al_2007.txt']
    C = np.empty([3,29])
    C[:] = np.nan
    U = np.copy(C)
    # extract the data
    for i,data in enumerate(datasets):
        filename=join(dataIn,data)
        df=pd.read_csv(filename,
                       header = None, engine='python')
        df.columns = ['U10','CD']
        u10=df['U10'][:]
        U[i,0:len(u10)]=u10.to_numpy(dtype=float) 
        C[i,0:len(u10)]=df['CD'][:].to_numpy(dtype=float)
        
    return U,C