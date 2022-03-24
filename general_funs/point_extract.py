import os
import numpy as np
from scipy.spatial import KDTree
import netCDF4 as nc
import iris
import cf_units
import waveutil.stash as wstash
from collections import defaultdict
from datetime_util import date2num, num2date, cf_time_unit
import cartopy.crs as ccrs


# Variable scale factors for netCDF file (variables using SHORT datatype)
NCVAR_SCALE_FACTOR = defaultdict(lambda: 0.1,
    hs=0.002,
    wdir=0.1,
    mdir=0.1,
    depth=0.1
)

def get_site_indices(lsm, outfile, buoyfiles=['gbl_intercomparison_sites.txt','wavenet_sites.txt', 'uk_ports.txt']):

    if not isinstance(lsm, iris.cube.Cube):
        # assume is lsm is a filename and try loading cube:
        print("[INFO] Trying to load LSM cube from filename: %s" % str(lsm))
        lsm = iris.load_cube(str(lsm), wstash.LSM)

    print("[INFO] Building KDTree of sea points")
    print("[INFO]   - Getting model sea points")
    mlon = lsm.coord(axis='x').points
    mlon0 = mlon[0]
    mlonn = mlon[-1]
    mlat = lsm.coord(axis='y').points
    cellsize = np.sqrt( (mlon[1] - mlon[0])**2 + (mlat[1] - mlat[0])**2)
    print("cellsize: %.4f" % cellsize)
    mmlon, mmlat = np.meshgrid(mlon, mlat)
    isea = np.where(lsm.data != 0)



    print("[INFO]   - Building tree")
    kdtree = KDTree(np.c_[mmlon[isea], mmlat[isea]])

    print("[INFO] Loading observation locations")
    lat = []
    lon = []
    name = []
    buoyid = []
    period_type = []

    for buoyfile in buoyfiles:
    #for buoyfile in ['uk_ports.txt']:
        with open(buoyfile, 'r') as fid:
            for line in fid:
                if line.startswith('#'):
                    continue
                toks = line.split()
                buoyid.append(toks[0])
                lat.append(float(toks[1]))
                lon.append(float(toks[2]))
                period_type.append(int(toks[3]))
                name.append(' '.join(toks[4:]))
    ##--

    lat = np.array(lat)
    lon = np.array(lon)
    buoyid = np.array(buoyid)
    period_type = np.array(period_type)
    name = np.array(name)

    # if model is rotated pole, then rotate obs:
    cs = lsm.coord(axis='x').coord_system
    if isinstance(cs, iris.coord_systems.RotatedGeogCS):
        print("[INFO] Rotating obs to model grid pole.")
        lon, lat = iris.analysis.cartography.rotate_pole(lon, lat,
                cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

    # Aling longitudes to model grid range
    lon[lon < mlon0] += 360.
    lon[lon > mlonn] -= 360.

    ## Do KD Search:
    ll = np.c_[lon,lat]
    print("[INFO] Finding nearest neighbours.")
    dist, idx = kdtree.query(ll, distance_upper_bound=cellsize)

    m = np.where(dist != np.inf)
    print("These ports removed:")
    m2 = dist == np.inf
    print(buoyid[m2])

    buoyid = buoyid[m]
    period_type = period_type[m]
    name = name[m]
    dist = dist[m]
    idx = idx[m]
    lon = lon[m]
    lat = lat[m]



    # get indicies into flat seapoint array...
    midx = np.ravel_multi_index(isea, lsm.shape)[idx]

    # ... and the unravel into ix and iy values:
    midx = np.unravel_index(midx, lsm.shape)

    ## Write nearest neighbour points to output file
    with open(outfile, 'w') as fid:
        fid.write("#%9s %9s %9s %9s %9s %8s %8s\n" % (
            "BuoyId","ObLat","ObLon","ModLat","ModLon","ModIX","ModIY"))

        for i in np.arange(np.alen(idx)):
            fid.write("%10s %9.3f %9.3f %9.3f %9.3f %8d %8d\n" % (
                    buoyid[i], lat[i], lon[i], mlat[midx[0][i]], mlon[midx[1][i]],
                    midx[0][i], midx[1][i]))
    ## --

    return midx,lon,lat

def extract_point_timeseries_nc(infiles, ind, obs_id, ncout, obs_lat=None,
                                obs_lon=None, variables=None):
    """ Extracts point data from the netCDF file(s) specifcied in infiles.
        `ind` should be an array of pre-calculated model indices. This will
        be a 2D array (ilat, ilon) for regular grids or a 1D array (seapoint,)
        for SMC grids.

        `variables` is a list of netCDF variable names to to extract from the
        input files. Defaults to ['hs','tm','tz','mdir','wdir','wspd']
    """

    if variables is None:
        variables = ['hs', 'tm', 'tz', 'mdir', 'wdir', 'wspd']

    if isinstance(infiles, (str, bytearray)):
        infiles = [infiles]

    npnts = len(ind[0])
    assert(len(obs_id) == npnts)
    if obs_lon is not None:
        assert(len(obs_lon) == npnts)
    if obs_lat is not None:
        assert(len(obs_lat) == npnts)

    tunits = cf_time_unit()

    with nc.Dataset(ncout, mode='w', clobber=True, format="NETCDF4") as dout:
        ## setup output file
        dout.createDimension('time', None)
        dout.createDimension('point', size=npnts)

        v = dout.createVariable('time', np.double, dimensions=('time'))
        v.standard_name = 'time'
        v.long_name = 'time'
        v.units = tunits
        v.axis = 'T'

        v = dout.createVariable('point', np.str, dimensions=('point'))
        v[:] = np.array(obs_id)

        if obs_lat is not None:
            v = dout.createVariable('obs_lat', np.float32, dimensions=('point'))
            v.long_name = 'observation_latitude'
            v.standard_name = 'latitude'
            v.units = 'degrees'
            v[:] = obs_lat

        if obs_lon is not None:
            v = dout.createVariable('obs_lon', np.float32, dimensions=('point'))
            v.long_name = 'observation_longitude'
            v.standard_name = 'longitude'
            v.units = 'degrees'
            v[:] = obs_lon
            # TODO: align longitudes

        v = dout.createVariable('mod_lat', np.float32, dimensions=('point'))
        v.long_name = 'mod_latitude'
        v.standard_name = 'latitude'
        v.units = 'degrees'

        v = dout.createVariable('mod_lon', np.float32, dimensions=('point'))
        v.long_name = 'mod_longitude'
        v.standard_name = 'longitude'
        v.units = 'degrees'

        for varname in variables:
            # depth is currently time invariant - just store single value
            dims = ('point',) if varname is 'depth' else ('time', 'point')

            if varname in ['hs', 'tm', 'tz', 'mdir', 'wdir', 'wspd']:
               # values stored as short integers with scaling factor
               # for efficiency:
               v = dout.createVariable(varname, np.short, dimensions=dims)
               v.scale_factor = NCVAR_SCALE_FACTOR[varname]
               # TODO: CF meta-data
            else:
               v              = dout.createVariable(varname, np.float, dimensions=dims)
               v.scale_factor = 1.0

        # Loop over files, extracting point data:
        n = 0
        for ifile, ncfn in enumerate(infiles):
            with nc.Dataset(ncfn, mode='r') as d:
                if ifile == 0:
                    # populate model lat/lons
                    if len(ind) > 1:
                        mlat = d.variables['latitude'][ind[0]]
                        mlon = d.variables['longitude'][ind[1]]
                    else:
                        mlat = d.variables['latitude'][ind[0]]
                        mlon = d.variables['longitude'][ind[0]]

                    # ToDo - derotate, if needed:
                    if 'rotated_pole' in d.variables:
                        plon = d.variables['rotated_pole'].grid_north_pole_longitude
                        plat = d.variables['rotated_pole'].grid_north_pole_latitude
                        stdll = ccrs.PlateCarree().transform_points(
                                ccrs.RotatedPole(pole_latitude=plat, pole_longitude=plon),
                                mlon, mlat)
                        mlon = stdll[:,0]
                        mlat = stdll[:,1]
                    ##--

                    dout.variables['mod_lat'][:] = mlat
                    dout.variables['mod_lon'][:] = mlon
                #-- ifile == 0

                # add model time
                t = d.variables['time']
                if t.units != tunits:
                    # convert time units, if required:
                    t = nc.num2date(t[:], t.units)
                    t = date2num(t)
                else:
                    t = t[:]
                ntimes = len(t)
                dout.variables['time'][n:n+ntimes] = np.rint(t)  # to nearest second

                def _extract(v, ind, itime):
                    if len(ind) > 1:
                        # Can't do fancy indexing on NC variables:
                        return v[itime][ind]
                    # Can extract SMC points with orthogonal NC indexing
                    return v[itime,ind[0]]


                # extract point data for variables:
                for varname in variables:
                    print(f"  {varname}")

                    if varname != 'wspd' and varname not in d.variables:
                        print(f"[WARN] Model file does not contain variable {varname}")
                        continue

                    # special handling of depth (just take first time):
                    if varname == 'depth':
                        dout.variables[varname][:] = _extract(
                            d.variables[varname], ind, 0)
                        continue

                    # for all other variables, extract for each time:
                    for itime in range(ntimes):
                        #print(f"    - {itime}")

                        if varname == 'wspd' and 'wspd' not in d.variables:
                            # special case for wind speed: derive from uwnd and vwnd
                            uwnd = _extract(d.variables['uwnd'], ind, itime)
                            vwnd = _extract(d.variables['vwnd'], ind, itime)
                            dout.variables[varname][n+itime,:] = np.sqrt(
                                uwnd**2 + vwnd**2)
                        else:
                            dout.variables[varname][n+itime,:] = _extract(
                                d.variables[varname], ind, itime)
                        #
                    #--itime
                #--varname
            ##-- d
            n += ntimes
        ##-- ncfn
    #-- dout


def extract_point_timeseries(infiles, pointsfile, fnout):

    import subprocess
    import shlex

    if isinstance(infiles, str):
        infiles = [infiles]

    flds = {
        'hs': wstash.HSIG,
        'tm': wstash.TM01,
        'tz': wstash.TM02,
        'mdir': wstash.MDIR,
        'wspd': wstash.WSPD,
        'wdir': wstash.WDIR
    }

    scale_factor = defaultdict(lambda: 0.01)
    scale_factor['hs'] = 0.002
    scale_factor['wdir'] = 0.1
    scale_factor['mdir'] = 0.1
    scale_factor['depth'] = 0.1


    # read points file
    print("[INFO] Reading points file")
    pnts = np.genfromtxt(pointsfile, comments='#',
            dtype=[('id','S10'),('olat','f4'),('olon','f4'),
            ('mlat','f4'),('mlon','f4'),('iy','i4'),('ix','i4')])
    npnts = np.alen(pnts)

    # Determine file type from header magic:
    with open(infiles[0], 'rb') as fd:
        hdr = fd.read(10)
    if hdr[:4] == b'\211HDF':
        fileformat = 'nc'
    else:
        fileformat = 'pp'


    # open input model file:
    if fileformat == 'pp':
        print("[INFO] Loading raw PP cubes...")
        # TODO: Check loading "raw" then extracting is acutally quicker.....
        c = iris.load_raw(infiles, callback=wstash.irisWaveCubeLoadCallback)

        print("[INFO] Extracting Depth")
        c_dep = c.extract(wstash.DEP)[0]   # just first field

        print("[INFO] Extracting Hs")
        hs = c.extract(wstash.HSIG).merge_cube()
        print("[INFO] Extracting Tm")
        tm = c.extract(wstash.TM01).merge_cube()
        print("[INFO] Extracting Tz")
        tz = c.extract(wstash.TM02).merge_cube()
        print("[INFO] Extracting mean dir")
        mdir = c.extract(wstash.MDIR).merge_cube()
        print("[INFO] Extracting wspd")
        wspd = c.extract(wstash.WSPD).merge_cube()
        print("[INFO] Extracting wdir")
        wdir = c.extract(wstash.WDIR).merge_cube()

        winduv = False
    else:
        # netCDF
        hs = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 'hs'))
        tm = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 't01'))
        tz = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 't02'))
        mdir = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 'dir'))
        uwnd = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 'uwnd'))
        vwnd = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 'vwnd'))

        winduv = True

#       # We might not have depth:
        try:
            c_dep = iris.load_cube(infiles, iris.Constraint(cube_func=lambda c: c.var_name == 'dpt'))[0]
        except:
            print("[WARN] No depth field in netCDF file")
            c_dep = None
    ##

    ntimes = hs.coord('time').points.size

    # output timeseries arrays:
    hs_ts = np.ndarray((ntimes,npnts), dtype=np.float16)
    tz_ts = np.ndarray((ntimes,npnts), dtype=np.float16)
    tm_ts = np.ndarray((ntimes,npnts), dtype=np.float16)
    mdir_ts = np.ndarray((ntimes,npnts), dtype=np.float16)
    wspd_ts = np.ndarray((ntimes,npnts), dtype=np.float16)
    wdir_ts = np.ndarray((ntimes,npnts), dtype=np.float16)

    uwnd_ts = np.ndarray((ntimes,npnts), dtype=np.float32)
    vwnd_ts = np.ndarray((ntimes,npnts), dtype=np.float32)

    for infld, tsfld in zip([hs,tm,tz,mdir],[hs_ts,tm_ts,tz_ts,mdir_ts]):
        print("[INFO] Processing: field %s" % infld.long_name)
        for itime in np.arange(ntimes):
            dat = infld[itime].data
            tsfld[itime,:] = dat[pnts['iy'],pnts['ix']]

    # derived wspd?
    if winduv:
        for infld, tsfld in zip([uwnd,vwnd],[uwnd_ts,vwnd_ts]):
            print("[INFO] Processing: field %s" % infld.long_name)
            for itime in np.arange(ntimes):
                dat = infld[itime].data
                tsfld[itime,:] = dat[pnts['iy'],pnts['ix']]
        wspd_ts = np.sqrt(uwnd_ts**2 + vwnd_ts**2)
        wdir_ts = np.ma.mod(630.0 - np.rad2deg(np.ma.arctan2(vwnd_ts,uwnd_ts)), 360)
    else:
        for infld, tsfld in zip([wspd,wdir],[wspd_ts,wdir_ts]):
            print("[INFO] Processing: field %s" % infld.long_name)
            for itime in np.arange(ntimes):
                dat = infld[itime].data
                tsfld[itime,:] = dat[pnts['iy'],pnts['ix']]
    ##--

    # extract depth:
    if c_dep is not None:
        dep = c_dep.data[pnts['iy'],pnts['ix']]
    else:
        dep = 0

    # write data to netCDF file:
    print("[INFO] Creating output file")
    if os.path.exists(fnout):
        os.unlink(fnout)

    t = hs.coord('time')
    tunits = 'days since 1990-01-01T00:00:00Z'
    if t.units.name != tunits:
        # naturalise units
        tmp = nc.num2date(t.points[:], t.units.name)
        t = nc.date2num(tmp, tunits)
    else:
        t = t.points[:]

    dout = nc.Dataset(fnout, mode='w', clobber=True, format="NETCDF4")
    dout.createDimension('time', size=ntimes)
    dout.createDimension('point', size=npnts)

    var = dout.createVariable('time', np.double, dimensions=('time'))
    var.standard_name = 'time'
    var.units = tunits
    var.conventions = 'relative julian days with decimal part (as parts of the day)'
    var.long_name = 'julian day (UT)'
    var.axis = 'T'
    var[:] = t

    var = dout.createVariable('point', np.str, dimensions=('point'))
    var[:] = pnts['id'].astype(np.object)

    var = dout.createVariable('obs_lat', np.float32, dimensions=('point'))
    var[:] = pnts['olat']

    var = dout.createVariable('obs_lon', np.float32, dimensions=('point'))
    var[:] = pnts['olon']

    var = dout.createVariable('mod_lat', np.float32, dimensions=('point'))
    var[:] = pnts['mlat']

    var = dout.createVariable('mod_lon', np.float32, dimensions=('point'))
    var[:] = pnts['mlon']

    # populate diagnostic variables:
    for fld in flds:
        print("[INFO] Saving %s" % fld)
        v = dout.createVariable(fld, np.short, dimensions=('time','point'))
        v.scale_factor = scale_factor[fld]
        v[:] = locals()["%s_ts" % fld]

    # add model depths:
    v = dout.createVariable('mod_depth', np.float32, dimensions=('point'))
    v[:] = dep


    dout.close()


def write_kml_pointfile(pointsfile, model, kmlfile, smc=False, only_ports=True):

    place_kml = """
    <Placemark>
        <name>{0} ({1})</name>
        <styleUrl>#default</styleUrl>
        <description>{2}\n{1} config</description>
        <Point>
            <coordinates>{3:.4f},{4:.4f},0</coordinates>
        </Point>
    </Placemark>
    """

    if(model.lower() == 'euro'):
        pole = (177.5, 37.55)
        marker_color = 'green'
    elif(model.lower() == 'uk'):
        pole = (177.5, 37.5)
        marker_color = 'red'
    else:
        pole = None
        marker_color = 'cyan'

    # read points file
    if smc:
        pnts = np.genfromtxt(pointsfile, dtype=[('id','S10'),('olat','f4'),('olon','f4'),
            ('mlat','f4'),('mlon','f4'),('seapoint','i4')])
    else:
        pnts = np.genfromtxt(pointsfile, comments='#',
            dtype=[('id','S10'),('olat','f4'),('olon','f4'),
            ('mlat','f4'),('mlon','f4'),('iy','i4'),('ix','i4')])

    if only_ports:
        vecfun = np.vectorize(lambda x: (x.isalpha() and len(x) == 4) or x == "DUN1")
        m = vecfun(pnts['id'])
        pnts = pnts[m]

    if pole is not None:
        from cartopy import crs
        rtd = crs.RotatedPole(pole_latitude=pole[1], pole_longitude=pole[0])
        ll = crs.PlateCarree().transform_points(rtd, pnts['mlon'], pnts['mlat'])
        pnts['mlon'] = ll[:,0]
        pnts['mlat'] = ll[:,1]


    with open(kmlfile, 'w') as fkml:
        fkml.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fkml.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        fkml.write('<Document>\n')
        fkml.write('<Style id="default">\n')
        fkml.write('    <IconStyle>\n')
        fkml.write('        <scale>0.5</scale>\n')
        fkml.write('        <Icon>\n')
        #fkml.write('            <href>http://maps.google.com/mapfiles/kml/pal4/icon49.png</href>\n')
        fkml.write('            <href>http://www-nwp/~frey/img/markers/marker_circ_%s.png</href>\n' % marker_color)
        fkml.write('        </Icon>\n')
        fkml.write('    </IconStyle>\n')
        fkml.write('</Style>\n')

        for pnt in pnts:
            print(pnt['id'])
            fkml.write(place_kml.format(pnt['id'], model, pnt['id'], pnt['mlon'], pnt['mlat']))

        fkml.write('</Document>\n')
        fkml.write('</kml>\n')
