import numpy as np
from os.path import join
from collections import OrderedDict
import datetime as dt

"""
Created on Fri Feb 12 16:42:27 2021

Module with functions to read JCOMM-WFVS observations (i.e., from LCWFV_collocations) 

@author: nvalient
"""

def genMonthstats(centres):
    """Generate the monthly stats from centres data"""

    monthstats = {}
    for centre, cvals in centres.items():
        if centre == "Obs":
            datacount = np.float(np.ma.count(cvals["data"]))
            datamean = np.ma.mean(cvals["data"])
            datastd = np.ma.std(cvals["data"])
            monthstats[centre] = {
                "count": datacount,
                "mean": datamean,
                "std": datastd,
            }
        else:
            errors = cvals["data"] - centres["Obs"]["data"]
            ## diagnostic scatter plot to check data
            # if centre == 'MetOffice':
            #    plt.scatter(cvals['data'], centres['Obs']['data'])
            #    plt.show()
            # apply a qc threshold to the error data
            errors = np.ma.masked_where(np.abs(errors) > 15.0, errors)
            datacount = np.float(np.ma.count(errors))
            if datacount > 100.0:
                datamean = np.ma.mean(errors)
                datastd = np.ma.std(errors)
                monthstats[centre] = {
                    "count": datacount,
                    "bias": datamean,
                    "stdd": datastd,
                }
            else:
                monthstats[centre] = {
                    "count": -999.9,
                    "bias": -999.9,
                    "stdd": -999.9,
                }

    return monthstats

def LCread(filein, leadtime: int = 0):
    """Read data from LCWFV files - without ID"""

    centres = {
        "Time": {"col": "valid_date", "data": []},
        "Lat": {"col": "latitude", "data": []},
        "Lon": {"col": "longitude", "data": []},
        "Obs": {"col": "obs", "data": []},
        "MetOffice": {"col": "egrr", "data": []},
        "ECMWF": {"col": "ecmf", "data": []},
        "NCEP": {"col": "kwbc", "data": []},
        "MeteoFrance": {"col": "lfpw", "data": []},
        "DWD": {"col": "edzw", "data": []},
        "BoM": {"col": "ammc", "data": []},
        "JMA": {"col": "rjtd", "data": []},
        "KMA": {"col": "rksl", "data": []}
    }
#"Station_id": {"col": "station:id", "data": []}
    print(f"Reading JCOMM data from {filein}")
    with open(filein, "r") as reader:
        rddata = reader.readlines()
    #        inp.close()

    # leadtime string
    ltstr = f"{leadtime:03d}"

    # get indices for columns
    hdrs = rddata[0].split()
    for centre in centres.keys():
        if centres[centre]["col"] in hdrs:
            centres[centre]["idx"] = hdrs.index(centres[centre]["col"])
        else:
            centres[centre]["idx"] = None

    # read the file data
    for idata in range(1, len(rddata)):
        tdata = rddata[idata].split()
        if tdata[1] != ltstr:
            continue
        # containers for data checks
        for centre, cvals in centres.items():
            if cvals["idx"] is not None:
                cvals["data"].append(np.float(tdata[cvals["idx"]]))
            else:
                cvals["data"].append(99.99)
    # masked
    for centre, cvals in centres.items():
        cvals["data"] = np.array(cvals["data"])
        cvals["data"] = np.ma.masked_where(cvals["data"] > 99.0, cvals["data"])

    return centres

def LCreadAll(File):
    """Read data from LCWFV files"""
    
    matrix = [] 
    n = 0
    for line in open(File):
        List = line.strip().split(' ')
        List = list(OrderedDict.fromkeys(List))
        #print List
        row = []
        try:
            if int(List[2]) == 0:
               n += 1
               row.append(int(List[0]))   # date
               row.append(float(List[3])) # lat
               row.append(float(List[4])) # lon
               row.append(float(List[5])) # obs
               row.append(str(List[-1]))  # station ID    
        except:
            ValueError
        if len(row) > 0: 
           #print row
           matrix.append(row)
    matrix = np.asarray(matrix)        
    #print(n)

    ids   = np.unique(matrix[:,-1])
    TIMES  = np.unique(matrix[:,0])
    
    OBS = []
    for ID in ids:
        STATION = matrix[np.where(matrix==ID)[0],:]
        OBS.append(STATION)

    return OBS, TIMES

def read_wfvsNew(dirIn, filename, fileVarN, addmodel=None, verbose=False):
    """convert variables in centres into arrays"""
    
    stn_fileName = '_data_'
    
    month = filename[-10:-4]       
    
    hsarr = []
    vtarr = []
    latsarr = []
    lonsarr = []
    idsarr = []

    filein = join(dirIn,fileVarN+stn_fileName+month+'.txt')
    #leadtimes = [0]
    OBS, TIMES = LCreadAll(filein)
    fillarr    = np.ones(len(TIMES))*-99.99
    # Matrix in OBS is [date,lat,lon,obs,station ID]      
    for iID in range(np.shape(OBS)[0]): # loop over number of locations
        data = OBS[iID][:]
        tdata = data[:,0] 
        latsarr.append(data[0,1])
        lonsarr.append(data[0,2])
        idsarr.append(data[0,4])
        # if timeseries is complete
        if len(tdata) == len(TIMES):
            vtarr = TIMES # array of time
            if iID == 0:                
                hsarr = data[:,3]
            else:
                hsarr = np.vstack((hsarr,data[:,3])) # (locations,time)                  

        # condition for those locations that do not record all the time
        elif len(tdata) != len(TIMES):
            hs0 = fillarr              
            hstemp = data[:,3]
            # loop over time as some locations do not record all the time
            cont = 0
            for it in np.shape(len(tdata)):
                mdata = np.where(TIMES == tdata[it])
                hs0[mdata] = hstemp[cont]
                cont += 1
            if iID == 0:
                hsarr = hs0
            else:
                hsarr = np.vstack((hsarr,hs0))          
                    
    return idsarr, latsarr, lonsarr, vtarr, hsarr

def read_wfvsAllVar(dirIn, filename):
    """Final function that returns in-situ obs as arrays with dimension: [locations,times]"""
    # File name - individual files per variable
    #fileVarN = ['SWH','PP1D','V10']
    hsidsarr, hslatsarr, hslonsarr, vtarr, hsarr = read_wfvsNew(dirIn, filename, 'SWH')
    tpidsarr, tplatsarr, tplonsarr, tpvtarr, tparr0 = read_wfvsNew(dirIn, filename, 'PP1D')
    wsidsarr, wslatsarr, wslonsarr, wsvtarr, wsarr0 = read_wfvsNew(dirIn, filename, 'V10')
    lat = np.array(hslatsarr)
    lon = np.array(hslonsarr)
    tpidsarr = np.array(tpidsarr)
    wsidsarr = np.array(wsidsarr)    
    
    # Check if id coincide --> start with Hs; when no data fill with -99.99
    tparr1 = np.ones((len(hsidsarr),len(vtarr)))*-99.99
    wsparr1 = np.ones((len(hsidsarr),len(vtarr)))*-99.99
    ia = 0
    for ihs in hsidsarr:
        mmtp = tpidsarr == ihs 
        mmwsp = wsidsarr == ihs
        if True in mmtp:
            tparr1[ia,:] = tparr0[np.where(tpidsarr==ihs),:]
        elif True in mmwsp:
            wsparr1[ia,:] = wsarr0[np.where(wsidsarr==ihs),:]
        else:
            continue
        ia += 1
    
    idsarr = np.array(hsidsarr)   
    tparr  = tparr1
    wsparr = wsparr1
    dates   = vtarr
    vtarr = [dt.datetime.strptime(date, '%Y%m%d%H') for date in dates]
    
    return idsarr, lat, lon, vtarr, hsarr, tparr, wsparr


# def read_wfvsNewNO(dirIn, filename, addmodel=None, verbose=False):
#     """convert variables in centres into arrays - BROKEN  FC!!!"""
    
#     stn_fileName = '_data_'
    
#     month = filename[-10:-4]    
#     # File name - individual files pero variable
#     fileVarN = ['SWH','PP1D','V10']
    
#     parameters = {"Hs": "SWH", "Tp": "PPD1", "Ws": "V10"}
#     for ii in range(len(parameters)):
#         filein = join(dirIn,fileVarN[ii]+stn_fileName+month+'.txt')
#         #leadtimes = [0]
#         centres = LCread(filein, leadtime=0)
#         if ii == 0:
#             t = centres["Time"]["data"]
#             lat = centres["Lat"]["data"]
#             lon = centres["Lon"]["data"]
#             hs = centres["Obs"]["data"]
#             ID = centres["Station_id"]["data"]
#         elif ii == 1:
#             Tp = centres["Obs"]["data"]
#         elif ii == 2:
#             wspd = centres["Obs"]["data"]
    
#     # Masked by ID
#     id_counter = np.unique(ID)
#     idsarr     = id_counter
#     hsarr = []
#     tparr = []
#     wsarr = []
#     vtarr = []
#     latsarr = []
#     lonsarr = []
#     for id_c in id_counter:
#         hsarr.append(hs[id_c==ID])
#         tparr.append(Tp[id_c==ID])
#         wsarr.append(wspd[id_c==ID])
#         vtarr.append(t[id_c==ID])
#         latsarr.append(lat[id_c==ID])
#         lonsarr.append(lon[id_c==ID])
    
#     return idsarr, latsarr, lonsarr, vtarr, wsarr, hsarr, tparr