import numpy as np
from collections import OrderedDict
import os
import csv

def obs_reader(File):
  
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

    OBS = []
    for ID in ids:
        STATION = matrix[np.where(matrix==ID)[0],:]
        OBS.append(STATION)

    return OBS

def get_mean_stats(stat):
    return np.nanmean(np.array(stat))

def read_summary_csv(fileIn):
    """Function to read the .CSV files with the summary stats;
    Stats are per location (IDs)
    
    Average stats are printed"""
    
    if os.path.exists(fileIn) == True:     
        lat = []
        lon = []
        mean = []
        b   = []
        rmse = []
        r   = []
        std = []
              
        # Read SumamrStats.CSV file per variable
        header = ['ID','Lat','Lon','Samples','Model Mean','Model Std','Ob Mean','Ob Std',
                  'Bias','RMSD','Error Std','SI','Sym Slope','R value','Slope','Offset'] 
                  
        with open(fileIn, mode='r') as csvfile:
            print('[Info] Reading '+fileIn)
            next(csvfile) # skip first line
            csvreader = csv.DictReader(csvfile, delimiter=',')

            for lines in csvreader:
                mean.append(lines['Ob Mean'])
                b.append(lines['Bias'])
                rmse.append(lines['RMSD'])
                lat.append(lines['Lat'])
                lon.append(lines['Lon'])
                r.append(lines['R value'])
                std.append(lines['Error Std'])
        
        for i in range(len(b)):
            if b[i] == 'nan':
                b[i] = np.nan
                mean[i] = np.nan
                r[i] = np.nan
                rmse[i] = np.nan
                std[i] = np.nan
            else:
                mean[i] = float(mean[i])
                b[i] = float(b[i])
                r[i] = float(r[i])
                rmse[i] = float(rmse[i])
                std[i] = float(std[i])

        lat = list(map(float,lat))
        lon = list(map(float,lon))
        
        #print('[Debug] bias = ' + str(b1))            
        Mstats = [get_mean_stats(b),get_mean_stats(r),get_mean_stats(rmse),get_mean_stats(std)] 
        print('Mean stats '+' - [bias,r,rmsd,std] = ')
        print(str(Mstats))
                
    return lat, lon, mean, b, rmse, std, r 
