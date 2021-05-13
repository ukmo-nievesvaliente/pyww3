import numpy as np
from collections import OrderedDict

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
