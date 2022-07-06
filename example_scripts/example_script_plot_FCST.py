#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 10:03:35 2022

EXAMPLE SCRIPT TO PLOT RESULTS FROM FORECAST VERIFICATION 

@author: nvalient
"""

import sys
import os
from os.path import join
sys.path.insert(0,'/net/home/h01/nvalient/nvalient-python/')
# import local library
import pyww3.plot.plot_obs_mean_stats as pltmeanstats


DATADIR = '/data/users/nvalient/verification/verification_PS45-FCST/'
OUT_DIR = '/scratch/nvalient/img_testing'
RUN_NAME = ['ps45-ukS','ps45-ukW']
RUN_NAME_CTRL = ['ps45CTRL-ukS','ps45CTRL-ukW']
date_strS = ['20190720_20190814','20191204_20200124'] # Add date string with _ in between
FCST_LEN = ['T+24','T+48','T+72','T+96','T+120','T+144']
obs_type = 'WFVS' # choose between WFVS, MA_SUP03, WAVENET, SHPSYN // Depending on the observation type used, verified plotted variables will differ 

# Call the function for plotting a verification period 
for i,dateS in enumerate(date_strS):
    # Create folders for independent verification periods; otherwise the .png will be overwritten
    # check if dirout exists
    if not os.path.exists(join(OUT_DIR,dateS)):
        os.makedirs(join(OUT_DIR,dateS))
    inFolder=join(DATADIR,RUN_NAME[i])
    inFolderCTRL = join(DATADIR,RUN_NAME_CTRL[i])
    pltmeanstats.plot_FCST(inFolder,inFolderCTRL,dateS,FCST_LEN,obs_type,join(OUT_DIR,dateS))