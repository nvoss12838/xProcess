"""
Class for holding newtork of GPS sites

Nick Voss USF Geodesy

borrowed heavily from obspys station
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from future.utils import native_str

import glob
import io
import copy
import os
import warnings

import numpy as np
from timeseries import TimeSeries
class Station(object):
    """
    From the StationXML definition:
        This type represents a Station epoch. It is common to only have a
        single station epoch with the station's creation and termination dates
        as the epoch start and end dates.
    """
    def __init__(self,name,latitude, longitude, elevation, data=None,
                 monument=None,setting=None, antenna=None,jump_list=None,
                 operators=None,ts=None):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.data = data or []
        self.monument = monument
        self.setting = setting
        self.operators = operators or []
        self.ts = ts

    def plot_map(self):
        '''
        plot a map of the station
        '''
        return
    def process(self,tree=None,start_date=None,end_date=None):
        '''
        Process each day between start date and end date if it exists.
        tree is a tree object
        '''
        import glob
        import pandas as pd
        #get all files in directory *data*
        files = glob.glob(self.data)
        files = files[0:200]
        if self.ts is None:
            self.ts = TimeSeries()
        for f in files:
            #figure out what day it is
            #? read the rinex file check if day is already in TS?
            #process each day
            if self.ts is None:
                self.ts = TimeSeries()
            os.system('tropNominal.py -m VMF1 -b %f -e %e -stns %s -append -o tdpIn.tdp'%(start,end,self.station))
            # change tree to add in correct second order ionosphere
            
            if tree is None:
                os.system('gd2e.py -runType=PPP -rnxFile %s -gdCov -nProcessors=4 -GNSSproducts /home/nvoss/orbits/sideshow.jpl.nasa.gov/pub/JPL_GPS_Products/Final'%(f))
            else:
                os.system('gd2e.py -runType=PPP -treeS Trees -rnxFile %s -gdCov -GNSSproducts /home/nvoss/orbits/sideshow.jpl.nasa.gov/pub/JPL_GPS_Products/Final'%(f))
            #grep the summary file

            with open('smoothFinal.gdcov') as fil:
                content = fil.readlines() #read the files line by line
                X = content[1]
                Y = content[2]
                Z = content[3]
                time,x,ux = splitline(X)
                time,y,uy = splitline(Y)
                time,z,uz = splitline(Z)
                time1 = pd.datetime(2000,1,1,11,59,47)
                time2 = pd.to_datetime(0)
                timedif = time1-time2
                time = pd.to_datetime(time,unit='s')+timedif
                df = pd.DataFrame([[time,x,y,z,ux,uy,uz]],columns=['Time','X','Y','Z','UX','UY','UZ'])
                self.ts.add_data(df)
                os.system('mkdir %s'%(f.split('/')[7]))
                os.system('mkdir %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                os.system('cp runAgain %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                os.system('cp rtgx_ppp_0.tree.err0_0 %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                os.system('cp rtgx_ppp_0.tree.log0_0 %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))

            #append position and uncertainty to TS
def splitline(line):
  index,sta,time,position,unc = line.split(' ')
  return int(time),float(position),float(unc)
