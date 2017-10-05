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
import fileinput

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
        import fileinput
        import subprocess

        #get all files in directory *data*
        files = glob.glob(self.data)
        #files = files[0:2]
        if self.ts is None:
            self.ts = TimeSeries()
        for f in files:
            #figure out what day it is
            #? read the rinex file check if day is already in TS?
            #process each day
            if self.ts is None:
                self.ts = TimeSeries()

            # change tree to add in correct second order ionosphere
            #replace the ionex file summary file
            year = f.split('/')[7]
            print(f.split('/')[8])
            day = f.split('/')[8][4:7]
            print(day,year)
            ionFile = ('/home/nvoss/goa-var/cddis.gsfc.nasa.gov/gps/products/ionex/%s/%s/jplg%s0.%si'%(year,day,day,year[2:]))
            print(ionFile)

            date = subprocess.check_output('doy2date %s %s'%(day,year),shell=True)
            if f==files[0]:
                startDate = date
            start = subprocess.check_output('date2sec %s 00:00:00'%(date[0:10]),shell=True)[0:11]
            end = subprocess.check_output('date2sec %s 23:59:59'%(date[0:10]),shell=True)[0:11]
            if not os.path.exists('%s.gdcov'%(date[:-1]) and year!='2006'):
                #print('tropNominal.py -m VMF1 -b %s -e %s -stns %s -append -o tdpIn.tdp'%(start,end,self.name))
                #os.system('tropNominal.py -m VMF1 -b %s -e %s -stns %s -append -o tdpIn.tdp'%(start,end,self.name))
                #x = fileinput.input(files='Trees/Nick_0.tree',inplace=1)
                #for line in x:
                    #if 'IONEXFILE ==*' in line:
                        #line = line.replace("IONEXFILE ==*",'IONEXFILE == "%s"'%(ionFile))
                replace('Trees/Nick_0.tree',"GLOBAL_EPOCH ==",'GLOBAL_EPOCH == %s'%(start))
                replace('Trees/Nick_0.tree',"IONEXFILE ==",'IONEXFILE == %s'%(ionFile))
                if tree is None:
                    os.system('gd2e.py -runType=PPP -rnxFile %s -gdCov -nProcessors=4 -GNSSproducts /home/nvoss/orbits/sideshow.jpl.nasa.gov/pub/JPL_GPS_Products/Final'%(f))
                    os.system('gd2e.py -runType=PPP -treeS Trees -rnxFile %s -nProcessors=4 -gdCov -GNSSproducts /home/nvoss/orbits/sideshow.jpl.nasa.gov/pub/JPL_GPS_Products/Final'%(f))
                os.system('cp smoothFinal.gdcov %s.gdcov'%(date[:-1]))
                print('sAVED')
                with open('smoothFinal.gdcov') as fil:
                    content = fil.readlines() #read the files line by line
                    X = content[1]
                    Y = content[2]
                    Z = content[3]
                    #time,x,ux = splitline(X)
                    #time,y,uy = splitline(Y)
                    #time,z,uz = splitline(Z)
                    #time1 = pd.datetime(2000,1,1,11,59,47)
                    #time2 = pd.to_datetime(0)
                    #timedif = time1-time2
                    #time = pd.to_datetime(time,unit='s')+timedif
                    #df = pd.DataFrame([[time,x,y,z,ux,uy,uz]],columns=['Time','X','Y','Z','UX','UY','UZ'])
                    #self.ts.add_data(df)
                    os.system('mkdir %s'%(f.split('/')[7]))
                    os.system('mkdir %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                    os.system('cp runAgain %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                    os.system('cp rtgx_Nick_0.tree.err0_0 %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))
                    os.system('cp rtgx_Nick_0.tree.log0_0 %s/%s'%(f.split('/')[7],f.split('/')[-1][4:7]))

        os.system('netSeries.py -r %s.gdcov -i *.gdcov'%(startDate[:-1]))
            #append position and uncertainty to TS
def splitline(line):
  index,sta,time,position,unc = line.split(' ')
  return int(time),float(position),float(unc)

import re

import fileinput
import sys
def replace(fil,inp,output):
    for line in fileinput.input(fil,inplace=True):
        # Whatever is written to stdout or with print replaces
        # the current line
        if line.startswith(inp):
            print(output)
        else:
            sys.stdout.write(line)
