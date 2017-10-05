"""
Class for holding newtork of GPS sites

Nick Voss USF Geodesy

borrowed heavily from obspys catalog
"""

import glob
import io
import copy
import os
import warnings

import numpy as np
import pandas as pd
from .station import Station


class Network(object):
    """
    Holds a collection of GPS stations as well as metadata on the Network
    with a few useful functions
    """
    def __init__(self, stations=None, **kwargs):
        if not stations:
            self.stations = []
        else:
            self.stations = stations
        self.comments = kwargs.get("comments", [])
        self._set_resource_id(kwargs.get("resource_id", None))
        self.description = kwargs.get("description", "")
        self._set_creation_info(kwargs.get("creation_info", None))
        def _get_resource_id(self):
        return self.__dict__['resource_id']

    def addStations(self,stationList):
        '''
        read a text file containtining:
        Station Lon Lat Height
        add to stations
        '''
        df = pd.read_fwf(stationList,names=['name','Lon','Lat','Height'])
        stations = []
        for i,name in enumerate(df.name):
            stations.append(Station(name.df.Lat[i],df.Lon[i],df.Height[i]))
        self.stations = stations
    def process(self):
        for station in stations:
            #make a directory to store processed station
            #change to the directory
            #copy the tree overview
            #process the station
            station.process(tree='Trees')
