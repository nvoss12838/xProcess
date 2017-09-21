"""
Class for holding newtork of GPS sites

Nick Voss USF Geodesy

borrowed heavily from obspys catalog
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

from .station import Station

EVENT_ENTRY_POINTS = ENTRY_POINTS['event']
EVENT_ENTRY_POINTS_WRITE = ENTRY_POINTS['event_write']

class Network(object):
    "
    Holds a collection of GPS stations as well as metadata on the Network
    with a few useful functions
    "
    def __init__(self, stations=None, **kwargs):
        if not events:
            self.stations = []
        else:
            self.stations = stations
        self.comments = kwargs.get("comments", [])
        self._set_resource_id(kwargs.get("resource_id", None))
        self.description = kwargs.get("description", "")
        self._set_creation_info(kwargs.get("creation_info", None))
        def _get_resource_id(self):
        return self.__dict__['resource_id']

    def _set_resource_id(self, value):
        if type(value) == dict:
            value = ResourceIdentifier(**value)
        elif type(value) != ResourceIdentifier:
            value = ResourceIdentifier(value)
        self.__dict__['resource_id'] = value

    resource_id = property(_get_resource_id, _set_resource_id)

    def _get_creation_info(self):
        return self.__dict__['creation_info']

    def _set_creation_info(self, value):
        if type(value) == dict:
            value = CreationInfo(**value)
        elif type(value) != CreationInfo:
            value = CreationInfo(value)
        self.__dict__['creation_info'] = value

    creation_info = property(_get_creation_info, _set_creation_info)

    def __add__(self, other):
        """
        Method to add two catalogs.
        """
        if isinstance(other, Event):
            other = Network([other])
        if not isinstance(other, Network):
            raise TypeError
        events = self.stations + other.stations
        return self.__class__(stations=stations)

    def __delitem__(self, index):
        """
        Passes on the __delitem__ method to the underlying list of traces.
        """
        return self.events.__delitem__(index)

    def __eq__(self, other):
        """
        __eq__ method of the Catalog object.
        :type other: :class:`~obspy.core.event.Catalog`
        :param other: Catalog object for comparison.
        :rtype: bool
        :return: ``True`` if both Catalogs contain the same events.
        .. rubric:: Example
        >>> from obspy.core.event import read_events
        >>> cat = read_events()
        >>> cat2 = cat.copy()
        >>> cat is cat2
        False
        >>> cat == cat2
        True
        """
        if not isinstance(other, Network):
            return False
        if self.stations != other.stations:
            return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, index):
        """
        __getitem__ method of the Catalog object.
        :return: Event objects
        """
        if index == "extra":
            return self.__dict__[index]
        if isinstance(index, slice):
            return self.__class__(station=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def __getslice__(self, i, j, k=1):
        """
        __getslice__ method of the Catalog object.
        :return: Catalog object
        """
        # see also http://docs.python.org/reference/datamodel.html
        return self.__class__(events=self.stations[max(0, i):max(0, j):k])

    def __iadd__(self, other):
        """
        Method to add two catalog with self += other.
        It will extend the current Catalog object with the events of the given
        Catalog. Events will not be copied but references to the original
        events will be appended.
        :type other: :class:`~obspy.core.event.Catalog` or
            :class:`~obspy.core.event.event.Event`
        :param other: Catalog or Event object to add.
        """
        if isinstance(other, Station):
            other = Network(stations=[other])
        if not isinstance(other, Network):
            raise TypeError
        self.extend(other.stations)
        return self

    def __iter__(self):
        """
        Return a robust iterator for Events of current Catalog.
        Doing this it is safe to remove events from catalogs inside of
        for-loops using catalog's :meth:`~obspy.core.event.Catalog.remove`
        method. Actually this creates a new iterator every time a event is
        removed inside the for-loop.
        """
        return list(self.stations).__iter__()

    def __len__(self):
        """
        Returns the number of Events in the Catalog object.
        """
        return len(self.stations)

    count = __len__

    def __setitem__(self, index, station):
        """
        __setitem__ method of the Catalog object.
        """
        if not isinstance(index, (str, native_str)):
            self.events.__setitem__(index, station)
        else:
            super(Catalog, self).__setitem__(index, station)

    def __str__(self, print_all=False):
        """
        Returns short summary string of the current catalog.
        It will contain the number of Events in the Catalog and the return
        value of each Event's :meth:`~obspy.core.event.event.Event.__str__`
        method.
        :type print_all: bool, optional
        :param print_all: If True, all events will be printed, otherwise a
            maximum of ten event will be printed.
            Defaults to False.
        """
        out = str(len(self.events)) + ' Stations(s) in Network:\n'
        if len(self) <= 10 or print_all is True:
            out += "\n".join([ev.short_str() for ev in self])
        else:
            out += "\n".join([ev.short_str() for ev in self[:2]])
            out += "\n...\n"
            out += "\n".join([ev.short_str() for ev in self[-2:]])
            out += "\nTo see all events call " + \
                   "'print(CatalogObject.__str__(print_all=True))'"
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(print_all=p.verbose))

    def append(self, station):
        """
        Appends a single Event object to the current Catalog object.
        """
        if isinstance(station, Station):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single Event object as an argument.'
            raise TypeError(msg)

    def clear(self):
        """
        Clears event list (convenient method).
        .. rubric:: Example
        >>> from obspy.core.event import read_events
        >>> cat = read_events()
        >>> len(cat)
        3
        >>> cat.clear()
        >>> cat.events
        []
        """
        self.stations = []

    def copy(self):
        """
        Returns a deepcopy of the Catalog object.
        :rtype: :class:`~obspy.core.stream.Catalog`
        :return: Copy of current catalog.
        .. rubric:: Examples
        1. Create a Catalog and copy it
            >>> from obspy.core.event import read_events
            >>> cat = read_events()
            >>> cat2 = cat.copy()
           The two objects are not the same:
            >>> cat is cat2
            False
           But they have equal data:
            >>> cat == cat2
            True
        2. The following example shows how to make an alias but not copy the
           data. Any changes on ``cat3`` would also change the contents of
           ``cat``.
            >>> cat3 = cat
            >>> cat is cat3
            True
            >>> cat == cat3
            True
        """
        return copy.deepcopy(self)

    def extend(self, station_list):
        """
        Extends the current Catalog object with a list of Event objects.
        """
        if isinstance(station_list, list):
            for _i in station_list:
                # Make sure each item in the list is a event.
                if not isinstance(_i, Station):
                    msg = 'Extend only accepts a list of Event objects.'
                    raise TypeError(msg)
            self.events.extend(station_list)
        elif isinstance(event_list, Network):
            self.events.extend(station_list.events)
        else:
            msg = 'Extend only supports a list of Event objects as argument.'
            raise TypeError(msg)

    def write(self, filename, format, **kwargs):
        """
        not implemented
        will use pickle or sql database?
        """
        return

    def plot(self, projection='global', resolution='l',
             continent_fill_color='0.9', water_fill_color='1.0',
             label='station', color=k, colormap=None, show=True,
             outfile=None, method=None, fig=None, title=None,
             **kwargs):  # @UnusedVariable
        """
        Creates preview map of all events in current Catalog object.
        :type projection: str, optional
        :param projection: The map projection. Currently supported are:
            * ``"global"`` (Will plot the whole world.)
            * ``"ortho"`` (Will center around the mean lat/long.)
            * ``"local"`` (Will plot around local events)
            Defaults to "global"
        :type resolution: str, optional
        :param resolution: Resolution of the boundary database to use. Will be
        :param projection: The map projection. Currently supported are:
            * ``"global"`` (Will plot the whole world.)
            * ``"l"`` (low)
            * ``"i"`` (intermediate)
            * ``"h"`` (high)
            * ``"f"`` (full)
            Defaults to ``"l"``
        :type continent_fill_color: Valid matplotlib color, optional
        :param continent_fill_color:  Color of the continents. Defaults to
            ``"0.9"`` which is a light gray.
        :type water_fill_color: Valid matplotlib color, optional
        :param water_fill_color: Color of all water bodies.
            Defaults to ``"white"``.
        :type label: str, optional
        :param label: Events will be labelled based on the chosen property.
            Possible values are:
            * ``"magnitude"``
            * ``None``
            Defaults to ``"magnitude"``
        :type color: str, optional
        :param color: The events will be color-coded based on the chosen
            property. Possible values are:
            * ``"date"``
            * ``"depth"``
            Defaults to ``"depth"``
        :type colormap: str, any matplotlib colormap, optional
        :param colormap: The colormap for color-coding the events.
            The event with the smallest property will have the
            color of one end of the colormap and the event with the biggest
            property the color of the other end with all other events in
            between.
            Defaults to None which will use the default colormap for the date
            encoding and a colormap going from green over yellow to red for the
            depth encoding.
        :type show: bool
        :param show: Whether to show the figure after plotting or not. Can be
            used to do further customization of the plot before
            showing it. Has no effect if `outfile` is specified.
        :type outfile: str
        :param outfile: Output file path to directly save the resulting image
            (e.g. ``"/tmp/image.png"``). Overrides the ``show`` option, image
            will not be displayed interactively. The given path/filename is
            also used to automatically determine the output format. Supported
            file formats depend on your matplotlib backend.  Most backends
            support png, pdf, ps, eps and svg. Defaults to ``None``.
            The figure is closed after saving it to file.
        :type method: str
        :param method: Method to use for plotting. Possible values are:
            * ``'basemap'`` to use the Basemap library
            * ``'cartopy'`` to use the Cartopy library
            * ``None`` to pick the best available library
            Defaults to ``None``.
        :type fig: :class:`matplotlib.figure.Figure` (or
            :class:`matplotlib.axes.Axes`)
        :param fig: Figure instance to reuse, returned from a previous
            inventory/catalog plot call with `method=basemap`.
            If a previous basemap plot is reused, any kwargs regarding the
            basemap plot setup will be ignored (i.e.  `projection`,
            `resolution`, `continent_fill_color`, `water_fill_color`). Note
            that multiple plots using colorbars likely are problematic, but
            e.g. one station plot (without colorbar) and one event plot (with
            colorbar) together should work well.
            If an :class:`~matplotlib.axes.Axes` is supplied, the given axes is
            used to plot into and no colorbar will be produced.
        :type title: str
        :param title: Title above plot. If left ``None``, an automatic title
            will be generated. Set to ``""`` for no title.
        :returns: Figure instance with the plot.
        .. rubric:: Examples
        Mollweide projection for global overview:
        >>> from obspy import read_events
        >>> cat = read_events()
        >>> cat.plot()  # doctest:+SKIP
        .. plot::
            from obspy import read_events
            cat = read_events()
            cat.plot()
        Orthographic projection:
        >>> cat.plot(projection="ortho")  # doctest:+SKIP
        .. plot::
            from obspy import read_events
            cat = read_events()
            cat.plot(projection="ortho")
        Local (Albers equal area) projection:
        >>> cat.plot(projection="local")  # doctest:+SKIP
        .. plot::
            from obspy import read_events
            cat = read_events()
            cat.plot(projection="local")
        Combining a station and event plot (uses basemap):
        >>> from obspy import read_inventory, read_events
        >>> inv = read_inventory()
        >>> cat = read_events()
        >>> fig = inv.plot(method=basemap, show=False)  # doctest:+SKIP
        >>> cat.plot(method=basemap, fig=fig)  # doctest:+SKIP
        .. plot::
            from obspy import read_inventory, read_events
            inv = read_inventory()
            cat = read_events()
            fig = inv.plot(show=False)
            cat.plot(fig=fig)
        """
        from obspy.imaging.maps import plot_map, _plot_basemap_into_axes
        import matplotlib
        import matplotlib.pyplot as plt

        if color not in ('date', 'depth'):
            raise ValueError('Events can be color coded by date or depth. '
                             "'%s' is not supported." % (color,))
        if label not in (None, 'magnitude', 'depth'):
            raise ValueError('Events can be labeled by magnitude or events can'
                             ' not be labeled. '
                             "'%s' is not supported." % (label,))

        # lat/lon coordinates, magnitudes, dates
        lats = []
        lons = []
        labels = []
        mags = []
        colors = []
        times = []
        for event in self:
            if not event.origins:
                msg = ("Event '%s' does not have an origin and will not be "
                       "plotted." % str(event.resource_id))
                warnings.warn(msg)
                continue
            if not event.magnitudes:
                msg = ("Event '%s' does not have a magnitude and will not be "
                       "plotted." % str(event.resource_id))
                warnings.warn(msg)
                continue
            origin = event.preferred_origin() or event.origins[0]
            lats.append(origin.latitude)
            lons.append(origin.longitude)
            times.append(origin.time)
            magnitude = event.preferred_magnitude() or event.magnitudes[0]
            mag = magnitude.mag
            mags.append(mag)
            labels.append(('  %.1f' % mag) if mag and label == 'magnitude'
                          else '')
            if color == 'date':
                c_ = origin.get('time') or np.nan
            else:
                c_ = (origin.get('depth') or np.nan) / 1e3
            colors.append(c_)

        # Create the colormap for date based plotting.
        if colormap is None:
            colormap = obspy_sequential

        if title is None:
            if len(lons) > 1:
                # if we have a `None` in the origin time list it likely ends up
                # as min and/or max and causes problems..
                times_ = np.ma.masked_equal(times, None).compressed()
                min_time = times_.min()
                max_time = times_.max()
                title = (
                    "{event_count} events ({start} to {end}) "
                    "- Color codes {colorcode}, size the magnitude".format(
                        event_count=len(self.events),
                        start=min_time.strftime("%Y-%m-%d"),
                        end=max_time.strftime("%Y-%m-%d"),
                        colorcode="origin time"
                                  if color == "date"
                                  else "depth"))
            else:
                title = "Event at %s" % times[0].strftime("%Y-%m-%d")

        if color not in ("date", "depth"):
            msg = "Invalid option for 'color' parameter (%s)." % color
            raise ValueError(msg)

        min_size = 2
        max_size = 30
        min_size_ = min(mags) - 1
        max_size_ = max(mags) + 1
        if len(lons) > 1:
            frac = [(0.2 + (_i - min_size_)) / (max_size_ - min_size_)
                    for _i in mags]
            size_plot = [(_i * (max_size - min_size)) ** 2 for _i in frac]
        else:
            size_plot = 15.0 ** 2

        if isinstance(fig, matplotlib.axes.Axes):
            if method is not None and method != "basemap":
                msg = ("Plotting into an matplotlib.axes.Axes instance "
                       "currently only implemented for `method='basemap'`.")
                raise NotImplementedError(msg)
            ax = fig
            fig = ax.figure
            _plot_basemap_into_axes(
                ax=ax, lons=lons, lats=lats, size=size_plot,
                color=colors, bmap=None, labels=labels,
                projection=projection, resolution=resolution,
                continent_fill_color=continent_fill_color,
                water_fill_color=water_fill_color,
                colormap=colormap, marker="o", title=title,
                show=False, **kwargs)
        else:
            fig = plot_map(method, lons, lats, size_plot, colors, labels,
                           projection=projection, resolution=resolution,
                           continent_fill_color=continent_fill_color,
                           water_fill_color=water_fill_color,
                           colormap=colormap, marker="o", title=title,
                           show=False, fig=fig, **kwargs)

        if outfile:
            fig.savefig(outfile)
            plt.close(fig)
        else:
            if show:
                plt.show()

        return fig
