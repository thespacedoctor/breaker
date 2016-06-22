#!/usr/local/bin/python
# encoding: utf-8
"""
*Continually listen to gracedb and await for new waves and new wave maps*

:Author:
    David Young

:Date Created:
    June 22, 2016
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import numpy as np
from ligo.gracedb.rest import GraceDb, HTTPError
from ligo.gracedb.rest import GraceDbBasic
from astrocalc.times import conversions


class listen():
    """
    *The worker class for the listen module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``label`` -- filter wave event by label. Default *EM_READY*
        - ``farThreshold`` -- the false alarm rate threshold. Default *1e-7*
        - ``startMJD`` -- startMJD. Default *56658.0 (2014-01-01)*
        - ``endMJD`` -- endMJD. Default *69807.0 (2050-01-01)*
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            label="EM_READY",
            farThreshold=1e-7,
            startMJD=56658.0,
            endMJD=69807.0
    ):
        self.log = log
        log.debug("instansiating a new 'listen' object")
        self.settings = settings
        self.label = label
        self.farThreshold = farThreshold
        self.startMJD = startMJD
        self.endMJD = startMJD

        # xt-self-arg-tmpx

        # Initial Actions
        self.mapDirectory = self.settings["gw maps directory"]
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.mapDirectory):
            os.makedirs(self.mapDirectory)

        # CONVERT MJD TO GPS
        converter = conversions(
            log=self.log
        )
        gpsZeroMjd = converter.ut_datetime_to_mjd(
            utDatetime="1980-01-06 00:00:00")
        startGPS = (startMJD - float(gpsZeroMjd)) * 60 * 60 * 24.
        endGPS = (endMJD - float(gpsZeroMjd)) * 60 * 60 * 24.

        # INSTANTIATE GRACEDB CLIENT WITH A QUERY STRING
        self.client = GraceDbBasic()
        eventString = '%(label)s far <%(farThreshold)s %(startGPS)s  .. %(endGPS)s' % locals(
        )

        # REST API RETURNS AN ITERATOR
        self.events = self.client.events(
            query=eventString,
            orderby=None,
            count=None,
            columns=None
        )

        return None

    def get_maps(
            self):
        """*download the maps for all events in our time range*

        **Usage:**

            .. code-block:: python 

                from breaker.gracedb import listen
                downloader = listen(
                    log=log,
                    settings=settings,
                    label="EM_READY",
                    farThreshold=1e-7,
                    startMJD=56658.0,
                    endMJD=69807.0
                )
                downloader.get_maps()
        """
        self.log.info('starting the ``get_maps`` method')

        eventKeys = ['graceid', 'gpstime', 'group', 'links', 'created',
                     'far', 'instruments', 'labels', 'nevents', 'submitter', 'search', 'likelihood']
        fileorder = ['LALInference_skymap.fits.gz',
                     'bayestar.fits.gz', 'LIB_skymap.fits.gz', 'skymap.fits.gz']

        for event in self.events:
            eventinfo = {}
            for key in eventKeys:
                if not key in event:
                    self.log.warning(
                        "`%(key)s` not in event %(event)s" % locals())
                    continue
                eventinfo[key] = event[key]
            eventinfo['gpstime'] = float(eventinfo['gpstime'])

            if eventinfo['far'] > self.farThreshold:
                far = eventinfo['far']
                farthres = self.farThreshold
                self.info.warning(
                    "event %(event)s does not pass FAR threshold of %(farthres)s. (FAR = %(far)s)" % locals())
                continue

            allMaps = []

            for lvfile in fileorder:
                try:
                    aMap = self.client.files(eventinfo['graceid'], lvfile)
                    allMaps.append(aMap)
                    self._write_map_to_disk(
                        sMap=aMap,
                        mapName=lvfile,
                        waveId=eventinfo['graceid']
                    )
                except:
                    eventId = eventinfo['graceid']
                    self.log.info(
                        "The %(lvfile)s path for %(eventId)s does not seem to exist yet" % locals())

            if len(allMaps) == 0:
                eventId = eventinfo['graceid']
                self.log.warning(
                    'cound not download skymaps for event %(eventId)s' % locals())

        self.log.info('completed the ``get_maps`` method')
        return None

    def _write_map_to_disk(
            self,
            sMap,
            mapName,
            waveId):
        """*Given a skymap, write it to disk*

        **Key Arguments:**
            - ``sMap`` -- the skymap from gracedb.
            - ``mapName`` -- the map flavour.
            - ``waveId`` -- the graceDB id for the wave.
        """
        self.log.info('starting the ``_write_map_to_disk`` method')

        outputDir = self.mapDirectory + "/" + waveId.upper()
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        print "Downloading %(mapName)s for GW event %(waveId)s " % locals()

        mapPath = outputDir + "/" + mapName

        skymapfile = open(mapPath, 'w')
        skymapfile.write(sMap.read())
        skymapfile.close()

        self.log.info('completed the ``_write_map_to_disk`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
