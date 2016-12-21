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
import yaml
import codecs
import json
import collections
import numpy as np
from ligo.gracedb.rest import GraceDb, HTTPError
from ligo.gracedb.rest import GraceDbBasic
from astrocalc.times import conversions
from astrocalc.times import now as mjdNow
import time
from astropy.time import Time


class listen():
    """
    *The listen object; connects to GraceDB and 'listens' for new wave events and new skymaps*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``label`` -- filter wave event by label. Default *ADVOK & EM_READY*
        - ``farThreshold`` -- the false alarm rate threshold. Default *1e-7*
        - ``startMJD`` -- startMJD. Default *57266.0 (2015-09-01)*
        - ``endMJD`` -- endMJD. Default *69807.0 (2050-01-01)*
        - ``daemon`` -- run in daemon mode. Set to 'True' to run every 30 sec, or pass in an integer to set a custom time frequency
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            label="ADVOK & EM_READY",
            farThreshold=1e-7,
            startMJD=57266.0,
            endMJD=False,
            daemon=False
    ):
        self.log = log
        log.debug("instansiating a new 'listen' object")
        self.settings = settings
        self.label = label
        self.farThreshold = farThreshold
        self.startMJD = startMJD
        self.endMJD = endMJD
        self.daemon = daemon

        # SET DEFAULT FREQUENCY TO 30 SEC
        if self.daemon == True:
            self.daemon = 30
        # xt-self-arg-tmpx

        if not self.endMJD:
            self.endMJD = Time.now().mjd + 20. / (60. * 24.)

        # Initial Actions
        self.mapDirectory = self.settings["gw maps directory"]
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.mapDirectory):
            os.makedirs(self.mapDirectory)

        # INSTANTIATE GRACEDB CLIENT WITH A QUERY STRING
        u = None
        p = None
        # CHECK FOR ROBOT CREDENTIALS IN BREAKER SETTINGS ELSE RELY ON .netrc
        # FILE
        if "graceDB robot credentials" in self.settings:
            u = self.settings["graceDB robot credentials"]["username"]
            p = self.settings["graceDB robot credentials"]["password"]
        self.client = GraceDbBasic(username=u, password=p)

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
                    label="ADVOK & EM_READY",
                    farThreshold=1e-7,
                    startMJD=56658.0,
                    endMJD=False,
                    daemon=False
                )
                downloader.get_maps()
        """
        self.log.info('starting the ``get_maps`` method')

        # VARIABLES
        fileorder = ['LALInference_skymap.fits.gz',
                     'bayestar.fits.gz', 'LIB_skymap.fits.gz', 'skymap.fits.gz', 'LALInference3d.fits.gz', 'bayestar3d.fits.gz']
        stop = False

        # INPUT TIME-VALUES CAN BE SCALAR OR AN ARRAY
        # GET TIME FOR THE VERY START OF LV OPERATIONS
        startOfLV = Time(
            "2015-09-01T00:00:00",
            format='isot',
            scale='utc'
        )

        while stop == False:
            # DAEMON MODE - LISTEN FROM START OF LV OPERATIONS UNTIL NOW + 20
            # MIN
            if self.daemon:
                now = Time.now()
                startGPS = startOfLV.gps
                startUTC = startOfLV.value
                # 20 MINS FROM NOW
                endGPS = now.gps + 1200.
                endUTC = now.isot
            else:
                # NON-DAEMON MODE - DOWNLOAD EVENT MAPS WITHIN THE GIVE
                # TIME-RANGE
                stop = True
                times = [self.startMJD, self.endMJD]
                t = Time(
                    times,
                    format='mjd',
                    scale='utc'
                )
                startGPS = t[0].gps
                endGPS = t[1].gps
                startUTC = t[0].isot
                endUTC = t[1].isot
            self.log.info(
                "checking for events detected between GPS times %(startGPS)s (%(startUTC)s UTC) and %(endGPS)s (%(endUTC)s UTC)" % locals())

            # BUILD THE EVENT QUERY STRING
            label = self.label
            farThreshold = self.farThreshold
            self.eventString = '%(label)s far < %(farThreshold)s %(startGPS)s  .. %(endGPS)s' % locals(
            )
            # QUERY GRACEDB REST API - RETURNS AN ITERATOR
            self.events = self.client.events(
                query=self.eventString,
                orderby=None,
                count=None,
                columns=None
            )

            # ITERATE OVER RESULTS
            oldEvents = 0
            newEvents = 0
            for event in self.events:
                waveId = event['graceid']

                # GET LATEST METADATA FOR THE EVENT FROM GRACEDB
                meta = self._get_event_meta_data(event=event)

                # IF THERE IS NO META IT IS BECOME EVENT ABOVE FAR OR HAD
                # INCORRECT LABELS
                if not meta:
                    continue

                # DETERMINE IF THIS SYSTEM HAS SEEN THE EVENT BEFORE
                newEvent = self._is_the_a_new_event(waveId=event['graceid'])

                if newEvent:
                    newEvents += 1
                    print """NEW GRAVITATIONAL WAVE EVENT FOUND ...
    GraceDB ID: %(waveId)s""" % locals()
                else:
                    oldEvents += 1

                allMaps = []

                # CHECK FOR NEW EVENT SKYMAPS
                maps = {}
                for lvfile in fileorder:
                    matchedFiles = {}
                    files = self.client.files(event['graceid'])
                    for k, v in files.json().iteritems():
                        if lvfile == k.split(",")[0]:
                            matchedFiles[k] = v

                    omatchedFiles = collections.OrderedDict(
                        sorted(matchedFiles.items(), reverse=True))
                    maps[lvfile] = False

                    count = 1

                    for k, v in omatchedFiles.iteritems():
                        count += 1
                        if maps[lvfile] == False:
                            try:
                                aMap = self.client.files(
                                    event['graceid'], k)
                                self._write_map_to_disk(
                                    sMap=aMap,
                                    mapName=lvfile,
                                    waveId=event['graceid']
                                )
                                allMaps.append(aMap)
                                maps[lvfile] = True
                            except:
                                eventId = event['graceid']
                                self.log.error(
                                    "The %(lvfile)s path for %(eventId)s does not seem to exist yet" % locals())

                if len(allMaps) == 0:
                    eventId = event['graceid']
                    self.log.warning(
                        'cound not download skymaps for event %(eventId)s' % locals())

                # DUMP THE KNOWN EVENT METADATA BESIDE MAPS
                meta["Maps"] = maps
                fileName = self.mapDirectory + "/" + waveId + "/meta.yaml"
                stream = file(fileName, 'w')
                yaml.dump(meta, stream, default_flow_style=False)
                stream.close()

                # PRINT METADATA TO SCREEN IF THIS IS THE FIRST TIME THIS
                # SYSTEM HAS SEEN THE EVENT
                if newEvent:
                    try:
                        self.log.debug(
                            "attempting to open the file %s" % (fileName,))
                        readFile = codecs.open(
                            fileName, encoding='utf-8', mode='r')
                        thisData = readFile.read()
                        readFile.close()
                    except IOError, e:
                        message = 'could not open the file %s' % (
                            fileName,)
                        self.log.critical(message)
                        raise IOError(message)

                    print "\nMETADATA FOR %(waveId)s ..." % locals()
                    print thisData
                    readFile.close()

            if stop == False:
                freq = self.daemon
                print "%(oldEvents)s archived and %(newEvents)s events found, will try again in %(freq)s secs" % locals()
                time.sleep(freq)

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
        mapPath = outputDir + "/" + mapName

        if not os.path.exists(mapPath):
            print "NEW MAP FOUND FOR GW EVENT %(waveId)s ... " % locals()
            print "    Downloading %(mapName)s" % locals()
            skymapfile = open(mapPath, 'w')
            skymapfile.write(sMap.read())
            skymapfile.close()
        else:
            self.log.info("%(mapName)s has already been downloaded" % locals())

        self.log.info('completed the ``_write_map_to_disk`` method')
        return None

    def _get_event_meta_data(
            self,
            event):
        """*query graceDB and parse the event metadata into dictionary*

        **Key Arguments:**
            - ``event`` -- the event data from graceDB

        **Return:**
            - ``None`` or ``meta`` -- a dictionary of event metadata pulled from graceDB, or none if event doesn't match our requirements
        """
        self.log.info('starting the ``_get_event_meta_data`` method')

        eventKeys = ['graceid', 'gpstime', 'group', 'links', 'created',
                     'far', 'instruments', 'labels', 'nevents', 'submitter', 'search', 'likelihood']
        eventinfo = {}
        mjds = [-1, -1]
        timediff = -1

        # CHECK ALL THE EVENT KEYS EXIST - ELSE PASS ON THIS EVENT
        for key in eventKeys:
            if not key in event:
                self.log.info(
                    "`%(key)s` not in event %(event)s" % locals())
                return None
            eventinfo[key] = event[key]
        eventinfo['gpstime'] = float(eventinfo['gpstime'])

        # INJ : EVENT RESULTS FROM AN INJECTION (FAKE) - PASS
        if "INJ" in eventinfo['labels']:
            self.log.info(
                "event %(event)s has an INJ label" % locals())
            return None

        # QUERY THE FALSE ALARM RATE. FAR TOO HIGH == PASS
        if eventinfo['far'] > self.farThreshold:
            far = eventinfo['far']
            farthres = self.farThreshold
            self.log.info(
                "event %(event)s does not pass FAR threshold of %(farthres)s. (FAR = %(far)s)" % locals())
            return None

        self.log.info("Getting info for %s" % event["graceid"])

        if 'extra_attributes' in event:
            # LOOK UP THE EXTRA-ATTRIBUTES VALUE (NEED TO BE LV-MEMBER TO VIEW)
            if 'SingleInspiral' in event['extra_attributes']:
                eventinfo['singles'] = {}
                for single in event['extra_attributes']['SingleInspiral']:
                    eventinfo['singles'][single['ifo']] = single
                    eventinfo['singles'][single['ifo']]['gpstime'] = single[
                        'end_time'] + 10**-9 * single['end_time_ns']

                if ("H1" in eventinfo['singles']) and ("L1" in eventinfo['singles']):
                    eventinfo["H1_L1_difference"] = eventinfo['singles']['H1'][
                        "gpstime"] - eventinfo['singles']['L1']["gpstime"]
                    t = Time([eventinfo['singles']['H1']["gpstime"], eventinfo[
                             'singles']['L1']["gpstime"]], format='gps', scale='utc')
                    mjds = t.mjd
                    timediff = eventinfo["H1_L1_difference"]

            try:
                self.log.debug("Looking for cWB file for the event %s" %
                               (eventinfo['graceid'],))
                # READ THIS TRIGGER FILE FROM GRACEDB
                r = self.client.files(eventinfo['graceid'],
                                      "trigger_%.4f.txt" % eventinfo['gpstime'])
                exists = True
            except Exception, e:
                self.log.info("No cWB file found for the event %s" %
                              (eventinfo['graceid'],))
                exists = False

            if exists:
                cwbfile = open('/tmp/trigger.txt', 'w')
                cwbfile.write(r.read())
                cwbfile.close()

                eventinfo['burst'] = {}
                lines = [line.rstrip('\n')
                         for line in open('/tmp/trigger.txt')]
                for line in lines:
                    lineSplit = line.split(":")
                    if len(lineSplit) < 2:
                        continue
                    key = lineSplit[0]
                    value = filter(None, lineSplit[1].split(" "))
                    eventinfo['burst'][lineSplit[0]] = value

                ifo1 = eventinfo['burst']['ifo'][0]
                gps1 = float(eventinfo['burst']['time'][0])

                ifo2 = eventinfo['burst']['ifo'][1]
                gps2 = float(eventinfo['burst']['time'][1])

                eventinfo['burst'][ifo1] = {}
                eventinfo['burst'][ifo1]['gpstime'] = gps1

                eventinfo['burst'][ifo2] = {}
                eventinfo['burst'][ifo2]['gpstime'] = gps2

                if ("H1" in eventinfo['burst']) and ("L1" in eventinfo['burst']):
                    eventinfo["H1_L1_difference"] = eventinfo['burst']['H1'][
                        "gpstime"] - eventinfo['burst']['L1']["gpstime"]
                    t = Time([eventinfo['burst']['H1']["gpstime"], eventinfo[
                        'burst']['L1']["gpstime"]], format='gps', scale='utc')
                    mjds = t.mjd
                    timediff = eventinfo["H1_L1_difference"]

        # FILL IN META DICTIONARY
        meta = {}
        meta["GraceDB ID"] = str(event['graceid'])
        meta["GPS Event Time"] = event["gpstime"]
        meta["Discovery Group"] = str(event["group"])
        meta["Detection Pipeline"] = str(event["pipeline"])
        meta["Discovery Search Type"] = str(event["search"])
        # meta["Links"] = event["links"]
        meta["Date Added to GraceDB"] = str(event["created"])
        meta["False Alarm Rate"] = str(event["far"]) + " Hz"
        meta["Event Submitter"] = str(event["submitter"])
        meta["Detection Interferometers"] = str(event["instruments"])
        if 'extra_attributes' in event:
            meta["Hanford MJD"] = float("%.10f" % (mjds[0],))
            meta["Livingston MJD"] = float("%.10f" % (mjds[1],))
            meta["MJD Difference Seconds"] = float("%.10f" % (timediff,))

        self.log.info('completed the ``_get_event_meta_data`` method')
        return meta

    def _is_the_a_new_event(
            self,
            waveId):
        """*determine if this system has seen this event yet*

        **Key Arguments:**
            - ``waveId`` -- the graceDB id for the wave.

        **Return:**
            - ``newEvent`` -- is this a new event (Boolean)
        """
        self.log.info('starting the ``_is_the_a_new_event`` method')

        newEvent = False
        outputDir = self.mapDirectory + "/" + waveId.upper()
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outputDir):
            newEvent = True
            os.makedirs(outputDir)

        self.log.info('completed the ``_is_the_a_new_event`` method')
        return newEvent

    # use the tab-trigger below for new method
    # xt-class-method
