#!/usr/local/bin/python
# encoding: utf-8
"""
*Plot the maps showing the timeline of observations of the gravitation wave sky-locations*

:Author:
    David Young

:Date Created:
    November  5, 2015
"""
################# GLOBAL IMPORTS ####################
# SUPPRESS MATPLOTLIB WARNINGS
from __future__ import unicode_literals
import warnings
warnings.filterwarnings("ignore")
import sys
import os
os.environ['TERM'] = 'vt100'
import healpy as hp
import numpy as np
import math
from datetime import datetime
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
# from shapely.geometry import Point
# from shapely.ops import cascaded_union
from matplotlib.path import Path
from matplotlib.pyplot import savefig
import matplotlib.patches as patches
import matplotlib.path as mpath
from matplotlib.projections.geo import GeoAxes
from astropy import wcs as awcs
from astropy.io import fits
from fundamentals import tools, times
from fundamentals.mysql import readquery
from crowdedText import adjust_text
from astrocalc.times import now


class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""

    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2 * np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)


class plot_wave_observational_timelines():
    """
    *TPlot the maps showing the timeline of observations of the gravitation wave sky-locations*

    You can plot either the history (looking back from now) or timeline (looking forward from date of GW detection) of the survey.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``plotType`` -- history (looking back from now) or timeline (looking forward from date of GW detection)
        - ``gwid`` -- a given graviational wave ID. If given only maps for this wave shall be plotted. Default *False* (i.e. plot all waves)
        - ``projection`` -- projection for the plot. Default *mercator* [mercator|gnomonic|mollweide]
        - ``probabilityCut`` -- remove footprints where probability assigned to the healpix pixel found at the center of the exposure is ~0.0. Default *False*
        - ``databaseConnRequired`` -- are the database connections going to be required? Default *True*
        - ``allPlots`` -- plot all timeline plot (including the CPU intensive -21-0 days and all transients/footprints plots). Default *False*
        - ``telescope`` -- select an individual telescope. Default *False*. [ps1|atlas]
        - ``timestamp`` -- add a timestamp to the plot to show when it was created. Default *True*
        - ``filters`` -- only plot certain filters. Default *False*

    **Usage:**

        To plot a the history of a specific wave:

        .. code-block:: python

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="history",
                   gwid="G184098",
                   projection="mercator"
            )
            plotter.get()

        or to plot all waves in the settings file with a mollweide projection:

        .. code-block:: python

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="history",
                   projection="mollweide"
            )
            plotter.get()

        To plot the timeline of the survey, change ``plotType="timeline"``:

        .. code-block:: python

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="timeline"
            )
            plotter.get()
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            plotType=False,
            gwid=False,
            projection="mercator",
            probabilityCut=False,
            databaseConnRequired=True,
            allPlots=False,
            telescope=False,
            timestamp=True,
            filters=False
    ):
        self.log = log
        log.debug("instantiating a new 'plot_wave_observational_timelines' object")
        self.settings = settings
        self.plotType = plotType
        self.gwid = gwid
        self.projection = projection
        self.probabilityCut = probabilityCut
        self.allPlots = allPlots
        self.telescope = telescope
        self.timestamp = timestamp
        self.filters = filters

        # xt-self-arg-tmpx

        # Initial Actions

        if self.settings and databaseConnRequired:
            # CONNECT TO THE VARIOUS DATABASES REQUIRED
            from breaker import database
            db = database(
                log=self.log,
                settings=self.settings
            )
            self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn, self.atlasDbConn, self.ps13piDbConn = db.get()
        else:
            self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = False, False, False

        self.log.debug(
            'connected to databases')

        return None

    def get(self):
        """
        *Generate the plots*
        """
        self.log.debug('starting the ``get`` method')

        if self.plotType == "history":
            self.get_history_plots()
        elif self.plotType == "timeline":
            self.get_timeline_plots()

        self.log.debug('completed the ``get`` method')
        return None

    def get_gw_parameters_from_settings(
            self,
            gwid,
            inPastDays=False,
            inFirstDays=False,
            maxProbCoordinate=[0, 0],
            stackOnly=False):
        """
        *Query the settings file and database for PS1 Pointings, PS1 discovered transients and plot parameters relatiing to the given gravitational wave (``gwid``)*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection). A tuple (start day, end day).
            - ``maxProbCoordinate`` -- the sky-coordiante of the pixel containing the highest likelihood (calculated from map).

        **Return:**
            - ``plotParameters`` -- the parameters used for the plots
            - ``ps1Transients`` -- the transients to add to the plot
            - ``ps1Pointings`` -- the pointings to place on the plot

        **Usage:**

            .. code-block:: python

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                       settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
                    gwid="G211117"
                )
                print plotParameters

                # OUT: {'raRange': 90.0, 'centralCoordinate': [55.0, 27.5],
                # 'decRange': 85.0}

                print ps1Transients

                # OUT: ({'local_designation': u'5L3Gbaj', 'ra_psf':
                # 39.29767123836419, 'ps1_designation': u'PS15don', 'dec_psf':
                # 19.055638423458053}, {'local_designation': u'5L3Gcbu',
                # 'ra_psf': 40.06271352712189, 'ps1_designation': u'PS15dox',
                # 'dec_psf': 22.536709765810823}, {'local_designation':
                # u'5L3Gcca', 'ra_psf': 41.97569854977185, 'ps1_designation':
                # u'PS15doy', 'dec_psf': 21.773344501616435},
                # {'local_designation': u'5L3Gbla', 'ra_psf':
                # 50.732664347994714, 'ps1_designation': u'PS15dcq', 'dec_psf':
                # 34.98988923347591}, {'local_designation': u'6A3Gcvu',
                # 'ra_psf': 34.77565307934415, 'ps1_designation': u'PS16ku',
                # 'dec_psf': 10.629310832257824}, {'local_designation':
                # u'5L3Gcel', 'ra_psf': 38.24898916543392, 'ps1_designation':
                # u'PS15dpn', 'dec_psf': 18.63530332013424},
                # {'local_designation': u'5L3Gcvk', 'ra_psf':
                # 40.13754684778398, 'ps1_designation': u'PS15dpz', 'dec_psf':
                # 23.003023065333267}, ....

                print ps1Pointings

                # OUT: ({'raDeg': 37.1814041667, 'mjd': 57388.2124067,
                # 'decDeg': 18.9258969444}, {'raDeg': 37.1813666667, 'mjd':
                # 57388.2140101, 'decDeg': 18.9259066667}, ...

            It can also be useful to give time-limits for the request to get the observations and discoveries from the past few days (``inPastDays``), or for the first few days after wave detection (``inFirstDays``). So for the past week:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inPastDays=7
                )

            Or the first 3 days since wave detection:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inFirstDays=(0,3)
                )
        """
        self.log.debug(
            'completed the ````get_gw_parameters_from_settings`` method')

        plotParameters = self.settings["gravitational waves"][gwid]["plot"]
        if not plotParameters:
            plotParameters = {}

        if "centralCoordinate" not in plotParameters:
            plotParameters["centralCoordinate"] = maxProbCoordinate

        if "raRange" not in plotParameters:
            plotParameters["centralCoordinate"] = list(
                plotParameters["centralCoordinate"])
            plotParameters["centralCoordinate"][1] = 0.

        # GRAB PS1 TRANSIENTS FROM THE DATABASE
        ps1Transients, atlasTransients = self._get_ps1_atlas_transient_candidates(
            gwid=gwid,
            mjdStart=self.settings["gravitational waves"][
                gwid]["mjd"],
            mjdEnd=self.settings["gravitational waves"][
                gwid]["mjd"] + 31.,
            plotParameters=plotParameters,
            inPastDays=inPastDays,
            inFirstDays=inFirstDays,
            maxProbCoordinate=maxProbCoordinate
        )

        self.log.debug(
            'finished getting the PS1 transients')

        # GRAB PS1 & ATLAS POINTINGS FROM THE DATABASE
        ps1Pointings = self._get_ps1_pointings(
            gwid, inPastDays, inFirstDays, stackOnly=stackOnly)
        atlasPointings = self._get_atlas_pointings(
            gwid, inPastDays, inFirstDays)

        self.log.debug(
            'completed the ``get_gw_parameters_from_settings`` method')

        if self.telescope == "ps1":
            atlasPointings = []
            atlasTransients = []
        elif self.telescope == "atlas":
            ps1Transients = []
            ps1Pointings = []

        return plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients

    def _get_ps1_atlas_transient_candidates(
            self,
            gwid,
            mjdStart,
            mjdEnd,
            plotParameters,
            inPastDays,
            inFirstDays,
            maxProbCoordinate):
        """
        *get ps1 transient candidates*

        **Key Arguments:**
            - ``gwid`` -- the gravitational wave id
            - ``mjdStart`` -- earliest mjd of discovery
            - ``mjdEnd`` -- latest mjd of discovery
            - ``plotParameters`` -- the parameters of the plot (for spatial & temporal parameters etc)
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection).  A tuple (start day, end day).

        **Return:**
            - ``ps1Transients`` -- the transients to add to the plot
        """
        self.log.debug(
            'completed the ````_get_ps1_atlas_transient_candidates`` method')

        # UNPACK THE PLOT PARAMETERS
        if "centralCoordinate" in plotParameters:
            centralCoordinate = plotParameters["centralCoordinate"]
        else:
            centralCoordinate = maxProbCoordinate

        if "raRange" not in plotParameters:
            raRange = 360.
            decRange = 180.
            raMax = 360.
            raMin = 0.
            decMax = 90.
            decMin = -90.
        else:
            raRange = plotParameters["raRange"]
            decRange = plotParameters["decRange"]
            raMax = centralCoordinate[0] + raRange / 2.
            raMin = centralCoordinate[0] - raRange / 2.
            decMax = centralCoordinate[1] + decRange / 2.
            decMin = centralCoordinate[1] - decRange / 2.

        if inPastDays:
            nowMjd = now(
                log=self.log
            ).get_mjd()
            mjdStart = nowMjd - inPastDays
            mjdEnd = 1000000000000000

        if inFirstDays:
            mjdStart = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[1]
            if inFirstDays[1] == 0 and inFirstDays[0] == 0:
                mjdEnd = 10000000000

        if raMin >= 0 and raMax <= 360.:
            sqlQuery = u"""
                SELECT ps1_designation, local_designation, ra_psf, dec_psf FROM tcs_transient_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and (ra_psf between %(raMin)s and %(raMax)s) and (`dec_psf` between %(decMin)s and %(decMax)s) ;
            """ % locals()
        elif raMin < 0.:
            praMin = 360. + raMin
            sqlQuery = u"""
                SELECT ps1_designation, local_designation, ra_psf, dec_psf FROM tcs_transient_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and ((ra_psf between %(praMin)s and 360.) or (ra_psf between 0. and %(raMax)s)) and (`dec_psf` between %(decMin)s and %(decMax)s) ;
            """ % locals()
        elif raMax > 360.:
            praMax = raMax - 360.
            sqlQuery = u"""
                SELECT ps1_designation, local_designation, ra_psf, dec_psf FROM tcs_transient_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and ((ra_psf between %(raMin)s and 360.) or (ra_psf between 0. and %(praMax)s)) and (`dec_psf` between %(decMin)s and %(decMax)s) ;
            """ % locals()

        ps1Transients = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ps1gwDbConn
        )

        if raMin >= 0 and raMax <= 360.:
            sqlQuery = u"""
                SELECT atlas_designation, ra, `dec` FROM atlas_diff_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and (ra between %(raMin)s and %(raMax)s) and (`dec` between %(decMin)s and %(decMax)s) ;
            """ % locals()
        elif raMin < 0.:
            araMin = 360. + raMin
            sqlQuery = u"""
                SELECT atlas_designation, ra, `dec` FROM atlas_diff_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and ((ra between %(araMin)s and 360.) or (ra between 0. and %(raMax)s)) and (`dec` between %(decMin)s and %(decMax)s) ;
            """ % locals()
        elif raMax > 360.:
            araMax = raMax - 360.
            sqlQuery = u"""
                SELECT atlas_designation, ra, `dec` FROM atlas_diff_objects o, tcs_latest_object_stats s where o.detection_list_id in (1,2) and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and ((ra between %(raMin)s and 360.) or (ra between 0. and %(araMax)s)) and (`dec` between %(decMin)s and %(decMax)s) ;
            """ % locals()

        atlasTransients = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.atlasDbConn
        )

        self.log.debug(
            'completed the ``_get_ps1_atlas_transient_candidates`` method')
        return ps1Transients, atlasTransients

    def _get_ps1_pointings(
            self,
            gwid,
            inPastDays,
            inFirstDays,
            stackOnly=False):
        """
        *get ps1 pointings to add to the plot*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection). A tuple (start day, end day).

        **Return:**
            - ``ps1Pointings`` -- the pointings to place on the plot
        """
        self.log.debug('starting the ``_get_ps1_pointings`` method')

        # DETERMINE THE TEMPORAL CONSTRAINTS FOR MYSQL QUERY
        if inPastDays != False or inPastDays == 0:
            nowMjd = now(
                log=self.log
            ).get_mjd()
            mjdStart = nowMjd - inPastDays
            mjdEnd = 10000000000
            if inPastDays == 0:
                mjdStart = 0.0

        if inFirstDays:
            mjdStart = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[1]
            if inFirstDays[1] == 0 and inFirstDays[0] == 0:
                mjdEnd = 10000000000

        if inPastDays == False and inFirstDays == False:
            mjdStart = self.settings["gravitational waves"][gwid]["mjd"]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["mjd"] + 31.

        sqlQuery = u"""
            SELECT raDeg, decDeg, mjd, exp_time, filter, limiting_mag FROM ps1_pointings where gw_id like "%%%(gwid)s%%" and mjd between %(mjdStart)s and %(mjdEnd)s
        """ % locals()

        sqlQuery = u"""
            SELECT raDeg, decDeg, mjd_registered as mjd, etime as exp_time, f as filter FROM ps1_nightlogs where gw_id like "%%%(gwid)s%%" and mjd_registered between %(mjdStart)s and %(mjdEnd)s
        """ % locals()

        if self.filters:
            filters = ' and filter in ("' + ('", "').join(self.filters) + '")'
        else:
            filters = ""

        if not stackOnly:

            sqlQuery = u"""
                SELECT distinct * from (
    SELECT distinct raDeg, decDeg, mjd, filter, p.skycell_id as exp_id, exp_time, limiting_mag as limiting_magnitude FROM ps1_stack_stack_diff_skycells p, ps1_skycell_gravity_event_annotations s, ps1_skycell_map m where mjd between %(mjdStart)s and %(mjdEnd)s and p.skycell_id=s.skycell_id and p.skycell_id=m.skycell_id and gracedb_id = "%(gwid)s" %(filters)s and s.prob_coverage > 1e-6
    UNION
    SELECT distinct raDeg, decDeg, mjd, filter, p.skycell_id as exp_id, exp_time, limiting_mag as limiting_magnitude FROM ps1_warp_stack_diff_skycells p, ps1_skycell_gravity_event_annotations s, ps1_skycell_map m where mjd between %(mjdStart)s and %(mjdEnd)s and p.skycell_id=s.skycell_id and p.skycell_id=m.skycell_id and gracedb_id = "%(gwid)s" %(filters)s and s.prob_coverage > 1e-6) p order by mjd;
            """ % locals()

        else:

            sqlQuery = u"""
                SELECT distinct raDeg, decDeg, mjd, filter, p.skycell_id as exp_id, exp_time, limiting_mag as limiting_magnitude FROM ps1_stack_stack_diff_skycells p, ps1_skycell_gravity_event_annotations s, ps1_skycell_map m where mjd between %(mjdStart)s and %(mjdEnd)s and p.skycell_id=s.skycell_id and p.skycell_id=m.skycell_id and gracedb_id = "%(gwid)s" %(filters)s and s.prob_coverage > 1e-6
            """ % locals()

        ps1Pointings = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.debug('completed the ``_get_ps1_pointings`` method')
        return ps1Pointings

    def _get_atlas_pointings(
            self,
            gwid,
            inPastDays,
            inFirstDays):
        """
        *get atlas pointings to add to the plot*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection).  A tuple (start day, end day).

        **Return:**
            - ``atlasPointings`` -- the pointings to place on the plot
        """
        self.log.debug('starting the ``_get_atlas_pointings`` method')

        # DETERMINE THE TEMPORAL CONSTRAINTS FOR MYSQL QUERY
        if inPastDays != False or inPastDays == 0:
            nowMjd = now(
                log=self.log
            ).get_mjd()
            mjdStart = nowMjd - inPastDays
            mjdEnd = 10000000000
            if inPastDays == 0:
                mjdStart = 0.0

        if inFirstDays:
            mjdStart = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["mjd"] + inFirstDays[1]
            if inFirstDays[1] == 0 and inFirstDays[0] == 0:
                mjdEnd = 10000000000

        if inPastDays == False and inFirstDays == False:
            mjdStart = self.settings["gravitational waves"][gwid]["mjd"]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["mjd"] + 31.

        sqlQuery = u"""
            SELECT p.atlas_object_id as exp_id, p.raDeg, p.decDeg, mjd, exp_time, filter, limiting_magnitude FROM atlas_pointings p, atlas_exposure_gravity_event_annotations a where a.prob_coverage > 1e-3 and a.gracedb_id = "%(gwid)s" and a.atlas_object_id=p.atlas_object_id and gw_id like "%%%(gwid)s%%" and mjd between %(mjdStart)s and %(mjdEnd)s order by mjd;
        """ % locals()

        atlasPointings = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.debug('completed the ``_get_atlas_pointings`` method')
        return atlasPointings

    def generate_probability_plot(
            self,
            pathToProbMap,
            gwid,
            mjdStart=False,
            timeLimitLabel=False,
            timeLimitDay=False,
            fileFormats=["pdf"],
            folderName="",
            plotType="timeline",
            plotParameters=False,
            ps1Transients=[],
            atlasTransients=[],
            ps1Pointings=[],
            atlasPointings=[],
            projection="mercator",
            raLimit=False,
            probabilityCut=False,
            outputDirectory=False,
            fitsImage=False,
            allSky=False,
            center=False,
            symLink=True):
        """
        *Generate a single probability map plot for a given gravitational wave and save it to file*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``plotParameters`` -- the parameters of the plot (for spatial & temporal parameters etc).
            - ``ps1Transients`` -- the PS1 transients to add to the plot. Default **[]**
            - ``atlasTransients`` -- the ATLAS transients to add to the plot. Default **[]**
            - ``ps1Pointings`` -- the PS1 pointings to place on the plot. Default **[]**
            - ``atlasPointings`` -- the atlas pointings to add to the plot. Default **[]**
            - ``pathToProbMap`` -- path to the FITS file containing the probability map of the wave
            - ``mjdStart`` -- earliest mjd of discovery
            - ``timeLimitLabel`` -- the labels of the time contraints (for titles)
            - ``timeLimitDay`` -- the time limits (in ints)
            - ``raLimit`` -- ra limit at twilight
            - ``fileFormats`` -- the format(s) to output the plots in (list of strings) Default **["pdf"]**
            - ``folderName`` -- the name of the folder to add the plots to
            - ``plotType`` -- history (looking back from now) or timeline (looking forward from date of GW detection). Default **timeline**
            - ``projection`` -- projection for the plot. Default *mercator*. [mercator|mollweide|gnomonic]
            - ``probabilityCut`` -- remove footprints where probability assigned to the healpix pixel found at the center of the exposure is ~0.0. Default *False*
            - ``outputDirectory`` -- can be used to override the output destination in the settings file
            - ``fitsImage`` -- generate a FITS image file of map
            - ``allSky`` -- generate an all-sky map (do not use the ra, dec window in the breaker settings file). Default *False*
            - ``center`` -- central longitude in degrees. Default *0*.


        **Return:**
            - None

        **Usage:**

            First you neeed to collect your data and a few plot parameters:

            .. code-block:: python

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inFirstDays=(0,7)
                )

            Then you can pass in these parameter to generate a plot:

            .. code-block:: python

                plotter.generate_probability_plot(
                    gwid="G211117",
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    atlasTransient=atlasTransient,
                    ps1Pointings=ps1Pointings,
                    atlasPointings=altasPointings,
                    pathToProbMap="/Users/Dave/config/breaker/maps/G211117/LALInference_skymap.fits",
                    mjdStart=57382.,
                    timeLimitLabel="",
                    timeLimitDay=(0, 5),
                    raLimit=False,
                    fileFormats=["pdf"],
                    folderName="survey_timeline_plots",
                    projection="mercator",
                    plotType="timeline",
                    probabilityCut=True,
                    outputDirectory=False
                )

        """
        self.log.debug('starting the ``generate_probability_plot`` method')

        import matplotlib.pyplot as plt
        import healpy as hp
        from matplotlib.font_manager import FontProperties
        font = FontProperties()

        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        # VARIABLES
        font.set_family("Arial")
        pixelSizeDeg = 0.066667
        unit = "likelihood"
        cmap = "YlOrRd"
        colorBar = False

        # DETERMINE PATH TO MAP AND IF IT IS THE PREFERRED MAP AT THIS TIME
        bestMap = False
        if pathToProbMap is None or not pathToProbMap:
            pathToProbMap = self.settings[
                "gravitational waves"][gwid]["mapPath"]

        pathToProbMap = os.path.abspath(pathToProbMap)

        if gwid in self.settings["gravitational waves"]:
            bestMapPath = self.settings[
                "gravitational waves"][gwid]["mapPath"]
            bestMapPath = os.path.abspath(bestMapPath)
            if pathToProbMap == bestMapPath:
                bestMap = True

        mapBasename = os.path.basename(pathToProbMap)
        mapBasename = os.path.splitext(mapBasename)[0]
        mapBasename = os.path.splitext(mapBasename)[0]
        mapBasename = os.path.splitext(mapBasename)[0]

        # INITIALISE FIGURE
        fig = plt.figure()

        # READ HEALPIX MAPS FROM FITS FILE
        # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
        # ARRAY OF PROBABILITIES (3,072 ROWS)
        # READ IN THE HEALPIX FITS FILE
        aMap, mapHeader = hp.read_map(pathToProbMap, 0, h=True, verbose=False)
        # DETERMINE THE SIZE OF THE HEALPIXELS
        nside = hp.npix2nside(len(aMap))

        # MIN-MAX PROB VALUES TO ADJUST MAP CONTRAST
        vmin = min(aMap)
        vmax = max(aMap) * 0.9

        totalProb = sum(aMap)
        # print "Total Probability for the entire sky is %(totalProb)s" %
        # locals()

        # UNPACK THE PLOT PARAMETERS
        # FIND THE COORDINATES OF THE CORE LIKEIHOOD
        maxProbHealpix = aMap.argmax()
        maxCoordinate = hp.pix2ang(nside, maxProbHealpix, lonlat=True)
        print "The %(gwid)s %(mapBasename)s map's maximum likelihood is centered at %(maxCoordinate)s" % locals()

        if center == False:
            center = maxCoordinate[0]

        if plotParameters:
            centralCoordinate = plotParameters["centralCoordinate"]
        else:
            centralCoordinate = [center, 0]

        # CREATE A NEW WCS OBJECT
        w = awcs.WCS(naxis=2)
        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = np.array([pixelSizeDeg, pixelSizeDeg])
        # WORLD COORDINATES AT REFERENCE PIXEL
        w.wcs.crval = centralCoordinate

        projectionDict = {
            "mollweide": "MOL",
            "gnomonic": "MER",
            "mercator": "MER",
            "cartesian": "CAR"
        }

        if projection in ["mollweide"]:
            # MAP VISULISATION RATIO IS ALWAYS 1/2
            xRange = 2000
            yRange = xRange / 2.

            # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
            w.wcs.crpix = [xRange / 2., yRange / 2.]

            # FULL-SKY MAP SO PLOT FULL RA AND DEC RANGES
            # DEC FROM 180 to 0
            theta = np.linspace(np.pi, 0, yRange)

            latitude = np.radians(np.linspace(-90, 90, yRange))
            # RA FROM -180 to +180
            phi = np.linspace(-np.pi, np.pi, xRange)
            longitude = np.radians(
                np.linspace(-180, 180, xRange))
            X, Y = np.meshgrid(longitude, latitude)

            # PROJECT THE MAP TO A RECTANGULAR MATRIX xRange X yRange
            PHI, THETA = np.meshgrid(phi, theta)
            healpixIds = hp.ang2pix(nside, THETA, PHI)
            probs = aMap[healpixIds]

            # healpixIds = np.reshape(healpixIds, (1, -1))[0]

            # CTYPE FOR THE FITS HEADER
            thisctype = projectionDict[projection]
            w.wcs.ctype = ["RA---%(thisctype)s" %
                           locals(), "DEC--%(thisctype)s" % locals()]

            # ALL PROJECTIONS IN FITS SEEM TO BE MER
            w.wcs.ctype = ["RA---MER" %
                           locals(), "DEC--MER" % locals()]

            stampProb = np.sum(aMap)
            print "Probability for the plot stamp is %(stampProb)s" % locals()

            # MATPLOTLIB IS DOING THE PROJECTION
            ax = fig.add_subplot(111, projection=projection)

            # RASTERIZED MAKES THE MAP BITMAP WHILE THE LABELS REMAIN VECTORIAL
            # FLIP LONGITUDE TO THE ASTRO CONVENTION
            image = ax.pcolormesh(longitude[
                ::-1], latitude, probs, rasterized=True, cmap=cmap)

            # GRATICULE
            ax.set_longitude_grid(30)
            ax.set_latitude_grid(15)
            ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(30))
            ax.set_longitude_grid_ends(90)

            # CONTOURS - NEED TO ADD THE CUMMULATIVE PROBABILITY
            i = np.flipud(np.argsort(aMap))
            cumsum = np.cumsum(aMap[i])
            cls = np.empty_like(aMap)
            cls[i] = cumsum * 99.99999999 * stampProb

            # EXTRACT CONTOUR VALUES AT HEALPIX INDICES
            contours = []
            contours[:] = [cls[i] for i in healpixIds]
            # contours = np.reshape(np.array(contours), (yRange, xRange))

            CS = ax.contour(longitude[::-1], latitude,
                            contours, linewidths=.5, alpha=0.7, zorder=2)

            CS.set_alpha(0.5)
            CS.clabel(fontsize=10, inline=True,
                      fmt='%2.1f', fontproperties=font, alpha=0.0)

            # COLORBAR
            if colorBar:
                cb = fig.colorbar(image, orientation='horizontal',
                                  shrink=.6, pad=0.05, ticks=[0, 1])
                cb.ax.xaxis.set_label_text("likelihood")
                cb.ax.xaxis.labelpad = -8
                # WORKAROUND FOR ISSUE WITH VIEWERS, SEE COLORBAR DOCSTRING
                cb.solids.set_edgecolor("face")

            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)
            # lon.set_ticks_position('bt')
            # lon.set_ticklabel_position('b')
            # lon.set_ticklabel(size=20)
            # lat.set_ticklabel(size=20)
            # lon.set_axislabel_position('b')
            # lat.set_ticks_position('lr')
            # lat.set_ticklabel_position('l')
            # lat.set_axislabel_position('l')

            # # REMOVE TICK LABELS
            # ax.xaxis.set_ticklabels([])
            # ax.yaxis.set_ticklabels([])
            # # REMOVE GRID
            # ax.xaxis.set_ticks([])
            # ax.yaxis.set_ticks([])

            # REMOVE WHITE SPACE AROUND FIGURE
            spacing = 0.01
            plt.subplots_adjust(bottom=spacing, top=1 - spacing,
                                left=spacing, right=1 - spacing)

            plt.grid(True)

        elif projection in ["mercator", "cartesian"]:

            if allSky or "raRange" not in plotParameters:
                centralCoordinate = list(centralCoordinate)
                centralCoordinate[1] = 0.
                raRange = 360.
                decRange = 180.
            else:
                # UNPACK THE PLOT PARAMETERS
                raRange = plotParameters["raRange"]
                decRange = plotParameters["decRange"]

            raMax = centralCoordinate[0] + raRange / 2.
            raMin = centralCoordinate[0] - raRange / 2.
            decMax = centralCoordinate[1] + decRange / 2.
            decMin = centralCoordinate[1] - decRange / 2.

            # DETERMINE THE PIXEL GRID X,Y RANGES
            xRange = int(raRange / pixelSizeDeg)
            yRange = int(decRange / pixelSizeDeg)
            if projection == "mercator" and allSky:
                # yRange = yRange
                yRange = yRange * 2
            largest = max(xRange, yRange)
            # xRange = largest
            # yRange = largest

            # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
            w.wcs.crpix = [xRange / 2., yRange / 2.]

            # FOR AN ORTHOGONAL GRID THE CRPIX2 VALUE MUST BE ZERO AND CRPIX2
            # MUST REFLECT THIS
            w.wcs.crpix[1] -= w.wcs.crval[1] / w.wcs.cdelt[1]
            w.wcs.crval[1] = 0

            # USE THE "GNOMONIC" PROJECTION ("COORDINATESYS---PROJECTION")
            ctype = projectionDict[projection]
            w.wcs.ctype = ["RA---" + ctype, "DEC--" + ctype]

            # CREATE A PIXEL GRID - 2 ARRAYS OF X, Y
            columns = []
            px = np.tile(np.arange(0, xRange), yRange)
            py = np.repeat(np.arange(0, yRange), xRange)

            # CONVERT THE PIXELS TO WORLD COORDINATES
            wr, wd = w.wcs_pix2world(px, py, 1)

            # MAKE SURE RA IS +VE
            nr = []
            nr[:] = [r if r > 0 else r + 360. for r in wr]
            wr = np.array(nr)

            # THETA: IS THE POLAR ANGLE, RANGING FROM 0 AT THE NORTH POLE TO PI AT THE SOUTH POLE.
            # PHI: THE AZIMUTHAL ANGLE ON THE SPHERE FROM 0 TO 2PI
            # CONVERT DEC TO THE REQUIRED HEALPIX FORMAT
            nd = -wd + 90.

            # CONVERT WORLD TO HEALPIX INDICES
            healpixIds = hp.ang2pix(nside, theta=nd * DEG_TO_RAD_FACTOR,
                                    phi=wr * DEG_TO_RAD_FACTOR)

            # NOW READ THE VALUES OF THE MAP AT THESE HEALPIX INDICES
            uniqueHealpixIds = np.unique(healpixIds)
            probs = []
            probs[:] = [aMap[i] for i in healpixIds]

            uniProb = []
            uniProb[:] = [aMap[i] for i in uniqueHealpixIds]

            stampProb = np.sum(uniProb)
            print "Probability for the plot stamp is %(stampProb)s" % locals()

            # RESHAPE THE ARRAY AS BITMAP
            probs = np.reshape(np.array(probs), (yRange, xRange))

            # CREATE THE FITS HEADER WITH WCS
            header = w.to_header()
            # CREATE THE FITS FILE
            hdu = fits.PrimaryHDU(header=header, data=probs)

            # GRAB THE WCS FROM HEADER GENERATED EARLIER
            from astropy.wcs import WCS
            from astropy.visualization.wcsaxes import WCSAxes

            wcs = WCS(hdu.header)
            # USE WCS AS THE PROJECTION
            ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=wcs)
            # note that the axes have to be explicitly added to the figure
            ax = fig.add_axes(ax)

            # PLOT MAP WITH PROJECTION IN HEADER
            im = ax.imshow(probs,
                           cmap=cmap, origin='lower', alpha=0.7, zorder=1, vmin=vmin, vmax=vmax, aspect='auto')

            # CONTOURS - NEED TO ADD THE CUMMULATIVE PROBABILITY
            i = np.flipud(np.argsort(aMap))
            cumsum = np.cumsum(aMap[i])
            cls = np.empty_like(aMap)
            cls[i] = cumsum * 100 * stampProb
            cls[i] = cumsum * 100

            # EXTRACT CONTOUR VALUES AT HEALPIX INDICES
            contours = []
            contours[:] = [cls[i] for i in healpixIds]
            contours = np.reshape(np.array(contours), (yRange, xRange))

            # PLOT THE CONTOURS ON THE SAME PLOT
            CS = plt.contour(contours, linewidths=1,
                             alpha=0.3, zorder=3)
            plt.clabel(CS, fontsize=12, inline=1,
                       fmt='%2.1f', fontproperties=font)

            # RESET THE AXES TO THE FRAME OF THE FITS FILE
            ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
            ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

            # THE COORDINATES USED IN THE PLOT CAN BE ACCESSED USING THE COORDS
            # ATTRIBUTE (NOT X AND Y)
            lon = ax.coords[0]
            lat = ax.coords[1]

            lon.set_axislabel('RA (deg)', minpad=0.87,
                              size=20)
            lat.set_axislabel('DEC (deg)', minpad=0.87,
                              size=20)
            lon.set_major_formatter('d')
            lat.set_major_formatter('d')

            # THE SEPARATORS FOR ANGULAR COORDINATE TICK LABELS CAN ALSO BE SET BY
            # SPECIFYING A STRING
            # lat.set_separator(':-s')
            # SET THE APPROXIMATE NUMBER OF TICKS, WITH COLOR & PREVENT OVERLAPPING
            # TICK LABELS FROM BEING DISPLAYED.
            lon.set_ticks(number=6, color='#657b83',
                          exclude_overlapping=True, size=10)
            lat.set_ticks(number=10, color='#657b83',
                          exclude_overlapping=True, size=10)

            # MINOR TICKS NOT SHOWN BY DEFAULT
            lon.display_minor_ticks(True)
            lat.display_minor_ticks(True)
            # lat.set_minor_frequency(2)

            # CUSTOMISE TICK POSITIONS (l, b, r, t == left, bottom, right, or
            # top)
            lon.set_ticks_position('bt')
            lon.set_ticklabel_position('b')
            lon.set_ticklabel(size=20)
            lat.set_ticklabel(size=20)
            lon.set_axislabel_position('b')

            # HIDE AXES
            # lon.set_ticklabel_position('')
            # lat.set_ticklabel_position('')
            # lon.set_axislabel('', minpad=0.5, fontsize=12)
            # lat.set_axislabel('', minpad=0.5, fontsize=12)

            # ADD A GRID
            ax.coords.grid(color='#657b83', alpha=0.5, linestyle='dashed')

            lat.set_ticks_position('l')
            lat.set_ticklabel_position('l')
            lat.set_axislabel_position('l')

            lat.set_ticks_position('r')
            lat.set_ticklabel_position('r')
            lat.set_axislabel_position('r')

            plt.gca().invert_xaxis()

            # lon.set_ticks(number=20)
            # lat.set_ticks(number=3)

        elif projection == "gnomonic":
            # UNPACK THE PLOT PARAMETERS

            if allSky and "raRange" not in plotParameters:
                centralCoordinate = list(centralCoordinate)
                centralCoordinate[1] = 0.
                raRange = 360.
                decRange = 180.
            else:
                raRange = plotParameters["raRange"]
                decRange = plotParameters["decRange"]

            raMax = centralCoordinate[0] + raRange / 2.
            raMin = centralCoordinate[0] - raRange / 2.
            decMax = centralCoordinate[1] + decRange / 2.
            decMin = centralCoordinate[1] - decRange / 2.

            # DETERMINE THE PIXEL GRID X,Y RANGES
            xRange = int(raRange / pixelSizeDeg)
            yRange = int(decRange / pixelSizeDeg)
            largest = max(xRange, yRange)
            # xRange = largest
            # yRange = largest

            # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
            w.wcs.crpix = [xRange / 2., yRange / 2.]

            # USE THE "GNOMONIC" PROJECTION ("COORDINATESYS---PROJECTION")
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

            # CREATE A PIXEL GRID - 2 ARRAYS OF X, Y
            columns = []
            px = np.tile(np.arange(0, xRange), yRange)
            py = np.repeat(np.arange(0, yRange), xRange)

            # CONVERT THE PIXELS TO WORLD COORDINATES
            wr, wd = w.wcs_pix2world(px, py, 1)

            # MAKE SURE RA IS +VE
            nr = []
            nr[:] = [r if r > 0 else r + 360. for r in wr]
            wr = np.array(nr)

            # THETA: IS THE POLAR ANGLE, RANGING FROM 0 AT THE NORTH POLE TO PI AT THE SOUTH POLE.
            # PHI: THE AZIMUTHAL ANGLE ON THE SPHERE FROM 0 TO 2PI
            # CONVERT DEC TO THE REQUIRED HEALPIX FORMAT
            nd = -wd + 90.

            # CONVERT WORLD TO HEALPIX INDICES
            healpixIds = hp.ang2pix(nside, theta=nd * DEG_TO_RAD_FACTOR,
                                    phi=wr * DEG_TO_RAD_FACTOR)

            # NOW READ THE VALUES OF THE MAP AT THESE HEALPIX INDICES
            uniqueHealpixIds = np.unique(healpixIds)
            probs = []
            probs[:] = [aMap[i] for i in healpixIds]

            uniProb = []
            uniProb[:] = [aMap[i] for i in uniqueHealpixIds]

            stampProb = np.sum(uniProb)
            print "Probability for the plot stamp is %(stampProb)s" % locals()

            # RESHAPE THE ARRAY AS BITMAP
            probs = np.reshape(np.array(probs), (yRange, xRange))

            # CREATE THE FITS HEADER WITH WCS
            header = w.to_header()
            # CREATE THE FITS FILE
            hdu = fits.PrimaryHDU(header=header, data=probs)

            # GRAB THE WCS FROM HEADER GENERATED EARLIER
            from astropy.wcs import WCS
            from astropy.visualization.wcsaxes import WCSAxes

            wcs = WCS(hdu.header)
            # USE WCS AS THE PROJECTION
            ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=wcs)
            # note that the axes have to be explicitly added to the figure
            ax = fig.add_axes(ax)

            # PLOT MAP WITH PROJECTION IN HEADER
            im = ax.imshow(probs,
                           cmap=cmap, origin='lower', alpha=0.7, zorder=1, vmin=vmin, vmax=vmax, aspect='auto')

            # CONTOURS - NEED TO ADD THE CUMMULATIVE PROBABILITY
            i = np.flipud(np.argsort(aMap))
            cumsum = np.cumsum(aMap[i])
            cls = np.empty_like(aMap)
            cls[i] = cumsum * 100 * stampProb

            # EXTRACT CONTOUR VALUES AT HEALPIX INDICES
            contours = []
            contours[:] = [cls[i] for i in healpixIds]
            contours = np.reshape(np.array(contours), (yRange, xRange))

            # PLOT THE CONTOURS ON THE SAME PLOT
            CS = plt.contour(contours, linewidths=1,
                             alpha=0.3, zorder=3)
            plt.clabel(CS, fontsize=12, inline=1,
                       fmt='%2.1f', fontproperties=font)

            # RESET THE AXES TO THE FRAME OF THE FITS FILE
            ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
            ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

            # THE COORDINATES USED IN THE PLOT CAN BE ACCESSED USING THE COORDS
            # ATTRIBUTE (NOT X AND Y)
            lon = ax.coords[0]
            lat = ax.coords[1]

            lon.set_axislabel('RA (deg)', minpad=0.87,
                              size=20)
            lat.set_axislabel('DEC (deg)', minpad=0.87,
                              size=20)
            lon.set_major_formatter('d')
            lat.set_major_formatter('d')

            # THE SEPARATORS FOR ANGULAR COORDINATE TICK LABELS CAN ALSO BE SET BY
            # SPECIFYING A STRING
            lat.set_separator(':-s')
            # SET THE APPROXIMATE NUMBER OF TICKS, WITH COLOR & PREVENT OVERLAPPING
            # TICK LABELS FROM BEING DISPLAYED.
            lon.set_ticks(number=4, color='#657b83',
                          exclude_overlapping=True, size=10)
            lat.set_ticks(number=10, color='#657b83',
                          exclude_overlapping=True, size=10)

            # MINOR TICKS NOT SHOWN BY DEFAULT
            lon.display_minor_ticks(True)
            lat.display_minor_ticks(True)
            lat.set_minor_frequency(2)

            # CUSTOMISE TICK POSITIONS (l, b, r, t == left, bottom, right, or
            # top)
            lon.set_ticks_position('bt')
            lon.set_ticklabel_position('b')
            lon.set_ticklabel(size=20)
            lat.set_ticklabel(size=20)
            lon.set_axislabel_position('b')
            lat.set_ticks_position('lr')
            lat.set_ticklabel_position('l')
            lat.set_axislabel_position('l')

            # HIDE AXES
            # lon.set_ticklabel_position('')
            # lat.set_ticklabel_position('')
            # lon.set_axislabel('', minpad=0.5, fontsize=12)
            # lat.set_axislabel('', minpad=0.5, fontsize=12)

            # ADD A GRID
            ax.coords.grid(color='#657b83', alpha=0.5, linestyle='dashed')
            plt.gca().invert_xaxis()
            lon.set_ticks(number=20)

        else:
            self.log.error(
                'please give a valid projection. The projection given was `%(projection)s`.' % locals())

        header = w.to_header()
        # CREATE THE FITS FILE
        hdu = fits.PrimaryHDU(header=header, data=probs)

        # ADD RA LIMIT
        if raLimit:
            x = np.ones(100) * raLimit
            y = np.linspace(decMin, decMax, 100)
            ax.plot(x, y, 'b--', transform=ax.get_transform('fk5'))

        # SETUP TITLE OF PLOT
        from datetime import datetime, date, time
        now = datetime.now()
        now = now.strftime("%Y-%m-%d")
        timeRangeLabel = "NULL"
        if plotType == "timeline" and timeLimitDay:
            start = timeLimitDay[0]
            end = timeLimitDay[1]
            if self.filters:
                f = ("").join(self.filters)
            else:
                f = ""
            if self.telescope:
                t = self.telescope
            else:
                t = ""
            plotTitle = "%(gwid)s %(timeLimitLabel)s %(projection)s %(t)s %(f)s" % locals(
            )
            timeRangeLabel = timeLimitLabel.lower().replace(
                "in", "").replace("between", "").strip()
            if timeLimitLabel == "no limit":
                plotTitle = "%(gwid)s" % locals(
                )
                timeRangeLabel = "all transients"
        elif plotType == "timeline":

            plotTitle = "%(gwid)s %(mapBasename)s skymap %(projection)s" % locals(
            )
            timeRangeLabel = ""

        else:
            timeRangeLabel = ""

        subTitle = "(updated %(now)s)" % locals()
        if timeLimitDay == 0 or plotType == "timeline":
            subTitle = ""

        # ax.set_title(plotTitle + "\n", fontsize=10)
        # GRAB PS1 POINTINGS
        pointingArray = []

        from matplotlib.patches import Ellipse
        from matplotlib.patches import Circle
        from matplotlib.patches import Rectangle

        # ps1Pointings = []
        # PS1 POINTINGS (NOT SKYCELLS)
        if len(ps1Pointings) and "skycell" not in ps1Pointings[0]["exp_id"]:

            for psp in ps1Pointings:
                raDeg = psp["raDeg"]
                decDeg = psp["decDeg"]

                # REMOVE LOWER PROBABILITY FOOTPRINTS
                phi = raDeg
                if phi > 180.:
                    phi = phi - 360.
                theta = -decDeg + 90.
                healpixId = hp.ang2pix(
                    nside, theta * DEG_TO_RAD_FACTOR, phi * DEG_TO_RAD_FACTOR)
                probs = aMap[healpixId]
                probs = float("%0.*f" % (7, probs))
                if probabilityCut and probs == 0.:
                    continue

                # height = 2.8
                height = 0.2
                width = height / math.cos(decDeg * DEG_TO_RAD_FACTOR)

                # MULTIPLE CIRCLES
                if projection in ["mercator", "gnomonic", "cartesian"]:
                    circ = Ellipse(
                        (raDeg, decDeg), width=width, height=height, alpha=0.2, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
                else:
                    if raDeg > 180.:
                        raDeg = raDeg - 360.
                    circ = Ellipse(
                        (-raDeg * DEG_TO_RAD_FACTOR, decDeg * DEG_TO_RAD_FACTOR), width=width * DEG_TO_RAD_FACTOR, height=height * DEG_TO_RAD_FACTOR, alpha=0.2, color='#859900', fill=True, zorder=3)

                ax.add_patch(circ)
        else:
            plotted = []
            patches = []
            for psp in ps1Pointings:

                if psp["exp_id"] in plotted:
                    continue
                else:
                    plotted.append(psp["exp_id"])

                patch = add_square_fov(
                    log=self.log,
                    raDeg=psp["raDeg"],
                    decDeg=psp["decDeg"],
                    nside=nside,
                    aMap=aMap,
                    fovSide=0.4,
                    axes=ax,
                    projection=projection)

                if patch:
                    patches.append(patch)

            if projection in ["mercator", "gnomonic", "cartesian"]:
                ax.add_collection(PatchCollection(patches, alpha=0.2,
                                                  color='#859900', zorder=2, transform=ax.get_transform('fk5')))
            else:
                ax.add_collection(PatchCollection(patches, alpha=0.2,
                                                  color='#859900', zorder=2))

        patches = []
        for atp in atlasPointings:

            if atp["exp_id"] in plotted:
                continue
            else:
                plotted.append(atp["exp_id"])

            patch = add_square_fov(
                log=self.log,
                raDeg=atp["raDeg"],
                decDeg=atp["decDeg"],
                nside=nside,
                aMap=aMap,
                fovSide=5.46,
                axes=ax,
                projection=projection)

            if patch:
                patches.append(patch)

        if projection in ["mercator", "gnomonic", "cartesian"]:
            ax.add_collection(PatchCollection(patches, alpha=0.2,
                                              color="#6c71c4", zorder=2, transform=ax.get_transform('fk5')))
        else:
            ax.add_collection(PatchCollection(patches, alpha=0.2,
                                              color="#6c71c4", zorder=2))

        # ADD DATA POINTS FOR TRANSIENTS
        names = []
        ra = []
        dec = []
        raRad = []
        decRad = []
        texts = []

        # ps1Transients = []
        for trans in ps1Transients:
            # if trans["ps1_designation"] in ["PS15dpg", "PS15dpp", "PS15dpq", "PS15don", "PS15dpa", "PS15dom"]:
            #     continue

            name = trans["ps1_designation"]
            if name is None:
                name = trans["local_designation"]
            names.append(name)
            raDeg = trans["ra_psf"]
            decDeg = trans["dec_psf"]
            print name, raDeg, decDeg
            ra.append(raDeg)
            dec.append(decDeg)
            raRad.append(-raDeg * DEG_TO_RAD_FACTOR)
            decRad.append(decDeg * DEG_TO_RAD_FACTOR)

        if len(ra) > 0:
            # MULTIPLE CIRCLES
            if projection in ["mercator", "gnomonic", "cartesian"]:
                ax.scatter(
                    x=np.array(ra),
                    y=np.array(dec),
                    transform=ax.get_transform('fk5'),
                    s=6,
                    c='#036d09',
                    edgecolor='#036d09',
                    alpha=1,
                    zorder=4
                )
                xx, yy = w.wcs_world2pix(np.array(ra), np.array(dec), 0)
                # ADD TRANSIENT LABELS
                for r, d, n in zip(xx, yy, names):
                    texts.append(ax.text(
                        r,
                        d,
                        n,
                        color='#036d09',
                        fontsize=10,
                        zorder=4,
                        family='monospace'
                    ))

                if len(texts):
                    adjust_text(
                        xx,
                        yy,
                        texts,
                        expand_text=(1.2, 1.6),
                        expand_points=(1.2, 3.2),
                        va='center',
                        ha='center',
                        force_text=2.0,
                        force_points=0.5,
                        lim=1000,
                        precision=0,
                        only_move={},
                        text_from_text=True,
                        text_from_points=True,
                        save_steps=False,
                        save_prefix='',
                        save_format='png',
                        add_step_numbers=True,
                        min_arrow_sep=50.0,
                        draggable=True,
                        arrowprops=dict(arrowstyle="-", color='#036d09', lw=1.2,
                                        patchB=None, shrinkB=0, connectionstyle="arc3,rad=0.1", zorder=3, alpha=0.5),
                        fontsize=10,
                        family='monospace'
                    )
            else:
                ax.scatter(
                    x=np.array(raRad),
                    y=np.array(decRad),
                    s=6,
                    c='#dc322f',
                    edgecolor='#dc322f',
                    alpha=1,
                    zorder=4
                )

        # ADD DATA POINTS FOR TRANSIENTS
        names = []
        ra = []
        dec = []
        raRad = []
        decRad = []
        texts = []
        # atlasTransients = []
        for trans in atlasTransients:
            # if trans["ps1_designation"] in ["PS15dpg", "PS15dpp", "PS15dpq", "PS15don", "PS15dpa", "PS15dom"]:
            #     continue

            name = trans["atlas_designation"]
            names.append(name)
            raDeg = trans["ra"]
            decDeg = trans["dec"]
            ra.append(raDeg)
            dec.append(decDeg)
            raRad.append(-raDeg * DEG_TO_RAD_FACTOR)
            decRad.append(decDeg * DEG_TO_RAD_FACTOR)

            print name, raDeg, decDeg

        if len(ra) > 0:
            # MULTIPLE CIRCLES
            if projection in ["mercator", "gnomonic", "cartesian"]:
                ax.scatter(
                    x=np.array(ra),
                    y=np.array(dec),
                    transform=ax.get_transform('fk5'),
                    s=6,
                    c='#5c0bb0',
                    edgecolor='#5c0bb0',
                    alpha=1,
                    zorder=4
                )
                xx, yy = w.wcs_world2pix(np.array(ra), np.array(dec), 0)
                # ADD TRANSIENT LABELS
                for r, d, n in zip(xx, yy, names):
                    texts.append(ax.text(
                        r,
                        d,
                        n,
                        fontsize=10,
                        color='#5c0bb0',
                        zorder=4,
                        family='monospace'
                    ))

                if len(texts):
                    adjust_text(
                        xx,
                        yy,
                        texts,
                        expand_text=(1.2, 1.6),
                        expand_points=(1.2, 3.2),
                        va='center',
                        ha='center',
                        force_text=2.0,
                        force_points=0.5,
                        lim=1000,
                        precision=0,
                        only_move={},
                        text_from_text=True,
                        text_from_points=True,
                        save_steps=False,
                        save_prefix='',
                        save_format='png',
                        add_step_numbers=True,
                        min_arrow_sep=50.0,
                        draggable=True,
                        arrowprops=dict(arrowstyle="-", color='#5c0bb0', lw=1.2,
                                        patchB=None, shrinkB=0, connectionstyle="arc3,rad=0.1", zorder=3, alpha=0.5),
                        fontsize=10,
                        family='monospace'
                    )
            else:
                ax.scatter(
                    x=np.array(raRad),
                    y=np.array(decRad),
                    s=6,
                    c='#dc322f',
                    edgecolor='#dc322f',
                    alpha=1,
                    zorder=4
                )

        # TIME-RANGE LABEL
        fig = plt.gcf()
        fWidth, fHeight = fig.get_size_inches()

        if projection in ["mercator", "gnomonic", "cartesian"]:
            fig.set_size_inches(8.0, 8.0)
            ax.text(0.95, 0.95, timeRangeLabel,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    color="#dc322f",
                    fontproperties=font,
                    fontsize=16,
                    zorder=4)

        else:
            ax.text(0.95, 0.95, timeRangeLabel,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    color="#dc322f",
                    fontproperties=font,
                    fontsize=16,
                    zorder=4)

        if self.timestamp:
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%d %H:%M.%S UTC")
            ax.text(0, 1.02, utcnow,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    color="#657b83",
                    fontproperties=font,
                    fontsize=6)

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if self.settings and not outputDirectory:
            plotDir = self.settings["output directory"] + "/" + gwid
        elif outputDirectory:
            plotDir = outputDirectory

        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = plotTitle.replace("  ", " ").replace("  ", " ").replace("  ", " ").replace(" ", "_").replace("-", "_").replace(
            "<", "lt").replace(">", "gt").replace(",", "").replace("\n", "_").replace("&", "and")
        figureName = """%(plotTitle)s""" % locals(
        )
        if timeLimitDay == 0:
            figureName = """%(plotTitle)s""" % locals(
            )
        if allSky and plotDir == ".":
            figureName = figureName + "_" + projection.title()
        if plotDir != ".":
            for f in fileFormats:
                if not os.path.exists("%(plotDir)s/%(folderName)s/%(f)s" % locals()):
                    os.makedirs("%(plotDir)s/%(folderName)s/%(f)s" % locals())
                figurePath = "%(plotDir)s/%(folderName)s/%(f)s/%(figureName)s.%(f)s" % locals()
                figurePath = figurePath.replace("_.", ".")

                if f == "pdf":
                    matplotlib.use('PDF')
                elif f == "png":
                    matplotlib.use('TkAgg')

                savefig(figurePath, bbox_inches='tight', dpi=300)
                # savefig(figurePath, dpi=300)

                if bestMap and allSky and symLink:
                    linkName = "%(plotDir)s/%(folderName)s/%(f)s/%(gwid)s_preferred_skymap_%(projection)s.%(f)s" % locals()
                    try:
                        os.remove(linkName)
                    except:
                        pass
                    os.symlink(figurePath, linkName)

            # if not os.path.exists("%(plotDir)s/%(folderName)s/fits" % locals()):
            #     os.makedirs("%(plotDir)s/%(folderName)s/fits" % locals())
            # pathToExportFits = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_map_%(projection)s.fits" % locals()
            # try:
            #     os.remove(pathToExportFits)
            # except:
            #     pass
            # hdu.writeto(pathToExportFits)
        else:
            for f in fileFormats:

                if f == "pdf":
                    matplotlib.use('PDF')
                elif f == "png":
                    matplotlib.use('TkAgg')

                figurePath = "%(plotDir)s/%(figureName)s.%(f)s" % locals()
                figurePath = figurePath.replace("_.", ".")
                savefig(figurePath, bbox_inches='tight', dpi=300)
                # savefig(figurePath, dpi=300)

            # pathToExportFits = "%(plotDir)s/%(gwid)s_skymap.fits" % locals()
            # try:
            #     os.remove(pathToExportFits)
            # except:
            #     pass
            # hdu.writeto(pathToExportFits)

        if fitsImage:
            self.generate_fits_image_map(
                gwid=gwid,
                pathToProbMap=pathToProbMap,
                folderName=folderName,
                outputDirectory=outputDirectory,
                center=center,
                bestMap=bestMap
            )

        self.log.debug('completed the ``generate_probability_plot`` method')
        return None

    def get_history_plots(
            self):
        """
        *plot the history plots*

        **Return:**
            - None

        **Usage:**

            .. code-block:: python

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                       settings=settings,
                       plotType="history",
                       gwid="G184098",
                       projection="mercator"
                )
                plotter.get()
        """
        self.log.debug('starting the ``get_history_plots`` method')

        timeLimitLabels = ["day", "2 days", "3 days", "4 days", "5 days", "6 days",
                           "7 days", "2 weeks", "3 weeks", "1 month", "2 months", "3 months", "no limit"]
        timeLimitDays = [1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 0]

        if self.gwid:
            theseIds = [self.gwid]
        else:
            theseIds = self.settings["gravitational waves"]

        for gwid in theseIds:
            for tday, tlabel in zip(timeLimitDays, timeLimitLabels):

                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = self.get_gw_parameters_from_settings(
                    gwid=gwid,
                    inPastDays=tday,
                    inFirstDays=False)

                pathToProbMap = self.settings[
                    "gravitational waves"][gwid]["mapPath"]
                if not os.path.exists(pathToProbMap):
                    message = "the path to the map %s does not exist on this machine" % (
                        pathToProbMap,)
                    self.log.critical(message)
                    raise IOError(message)

                mjdStart = self.settings["gravitational waves"][
                    gwid]["mjd"]

                self.generate_probability_plot(
                    gwid=gwid,
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    atlasTransients=atlasTransients,
                    ps1Pointings=ps1Pointings,
                    atlasPointings=atlasPointings,
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    timeLimitLabel=tlabel,
                    timeLimitDay=tday,
                    raLimit=False,
                    fileFormats=["png", "pdf"],
                    folderName="survey_history_plots",
                    plotType=self.plotType,
                    projection=self.projection,
                    probabilityCut=self.probabilityCut)

        self.log.debug('completed the ``get_history_plots`` method')
        return None

    def get_timeline_plots(
            self):
        """
        *plot the history plots*

        **Return:**
            - None

        **Usage:**

            .. code-block:: python

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                       settings=settings,
                       plotType="timeline",
                       gwid="G184098",
                       projection="mercator"
                )
                plotter.get()
        """
        self.log.debug('starting the ``get_timeline_plots`` method')

        if self.allPlots:
            timeLimitLabels = ["21 days pre-detection", "<1d", "1-2d",
                               "2-3d", "3-4d", "4-5d", "5-10d", "10-17d", "17-24d", "24-31d"]
            timeLimitDays = [(-21, 0), (0, 1), (1, 2), (2, 3), (3, 4),
                             (4, 5), (5, 10), (10, 17), (17, 24), (24, 31)]
        else:
            timeLimitLabels = ["0-1d", "1-2d", "2-3d", "3-4d",
                               "4-5d", "5-10d", "10-17d", "17-24d", "24-31d"]
            # timeLimitLabels = ["0-1d"]
            timeLimitDays = [(0, 1), (1, 2), (2, 3), (3, 4),
                             (4, 5), (5, 10), (10, 17), (17, 24), (24, 31)]

        raLimits = [134.25, 144.75, 152.25, 159.50, 167.0, 174.5]
        raLimits = [False, False, False, False,
                    False, False, False, False, False, False, False, False]

        if self.gwid:
            theseIds = [self.gwid]
        else:
            theseIds = self.settings["gravitational waves"]

        for gwid in theseIds:
            for tday, tlabel, raLimit in zip(timeLimitDays, timeLimitLabels, raLimits):

                pathToProbMap = self.settings[
                    "gravitational waves"][gwid]["mapPath"]
                aMap, mapHeader = hp.read_map(
                    pathToProbMap, 0, h=True, verbose=False)

                mapBasename = os.path.basename(pathToProbMap)
                mapBasename = os.path.splitext(mapBasename)[0]
                mapBasename = os.path.splitext(mapBasename)[0]
                mapBasename = os.path.splitext(mapBasename)[0]

                # DETERMINE THE SIZE OF THE HEALPIXELS
                nside = hp.npix2nside(len(aMap))
                maxProbHealpix = aMap.argmax()
                maxProbCoordinate = hp.pix2ang(
                    nside, maxProbHealpix, lonlat=True)

                plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = self.get_gw_parameters_from_settings(
                    gwid=gwid,
                    inPastDays=False,
                    inFirstDays=tday,
                    maxProbCoordinate=maxProbCoordinate)

                if not os.path.exists(pathToProbMap):
                    message = "the path to the map %s does not exist on this machine" % (
                        pathToProbMap,)
                    self.log.critical(message)
                    raise IOError(message)

                mjdStart = self.settings["gravitational waves"][
                    gwid]["mjd"]

                if "raRange" not in plotParameters:
                    allSky = True
                    center = False
                else:
                    allSky = False
                    center = maxProbCoordinate

                self.generate_probability_plot(
                    gwid=gwid,
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    ps1Pointings=ps1Pointings,
                    atlasTransients=atlasTransients,
                    atlasPointings=atlasPointings,
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    timeLimitLabel=tlabel,
                    timeLimitDay=tday,
                    raLimit=raLimit,
                    fileFormats=["png"],
                    folderName="survey_timeline_plots",
                    plotType=self.plotType,
                    projection=self.projection,
                    probabilityCut=self.probabilityCut,
                    allSky=allSky,
                    center=center,
                    symLink=False)

        self.log.debug('completed the ``get_timeline_plots`` method')
        return None

    def generate_fits_image_map(
            self,
            gwid,
            pathToProbMap,
            folderName="",
            outputDirectory=False,
            rebin=True,
            center=False,
            bestMap=False):
        """*generate fits image map from the LV-skymap (FITS binary table)*

        **Key Arguments:**
            - ``pathToProbMap`` -- path to the FITS file containing the probability map of the wave
            - ``outputDirectory`` -- can be used to override the output destination in the settings file
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``folderName`` -- the name of the folder to add the plots to
            - ``rebin`` -- rebin the final image to reduce size
            - ``center`` -- central longitude in degrees. Default *0*.
            - ``bestMap`` -- is this the prefered skymap. If so, add symlink to placeholder name for prefered map.

        **Return:**
            - None

        **Usage:**

            To generate an all-sky image from the LV FITS binary table healpix map run the following code:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings,
                    databaseConnRequired=False
                )
                plotter.generate_fits_image_map(
                    gwid="G211117",
                    pathToProbMap="/path/to/LV/healpix_map.fits",
                    outputDirectory="/path/to/output",
                    rebin=True
                )

        The size of the final FITS image map is ~1.1GB so it's probably best to rebin the image (~80MB) unless you really need the resolution.
        """
        self.log.debug('starting the ``generate_fits_image_map`` method')

        import healpy as hp
        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        mapBasename = os.path.basename(pathToProbMap)
        mapBasename = os.path.splitext(mapBasename)[0]
        mapBasename = os.path.splitext(mapBasename)[0]
        mapBasename = os.path.splitext(mapBasename)[0]

        # X, Y PIXEL COORDINATE GRID
        xRange = 10000
        yRange = xRange / 2

        # PIXELSIZE AS MAPPED TO THE FULL SKY
        pixelSizeDeg = 360. / xRange

        # READ HEALPIX MAPS FROM FITS FILE
        # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
        # ARRAY OF PROBABILITIES (3,072 ROWS)
        # READ IN THE HEALPIX FITS FILE
        aMap, mapHeader = hp.read_map(pathToProbMap, 0, h=True, verbose=False)
        # DETERMINE THE SIZE OF THE HEALPIXELS
        nside = hp.npix2nside(len(aMap))

        # FROM THE PIXEL GRID (xRange, yRange), GENERATE A MAP TO LAT (pi to 0) AND LONG (-pi to pi) THAT CAN THEN MAPS TO HEALPIX SKYMAP
        # RA FROM -pi to pi

        # FITS convention for fractional pixels is that the center of the lower left pixel is at (1,1), the lower left corner of the
        # lower left pixel is at (0.5,0.5). [1-index pix] refers to this
        # convention. arange index starts at 0 so we need to fix for this

        # FULL-SKY MAP SO PLOT FULL RA AND DEC RANGES
        # DEC FROM 180 to 0
        theta = np.linspace(np.pi - pixelSizeDeg / 2, +
                            pixelSizeDeg / 2, yRange)

        latitude = np.radians(np.linspace(-90 + pixelSizeDeg, 90, yRange))

        # FIND THE COORDINATES OF THE CORE LIKEIHOOD
        maxProbHealpix = aMap.argmax()
        maxCoordinate = hp.pix2ang(nside, maxProbHealpix, lonlat=True)
        print "The %(gwid)s %(mapBasename)s map's maximum likelihood is centered at %(maxCoordinate)s" % locals()

        if center == False:
            center = maxCoordinate[0]

        # RA FROM -180 to +180
        centralRa = center
        centralRaRad = centralRa * DEG_TO_RAD_FACTOR
        phi = np.linspace(-np.pi + centralRaRad + pixelSizeDeg / 2,
                          np.pi + centralRaRad - pixelSizeDeg / 2, xRange)

        longitude = np.radians(np.linspace(-180 + pixelSizeDeg, 180, xRange))
        X, Y = np.meshgrid(longitude, latitude)

        # PROJECT THE MAP TO A RECTANGULAR MATRIX xRange X yRange
        PHI, THETA = np.meshgrid(phi, theta)
        healpixIds = hp.ang2pix(nside, THETA, PHI)
        # GIVEN A HIGH ENOUGH RESOLUTION IN THE PIXEL GRID WE WILL HAVE
        # DUPLICATES - LET COUNT THEM ADD DIVIDE PROBABILITY EQUALLY
        unique, counts = np.unique(healpixIds, return_counts=True)

        # countDict has healpixid as keys and count rates as values (e.g. the
        # ID XXXX occurs in 78 pixels)
        countDict = dict(zip(unique, counts))

        # PROB SHAPE IS (5000, 10000) .. indexing measured from top left in
        # matrix.
        probs = aMap[healpixIds]

        # NOTE i = -y-direction, shape[0] is the y-axis range
        # j = x-direction, shape[1] is the x-axis range
        weightedProb = np.array([
            [probs[i, j] / countDict[healpixIds[i, j]]
                for j in xrange(probs.shape[1])]
            for i in xrange(probs.shape[0])
        ])

        if rebin == True:
            rebinSize = 4
            # resize by getting rid of extra columns/rows
            # xedge and yedge find the pixels we need to trim is rebinSize does
            # divide evenly into image
            yedge = np.shape(weightedProb)[0] % rebinSize
            xedge = np.shape(weightedProb)[1] % rebinSize
            weightedProb = weightedProb[yedge:, xedge:]

            # put image array into arrays of rebinSize x rebinSize - so (1250,
            # 4, 2500, 4)
            weightedProb = np.reshape(weightedProb, (np.shape(weightedProb)[
                                      0] / rebinSize, rebinSize, np.shape(weightedProb)[1] / rebinSize, rebinSize))

            # average each rebinSize x rebinSize array
            weightedProb = np.mean(weightedProb, axis=3) * rebinSize
            weightedProb = np.mean(weightedProb, axis=1) * rebinSize

            pixelSizeDeg = pixelSizeDeg * rebinSize
            xRange = int(xRange / rebinSize)
            yRange = int(yRange / rebinSize)

        # CREATE A NEW WCS OBJECT
        w = awcs.WCS(naxis=2)
        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = np.array([pixelSizeDeg, pixelSizeDeg])
        # WORLD COORDINATES AT REFERENCE PIXEL
        centralCoordinate = [centralRa, 0]
        w.wcs.crval = centralCoordinate
        cx = xRange / 2. + 0.5
        cy = yRange / 2. + 0.5
        # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
        w.wcs.crpix = [cx, cy]

        unweightedImageProb = np.sum(probs)
        self.log.info(
            "The total unweighted probability flux in the FITS images added to %(unweightedImageProb)s" % locals())

        totalImageProb = np.sum(weightedProb)
        self.log.info(
            "The total probability flux in the FITS images added to %(totalImageProb)s" % locals())

        # CTYPE FOR THE FITS HEADER
        w.wcs.ctype = ["RA---CAR" %
                       locals(), "DEC--CAR" % locals()]

        header = w.to_header()
        # CREATE THE FITS FILE
        hdu = fits.PrimaryHDU(header=header, data=weightedProb)

        # Recursively create missing directories
        if self.settings and not outputDirectory:
            plotDir = self.settings["output directory"] + "/" + gwid
        elif outputDirectory:
            plotDir = outputDirectory

        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        if plotDir != ".":
            if not os.path.exists("%(plotDir)s/%(folderName)s/fits" % locals()):
                os.makedirs("%(plotDir)s/%(folderName)s/fits" % locals())
            pathToExportFits = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_%(mapBasename)s_breaker_skymap.fits" % locals()
            try:
                os.remove(pathToExportFits)
            except:
                pass
            hdu.writeto(pathToExportFits)

            if bestMap and symLink:
                linkName = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_preferred_breaker_skymap.fits" % locals()
                print "The %(gwid)s prefered likeihood map is symlinked at `%(linkName)s`" % locals()
                try:
                    os.remove(linkName)
                except:
                    pass
                os.symlink(pathToExportFits, linkName)
        else:
            pathToExportFits = "%(plotDir)s/%(gwid)s_%(mapBasename)s_breaker_skymap.fits" % locals()
            try:
                os.remove(pathToExportFits)
            except:
                pass
            hdu.writeto(pathToExportFits)

        print "The %(gwid)s %(mapBasename)s likeihood map can be found here `%(pathToExportFits)s`" % locals()

        self.log.debug('completed the ``generate_fits_image_map`` method')
        return None


def y2lat(
    y,
    xRange,
    yRange
):
    # R is the radius of the sphere at the scale of the map as drawn
    y = y - yRange / 2
    R = xRange / (2. * np.pi)
    return -((np.pi / 2.0) - 2.0 * np.arctan(np.e ** (-y / R))) + np.pi / 2


def x2long(
    x,
    xRange
):
    # R is the radius of the sphere at the scale of the map as drawn
    R = xRange / (2. * np.pi)
    return x / R - np.pi


def add_square_fov(
        log,
        raDeg,
        decDeg,
        nside,
        aMap,
        fovSide,
        axes,
        projection):
    """*summary of function*

    **Key Arguments:**
        - ``dbConn`` -- mysql database connection
        - ``log`` -- logger

    **Return:**
        - None

    **Usage:**
        .. todo::

            add usage info
            create a sublime snippet for usage

        .. code-block:: python 

            usage code            
    """
    log.debug('starting the ``add_square_fov`` function')

    import math
    pi = (4 * math.atan(1.0))
    DEG_TO_RAD_FACTOR = pi / 180.0
    RAD_TO_DEG_FACTOR = 180.0 / pi

    deltaDeg = fovSide / 2
    if decDeg < 0:
        deltaDeg = -deltaDeg

    if projection in ["mercator", "gnomonic", "cartesian"]:
        widthDegTop = fovSide / \
            math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
        widthDegBottom = fovSide / \
            math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
        heightDeg = fovSide
        llx = (raDeg - widthDegBottom / 2)
        lly = decDeg - (heightDeg / 2)
        ulx = (raDeg - widthDegTop / 2)
        uly = decDeg + (heightDeg / 2)
        urx = (raDeg + widthDegTop / 2)
        ury = uly
        lrx = (raDeg + widthDegBottom / 2)
        lry = lly

        Path = mpath.Path
        path_data = [
            (Path.MOVETO, [llx, lly]),
            (Path.LINETO, [ulx, uly]),
            (Path.LINETO, [urx, ury]),
            (Path.LINETO, [lrx, lry]),
            (Path.CLOSEPOLY, [llx, lly])
        ]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)

        # EXCLUDE FOOTPRINTS THAT CROSS 360.
        if lrx > 180. and llx < 180:
            return None

        patch = patches.PathPatch(path)
    else:
        if raDeg > 180.:
            raDeg = raDeg - 360.

        widthRadTop = fovSide * DEG_TO_RAD_FACTOR / \
            math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
        widthRadBottom = fovSide * DEG_TO_RAD_FACTOR / \
            math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
        heightRad = fovSide * DEG_TO_RAD_FACTOR
        llx = -(raDeg * DEG_TO_RAD_FACTOR - widthRadBottom / 2)
        lly = decDeg * DEG_TO_RAD_FACTOR - (heightRad / 2)
        ulx = -(raDeg * DEG_TO_RAD_FACTOR - widthRadTop / 2)
        uly = decDeg * DEG_TO_RAD_FACTOR + (heightRad / 2)
        urx = -(raDeg * DEG_TO_RAD_FACTOR + widthRadTop / 2)
        ury = uly
        lrx = -(raDeg * DEG_TO_RAD_FACTOR + widthRadBottom / 2)
        lry = lly
        Path = mpath.Path
        path_data = [
            (Path.MOVETO, [llx, lly]),
            (Path.LINETO, [ulx, uly]),
            (Path.LINETO, [urx, ury]),
            (Path.LINETO, [lrx, lry]),
            (Path.CLOSEPOLY, [llx, lly])
        ]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = patches.PathPatch(path)

    log.debug('completed the ``add_square_fov`` function')
    return patch

# use the tab-trigger below for new function
# xt-def-function
