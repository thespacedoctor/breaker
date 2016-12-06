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
import warnings
warnings.filterwarnings("ignore")
import sys
import os
os.environ['TERM'] = 'vt100'
import healpy as hp
import numpy as np
import math
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
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


from matplotlib.projections.geo import GeoAxes


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
        - ``projection`` -- projection for the plot. Default *tan*
        - ``probabilityCut`` -- remove footprints where probability assigned to the healpix pixel found at the center of the exposure is ~0.0. Default *False*
        - ``databaseConnRequired`` -- are the database connections going to be required? Default *True*

    **Usage:**

        To plot a the history of a specific wave:

        .. code-block:: python

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="history",
                   gwid="G184098",
                   projection="tan"
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
            projection="tan",
            probabilityCut=False,
            databaseConnRequired=True
    ):
        self.log = log
        log.debug("instantiating a new 'plot_wave_observational_timelines' object")
        self.settings = settings
        self.plotType = plotType
        self.gwid = gwid
        self.projection = projection
        self.probabilityCut = probabilityCut

        # xt-self-arg-tmpx

        # Initial Actions

        if self.settings and databaseConnRequired:
            # CONNECT TO THE VARIOUS DATABASES REQUIRED
            from breaker import database
            db = database(
                log=self.log,
                settings=self.settings
            )
            self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = db.get()
        else:
            self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = False, False, False

        self.log.debug(
            'connected to databases')

        return None

    def get(self):
        """
        *Generate the plots*
        """
        self.log.info('starting the ``get`` method')

        if self.plotType == "history":
            self.get_history_plots()
        elif self.plotType == "timeline":
            self.get_timeline_plots()

        self.log.info('completed the ``get`` method')
        return None

    def get_gw_parameters_from_settings(
            self,
            gwid,
            inPastDays=False,
            inFirstDays=False):
        """
        *Query the settings file and database for PS1 Pointings, PS1 discovered transients and plot parameters relatiing to the given gravitational wave (``gwid``)*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection). A tuple (start day, end day).

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
                plotParameters, ps1Transients, ps1Pointings, atlasPointings = plotter.get_gw_parameters_from_settings(
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
                plotParameters, ps1Transients, ps1Pointings, atlasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inPastDays=7
                )

            Or the first 3 days since wave detection:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, atlasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inFirstDays=(0,3)
                )
        """
        self.log.info(
            'starting the ``get_gw_parameters_from_settings`` method')

        plotParameters = self.settings["gravitational waves"][gwid]["plot"]

        # GRAB PS1 TRANSIENTS FROM THE DATABASE
        ps1Transients = self._get_ps1_transient_candidates(
            gwid=gwid,
            mjdStart=self.settings["gravitational waves"][
                gwid]["time"]["mjdStart"],
            mjdEnd=self.settings["gravitational waves"][
                gwid]["time"]["mjdEnd"],
            plotParameters=plotParameters,
            inPastDays=inPastDays,
            inFirstDays=inFirstDays
        )

        self.log.debug(
            'finished getting the PS1 transients')

        # GRAB PS1 & ATLAS POINTINGS FROM THE DATABASE
        ps1Pointings = self._get_ps1_pointings(gwid, inPastDays, inFirstDays)
        atlasPointings = self._get_atlas_pointings(
            gwid, inPastDays, inFirstDays)

        self.log.info(
            'completed the ``get_gw_parameters_from_settings`` method')
        return plotParameters, ps1Transients, ps1Pointings, atlasPointings

    def _get_ps1_transient_candidates(
            self,
            gwid,
            mjdStart,
            mjdEnd,
            plotParameters,
            inPastDays,
            inFirstDays):
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
        self.log.info('starting the ``_get_ps1_transient_candidates`` method')

        # UNPACK THE PLOT PARAMETERS
        centralCoordinate = plotParameters["centralCoordinate"]
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
            mjdStart = self.settings["gravitational waves"][gwid]["time"][
                "mjdStart"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["time"]["mjdStart"] + inFirstDays[1]
            if inFirstDays[1] == 0:
                mjdEnd = 10000000000

        sqlQuery = u"""
            SELECT ps1_designation, local_designation, ra_psf, dec_psf FROM tcs_transient_objects o, tcs_latest_object_stats s where o.detection_list_id = 2 and o.id=s.id and (s.earliest_mjd between %(mjdStart)s and %(mjdEnd)s) and (ra_psf between %(raMin)s and %(raMax)s) and (dec_psf between %(decMin)s and %(decMax)s) ;
        """ % locals()
        ps1Transients = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ps1gwDbConn
        )

        self.log.info('completed the ``_get_ps1_transient_candidates`` method')
        return ps1Transients

    def _get_ps1_pointings(
            self,
            gwid,
            inPastDays,
            inFirstDays):
        """
        *get ps1 pointings to add to the plot*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``inPastDays`` -- used for the `history` plots (looking back from today)
            - ``inFirstDays`` -- used in the `timeline` plots (looking forward from wave detection).  A tuple (start day, end day).

        **Return:**
            - ``ps1Pointings`` -- the pointings to place on the plot
        """
        self.log.info('starting the ``_get_ps1_pointings`` method')

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
            mjdStart = self.settings["gravitational waves"][gwid]["time"][
                "mjdStart"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["time"]["mjdStart"] + inFirstDays[1]
            if inFirstDays[1] == 0:
                mjdEnd = 10000000000

        sqlQuery = u"""
            SELECT raDeg, decDeg, mjd FROM ps1_pointings where gw_id = "%(gwid)s" and mjd between %(mjdStart)s and %(mjdEnd)s
        """ % locals()

        ps1Pointings = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.info('completed the ``_get_ps1_pointings`` method')
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
        self.log.info('starting the ``_get_atlas_pointings`` method')

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
            mjdStart = self.settings["gravitational waves"][gwid]["time"][
                "mjdStart"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["time"]["mjdStart"] + inFirstDays[1]
            if inFirstDays[1] == 0:
                mjdEnd = 10000000000

        sqlQuery = u"""
            SELECT atlas_object_id, raDeg, decDeg, mjd FROM atlas_pointings where gw_id = "%(gwid)s" and mjd between %(mjdStart)s and %(mjdEnd)s group by atlas_object_id;
        """ % locals()

        atlasPointings = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.info('completed the ``_get_atlas_pointings`` method')
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
            ps1Pointings=[],
            atlasPointings=[],
            projection="wcs",
            raLimit=False,
            probabilityCut=False,
            outputDirectory=False):
        """
        *Generate a single probability map plot for a given gravitational wave and save it to file*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``plotParameters`` -- the parameters of the plot (for spatial & temporal parameters etc).
            - ``ps1Transients`` -- the transients to add to the plot. Default **[]**
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
            - ``projection`` -- projection for the plot. Default *wcs*. [wcs|mollweide]
            - ``probabilityCut`` -- remove footprints where probability assigned to the healpix pixel found at the center of the exposure is ~0.0. Default *False*
            - ``outputDirectory`` -- can be used to override the output destination in the settings file


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
                plotParameters, ps1Transients, ps1Pointings, atlasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inFirstDays=(0,7)
                )

            Then you can pass in these parameter to generate a plot:

            .. code-block:: python

                plotter.generate_probability_plot(
                    gwid="G211117",
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    ps1Pointings=ps1Pointings,
                    atlasPointings=altasPointings,
                    pathToProbMap="/Users/Dave/config/breaker/maps/G211117/LALInference_skymap.fits",
                    mjdStart=57382.,
                    timeLimitLabel="",
                    timeLimitDay=(0, 5),
                    raLimit=False,
                    fileFormats=["pdf"],
                    folderName="survey_timeline_plots",
                    projection="tan",
                    plotType="timeline",
                    probabilityCut=True,
                    outputDirectory=False
                )

        """
        self.log.info('starting the ``generate_probability_plot`` method')

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
        if plotParameters:
            centralCoordinate = plotParameters["centralCoordinate"]
        else:
            centralCoordinate = [0, 0]

        # CREATE A NEW WCS OBJECT
        w = awcs.WCS(naxis=2)
        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = np.array([pixelSizeDeg, pixelSizeDeg])
        # WORLD COORDINATES AT REFERENCE PIXEL
        w.wcs.crval = centralCoordinate

        projectionDict = {
            "mollweide": "MOL",
            "aitoff": "AIT",
            "hammer": "AIT",
            "lambert": "ZEA",
            "polar": "TAN",
            "rectilinear": "MER"
        }

        if projection in ["mollweide", "aitoff", "hammer", "lambert", "polar", "rectilinear"]:
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
            longitude = np.radians(np.linspace(-180, 180, xRange))
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
            ax.set_longitude_grid(60)
            ax.set_latitude_grid(45)
            ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
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

            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
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

        elif projection in ["tan"]:

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
            from wcsaxes import datasets, WCS
            from astropy.wcs import WCS
            from wcsaxes import WCSAxes

            wcs = WCS(hdu.header)
            # USE WCS AS THE PROJECTION
            ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=wcs)
            # note that the axes have to be explicitly added to the figure
            ax = fig.add_axes(ax)

            # PLOT MAP WITH PROJECTION IN HEADER
            im = ax.imshow(probs,
                           cmap=cmap, origin='lower', alpha=0.7, zorder=1, vmin=vmin, vmax=vmax)

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
        if plotType == "history":
            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Transients\nDiscovered in Last %(timeLimitDay)02d Days" % locals()
            timeRangeLabel = "Past %(timeLimitDay)02d days (updated %(now)s)" % locals(
            )
            if timeLimitDay == 0:
                plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Transients\nDiscovered Since MJD %(mjdStart)s" % locals(
                )
                timeRangeLabel = "MJD > %(mjdStart)s" % locals()
        elif plotType == "timeline" and timeLimitDay:
            start = timeLimitDay[0]
            end = timeLimitDay[1]
            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Transients\nDiscovered %(timeLimitLabel)s of Wave Detection" % locals()
            timeRangeLabel = timeLimitLabel.lower().replace(
                "in", "").replace("between", "").strip()
            if timeLimitLabel == "no limit":
                plotTitle = "%(gwid)s Probability Map, PS1 Footprints\n& All Transients Discovered" % locals(
                )
                timeRangeLabel = "all transients"
        elif plotType == "timeline":
            plotTitle = "%(gwid)s Probability Map" % locals()
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

            height = 2.8
            width = height / math.cos(decDeg * DEG_TO_RAD_FACTOR)

            # MULTIPLE CIRCLES
            if projection in ["tan"]:
                circ = Ellipse(
                    (raDeg, decDeg), width=width, height=height, alpha=0.2, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
            else:
                if raDeg > 180.:
                    raDeg = raDeg - 360.
                circ = Ellipse(
                    (-raDeg * DEG_TO_RAD_FACTOR, decDeg * DEG_TO_RAD_FACTOR), width=width * DEG_TO_RAD_FACTOR, height=height * DEG_TO_RAD_FACTOR, alpha=0.2, color='#859900', fill=True, zorder=3)

            ax.add_patch(circ)

        # LEGEND FOR PS1
        if len(ps1Pointings):
            raDeg = 90.
            decDeg = 10.
            height = 3.5

            width = height / math.cos(decDeg * DEG_TO_RAD_FACTOR)
            if projection in ["tan"]:
                circ = Ellipse(
                    (raDeg, decDeg), width=width, height=height, alpha=0.2, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
                ax.text(
                    raDeg - 3,
                    decDeg - 1.,
                    "PS1",
                    fontsize=14,
                    zorder=4,
                    family='monospace',
                    transform=ax.get_transform('fk5')
                )
            else:
                height = 6.
                width = height / math.cos(decDeg * DEG_TO_RAD_FACTOR)
                raDeg = 285.
                decDeg = 38.
                if raDeg > 180.:
                    raDeg = raDeg - 360.
                circ = Ellipse(
                    (-raDeg * DEG_TO_RAD_FACTOR, decDeg * DEG_TO_RAD_FACTOR), width=width * DEG_TO_RAD_FACTOR, height=height * DEG_TO_RAD_FACTOR, alpha=0.2, color='#859900', fill=True, zorder=3)
                ax.text(
                    (-raDeg + 8.) * DEG_TO_RAD_FACTOR,
                    (decDeg - 3.) * DEG_TO_RAD_FACTOR,
                    "PS1",
                    fontsize=14,
                    zorder=4,
                    family='monospace'
                )
            ax.add_patch(circ)

        # ADD ATLAS POINTINGS
        atlasPointingSide = 5.46
        for atp in atlasPointings:
            # add a path patch
            atlasExpId = atp["atlas_object_id"]
            raDeg = atp["raDeg"]
            decDeg = atp["decDeg"]

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
            elif probabilityCut:
                pass
                # print atlasExpId

            deltaDeg = atlasPointingSide / 2
            if decDeg < 0:
                deltaDeg = -deltaDeg

            if projection in ["tan"]:
                widthDegTop = atlasPointingSide / \
                    math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
                widthDegBottom = atlasPointingSide / \
                    math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
                heightDeg = atlasPointingSide
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
                patch = patches.PathPatch(path, alpha=0.2,
                                          color='#6c71c4', fill=True, zorder=3, transform=ax.get_transform('fk5'))
            else:
                if raDeg > 180.:
                    raDeg = raDeg - 360.

                widthRadTop = atlasPointingSide * DEG_TO_RAD_FACTOR / \
                    math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
                widthRadBottom = atlasPointingSide * DEG_TO_RAD_FACTOR / \
                    math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
                heightRad = atlasPointingSide * DEG_TO_RAD_FACTOR
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
                patch = patches.PathPatch(path, alpha=0.2,
                                          color='#6c71c4', fill=True, zorder=3,)

            ax.add_patch(patch)

        # LEGEND FOR ATLAS
        if len(atlasPointings):
            if projection in ["tan"]:
                raDeg = 88.
                decDeg = 4.
                atlasPointingSide = 4.5

                widthDegTop = atlasPointingSide / \
                    math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
                widthDegBottom = atlasPointingSide / \
                    math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
                heightDeg = atlasPointingSide
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
                patch = patches.PathPatch(path, alpha=0.2,
                                          color='#6c71c4', fill=True, zorder=3, transform=ax.get_transform('fk5'))
                ax.text(
                    raDeg - 4,
                    decDeg - 1,
                    "ATLAS",
                    fontsize=14,
                    zorder=4,
                    family='monospace',
                    transform=ax.get_transform('fk5')
                )
            else:
                atlasPointingSide = 8.
                raDeg = 285.
                decDeg = 28.

                if raDeg > 180.:
                    raDeg = raDeg - 360.

                widthRadTop = atlasPointingSide * DEG_TO_RAD_FACTOR / \
                    math.cos((decDeg + deltaDeg) * DEG_TO_RAD_FACTOR)
                widthRadBottom = atlasPointingSide * DEG_TO_RAD_FACTOR / \
                    math.cos((decDeg - deltaDeg) * DEG_TO_RAD_FACTOR)
                heightRad = atlasPointingSide * DEG_TO_RAD_FACTOR
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
                patch = patches.PathPatch(path, alpha=0.2,
                                          color='#6c71c4', fill=True, zorder=3,)
                ax.text(
                    (-raDeg + 8.) * DEG_TO_RAD_FACTOR,
                    (decDeg - 3.) * DEG_TO_RAD_FACTOR,
                    "ATLAS",
                    fontsize=14,
                    zorder=4,
                    family='monospace'
                )

            ax.add_patch(patch)

        # ADD DATA POINTS FOR TRANSIENTS
        names = []
        ra = []
        dec = []
        raRad = []
        decRad = []
        texts = []

        for trans in ps1Transients:
            # if trans["ps1_designation"] in ["PS15dpg", "PS15dpp", "PS15dpq", "PS15don", "PS15dpa", "PS15dom"]:
            #     continue

            name = trans["ps1_designation"]
            if name is None:
                name = trans["local_designation"]
            names.append(name)
            raDeg = trans["ra_psf"]
            decDeg = trans["dec_psf"]
            ra.append(raDeg)
            dec.append(decDeg)
            raRad.append(-raDeg * DEG_TO_RAD_FACTOR)
            decRad.append(decDeg * DEG_TO_RAD_FACTOR)

        if len(ra) > 0:
            # MULTIPLE CIRCLES
            if projection in ["tan"]:
                ax.scatter(
                    x=np.array(ra),
                    y=np.array(dec),
                    transform=ax.get_transform('fk5'),
                    s=6,
                    c='black',
                    edgecolor='black',
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
                        arrowprops=dict(arrowstyle="-", color='black', lw=1.2,
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
        if projection == "tan":
            plt.text(
                xRange * 0.25,
                # xRange * 0.95,
                yRange * 0.93,
                timeRangeLabel,
                fontsize=16,
                zorder=4,
                color="#dc322f",
                fontproperties=font
            )
            plt.text(
                xRange * 0.1,
                # xRange * 0.95,
                yRange * 0.93,
                "",
                fontsize=20,
                zorder=4,
                color="black",
                fontproperties=font
            )
        else:
            plt.text(
                fWidth * 0.7,
                # xRange * 0.95,
                fHeight * 0.95,
                timeRangeLabel,
                fontsize=16,
                zorder=4,
                color="#dc322f",
                fontproperties=font
            )
            plt.text(
                fWidth * 0.87,
                # xRange * 0.95,
                fHeight * 0.315,
                "",
                fontsize=20,
                zorder=4,
                color="black",
                fontproperties=font
            )

        # Recursively create missing directories
        if self.settings and not outputDirectory:
            plotDir = self.settings["output directory"] + "/" + gwid
        elif outputDirectory:
            plotDir = outputDirectory

        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = plotTitle.replace(" ", "_").replace(
            "<", "lt").replace(">", "gt").replace(",", "").replace("\n", "_").replace("&", "and")
        figureName = """%(plotTitle)s""" % locals(
        )
        if timeLimitDay == 0:
            figureName = """%(plotTitle)s""" % locals(
            )
        if plotDir != ".":
            for f in fileFormats:
                if not os.path.exists("%(plotDir)s/%(folderName)s/%(f)s" % locals()):
                    os.makedirs("%(plotDir)s/%(folderName)s/%(f)s" % locals())
                figurePath = "%(plotDir)s/%(folderName)s/%(f)s/%(figureName)s_%(projection)s.%(f)s" % locals()
                savefig(figurePath, bbox_inches='tight', dpi=300)

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
                figurePath = "%(plotDir)s/%(figureName)s.%(f)s" % locals()
                savefig(figurePath, bbox_inches='tight', dpi=300)

            # pathToExportFits = "%(plotDir)s/%(gwid)s_skymap.fits" % locals()
            # try:
            #     os.remove(pathToExportFits)
            # except:
            #     pass
            # hdu.writeto(pathToExportFits)

        self.generate_fits_image_map(
            gwid=gwid,
            pathToProbMap=pathToProbMap,
            folderName=folderName,
            outputDirectory=outputDirectory
        )

        self.log.info('completed the ``generate_probability_plot`` method')
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
                       projection="tan"
                )
                plotter.get()
        """
        self.log.info('starting the ``get_history_plots`` method')

        timeLimitLabels = ["day", "2 days", "3 days", "4 days", "5 days", "6 days",
                           "7 days", "2 weeks", "3 weeks", "1 month", "2 months", "3 months", "no limit"]
        timeLimitDays = [1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 0]

        if self.gwid:
            theseIds = [self.gwid]
        else:
            theseIds = self.settings["gravitational waves"]

        for gwid in theseIds:
            for tday, tlabel in zip(timeLimitDays, timeLimitLabels):

                plotParameters, ps1Transients, ps1Pointings, atlasPointings = self.get_gw_parameters_from_settings(
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
                    gwid]["time"]["mjdStart"]

                self.generate_probability_plot(
                    gwid=gwid,
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    ps1Pointings=ps1Pointings,
                    atlasPointings=atlasPointings,
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    timeLimitLabel=tlabel,
                    timeLimitDay=tday,
                    raLimit=False,
                    fileFormats=["png"],
                    folderName="survey_history_plots",
                    plotType=self.plotType,
                    projection=self.projection,
                    probabilityCut=self.probabilityCut)

        self.log.info('completed the ``get_history_plots`` method')
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
                       projection="tan"
                )
                plotter.get()
        """
        self.log.info('starting the ``get_timeline_plots`` method')

        timeLimitLabels = ["in First 3 Days",
                           "Between 3-10 Days", "Between 10-17 Days", "Between 17-24 Days", "Between 24-31 Days", "> 31 Days", "no limit"]
        timeLimitDays = [(0, 3), (3, 10), (10, 17),
                         (17, 24), (24, 31), (31, 0), (0, 0)]
        raLimits = [134.25, 144.75, 152.25, 159.50, 167.0, 174.5]

        # timeLimitLabels = ["in First 3 Days"]
        # timeLimitDays = [(2, 5)]

        if self.gwid:
            theseIds = [self.gwid]
        else:
            theseIds = self.settings["gravitational waves"]

        for gwid in theseIds:
            for tday, tlabel, raLimit in zip(timeLimitDays, timeLimitLabels, raLimits):

                plotParameters, ps1Transients, ps1Pointings, atlasPointings = self.get_gw_parameters_from_settings(
                    gwid=gwid,
                    inPastDays=False,
                    inFirstDays=tday)

                pathToProbMap = self.settings[
                    "gravitational waves"][gwid]["mapPath"]
                if not os.path.exists(pathToProbMap):
                    message = "the path to the map %s does not exist on this machine" % (
                        pathToProbMap,)
                    self.log.critical(message)
                    raise IOError(message)

                mjdStart = self.settings["gravitational waves"][
                    gwid]["time"]["mjdStart"]

                self.generate_probability_plot(
                    gwid=gwid,
                    plotParameters=plotParameters,
                    ps1Transients=ps1Transients,
                    ps1Pointings=ps1Pointings,
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
                    probabilityCut=self.probabilityCut)

        self.log.info('completed the ``get_timeline_plots`` method')
        return None

    def generate_fits_image_map(
            self,
            gwid,
            pathToProbMap,
            folderName="",
            outputDirectory=False,
            rebin=True):
        """*generate fits image map from the LV-skymap (FITS binary table)*

        **Key Arguments:**
            - ``pathToProbMap`` -- path to the FITS file containing the probability map of the wave
            - ``outputDirectory`` -- can be used to override the output destination in the settings file
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``folderName`` -- the name of the folder to add the plots to
            - ``rebin`` -- rebin the final image to reduce size

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
        self.log.info('starting the ``generate_fits_image_map`` method')

        import healpy as hp
        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        # X, Y PIXEL COORDINATE GRID
        xRange = 10000
        yRange = xRange * 1.7

        # PIXELSIZE AS MAPPED TO THE FULL SKY
        pixelSizeDeg = 360. / xRange

        # READ HEALPIX MAPS FROM FITS FILE
        # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
        # ARRAY OF PROBABILITIES (3,072 ROWS)
        # READ IN THE HEALPIX FITS FILE
        aMap, mapHeader = hp.read_map(pathToProbMap, 0, h=True, verbose=False)
        # DETERMINE THE SIZE OF THE HEALPIXELS
        nside = hp.npix2nside(len(aMap))
        centralCoordinate = [0, 0]

        # FROM THE PIXEL GRID (xRange, yRange), GENERATE A MAP TO LAT (pi to 0) AND LONG (-pi to pi) THAT CAN THEN MAPS TO HEALPIX SKYMAP
        # RA FROM -pi to pi
        phi = x2long(np.arange(xRange), xRange)
        # DEC FROM pi to 0
        theta = y2lat(np.arange(yRange), xRange, yRange)

        # PROJECT THE MAP TO A RECTANGULAR MATRIX xRange X yRange
        PHI, THETA = np.meshgrid(phi, theta)
        healpixIds = hp.ang2pix(nside, THETA, PHI)
        # GIVEN A HIGH ENOUGH RESOLUTION IN THE PIXEL GRID WE WILL HAVE
        # DUPLICATES - LET COUNT THEM ADD DIVIDE PROBABILITY EQUALLY
        unique, counts = np.unique(healpixIds, return_counts=True)
        countDict = dict(zip(unique, counts))

        probs = aMap[healpixIds]
        weightedProb = np.array([[probs[i, j] / countDict[healpixIds[i, j]] for j in xrange(probs.shape[1])]
                                 for i in xrange(probs.shape[0])])

        if rebin == True:
            rebinSize = 4
            # resize by getting rid of extra columns/rows
            xedge = np.shape(weightedProb)[0] % rebinSize
            yedge = np.shape(weightedProb)[1] % rebinSize
            weightedProb = weightedProb[xedge:, yedge:]

            # put image array into arrays of rebinSize x rebinSize
            weightedProb = np.reshape(weightedProb, (np.shape(weightedProb)[
                                      0] / rebinSize, rebinSize, np.shape(weightedProb)[1] / rebinSize, rebinSize))

            # average each rebinSize x rebinSize array
            weightedProb = np.mean(weightedProb, axis=3)
            weightedProb = np.mean(weightedProb, axis=1)

            pixelSizeDeg = pixelSizeDeg * rebinSize
            xRange = int(xRange / rebinSize)
            yRange = int(yRange / rebinSize)

        # CREATE A NEW WCS OBJECT
        w = awcs.WCS(naxis=2)
        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = np.array([pixelSizeDeg, pixelSizeDeg])
        # WORLD COORDINATES AT REFERENCE PIXEL
        w.wcs.crval = centralCoordinate

        # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
        w.wcs.crpix = [xRange / 2., yRange / 2.]

        unweightedImageProb = np.sum(probs)
        self.log.info(
            "The total unweighted probability flux in the FITS images added to %(unweightedImageProb)s" % locals())

        totalImageProb = np.sum(weightedProb)
        self.log.info(
            "The total probability flux in the FITS images added to %(totalImageProb)s" % locals())

        # CTYPE FOR THE FITS HEADER
        w.wcs.ctype = ["RA---MER" %
                       locals(), "DEC--MER" % locals()]

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
            pathToExportFits = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_skymap.fits" % locals()
            try:
                os.remove(pathToExportFits)
            except:
                pass
            hdu.writeto(pathToExportFits)
        else:
            pathToExportFits = "%(plotDir)s/%(gwid)s_skymap.fits" % locals()
            try:
                os.remove(pathToExportFits)
            except:
                pass
            hdu.writeto(pathToExportFits)

        self.log.info('completed the ``generate_fits_image_map`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method


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
