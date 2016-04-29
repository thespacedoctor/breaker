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
from dryxPython import astrotools as dat
from fundamentals import tools, times
from dryxPython import mysql as dms
from adjustText import adjust_text


class plot_wave_observational_timelines():
    """
    *TPlot the maps showing the timeline of observations of the gravitation wave sky-locations*

    You can plot either the history (looking back from now) or timeline (looking forward from date of GW detection) of the survey.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``plotType`` -- history (looking back from now) or timeline (looking forward from date of GW detection)
        - ``gwid`` -- a given graviational wave ID. If given only maps for this wave shall be plotted. Default *False* (i.e. plot all waves)

    **Usage:**

        To plot a the history of a specific wave:

        .. code-block:: python 

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="history",
                   gwid="G184098"
            )
            plotter.get()

        or to plot all waves in the settings file:

        .. code-block:: python 

            from breaker.plots import plot_wave_observational_timelines
            plotter = plot_wave_observational_timelines(
                log=log,
                   settings=settings,
                   plotType="history"
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
            gwid=False
    ):
        self.log = log
        log.debug("instansiating a new 'plot_wave_observational_timelines' object")
        self.settings = settings
        self.plotType = plotType
        self.gwid = gwid
        # xt-self-arg-tmpx

        # Initial Actions
        # CONNECT TO THE VARIOUS DATABASES REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = db.get()

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
                plotParameters, ps1Transients, ps1Pointings, altasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117"
                )
                print plotParameters

                # OUT: {'raRange': 90.0, 'centralCoordinate': [55.0, 27.5], 'decRange': 85.0}

                print ps1Transients

                # OUT: ({'local_designation': u'5L3Gbaj', 'ra_psf': 39.29767123836419, 'ps1_designation': u'PS15don', 'dec_psf': 19.055638423458053}, {'local_designation': u'5L3Gcbu', 'ra_psf': 40.06271352712189, 'ps1_designation': u'PS15dox', 'dec_psf': 22.536709765810823}, {'local_designation': u'5L3Gcca', 'ra_psf': 41.97569854977185, 'ps1_designation': u'PS15doy', 'dec_psf': 21.773344501616435}, {'local_designation': u'5L3Gbla', 'ra_psf': 50.732664347994714, 'ps1_designation': u'PS15dcq', 'dec_psf': 34.98988923347591}, {'local_designation': u'6A3Gcvu', 'ra_psf': 34.77565307934415, 'ps1_designation': u'PS16ku', 'dec_psf': 10.629310832257824}, {'local_designation': u'5L3Gcel', 'ra_psf': 38.24898916543392, 'ps1_designation': u'PS15dpn', 'dec_psf': 18.63530332013424}, {'local_designation': u'5L3Gcvk', 'ra_psf': 40.13754684778398, 'ps1_designation': u'PS15dpz', 'dec_psf': 23.003023065333267}, ....

                print ps1Pointings

                # OUT: ({'raDeg': 37.1814041667, 'mjd': 57388.2124067, 'decDeg': 18.9258969444}, {'raDeg': 37.1813666667, 'mjd': 57388.2140101, 'decDeg': 18.9259066667}, ...

            It can also be useful to give time-limits for the request to get the observations and discoveries from the past few days (``inPastDays``), or for the first few days after wave detection (``inFirstDays``). So for the past week:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, altasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inPastDays=7
                )

            Or the first 3 days since wave detection:

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, altasPointings = plotter.get_gw_parameters_from_settings(
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

        # GRAB PS1 & ATLAS POINTINGS FROM THE DATABASE
        ps1Pointings = self._get_ps1_pointings(gwid, inPastDays, inFirstDays)
        altasPointings = self._get_atlas_pointings(
            gwid, inPastDays, inFirstDays)

        self.log.info(
            'completed the ``get_gw_parameters_from_settings`` method')
        return plotParameters, ps1Transients, ps1Pointings, altasPointings

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
            now = dat.getCurrentMJD()
            mjdStart = now - inPastDays
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

        ps1Transients = dms.execute_mysql_read_query(
            sqlQuery=sqlQuery,
            dbConn=self.ps1gwDbConn,
            log=self.log
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
            nowMjd = dat.getCurrentMJD()
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

        ps1Pointings = dms.execute_mysql_read_query(
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log
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
            nowMjd = dat.getCurrentMJD()
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
            SELECT raDeg, decDeg, mjd FROM atlas_pointings where gw_id = "%(gwid)s" and mjd between %(mjdStart)s and %(mjdEnd)s group by atlas_object_id;
        """ % locals()

        atlasPointings = dms.execute_mysql_read_query(
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log
        )

        self.log.info('completed the ``_get_atlas_pointings`` method')
        return atlasPointings

    def generate_probability_plot(
            self,
            gwid,
            plotParameters,
            ps1Transients,
            ps1Pointings,
            atlasPointings,
            pathToProbMap,
            mjdStart,
            timeLimitLabel,
            timeLimitDay,
            fileFormats,
            folderName,
            plotType,
            raLimit=False):
        """
        *Generate a single probability map plot for a given gravitational wave and save it to file*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``plotParameters`` -- the parameters of the plot (for spatial & temporal parameters etc)
            - ``ps1Transients`` -- the transients to add to the plot
            - ``ps1Pointings`` -- the PS1 pointings to place on the plot
            - ``atlasPointings`` -- the atlas pointings to add to the plot
            - ``pathToProbMap`` -- path to the FITS file containing the probability map of the wave
            - ``mjdStart`` -- earliest mjd of discovery
            - ``timeLimitLabel`` -- the labels of the time contraints (for titles)
            - ``timeLimitDay`` -- the time limits (in ints)
            - ``raLimit`` -- ra limit at twilight
            - ``fileFormats`` -- the format(s) to output the plots in (list of strings)
            - ``folderName`` -- the name of the folder to add the plots to
            - ``plotType`` -- history (looking back from now) or timeline (looking forward from date of GW detection)


        **Return:**
            - None

        **Usage:**
            .. todo::

                - add usage info
                - create a sublime snippet for usage

            First you neeed to collect your data and a few plot parameters:

            .. code-block:: python 

                from breaker.plots import plot_wave_observational_timelines
                plotter = plot_wave_observational_timelines(
                    log=log,
                    settings=settings
                )
                plotParameters, ps1Transients, ps1Pointings, altasPointings = plotter.get_gw_parameters_from_settings(
                    gwid="G211117",
                    inFirstDays=(0,7)
                )

        .. todo::

            - @review: when complete, clean methodName method
            - @review: when complete add logging
        """
        self.log.info('starting the ``generate_probability_plot`` method')

        from matplotlib.font_manager import FontProperties
        font = FontProperties()
        font.set_family("Arial")

        pixelSizeDeg = 0.066667

        # UNPACK THE PLOT PARAMETERS
        centralCoordinate = plotParameters["centralCoordinate"]
        raRange = plotParameters["raRange"]
        decRange = plotParameters["decRange"]

        raMax = centralCoordinate[0] + raRange / 2.
        raMin = centralCoordinate[0] - raRange / 2.
        decMax = centralCoordinate[1] + decRange / 2.
        decMin = centralCoordinate[1] - decRange / 2.

        # DETERMINE THE PIXEL GRID X,Y RANGES
        xRange = int(raRange / pixelSizeDeg)
        yRange = int(decRange / pixelSizeDeg)

        # CREATE A NEW WCS OBJECT
        import numpy
        from astropy import wcs as awcs
        from astropy.io import fits
        w = awcs.WCS(naxis=2)

        # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
        w.wcs.crpix = [xRange / 2., yRange / 2.]
        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = numpy.array([pixelSizeDeg, pixelSizeDeg])
        # WORLD COORDINATES AT REFERENCE PIXEL
        w.wcs.crval = centralCoordinate
        # USE THE "GNOMONIC" PROJECTION ("COORDINATESYS---PROJECTION")
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # CREATE THE FITS HEADER WITH WCS
        header = w.to_header()

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

        # READ HEALPIX MAPS FROM FITS FILE
        # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
        # ARRAY OF PROBABILITIES (3,072 ROWS)
        import healpy as hp
        aMap, mapHeader = hp.read_map(pathToProbMap, h=True)
        nside = hp.pixelfunc.get_nside(aMap)

        # import matplotlib.pyplot as plt
        # hp.mollview(aMap, title="mollview image RING", cmap="YlOrRd")
        # hp.graticule()
        # plt.show()

        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

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

        totalProb = sum(aMap)
        print "Total Probability for the entire sky is %(totalProb)s" % locals()
        stampProb = np.sum(uniProb)
        print "Probability for the plot stamp is %(stampProb)s" % locals()

        # RESHAPE THE ARRAY AS BITMAP
        probs = np.reshape(np.array(probs), (yRange, xRange))

        # CREATE THE FITS FILE
        hdu = fits.PrimaryHDU(header=header, data=probs)

        # CONTOURS - NEED TO ADD THE CUMMULATIVE PROBABILITY
        i = np.flipud(np.argsort(aMap))
        cumsum = np.cumsum(aMap[i])
        cls = np.empty_like(aMap)
        cls[i] = cumsum * 100 * stampProb

        # EXTRACT CONTOUR VALUES AT HEALPIX INDICES
        contours = []
        contours[:] = [cls[i] for i in healpixIds]
        contours = np.reshape(np.array(contours), (yRange, xRange))

        # GRAB THE WCS FROM HEADER GENERATED EARLIER
        from wcsaxes import datasets, WCS
        wcs = WCS(hdu.header)

        # PLOT MAP WITH PROJECTION IN HEADER
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
        im = ax.imshow(probs,
                       cmap="YlOrRd", origin='lower', alpha=0.7, zorder=1)

        # PLOT THE CONTOURS ON THE SAME PLOT
        CS = plt.contour(contours, linewidths=1, alpha=0.5, zorder=2)
        plt.clabel(CS, fontsize=12, inline=1, fmt='%2.1f', fontproperties=font)

        # RESET THE AXES TO THE FRAME OF THE FITS FILE
        ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
        ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

        # CLIP THE IMAGE TO THE FRAME
        # im.set_clip_path(ax.coords.frame.patch)

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

        # CUSTOMISE TICK POSITIONS (l, b, r, t == left, bottom, right, or top)
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

        # ADD RA LIMIT
        if raLimit:
            x = np.ones(100) * raLimit
            print x
            y = np.linspace(decMin, decMax, 100)
            print y
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
        elif plotType == "timeline":
            start = timeLimitDay[0]
            end = timeLimitDay[1]
            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Transients\nDiscovered %(timeLimitLabel)s of Wave Detection" % locals()
            timeRangeLabel = timeLimitLabel.lower().replace(
                "in", "").replace("between", "").strip()
            if timeLimitLabel == "no limit":
                plotTitle = "%(gwid)s Probability Map, PS1 Footprints\n& All Transients Discovered" % locals(
                )
                timeRangeLabel = "all transients"

        else:
            timeRangeLabel = ""

        subTitle = "(updated %(now)s)" % locals()
        if timeLimitDay == 0 or plotType == "timeline":
            subTitle = ""

        # ax.set_title(plotTitle + "\n", fontsize=10)

        # GRAB PS1 POINTINGS
        pointingArray = []

        from matplotlib.patches import Circle

        for psp in ps1Pointings:
            raDeg = psp["raDeg"]
            decDeg = psp["decDeg"]

            # MULTIPLE CIRCLES
            circ = Circle(
                (raDeg, decDeg), radius=1.4, alpha=0.2, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
            ax.add_patch(circ)

        # ADD ATLAS POINTINGS
        for atp in atlasPointings:
            raDeg = atp["raDeg"]
            decDeg = atp["decDeg"]

            # MULTIPLE CIRCLES
            circ = Circle(
                (raDeg, decDeg), radius=1.4, alpha=0.2, color='#6c71c4', fill=True, transform=ax.get_transform('fk5'), zorder=3)
            ax.add_patch(circ)

        # ADD DATA POINTS FOR TRANSIENTS
        names = []
        ra = []
        dec = []
        xLabelShift = []
        yLabelShift = []

        movedSourcesX = []
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

            #########
            import math
            pi = (4 * math.atan(1.0))
            DEG_TO_RAD_FACTOR = pi / 180.0
            raDecFactor = math.cos(centralCoordinate[1] * DEG_TO_RAD_FACTOR)
            shiftX = -0.01 * raRange * raDecFactor
            shiftY = -0.005 * decRange

            movedx = False
            movedy = False

            for otherTrans in ps1Transients:
                otherName = otherTrans["ps1_designation"]
                if otherName is None:
                    otherName = otherTrans["local_designation"]
                if name == otherName:
                    continue
                otherTranRa, otherTranDec = (
                    otherTrans["ra_psf"], otherTrans["dec_psf"])
                raDiff = (raDeg - otherTranRa) * raDecFactor
                decDiff = decDeg - otherTranDec

                if abs(raDiff) < 0.07 * raRange * raDecFactor and abs(decDiff) < 0.025 * decRange and otherName not in movedSourcesX:
                    if raDiff > 0 and movedx == False:
                        shiftX = shiftX + 0.09 * raRange * raDecFactor
                        movedx = True
                        movedSourcesX.append(name)
                    if decDiff > 0 and movedy == False:
                        shiftY = shiftY + 0.01 * decRange
                        movedy = True
                    elif decDiff < 0 and movedy == False:
                        shiftY = shiftY - 0.005 * decRange
                        movedy = True

            xLabelShift.append(shiftX)
            yLabelShift.append(shiftY)

        if len(ra) > 0:
            ax.scatter(
                x=np.array(ra),
                y=np.array(dec),
                transform=ax.get_transform('fk5'),
                s=6,
                c='#dc322f',
                edgecolor='#dc322f',
                alpha=1,
                zorder=4
            )

        texts = []
        if len(ra):
            xx, yy = w.wcs_world2pix(np.array(ra), np.array(dec), 1)

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
                    expand_points=(1.2, 1.2),
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
                    arrowprops=dict(arrowstyle="-", color='#dc322f', lw=0.6,
                                    patchB=None, shrinkB=0, connectionstyle="arc3,rad=0.1", zorder=3, alpha=0.5),
                    fontsize=10,
                    family='monospace'
                )

        # # TIME-RANGE LABEL
        ax.text(
            xRange * 0.2,
            # xRange * 0.95,
            yRange * 0.93,
            timeRangeLabel,
            fontsize=16,
            zorder=4,
            color="#dc322f",
            fontproperties=font
        )

        # Recursively create missing directories
        plotDir = self.settings["output directory"] + "/" + gwid
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = plotTitle.replace(" ", "_").replace(
            "<", "lt").replace(">", "gt").replace(",", "").replace("\n", "_").replace("&", "and")
        figureName = """%(plotTitle)s""" % locals(
        )
        if timeLimitDay == 0:
            figureName = """%(plotTitle)s""" % locals(
            )
        for f in fileFormats:
            if not os.path.exists("%(plotDir)s/%(folderName)s/%(f)s" % locals()):
                os.makedirs("%(plotDir)s/%(folderName)s/%(f)s" % locals())
            figurePath = "%(plotDir)s/%(folderName)s/%(f)s/%(figureName)s.%(f)s" % locals()
            savefig(figurePath, bbox_inches='tight', dpi=300)

        if not os.path.exists("%(plotDir)s/%(folderName)s/fits" % locals()):
            os.makedirs("%(plotDir)s/%(folderName)s/fits" % locals())
        pathToExportFits = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_stamp.fits" % locals()
        try:
            os.remove(pathToExportFits)
        except:
            pass
        hdu.writeto(pathToExportFits)

        self.log.info('completed the ``generate_probability_plot`` method')
        return None

    def get_history_plots(
            self):
        """
        *plot the history plots*

        **Key Arguments:**
            # -

        **Return:**
            - None

        **Usage:**
            .. todo::

                - add usage info
                - create a sublime snippet for usage

            .. code-block:: python 

                usage code 

        .. todo::

            - @review: when complete, clean methodName method
            - @review: when complete add logging
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

                plotParameters, ps1Transients, ps1Pointings, altasPointings = self.get_gw_parameters_from_settings(
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
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    timeLimitLabel=tlabel,
                    timeLimitDay=tday,
                    fileFormats=["png"],
                    folderName="survey_history_plots",
                    plotType=self.plotType)

        self.log.info('completed the ``get_history_plots`` method')
        return None

    def get_timeline_plots(
            self):
        """
        *plot the history plots*

        **Key Arguments:**
            # -

        **Return:**
            - None

        **Usage:**
            .. todo::

                - add usage info
                - create a sublime snippet for usage

            .. code-block:: python 

                usage code 

        .. todo::

            - @review: when complete, clean methodName method
            - @review: when complete add logging
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

                plotParameters, ps1Transients, ps1Pointings, altasPointings = self.get_gw_parameters_from_settings(
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
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    timeLimitLabel=tlabel,
                    timeLimitDay=tday,
                    raLimit=raLimit,
                    fileFormats=["png", "pdf"],
                    folderName="survey_timeline_plots",
                    plotType=self.plotType)

        self.log.info('completed the ``get_timeline_plots`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
