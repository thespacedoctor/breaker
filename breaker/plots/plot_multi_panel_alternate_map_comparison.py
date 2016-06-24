#!/usr/local/bin/python
# encoding: utf-8
"""
*Given a directory of maps, plot a multipanel plot with probabilities plotted on the same colorbar scale*

:Author:
    David Young

:Date Created:
    February  5, 2016
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
from matplotlib.path import Path
from matplotlib.pyplot import savefig
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from fundamentals import tools, times
from crowdedText import adjust_text
from breaker.plots.plot_wave_observational_timelines import plot_wave_observational_timelines


class plot_multi_panel_alternate_map_comparison():
    """
    *Given a directory of maps, plot a multipanel plot with probabilities plotted on the same colorbar scale*

    All files with the ``.fits`` extension in the map directory are assumed to be a healpix likelihood map.

    Settings for the maps ratios etc are lifted from the setttings file.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the ID of the gravitational wave (matched in settings file)
        - ``pathToMapDirectory`` -- path to a directory containing the maps

    **Usage:**

        .. code-block:: python 

            from breaker.plots import plot_multi_panel_alternate_map_comparison
            p = plot_multi_panel_alternate_map_comparison(
                log=log,
                settings=settings,
                gwid="G211117",
                pathToMapDirectory="/path/to/my/maps"
            ).get()

    The resulting PNG will look something like this for 2 maps:

     .. image:: https://i.imgur.com/y6WBrwa.png
        :width: 800px
        :alt: GW151226 Skymap Comparison

    More panels are added depending on the number of maps in ``pathToMapDirectory``. Here's an example of a 4 map comparison:


    """
    # Initialisation

    def __init__(
            self,
            log,
            gwid,
            pathToMapDirectory,
            settings=False,

    ):
        self.log = log
        log.debug(
            "instansiating a new 'plot_multi_panel_alternate_map_comparison' object")
        self.settings = settings
        self.gwid = gwid
        self.pathToMapDirectory = pathToMapDirectory
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
        *Generate the plot*
        """
        self.log.info('starting the ``get`` method')

        # GRAB THE SURVEY DETAILS FOR THE PLOTS
        p = plot_wave_observational_timelines(
            log=self.log,
            settings=self.settings,
            plotType="history"
        )

        plotParameters, ps1Transients, ps1Pointings, altasPointings = p.get_gw_parameters_from_settings(
            gwid=self.gwid,
            inPastDays=30000,
            inFirstDays=False
        )

        # AND GENERATE THE PLOT
        self._generate_map_comparison_plot(
            gwid=self.gwid,
            plotParameters=plotParameters,
            ps1Transients=False,
            ps1Pointings=ps1Pointings,
            pathToMapDirectory=self.pathToMapDirectory,
            fileFormats=["png", ],
            folderName="map_comparisons"
        )

        self.log.info('completed the ``get`` method')
        return plot_multi_panel_alternate_map_comparison

    def _generate_map_comparison_plot(
            self,
            gwid,
            plotParameters,
            ps1Transients,
            ps1Pointings,
            pathToMapDirectory,
            fileFormats,
            folderName):
        """
        *generate probability plot*

        **Key Arguments:**
            - ``gwid`` -- the unique ID of the gravitational wave to plot
            - ``plotParameters`` -- the parameters of the plot (for spatial & temporal parameters etc)
            - ``ps1Transients`` -- the transients to add to the plot
            - ``ps1Pointings`` -- the pointings to place on the plot
            - ``pathToMapDirectory`` -- path to a directory containing the maps
            - ``fileFormats`` -- the format(s) to output the plots in (list of strings)
            - ``folderName`` -- the name of the folder to add the plots to

        **Return:**
            - None
        """
        self.log.info('starting the ``_generate_map_comparison_plot`` method')

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

        basePath = pathToMapDirectory
        mapPaths = []
        for d in os.listdir(basePath):
            if os.path.isfile(os.path.join(basePath, d)) and ".fits" in d[-5:]:
                mapPaths.append(os.path.join(basePath, d))
                self.log.debug('found the map `%(d)s`' % locals())

        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        # FIND THE MAX PROB VALUE
        maxProb = 0.0
        for mapI in range(1, len(mapPaths) - 1):
            import healpy as hp
            pathToProbMap = mapPaths[mapI]
            aMap, mapHeader = hp.read_map(pathToProbMap, h=True)
            maxValue = np.max(aMap)
            if maxValue > maxProb:
                maxProb = maxValue

        # mapPaths = list(reversed(mapPaths))
        for mapI in range(len(mapPaths)):

            pathToProbMap = mapPaths[mapI]
            self.log.debug(
                'starting to work on  `%(pathToProbMap)s`' % locals())

            mapName = pathToProbMap.split("/")[-1]
            mapName = mapName.replace(".fits", "").replace("_", " ")

            # READ HEALPIX MAPS FROM FITS FILE
            # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
            # ARRAY OF PROBABILITIES (3,072 ROWS)
            import healpy as hp
            aMap, mapHeader = hp.read_map(pathToProbMap, h=True)
            nside = hp.pixelfunc.get_nside(aMap)

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

            if mapI == 0:
                # PLOT MAP WITH PROJECTION IN HEADER
                import matplotlib.pyplot as plt
                fig = plt.figure()
                #  [0.15, 0.1, 0.8, 0.8]

                fig, axes = plt.subplots(
                    figsize=(12, 8), nrows=1, ncols=len(mapPaths), subplot_kw=dict(projection=wcs))

            ax = axes[mapI]
            # YlOrRd  dRrOlY
            im = ax.imshow(probs,
                           cmap="YlOrRd", origin='lower', alpha=0.7, zorder=1,  norm=matplotlib.colors.PowerNorm(0.5, vmin=0.0, vmax=maxProb))

            # PLOT THE CONTOURS ON THE SAME PLOT
            CS = ax.contour(contours, linewidths=1, alpha=1.0, zorder=6)
            ax.clabel(CS, CS.levels[-1:], fontsize=8,
                      inline=1, fmt='%2.1f', zorder=7)

            # RESET THE AXES TO THE FRAME OF THE FITS FILE
            ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
            ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

            # CLIP THE IMAGE TO THE FRAME
            # im.set_clip_path(ax.coords.frame.patch)

            # THE COORDINATES USED IN THE PLOT CAN BE ACCESSED USING THE COORDS
            # ATTRIBUTE (NOT X AND Y)
            lon = ax.coords[0]
            lat = ax.coords[1]
            # lon.set_axislabel('RA (deg)', minpad=0.5, fontsize=12)
            lat.set_axislabel('DEC (deg)', minpad=0.5, fontsize=12)
            lon.set_major_formatter('d.d')
            lat.set_major_formatter('d.d')

            # THE SEPARATORS FOR ANGULAR COORDINATE TICK LABELS CAN ALSO BE SET BY
            # SPECIFYING A STRING
            lat.set_separator(':-s')
            # SET THE APPROXIMATE NUMBER OF TICKS, WITH COLOR & PREVENT OVERLAPPING
            # TICK LABELS FROM BEING DISPLAYED.
            lon.set_ticks(number=4, color='#657b83', exclude_overlapping=True)
            lat.set_ticks(number=10, color='#657b83', exclude_overlapping=True)

            # MINOR TICKS NOT SHOWN BY DEFAULT
            lon.display_minor_ticks(True)
            lat.display_minor_ticks(True)
            lat.set_minor_frequency(2)

            # CUSTOMISE TICK POSITIONS (l, b, r, t == left, bottom, right, or
            # top)
            lon.set_ticks_position('bt')
            lon.set_ticklabel_position('b')
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
            ax.invert_xaxis()

            # SETUP TITLE OF PLOT
            plotTitle = "%(gwid)s multi panel skymap comparison" % locals()

            # GRAB PS1 POINTINGS
            pointingArray = []

            for psp in ps1Pointings:
                raDeg = psp["raDeg"]
                decDeg = psp["decDeg"]

                # MULTIPLE CIRCLES
                from matplotlib.patches import Circle
                circ = Circle(
                    (raDeg, decDeg), radius=1.4, alpha=0.015, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
                ax.add_patch(circ)

            # # MAP NAME LABEL
            if len(mapPaths) < 3:
                textFontsize = 12
            elif len(mapPaths) == 3:
                textFontsize = 9
            else:
                textFontsize = 8
            ax.text(
                xRange * 0.07,
                # xRange * 0.95,
                yRange * 0.93,
                mapName,
                fontsize=textFontsize,
                zorder=4,
                color="#dc322f",
                horizontalalignment="right"
            )

        for mapI in range(1, len(mapPaths)):
            ax = axes[mapI]
            lat = ax.coords[1]
            lat.set_ticklabel_position('')
            lat.set_axislabel('', minpad=0.5, fontsize=12)

        middle = int(len(mapPaths) / 2.)
        ax = axes[middle]
        lon = ax.coords[0]
        lon.set_axislabel('RA (deg)', minpad=0.8, fontsize=12)

        # Recursively create missing directories
        plotDir = self.settings["output directory"] + "/" + gwid
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = plotTitle.replace(" ", "_").replace(
            "<", "lt").replace(">", "gt").replace(",", "").replace("\n", "_").replace("&", "and")
        figureName = """%(plotTitle)s""" % locals(
        )

        for f in fileFormats:
            if not os.path.exists("%(plotDir)s/%(folderName)s/%(f)s" % locals()):
                os.makedirs("%(plotDir)s/%(folderName)s/%(f)s" % locals())
            figurePath = "%(plotDir)s/%(folderName)s/%(f)s/%(figureName)s.%(f)s" % locals()
            savefig(figurePath, bbox_inches='tight', dpi=300)
            #savefig(figurePath, dpi=300)

        if not os.path.exists("%(plotDir)s/%(folderName)s/fits" % locals()):
            os.makedirs("%(plotDir)s/%(folderName)s/fits" % locals())
        pathToExportFits = "%(plotDir)s/%(folderName)s/fits/%(gwid)s_stamp.fits" % locals()
        try:
            os.remove(pathToExportFits)
        except:
            pass
        hdu.writeto(pathToExportFits)

    # use the tab-trigger below for new method
    # xt-class-method
