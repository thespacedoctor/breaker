#!/usr/local/bin/python
# encoding: utf-8
"""
*Plot the gravitation wave sky-location maps showing catalogued sources matched within the survey footprint*
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
import matplotlib.patches as patches
from HMpTy import htm
from fundamentals import tools, times
from fundamentals.mysql import readquery
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from astrocalc.times import now


class plot_wave_matched_source_maps():
    """
    *Plot the gravitation wave sky - location maps showing catalogued sources matched within the survey footprint*

    **Key Arguments: **
        - ``log`` - - logger
        - ``settings`` - - the settings dictionary
        - ``gwid`` - - the wave id

    **Usage: **

        .. code-block: : python

            from breaker import plot_wave_matched_source_maps
            p = plot_wave_matched_source_maps(
                log=log,
                settings=settings,
                gwid="G211117"
            )
            p.get_source_plots()
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            gwid=False
    ):
        self.log = log
        log.debug("instansiating a new 'plot_wave_matched_source_maps' object")
        self.settings = settings
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
        *get the plot_wave_matched_source_maps object - history or timeline*
        """
        self.log.info('starting the ``get`` method')

        self.get_source_plots()

        self.log.info('completed the ``get`` method')
        return None

    def _get_ps1_pointings(
            self,
            gwid,
            inPastDays,
            inFirstDays):
        """
        *get ps1 pointings to add to the plot*

        **Key Arguments: **
            - ``gwid`` - - the unique ID of the gravitational wave to plot
            - ``inPastDays`` - - used for the `history` plots(looking back from today)
            - ``inFirstDays`` - - used in the `timeline` plots(looking forward from wave detection)

        **Return: **
            - ``ps1Pointings`` - - the pointings to place on the plot
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
        else:
            print inPastDays

        if inFirstDays:
            mjdStart = self.settings["gravitational waves"][gwid]["time"][
                "mjdStart"] + inFirstDays[0]
            mjdEnd = self.settings["gravitational waves"][
                gwid]["time"]["mjdStart"] + inFirstDays[1]
            if inFirstDays[1] == 0:
                mjdEnd = 10000000000

        sqlQuery = u"""
            SELECT raDeg, decDeg FROM ps1_pointings where gw_id = "%(gwid)s" and mjd between %(mjdStart)s and %(mjdEnd)s
        """ % locals()
        ps1Pointings = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.info('completed the ``_get_ps1_pointings`` method')
        return ps1Pointings

    def _generate_probability_plot(
            self,
            gwid,
            plotParameters,
            ps1Pointings,
            pathToProbMap,
            mjdStart,
            raArray,
            decArray,
            fileFormats,
            folderName,
            plotTitle):
        """
        *generate probability plot*

        **Key Arguments: **
            - ``gwid`` - - the unique ID of the gravitational wave to plot
            - ``plotParameters`` - - the parameters of the plot(for spatial & temporal parameters etc)
            - ``ps1Pointings`` - - the pointings to place on the plot
            - ``pathToProbMap`` - - path to the FITS file containing the probability map of the wave
            - ``mjdStart`` - - earliest mjd of discovery
            - ``raArray`` - - the array of matched source RAs
            - ``decArray`` - - the array of matched source DECs
            - ``fileFormats`` - - the format(s) to output the plots in (list of strings)
            - ``folderName`` - - the name of the folder to add the plots to
            - ``plotTitle`` - - the plot title

        **Return: **
            - None
        """
        self.log.info('starting the ``_generate_probability_plot`` method')

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
        plt.clabel(CS, fontsize=7, inline=1, fmt='%2.1f')

        # RESET THE AXES TO THE FRAME OF THE FITS FILE
        ax.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
        ax.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

        # CLIP THE IMAGE TO THE FRAME
        # im.set_clip_path(ax.coords.frame.patch)

        # THE COORDINATES USED IN THE PLOT CAN BE ACCESSED USING THE COORDS
        # ATTRIBUTE (NOT X AND Y)
        lon = ax.coords[0]
        lat = ax.coords[1]
        lon.set_axislabel('RA (deg)', minpad=0.5, fontsize=12)
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

        # CUSTOMISE TICK POSITIONS (l, b, r, t == left, bottom, right, or top)
        lon.set_ticks_position('bt')
        lon.set_ticklabel_position('b')
        lon.set_axislabel_position('b')
        lat.set_ticks_position('lr')
        lat.set_ticklabel_position('l')
        lat.set_axislabel_position('l')

        # ADD A GRID
        ax.coords.grid(color='#657b83', alpha=0.5, linestyle='dashed')
        plt.gca().invert_xaxis()

        #######

        # GRAB PS1 POINTINGS
        pointingArray = []

        for psp in ps1Pointings:
            raDeg = psp["raDeg"]
            decDeg = psp["decDeg"]

            # MULTIPLE CIRCLES

            circ = Circle(
                (raDeg, decDeg), radius=1.4, alpha=0.08, color='#859900', fill=True, transform=ax.get_transform('fk5'), zorder=3)
            ax.add_patch(circ)

        # SOME EXTRA CIRCLES -- IF NEEDED
        otherCircles1 = [
            {"raDeg": 140.0,
             "decDeg": 6.0},
            {"raDeg": 153.0,
             "decDeg": -12.0},
        ]
        otherCircles2 = [
            {"raDeg": 149.0,
             "decDeg": 2.0},
            {"raDeg": 152.0,
             "decDeg": -7.0},
        ]
        otherCircles3 = [
            {"raDeg": 132.0,
             "decDeg": 4.0},
            {"raDeg": 154.0,
             "decDeg": -16.5},
        ]
        final = [
            {"raDeg": 140.0,
             "decDeg": 6.0},
            {"raDeg": 154.0,
             "decDeg": -16.5},
        ]
        otherCircles = final
        for circle in otherCircles:
            raDeg = circle["raDeg"]
            decDeg = circle["decDeg"]
            circ = Circle(
                (raDeg, decDeg), radius=2.5, alpha=1, color='#100983', fill=False, transform=ax.get_transform('fk5'), zorder=7)
            # ax.add_patch(circ)

        if len(raArray) > 1000000:
            alpha = 0.12
            s = 0.05
        else:
            s = 0.5
            alpha = 0.5

        print plotTitle
        print len(raArray)
        print alpha

        if len(raArray):
            ax.scatter(
                x=raArray,
                y=decArray,
                transform=ax.get_transform('fk5'),
                s=s,
                c="#dc322f",
                alpha=alpha,
                zorder=4,
                lw=0
            )

        # ANNOTATIONS
        if True == True:
            finRA, finDec = (150.5, 2.5)
            startRA, startDec = (159.0, 2.5)
            ax.annotate('COSMOS',
                        color='#100983',
                        xy=w.wcs_world2pix(finRA, finDec, 0),
                        xytext=w.wcs_world2pix(startRA, startDec, 0),
                        arrowprops=dict(arrowstyle='->',
                                        facecolor='#100983', ec='#100983',  zorder=10),
                        xycoords=ax.transData,
                        textcoords=ax.transData,
                        horizontalalignment='left',
                        verticalalignment='top',
                        zorder=10
                        )

            finRA, finDec = (151, -12)
            startRA, startDec = (138.5, -8.)
            ax.annotate('LCRS',
                        color='#fffdd5',
                        xy=w.wcs_world2pix(finRA, finDec, 0),
                        xytext=w.wcs_world2pix(startRA, startDec, 0),
                        arrowprops=dict(arrowstyle='->',
                                        facecolor='#100983', ec='#100983', zorder=10),
                        xycoords=ax.transData,
                        textcoords=ax.transData,
                        horizontalalignment='right',
                        verticalalignment='top',
                        zorder=10
                        )

            finRA, finDec = (147.5, -5)
            startRA, startDec = (138.0, -7)
            ax.annotate('LCRS',
                        color='#100983',
                        xy=w.wcs_world2pix(finRA, finDec, 0),
                        xytext=w.wcs_world2pix(startRA, startDec, 0),
                        arrowprops=dict(arrowstyle='->',
                                        facecolor='#100983', ec='#100983',  zorder=10),
                        xycoords=ax.transData,
                        textcoords=ax.transData,
                        horizontalalignment='right',
                        verticalalignment='top',
                        zorder=10
                        )

            finRA, finDec = (154.5, -10.5)
            startRA, startDec = (163, -7.5)
            ax.annotate('WINGS',
                        color='#fffdd5',
                        xy=w.wcs_world2pix(finRA, finDec, 0),
                        xytext=w.wcs_world2pix(startRA, startDec, 0),
                        arrowprops=dict(arrowstyle='->',
                                        facecolor='#100983', ec='#100983',  zorder=10),
                        xycoords=ax.transData,
                        textcoords=ax.transData,
                        horizontalalignment='left',
                        verticalalignment='top',
                        zorder=10
                        )

            finRA, finDec = (152, -6)
            startRA, startDec = (162, -7)
            ax.annotate('WINGS',
                        color='#100983',
                        xy=w.wcs_world2pix(finRA, finDec, 0),
                        xytext=w.wcs_world2pix(startRA, startDec, 0),
                        arrowprops=dict(arrowstyle='->',
                                        facecolor='#100983', ec='#100983',  zorder=10),
                        xycoords=ax.transData,
                        textcoords=ax.transData,
                        horizontalalignment='left',
                        verticalalignment='top',
                        zorder=10
                        )

            lc = (162, -3)
            rc = (122, 18)
            width = rc[0] - lc[0]
            height = rc[1] - lc[1]

            ax.add_patch(Rectangle((lc),
                                   width, height, fill=False, linestyle='dashed', ec='#100983', zorder=10, transform=ax.get_transform('fk5')))
            r, d = (161, 16)
            ax.text(r, d, r'SDSS DR6', color='#100983',
                    transform=ax.get_transform('fk5'))

        # Recursively create missing directories
        plotDir = self.settings["output directory"] + "/" + gwid
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = plotTitle.replace(" ", "_").replace(
            "<", "lt").replace(">", "gt").replace(",", "").replace("\n", "_").replace("&", "").replace("__", "_")
        figureName = """ %(plotTitle)s""" % locals(
        )

        for f in fileFormats:
            if not os.path.exists("%(plotDir)s/%(folderName)s/%(f)s" % locals()):
                os.makedirs("%(plotDir)s/%(folderName)s/%(f)s" % locals())
            figurePath = "%(plotDir)s/%(folderName)s/%(f)s/%(figureName)s.%(f)s" % locals()
            savefig(figurePath, bbox_inches='tight', dpi=300)

        self.log.info('completed the ``_generate_probability_plot`` method')
        return None

    def get_source_plots(
            self):
        """
        *plot the history plots*
        """
        self.log.info('starting the ``get_source_plots`` method')

        if self.gwid:
            theseWaves = [self.gwid]
        else:
            theseWaves = self.settings["gravitational waves"]

        for gwid in theseWaves:

            plotParameters = self.settings["gravitational waves"][gwid]["plot"]
            ps1Pointings = self._get_ps1_pointings(
                gwid,
                inPastDays=False,
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

            r = 0.15
            raArray, decArray = self._get_matched_sources(
                gwid,
                plotParameters,
                redshiftLimit=r)

            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Matched Catalogue Sources" % locals(
            )

            if r is not False:
                plotTitle += " Within z < %(r)s" % locals(
                )

            self._generate_probability_plot(
                gwid=gwid,
                plotParameters=plotParameters,
                ps1Pointings=ps1Pointings,
                pathToProbMap=pathToProbMap,
                mjdStart=mjdStart,
                raArray=raArray,
                decArray=decArray,
                fileFormats=["png", "pdf"],
                folderName="survey_matched_source_maps",
                plotTitle=plotTitle
            )

            return None

            # EMPTY PLOT - NO SOURCES
            plotTitle = "%(gwid)s Probability Map, PS1 Footprints" % locals(
            )

            self._generate_probability_plot(
                gwid=gwid,
                plotParameters=plotParameters,
                ps1Pointings=ps1Pointings,
                pathToProbMap=pathToProbMap,
                mjdStart=mjdStart,
                raArray=[],
                decArray=[],
                fileFormats=["png", "pdf"],
                folderName="survey_matched_source_maps",
                plotTitle=plotTitle
            )

            # PLOT ALL NED SOURCES
            raArray, decArray = self._get_matched_sources(
                gwid,
                plotParameters,
                redshiftLimit=False,
                allNed=True
            )

            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & All NED Sources" % locals(
            )

            self._generate_probability_plot(
                gwid=gwid,
                plotParameters=plotParameters,
                ps1Pointings=ps1Pointings,
                pathToProbMap=pathToProbMap,
                mjdStart=mjdStart,
                raArray=raArray,
                decArray=decArray,
                fileFormats=["png", "pdf"],
                folderName="survey_matched_source_maps",
                plotTitle=plotTitle
            )

            for r in [False, 0.15]:

                raArray, decArray = self._get_matched_sources(
                    gwid,
                    plotParameters,
                    redshiftLimit=r)

                plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Matched Catalogue Sources" % locals(
                )

                if r is not False:
                    plotTitle += " Within z < %(r)s" % locals(
                    )

                self._generate_probability_plot(
                    gwid=gwid,
                    plotParameters=plotParameters,
                    ps1Pointings=ps1Pointings,
                    pathToProbMap=pathToProbMap,
                    mjdStart=mjdStart,
                    raArray=raArray,
                    decArray=decArray,
                    fileFormats=["png", "pdf"],
                    folderName="survey_matched_source_maps",
                    plotTitle=plotTitle
                )

            # WITH 2MASS AND FAKERS
            raArray, decArray = self._get_matched_sources(
                gwid,
                plotParameters,
                redshiftLimit=r,
                match2mass=True
            )

            plotTitle = "%(gwid)s Probability Map, PS1 Footprints & Matched Sources\n with 2MASS Counterparts, Axis Measurements and within z < 0.15" % locals(
            )

            self._generate_probability_plot(
                gwid=gwid,
                plotParameters=plotParameters,
                ps1Pointings=ps1Pointings,
                pathToProbMap=pathToProbMap,
                mjdStart=mjdStart,
                raArray=raArray,
                decArray=decArray,
                fileFormats=["png", "pdf"],
                folderName="survey_matched_source_maps",
                plotTitle=plotTitle
            )

        self.log.info('completed the ``get_source_plots`` method')
        return None

    # use the tab-trigger below for new method
    def _get_matched_sources(
            self,
            gwid,
            plotParameters,
            redshiftLimit=False,
            allNed=False,
            match2mass=False):
        """
        *get matched sources*

        **Key Arguments: **
            - ``gwid`` -- gravitational wave ID
            - ``plotParameters`` -- plot parameters from settings
            - ``redshiftLimit`` -- limit in redshift for returned sources
            - ``allNed`` -- no limits on query
            - ``match2mass`` -- NED sources need to be 2MASS sources with semi-major axis measurement

        **Return: **
            - ``ra`` -- array of match NED source RAs
            - ``dec`` -- array of match NED source DECs
        """
        self.log.info('starting the ``_get_matched_sources`` method')

        # return self._sampled_area_only_points()

        if allNed == True:
            # UNPACK THE PLOT PARAMETERS
            centralCoordinate = plotParameters["centralCoordinate"]
            raRange = plotParameters["raRange"]
            decRange = plotParameters["decRange"]

            raMax = centralCoordinate[0] + raRange / 2.
            raMin = centralCoordinate[0] - raRange / 2.
            decMax = centralCoordinate[1] + decRange / 2.
            decMin = centralCoordinate[1] - decRange / 2.

            sqlQuery = u"""
                select raDeg, decDeg from tcs_cat_ned_stream where raDeg > %(raMin)s and raDeg < %(raMax)s and decDeg > %(decMin)s and decDeg < %(decMax)s
            """ % locals()
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.cataloguesDbConn
            )

        else:
            if redshiftLimit:
                redshiftClause = " and t.z is not null and t.z < %(redshiftLimit)s and (t.z_quality is null or t.z_quality not like 'PHOT%%') and (t.catalogue_object_subtype is null or t.catalogue_object_subtype not like '%%*%%')" % locals(
                )
            else:
                redshiftClause = ""

            if match2mass:
                match2massClause = " and t.2mass_id is not null and t.major_axis_arcsec is not null"
            else:
                match2massClause = ""

            tcs_cross_matches = "tcs_%(gwid)s_catalogued_sources" % locals()

            sqlQuery = u"""
                select t.catalogue_object_ra as raDeg, t.catalogue_object_dec as decDeg from ps1_pointings p, %(tcs_cross_matches)s t where p.ps1_exp_id=t.transient_object_id and gw_id = "%(gwid)s" %(redshiftClause)s %(match2massClause)s;
            """ % locals()
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn
            )

        print sqlQuery

        ra = []
        dec = []
        ra[:] = [row["raDeg"] for row in rows]
        dec[:] = [row["decDeg"] for row in rows]

        ra = np.array(ra)
        dec = np.array(dec)

        self.log.info('completed the ``_get_matched_sources`` method')
        return ra, dec

    def _sampled_area_only_points(
            self):
        """
        *sampled area only points*
        """
        self.log.info('starting the ``_sampled_area_only_points`` method')

        coords1 = [
            (140.0, 6.0),
            (153.0, -12.0)
        ]
        coords2 = [
            (149.0, 2.0),
            (152.0, -7.0)
        ]
        coords3 = [
            (132.0, 4.0),
            (154.0, -16.5)
        ]
        final = [
            (140.0, 6.0),
            (154.0, -16.5)
        ]
        coords = final

        # CREATE AN ARRAY OF RELEVANT HTMIDS AND FIND MAX AND MIN

        mesh16 = htm.HTM(16)
        theseArrays = []
        radius = 2.5
        ra = []
        dec = []
        for co in coords:
            ra1 = co[0]
            dec1 = co[1]
            thisArray = mesh16.intersect(
                ra1, dec1, radius, inclusive=True)
            hmax = thisArray.max()
            hmin = thisArray.min()

            ratio = float(hmax - hmin + 1) / float(thisArray.size)
            if ratio < 100 or thisArray.size > 2000:
                htmWhereClause = "where htm16ID between %(hmin)s and %(hmax)s" % locals(
                )
            else:
                s = StringIO()
                np.savetxt(s, thisArray, fmt='%d', newline=",")
                thesHtmIds = s.getvalue()[:-1]
                htmWhereClause = "where htm16ID in (%(thesHtmIds)s)" % locals()

            # FINALLY BUILD THE FULL QUERY
            sqlQuery = """select raDeg, decDeg, redshift, object_type from tcs_cat_ned_stream %(htmWhereClause)s and redshift is not null and redshift < 0.15 and (redshift_quality is null or redshift_quality not like 'PHOT%%') and (object_type is null or object_type not like "%%*%%") """ % locals(
            )
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.cataloguesDbConn
            )

            raList = []
            decList = []
            for row in rows:
                raList.append(row["raDeg"])
                decList.append(row["decDeg"])

            tRa = np.array([ra1])
            tDec = np.array([dec1])
            raList = np.array(raList)
            decList = np.array(decList)
            indexList1, indexList2, separation = mesh16.match(
                tRa, tDec, raList, decList, radius, maxmatch=0)
            redshiftList = []
            for i in range(indexList1.size):
                ra.append(rows[indexList2[i]]["raDeg"])
                dec.append(rows[indexList2[i]]["decDeg"])

        ra = np.array(ra)
        dec = np.array(dec)

        self.log.info('completed the ``_sampled_area_only_points`` method')
        return ra, dec

    # use the tab-trigger below for new method
    # xt-class-method
