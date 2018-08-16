#!/usr/local/bin/python
# encoding: utf-8
"""
*Generate stats for a gravitational wave survey campaign footprint*
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
from operator import itemgetter
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.pyplot import savefig
import matplotlib.patches as patches
import numpy
from astropy import wcs as awcs
from astropy.io import fits
from HMpTy import Matcher
from fundamentals import tools, times
from astrocalc.coords import separations
from breaker.plots import plot_wave_observational_timelines
from fundamentals.renderer import list_of_dictionaries
import collections


class survey_footprint():
    """
    *Generate stats for a gravitational wave survey campaign footprint*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the unique wave id
        - ``telescope`` -- select a single telescope. The default is to select all telescopes. Default *False* [False|ps1|atlas]

    **Usage:**

        To generate cummulative stats for the wave "G184098"

        .. code-block:: python

            from breaker.stats import survey_footprint
            stats = survey_footprint(
                log=log,
                settings=settings,
                gwid="G184098",
                telescope=False
            )
            stats.get()
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            gwid=False,
            telescope=False
    ):
        self.log = log
        log.debug("instansiating a new 'survey_footprint' object")
        self.settings = settings
        self.gwid = gwid
        self.telescope = telescope

        pi = (4 * math.atan(1.0))
        self.DEG_TO_RAD_FACTOR = pi / 180.0
        self.RAD_TO_DEG_FACTOR = 180.0 / pi

        # INITIAL ACTIONS - CONNECT TO THE DATABASE REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn, self.atlasDbConn, self.ps13piDbConn = db.get()

        return None

    def get(self):
        """
        *get the survey footprint stats and print to screen/file*

        **Return:**
            - ``None``
        """
        self.log.debug('starting the ``get`` method')

        # GRAB METADATA FROM THE DATABASES
        this = plot_wave_observational_timelines(
            log=self.log, settings=self.settings)
        plotParameters, ps1Transients, ps1Pointings, altasPointings, atlasTransients = this.get_gw_parameters_from_settings(
            gwid=self.gwid,
            stackOnly=False)

        if self.telescope == "atlas":
            pointings = altasPointings
            pointingSide = 5.46
        if self.telescope == "ps1":
            pointings = ps1Pointings
            pointingSide = 0.4
        telescope = self.telescope.upper()

        # SORT ALL POINTINGS VIA MJD
        pointings = sorted(list(pointings),
                           key=itemgetter('mjd'), reverse=False)

        nside, hpixArea, aMap, healpixIds, wr, wd = self._create_healpixid_coordinate_grid()

        print "EXPID, RA, DEC, MJD, EXPTIME, FILTER, LIM-MAG, EXP-AREA, EXP-LIKELIHOOD, CUM-AREA, CUM-LIKELIHOOD" % locals()

        allHealpixIds = np.array([])
        dictList = []
        iindex = 0
        count = len(pointings)
        cumArea = 0
        cumProb = 0
        for pti, pt in enumerate(pointings):
            pti = pti + 1

            if pti > 1:
                # Cursor up one line and clear line
                sys.stdout.write("\x1b[1A\x1b[2K")

            percent = (float(pti) / float(count)) * 100.
            print '%(pti)s/%(count)s (%(percent)1.1f%% done): summing total area and likelihood covered by %(telescope)s' % locals()

            thisDict = collections.OrderedDict(sorted({}.items()))

            pra = pt["raDeg"]
            pdec = pt["decDeg"]
            pmjd = pt["mjd"]
            pexpid = pt["exp_id"]
            pexptime = pt["exp_time"]
            pfilter = pt["filter"]
            plim = pt["limiting_magnitude"]

            # DETERMINE THE CORNERS FOR EACH ATLAS EXPOSURE AS MAPPED TO THE
            # SKY
            decCorners = (pdec - pointingSide / 2,
                          pdec + pointingSide / 2)
            corners = []
            for d in decCorners:
                if d > 90.:
                    d = 180. - d
                elif d < -90.:
                    d = -180 - d
                raCorners = (pra - (pointingSide / 2) / np.cos(d * self.DEG_TO_RAD_FACTOR),
                             pra + (pointingSide / 2) / np.cos(d * self.DEG_TO_RAD_FACTOR))
                for r in raCorners:
                    if r > 360.:
                        r = 720. - r
                    elif r < 0.:
                        r = 360. + r
                    corners.append(hp.ang2vec(r, d, lonlat=True))

            # FLIP CORNERS 3 & 4 SO HEALPY UNDERSTANDS POLYGON SHAPE
            corners = [corners[0], corners[1],
                       corners[3], corners[2]]

            # RETURN HEALPIXELS IN EXPOSURE AREA
            expPixels = hp.query_polygon(nside, np.array(
                corners))

            expProb = []
            expProb[:] = [aMap[i] for i in expPixels]
            expProb = sum(expProb)
            expArea = len(expPixels) * hpixArea
            if expProb / expArea < 2e-6:
                continue

            pindex = "%(iindex)05d" % locals()
            iindex += 1

            allHealpixIds = np.append(allHealpixIds, expPixels)
            allHealpixIds = np.unique(allHealpixIds)
            cumProb = []
            cumProb[:] = [aMap[int(i)] for i in allHealpixIds]
            cumProb = sum(cumProb)
            cumArea = len(allHealpixIds) * hpixArea
            thisDict["INDEX"] = pindex
            thisDict["EXPID"] = pexpid
            thisDict["RA"] = "%(pra)5.5f" % locals()
            thisDict["DEC"] = "%(pdec)5.5f" % locals()
            thisDict["MJD"] = "%(pmjd)6.6f" % locals()
            thisDict["EXPTIME"] = "%(pexptime)02.1f" % locals()
            thisDict["FILTER"] = pfilter
            try:
                thisDict["LIM-MAG"] = "%(plim)5.2f" % locals()
            except:
                thisDict["LIM-MAG"] = "NaN"
            # thisDict["EXP-AREA"] = expArea
            # thisDict["EXP-LIKELIHOOD"] = expProb
            thisDict["CUM-AREA"] = "%(cumArea)05.2f" % locals()
            thisDict["CUM-LIKELIHOOD"] = "%(cumProb)05.2f" % locals()
            dictList.append(thisDict)

        if not len(dictList):
            thisDict = {}
            thisDict["INDEX"] = "NULL"
            thisDict["EXPID"] = "NULL"
            thisDict["RA"] = "NULL"
            thisDict["DEC"] = "NULL"
            thisDict["MJD"] = "NULL"
            thisDict["EXPTIME"] = "NULL"
            thisDict["FILTER"] = "NULL"
            thisDict["LIM-MAG"] = "NULL"
            dictList.append(thisDict)

        print "AREA: %(cumArea)0.2f. PROB: %(cumProb)0.5f" % locals()

        printFile = self.settings["output directory"] + "/" + \
            self.gwid + "/" + self.gwid + "-" + self.telescope + "-coverage-stats.csv"

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.settings["output directory"] + "/" + self.gwid):
            os.makedirs(self.settings["output directory"] + "/" + self.gwid)

        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=dictList,
        )
        csvData = dataSet.csv(filepath=printFile)

        print "The coverage stats file was written to `%(printFile)s`" % locals()

        self.log.debug('completed the ``get`` method')
        return None

    def _create_healpixid_coordinate_grid(
            self):
        """*create healpixid coordinate grid*
        """
        self.log.debug('starting the ``_create_healpix_id_list`` method')

        # GET THE PROBABILITY MAP FOR THE GIVEN GWID
        pathToProbMap = self.settings[
            "gravitational waves"][self.gwid]["mapPath"]
        mapName = pathToProbMap.split("/")[-1].replace(".fits", "")
        if not os.path.exists(pathToProbMap):
            message = "the path to the map %s does not exist on this machine" % (
                pathToProbMap,)
            self.log.critical(message)
            raise IOError(message)

        # UNPACK THE PLOT PARAMETERS
        pixelSizeDeg = 0.05
        raRange = 360.
        decRange = 180.

        # DETERMINE THE PIXEL GRID X,Y RANGES
        xRange = int(raRange / pixelSizeDeg)
        yRange = int(decRange / pixelSizeDeg) * 2

        # CREATE A NEW WCS OBJECT
        w = awcs.WCS(naxis=2)

        # SET THE REQUIRED PIXEL SIZE
        w.wcs.cdelt = numpy.array([pixelSizeDeg, pixelSizeDeg])

        # SET THE REFERENCE PIXEL TO THE CENTRE PIXEL
        w.wcs.crpix = [xRange / 2., yRange / 2.]

        # FOR AN ORTHOGONAL GRID THE CRPIX2 VALUE MUST BE ZERO AND CRPIX2
        # MUST REFLECT THIS
        w.wcs.crpix[1] -= w.wcs.crval[1] / w.wcs.cdelt[1]

        # WORLD COORDINATES AT REFERENCE PIXEL
        w.wcs.crval = [0, 0]
        # USE THE "GNOMONIC" PROJECTION ("COORDINATESYS---PROJECTION")
        w.wcs.ctype = ["RA---MER", "DEC--MER"]

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

        aMap, mapHeader = hp.read_map(pathToProbMap, h=True)
        nside = hp.pixelfunc.get_nside(aMap)
        hpixArea = hp.nside2pixarea(nside, degrees=True)

        # import matplotlib.pyplot as plt
        # hp.mollview(aMap, title="mollview image RING", cmap="YlOrRd")
        # hp.graticule()
        # plt.show()

        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        pi = (4 * math.atan(1.0))
        self.DEG_TO_RAD_FACTOR = pi / 180.0
        self.RAD_TO_DEG_FACTOR = 180.0 / pi

        # THETA: IS THE POLAR ANGLE, RANGING FROM 0 AT THE NORTH POLE TO PI AT THE SOUTH POLE.
        # PHI: THE AZIMUTHAL ANGLE ON THE SPHERE FROM 0 TO 2PI
        # CONVERT DEC TO THE REQUIRED HEALPIX FORMAT
        nd = -wd + 90.

        # CONVERT WORLD TO HEALPIX INDICES (NON-UNIQUE IDS!)
        healpixIds = hp.ang2pix(nside, theta=nd * self.DEG_TO_RAD_FACTOR,
                                phi=wr * self.DEG_TO_RAD_FACTOR)

        self.log.debug('completed the ``_create_healpix_id_list`` method')
        return nside, hpixArea, aMap, healpixIds, wr, wd

    def annotate_exposures(
        self,
        exposures,
        pointingSide
    ):
        """
        *generate a the likeihood coverage of an exposure*

        **Key Arguments:**
            - ``exposures`` -- a dictionary of exposures with the unique exposures IDs as keys and (ra, dec) as tuple value.
            - ``squareSize`` -- the size of the FOV of the exposure/skycell

        **Return:**
            - ``exposureIDs`` -- a list of the exposureIDs  as they appear in the original input exposure dictionary
            - ``pointingSide`` -- a list of total likeihood coverage of the exposures

        **Usage:**

            See class docstring
        """
        self.log.debug('starting the ``annotate`` method')

        nside, hpixArea, aMap, healpixIds, wr, wd = self._create_healpixid_coordinate_grid()

        exposureIDs = []
        ra = []
        dec = []

        exposureIDs = []
        exposureIDs[:] = [t for t in exposures.keys()]
        ra = []
        dec = []
        ra[:] = [r[0] for r in exposures.values()]
        dec[:] = [d[1] for d in exposures.values()]

        probs = []
        for e, pra, pdec in zip(exposureIDs, ra, dec):
            # DETERMINE THE CORNERS FOR EACH ATLAS EXPOSURE AS MAPPED TO THE
            # SKY
            decCorners = (pdec - pointingSide / 2,
                          pdec + pointingSide / 2)
            corners = []
            for d in decCorners:
                if d > 90.:
                    d = 180. - d
                elif d < -90.:
                    d = -180 - d
                raCorners = (pra - (pointingSide / 2) / np.cos(d * self.DEG_TO_RAD_FACTOR),
                             pra + (pointingSide / 2) / np.cos(d * self.DEG_TO_RAD_FACTOR))
                for r in raCorners:
                    if r > 360.:
                        r = 720. - r
                    elif r < 0.:
                        r = 360. + r
                    corners.append(hp.ang2vec(r, d, lonlat=True))

            # FLIP CORNERS 3 & 4 SO HEALPY UNDERSTANDS POLYGON SHAPE
            corners = [corners[0], corners[1],
                       corners[3], corners[2]]

            # RETURN HEALPIXELS IN EXPOSURE AREA
            expPixels = hp.query_polygon(nside, np.array(
                corners))

            expProb = []
            expProb[:] = [aMap[i] for i in expPixels]
            expProb = sum(expProb)
            probs.append(expProb)

        self.log.debug('completed the ``annotate`` method')
        return exposureIDs, probs

    # use the tab-trigger below for new method
    # xt-class-method
