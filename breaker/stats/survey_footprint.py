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

        # INITIAL ACTIONS - CONNECT TO THE DATABASE REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn, self.atlasDbConn = db.get()

        return None

    def get(self):
        """
        *get the survey footprint stats and print to screen/file*

        **Return:**
            - ``None``
        """
        self.log.info('starting the ``get`` method')

        # GRAB METADATA FROM THE DATABASES
        this = plot_wave_observational_timelines(
            log=self.log, settings=self.settings)
        plotParameters, ps1Transients, ps1Pointings, altasPointings, atlasTransients = this.get_gw_parameters_from_settings(
            gwid=self.gwid)

        if self.telescope == "atlas":
            ps1Pointings = []
        if self.telescope == "ps1":
            altasPointings = []

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
        pixelSizeDeg = 0.066667
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
        aMap, mapHeader = hp.read_map(pathToProbMap, h=True)
        nside = hp.pixelfunc.get_nside(aMap)
        hpixArea = hp.nside2pixarea(nside, "degrees")

        # import matplotlib.pyplot as plt
        # hp.mollview(aMap, title="mollview image RING", cmap="YlOrRd")
        # hp.graticule()
        # plt.show()

        # HEALPY REQUIRES RA, DEC IN RADIANS AND AS TWO SEPERATE ARRAYS
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        # THETA: IS THE POLAR ANGLE, RANGING FROM 0 AT THE NORTH POLE TO PI AT THE SOUTH POLE.
        # PHI: THE AZIMUTHAL ANGLE ON THE SPHERE FROM 0 TO 2PI
        # CONVERT DEC TO THE REQUIRED HEALPIX FORMAT
        nd = -wd + 90.

        # CONVERT WORLD TO HEALPIX INDICES (NON-UNIQUE IDS!)
        healpixIds = hp.ang2pix(nside, theta=nd * DEG_TO_RAD_FACTOR,
                                phi=wr * DEG_TO_RAD_FACTOR)

        # NOW READ THE VALUES OF THE MAP AT THESE HEALPIX INDICES
        probs = []
        probs[:] = [aMap[i] for i in healpixIds]
        probs = np.array(probs)

        # CREATE THE FITS FILE
        hdu = fits.PrimaryHDU(header=header, data=probs)

        # SORT ALL POINTINGS VIA MJD

        ps1Pointings = sorted(list(ps1Pointings),
                              key=itemgetter('mjd'), reverse=False)
        altasPointings = sorted(list(altasPointings),
                                key=itemgetter('mjd'), reverse=False)

        # ONLY KEEP NON-ZERO PROB FOOTPRINTS
        tmpPointings = []
        ps1exposureRadius = 1.4
        moveBy = (ps1exposureRadius / 2)**0.5
        atlasPointingSide = 5.46
        fivePoints = [(0, 0), (1, 1), (-1, -1), (-1, 1), (1, -1)]
        for pt in ps1Pointings:
            pra = pt["raDeg"]
            pdec = pt["decDeg"]
            # REMOVE LOWER PROBABILITY FOOTPRINTS
            phi = pra
            if phi > 180.:
                phi = phi - 360.
            theta = -pdec + 90.

            # CONVERT WORLD TO HEALPIX INDICES (NON-UNIQUE IDS!)
            for m in fivePoints:
                healpixId = hp.ang2pix(nside, theta=(theta + m[0] * moveBy) * DEG_TO_RAD_FACTOR,
                                       phi=(phi + m[1] * moveBy) * DEG_TO_RAD_FACTOR)
                thisProb = aMap[healpixId]
                thisProb = float("%0.*f" % (7, thisProb))
                if thisProb != 0.:
                    tmpPointings.append(pt)
                    break
        ps1Pointings = tmpPointings

        tmpPointings = []
        moveBy = (atlasPointingSide / 2)**0.5
        for pt in altasPointings:
            pra = pt["raDeg"]
            pdec = pt["decDeg"]
            # REMOVE LOWER PROBABILITY FOOTPRINTS
            phi = pra
            if phi > 180.:
                phi = phi - 360.
            theta = -pdec + 90.

            # CONVERT WORLD TO HEALPIX INDICES (NON-UNIQUE IDS!)
            for m in fivePoints:
                healpixId = hp.ang2pix(nside, theta=(theta + m[0] * moveBy) * DEG_TO_RAD_FACTOR,
                                       phi=(phi + m[1] * moveBy) * DEG_TO_RAD_FACTOR)
                thisProb = aMap[healpixId]
                thisProb = float("%0.*f" % (7, thisProb))
                if thisProb != 0.:
                    tmpPointings.append(pt)
                    break
        altasPointings = tmpPointings

        healpixDictionary = {}
        cumCoveredArray = np.ones(healpixIds.size)
        finalMJDs = []
        finalAreas = []
        finalProbs = []

        # USE THIS COORDINATE SET MATCHER TO MATCH OTHER COORDINATE SETS
        # AGAINST
        coordinateSet = Matcher(
            log=self.log,
            ra=wr,
            dec=wd,
            depth=12
        )

        ps1RAs = []
        ps1RAs[:] = [pt["raDeg"] for pt in ps1Pointings]
        ps1DECs = []
        ps1DECs[:] = [pt["decDeg"] for pt in ps1Pointings]
        ps1mjds = []
        ps1mjds[:] = [pt["mjd"] for pt in ps1Pointings]

        # MATCH THE PS1 POINTINGS AGAINST EVERY PIXEL IN THE SKY
        matchIndices1, matchIndices2, seps = coordinateSet.match(
            ra=np.array(ps1RAs),
            dec=np.array(ps1DECs),
            radius=ps1exposureRadius,
            maxmatch=0  # 1 = closest, 0 = all
        )
        # WHAT ARE ALL OF THE UNIQUE HEALPIXELS COVERED BY THE POINTINGS
        healpixIDIndices = np.unique(matchIndices2)

        # GENERATE DICTIONARY OF PIXEL/PROBABILITIES
        for i in healpixIDIndices:
            if i not in healpixDictionary:
                healpixDictionary[healpixIds[i]] = aMap[healpixIds[i]]

        atlasRAs = []
        atlasRAs[:] = [pt["raDeg"] for pt in altasPointings]
        atlasDECs = []
        atlasDECs[:] = [pt["decDeg"] for pt in altasPointings]
        atlasmjds = []
        atlasmjds[:] = [pt["mjd"] for pt in altasPointings]

        allHealpixIds = np.array([])
        for pti, pt in enumerate(altasPointings):
            pra = pt["raDeg"]
            pdec = pt["decDeg"]
            pmjd = pt["mjd"]

            pointingHealpixIds = healpixIds[(((pdec - wd)**2)**0.5 < atlasPointingSide / 2) & (
                (((pra - wr) * np.cos(wd * DEG_TO_RAD_FACTOR))**2)**0.5 < atlasPointingSide / 2)]

            # pointingHealpixIds = []
            # pointingHealpixIds[:] = [h for r, d, h in zip(wr, wd, healpixIds) if (
            #     ((pdec - d)**2)**0.5 < atlasPointingSide / 2) and ((((pra - r) * math.cos(d * DEG_TO_RAD_FACTOR))**2)**0.5 < atlasPointingSide / 2)]
            allHealpixIds = np.append(allHealpixIds, pointingHealpixIds)

        allHealpixIds = np.unique(allHealpixIds)

        # GENERATE DICTIONARY OF PIXEL/PROBABILITIES
        for hid in allHealpixIds:
            if int(hid) not in healpixDictionary:
                healpixDictionary[int(hid)] = aMap[int(hid)]

        # NOW CALCULATE THE TOTAL AREA AND PROBABILITY COVERAGE
        cumArea = len(healpixDictionary) * hpixArea
        cumProb = sum(healpixDictionary.values())

        print "AREA: %(cumArea)0.2f. PROB: %(cumProb)0.5f" % locals()

        self.log.info('completed the ``get`` method')
        return None

    # xt-class-method
