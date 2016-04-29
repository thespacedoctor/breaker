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
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.pyplot import savefig
import matplotlib.patches as patches
from dryxPython import astrotools as dat
from dryxPython import mysql as dms
from fundamentals import tools, times


class survey_footprint():
    """
    *Generate stats for a gravitational wave survey campaign footprint*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the unique wave id

    **Usage:**

        .. todo::

            - add usage info
            - create a sublime snippet for usage

        .. code-block:: python 

            usage code 

    .. todo::

        - @review: when complete, clean _database class
        - @review: when complete add logging
        - @review: when complete, decide whether to abstract class to another module
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            gwid=False

    ):
        self.log = log
        log.debug("instansiating a new 'survey_footprint' object")
        self.settings = settings
        self.gwid = gwid

        # INITIAL ACTIONS - CONNECT TO THE DATABASE REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = db.get()

        return None

    # Method Attributes
    def get(self):
        """
        *get the survey_footprint stats and print to screen/file*

        **Return:**
            - ``survey_footprint``

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
        self.log.info('starting the ``get`` method')

        # GRAB METADATA FROM THE DATABASES
        from breaker.plots import plot_wave_observational_timelines
        this = plot_wave_observational_timelines(
            log=self.log, settings=self.settings)
        plotParameters, ps1Transients, ps1Pointings, altasPointings = this.get_gw_parameters_from_settings(
            gwid=self.gwid)

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
        hpixArea = hp.nside2pixarea(nside, "degrees")

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

        # CONVERT WORLD TO HEALPIX INDICES (NON-UNIQUE IDS!)
        healpixIds = hp.ang2pix(nside, theta=nd * DEG_TO_RAD_FACTOR,
                                phi=wr * DEG_TO_RAD_FACTOR)

        # NOW READ THE VALUES OF THE MAP AT THESE HEALPIX INDICES
        probs = []
        probs[:] = [aMap[i] for i in healpixIds]
        probs = np.array(probs)

        # CREATE THE FITS FILE
        hdu = fits.PrimaryHDU(header=header, data=probs)

        from operator import itemgetter
        ps1Pointings = list(ps1Pointings)
        ps1Pointings = sorted(
            ps1Pointings, key=itemgetter('mjd'), reverse=False)

        exposureRadius = 1.4
        healpixDictionary = {}
        totalPointings = len(ps1Pointings)
        cumCoveredArray = np.ones(healpixIds.size)
        finalMJDs = []
        finalAreas = []
        finalProbs = []

        for pti, pt in enumerate(ps1Pointings):

            pra = pt["raDeg"]
            pdec = pt["decDeg"]
            pmjd = pt["mjd"]
            # print pmjd

            decFactor = math.cos(pdec * DEG_TO_RAD_FACTOR)
            coveredPixels = []

            covered = -1
            for r, d, c in zip(wr, wd, cumCoveredArray):
                if c == 0:
                    # WE HAVE COVERED THIS PIXEL BEFORE -- MOVE ON
                    covered = 0
                elif (((pra - r) * decFactor)**2)**0.5 > exposureRadius:
                    covered = 0
                elif ((pdec - d)**2)**0.5 > exposureRadius:
                    covered = 0
                else:
                    angularSeparation, northSep, eastSep = dat.get_angular_separation(
                        log=self.log,
                        ra1=pra,
                        dec1=pdec,
                        ra2=r,
                        dec2=d
                    )
                    if angularSeparation < exposureRadius * 60 * 60:
                        covered = 1
                    else:
                        covered = 0
                coveredPixels.append(covered)

            # APPEND TO CUMMULATIVE COVERED ARRAY - SPEED THINGS UP!
            cumCoveredArray = cumCoveredArray * (1 - np.array(coveredPixels))

            # RESHAPE THE ARRAY AS BITMAP
            coveredPixelGrid = np.reshape(
                np.array(coveredPixels), (yRange, xRange))

            # GRAB THE WCS FROM HEADER GENERATED EARLIER
            from wcsaxes import datasets, WCS
            wcs = WCS(hdu.header)

            # PLOT MAP WITH PROJECTION IN HEADER
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
            im = ax.imshow(coveredPixelGrid,
                           cmap="YlOrRd", origin='lower', alpha=0.7, zorder=1)

            coveredHealPixels = coveredPixels * healpixIds
            coveredProbabilities = coveredPixels * probs

            for cx, cp in zip(coveredHealPixels, coveredProbabilities):
                if cx != 0:
                    healpixDictionary[cx] = cp

            cumArea = len(healpixDictionary) * hpixArea
            cumProb = sum(healpixDictionary.values())

            finalMJDs.append(pmjd)
            finalAreas.append(cumArea)
            finalProbs.append(cumProb)

            print "%(pti)s/%(totalPointings)s.  MJD: %(pmjd)s. AREA: %(cumArea)0.2f. PROB: %(cumProb)0.5f" % locals()

        gwid = self.gwid
        import codecs
        from datetime import datetime, date, time
        now = datetime.now()
        now = now.strftime("%Y%m%dt%H%M%S")
        pathToWriteFile = "/Users/Dave/Desktop/%(now)s_%(gwid)s_%(mapName)s.csv" % locals(
        )
        writeFile = codecs.open(pathToWriteFile, encoding='utf-8', mode='w')
        writeFile.write(
            "MJD,Cummulative Area (degrees),Cummulative Probability Covered\n" % locals())
        for m, a, p in zip(finalMJDs, finalAreas, finalProbs):
            writeFile.write("%(m)0.3f,%(a)0.2f,%(p)0.5f\n" % locals())
        writeFile.close()

        self.log.info('completed the ``get`` method')
        return survey_footprint

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx