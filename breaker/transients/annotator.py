#!/usr/local/bin/python
# encoding: utf-8
"""
*Annotate transients with supplimentary information resulting from the gravity events*

:Author:
    David Young

:Date Created:
    March 14, 2017
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import healpy as hp
import numpy as np
import math


class annotator():
    """
    *Annotate transients with supplimentary information resulting from the gravity events*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the ID of the gravity event to annotate the transients against (matched in settings file) 

    **Usage:**

        To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

        To initiate a annotator object, use the following:

        .. todo::

            - add a tutorial about ``annotator`` to documentation
            - create a blog post about what ``annotator`` does

        .. code-block:: python 

            from breaker.transients import annotator
            an = annotator(
                log=log,
                settings=settings,
                gwid="G211117"
            )
            transientNames, probs = an.annotate(transients) 
    """
    # INITIALISATION

    def __init__(
            self,
            log,
            gwid,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new 'annotator' object")
        self.settings = settings
        self.gwid = gwid

        return None

    def annotate(
        self,
        transients
    ):
        """
        *generate a list of the likihood contours the transients lie within*

        **Key Arguments:**
            - ``transients`` -- a dictionary of transients with the unique transient IDs as keys and (ra, dec) as tuple value.

        **Return:**
            - ``transientNames`` -- a list of the transient names as they appear in the original input transient dictionary
            - ``prosb`` -- a list of the nearest 10% likihood conntour the transient falls within (syncs with `transientNames` list)

        **Usage:**

            See class docstring
        """
        self.log.info('starting the ``annotate`` method')

        transientNames = []
        ra = []
        dec = []

        transientNames = []
        transientNames[:] = [t for t in transients.keys()]
        ra = []
        dec = []
        ra[:] = [r[0] for r in transients.values()]
        dec[:] = [d[1] for d in transients.values()]

        probs = self._generate_probability_distribution(ra, dec)

        self.log.info('completed the ``annotate`` method')
        return transientNames, probs

    def _generate_probability_distribution(
            self,
            ra,
            dec):
        """* generate probability distribution*

        **Key Arguments:**
            - ``ra`` -- a list of RA values for transients.
            - ``dec`` -- a list of DEC values for transients.

        **Return:**
            - ``prob`` -- a list of the probabilty contours the transients lie within (indices synced with RA and DEC lists)
        """
        self.log.info(
            'starting the ``_generate_probability_distribution`` method')

        pathToProbMap = self.settings[
            "gravitational waves"][self.gwid]["mapPath"]

        # READ HEALPIX MAPS FROM FITS FILE
        # THIS FILE IS A ONE COLUMN FITS BINARY, WITH EACH CELL CONTAINING AN
        # ARRAY OF PROBABILITIES
        # READ IN THE HEALPIX FITS FILE
        aMap, mapHeader = hp.read_map(pathToProbMap, 0, h=True, verbose=False)
        # DETERMINE THE SIZE OF THE HEALPIXELS
        nside = hp.npix2nside(len(aMap))
        totalProb = sum(aMap)

        # CONTOURS - NEED TO ADD THE CUMMULATIVE PROBABILITY
        # GET THE INDEXES ORDERED BY PROB OF PIXELS; HIGHEST TO LOWEST
        i = np.flipud(np.argsort(aMap))
        cumsum = np.cumsum(aMap[i])
        cls = np.empty_like(aMap)
        cls[i] = cumsum * 100

        healpixId = hp.pixelfunc.ang2pix(
            nside, ra, dec, nest=False, lonlat=True)

        vals = cls[healpixId]

        probs = []
        probs[:] = [roundup(p, 10) if p < 100. else 100 for p in vals]

        self.log.info(
            'completed the ``_generate_probability_distribution`` method')
        return probs

    # use the tab-trigger below for new method
    # xt-class-method


def roundup(val, resolution):
    return int(math.ceil(val / float(resolution))) * int(resolution)
