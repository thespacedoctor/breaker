#!/usr/local/bin/python
# encoding: utf-8
"""
*Plot the longitude coverage of a given wave over time for ATLAS and PS1*

:Author:
    David Young

:Date Created:
    March  1, 2017
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime


class longitude_coverage():
    """
    *The worker class for the longitude_coverage module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the gravitational wave id to plot the longitude coverage plot for
        - ``outputDirectory`` -- path to output the plots to (if false will read path from settings). Default *False*   
        - ``timestamp`` -- add a timestamp to the plot to show when it was created. Default *True* 

    **Usage:**

        To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

        To initiate a longitude_coverage object, use the following:

        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - add a tutorial about ``longitude_coverage`` to documentation
            - create a blog post about what ``longitude_coverage`` does

        .. code-block:: python

            usage code
    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            gwid,
            outputDirectory=False,
            settings=False,
            timestamp=True
    ):
        self.log = log
        log.debug("instansiating a new 'longitude_coverage' object")
        self.settings = settings
        self.gwid = gwid
        self.outputDirectory = outputDirectory
        self.timestamp = timestamp
        # xt-self-arg-tmpx

        from breaker.plots import plot_wave_observational_timelines
        from breaker import database

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # INITIAL ACTIONS
        # SETUP BREAKER DATABASE CONNECTIONS
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn, self.atlasDbConn, self.ps13piDbConn = db.get()

        plotter = plot_wave_observational_timelines(
            log=self.log,
            settings=self.settings
        )
        plotParameters, ps1Transients, self.ps1Pointings, self.atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
            gwid=self.gwid,
            inPastDays=False,  # int
            inFirstDays=(0, 7)  # (start, finish) tuple
        )

        # TIME OF WAVE DETECTION
        self.t0 = self.settings["gravitational waves"][
            gwid]["time"]["mjdStart"]

        return None

    def plot(self):
        """
        *get the longitude_coverage object*

        **Return:**
            - ``longitude_coverage``

        **Usage:**
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - update the package tutorial if needed

        .. code-block:: python 

            usage code 
        """
        self.log.info('starting the ``get`` method')

        gwid = self.gwid.upper()

        atlasTdelta = []
        atlasTdelta[:] = [(p["mjd"] - self.t0) *
                          24. for p in self.atlasPointings]
        atlasRa = []
        atlasRa[:] = [p["raDeg"] for p in self.atlasPointings]

        ps1Tdelta = []
        ps1Tdelta[:] = [(p["mjd"] - self.t0) * 24. for p in self.ps1Pointings]
        ps1Ra = []
        ps1Ra[:] = [p["raDeg"] for p in self.ps1Pointings]

        maxTD = max(max(ps1Tdelta), max(atlasTdelta)) * 1.2

        plt.scatter(
            x=np.array(ps1Tdelta),  # numpy array of x-points
            y=np.array(ps1Ra),  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=2,
            c="#859900",    # color or sequence of color, optional, default
            marker='o',
            cmap=None,
            norm=None,
            vmin=None,
            vmax=None,
            alpha=None,
            linewidths=None)

        plt.scatter(
            x=np.array(atlasTdelta),  # numpy array of x-points
            y=np.array(atlasRa),  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=2,
            c="#6c71c4",    # color or sequence of color, optional, default
            marker='o',
            cmap=None,
            norm=None,
            vmin=None,
            vmax=None,
            alpha=None,
            linewidths=None)

        ax = plt.gca()
        ax.set_ylabel("Right Ascension [deg.]")
        ax.set_xlabel(
            "Hrs. Since %(gwid)s  Detection [$log(t-t_0)$]" % locals())
        ax.set_xscale("log", nonposy='clip')
        ax.set_xlim(0.001, maxTD)

        if self.timestamp:
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%d %H:%M.%S UTC")
            ax.text(0, 1.02, utcnow,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    color="#657b83",
                    fontsize=6)

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if self.settings and not self.outputDirectory:
            plotDir = self.settings["output directory"] + "/" + gwid
        elif self.outputDirectory:
            plotDir = self.outputDirectory
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plotTitle = "%(gwid)s Sky Longitude Coverage Rate" % locals()
        figureName = """%(plotTitle)s""" % locals(
        )

        if plotDir != ".":
            for f in ["png", "pdf"]:
                if not os.path.exists("%(plotDir)s/sky_coverage_rate/%(f)s" % locals()):
                    os.makedirs(
                        "%(plotDir)s/sky_coverage_rate/%(f)s" % locals())
                figurePath = "%(plotDir)s/sky_coverage_rate/%(f)s/%(figureName)s.%(f)s" % locals()
                plt.savefig(figurePath, bbox_inches='tight', dpi=300)
        else:
            for f in ["png", "pdf"]:
                figurePath = "%(plotDir)s/%(figureName)s.%(f)s" % locals()
                plt.savefig(figurePath, bbox_inches='tight', dpi=300)

        self.log.info('completed the ``get`` method')
        return longitude_coverage

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
