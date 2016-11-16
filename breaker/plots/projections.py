#!/usr/local/bin/python
# encoding: utf-8
"""
*Generate various sky-projection maps given a Healpix likeihood map*

:Author:
    David Young

:Date Created:
    June 30, 2016
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from plot_wave_observational_timelines import plot_wave_observational_timelines


class projections():
    """
    *The worker class for the projections module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the map's gravitational wave ID 
        - ``healpixPath`` -- path to the LV Healpix likelihood map.
        - ``projection`` -- the desired projection to warp the map to. Default *mollweide* (nothing else added yet)
        - ``outputDirectory`` -- path to the directory to output the map to. 

    **Usage:**
        .. todo::

            - add usage info
            - create a sublime snippet for usage

        .. code-block:: python 

            usage code   

    .. todo::

        - @review: when complete, clean projections class
        - @review: when complete add logging
        - @review: when complete, decide whether to abstract class to another module
    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            gwid,
            healpixPath,
            outputDirectory,
            settings=False,
            projection="mollweide",
    ):
        self.log = log
        log.debug("instansiating a new 'projections' object")
        self.settings = settings
        self.healpixPath = healpixPath
        self.projection = projection
        self.outputDirectory = outputDirectory
        self.gwid = gwid

        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions
        plotter = plot_wave_observational_timelines(
            log=log,
            settings=settings
        )

        plotter.generate_probability_plot(
            gwid=gwid,
            pathToProbMap=healpixPath,
            fileFormats=["pdf"],
            outputDirectory=outputDirectory,
            projection=projection,
            plotType="timeline"
        )

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the projections object*

        **Return:**
            - ``projections``

        .. todo::

            - @review: when complete, clean get method
            - @review: when complete add logging
        """
        self.log.info('starting the ``get`` method')

        projections = None

        self.log.info('completed the ``get`` method')
        return projections

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
