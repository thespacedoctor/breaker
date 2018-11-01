#!/usr/local/bin/python
# encoding: utf-8
"""
*The CL tools for breaker*

:Author:
    David Young

:Date Created:
    October 29, 2015

Usage:
    breaker init
    breaker update [-naP] [-s <pathToSettingsFile>]
    breaker skymap [-oe] <gwid> [<pathToLVMap>] [-c <centerDeg>]
    breaker plot [-a] (timeline|history|sources) [-w <gwid>] [-t <telescope>] [-p <projection>] [-f <filters>] [-s <pathToSettingsFile>]
    breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
    breaker stats <gwid> [<telescope>] [-s <pathToSettingsFile>]
    breaker listen <far> (<mjdStart> <mjdEnd> | <inLastNMins>) [-s <pathToSettingsFile>]
    breaker listen -d [<far> <sec>] [-s <pathToSettingsFile>]
    breaker contour <gwid> <ra> <dec> 

    COMMANDS
    --------
    init                  setup the breaker settings file for the first time
    update                update the PS1 footprint table in breaker database and associate with GW-IDs. Optionally download overlapping NED source and also add to the database
    skymap                generate an all sky FITS & PDF image map given the path to the LV likeihood map (Meractor and Mollweide projections respectively)
    plot                  enter plotting mode
    timeline              plot from the epoch of the wave detection forward in time
    history               plot from now back in time over the last days, weeks and months
    stats                 generate some coverage stats for a given wave survey campaign
    sources               overplot map with NED sources found within the wave campaign footprint
    faker                 generate a catalogue of simulated transient sources in PS1 exposure ID footprint
    listen                connect to grace DB and download maps found within the given time range
    contour               determine within which likelihood contour a given transient location lies (nearest 10%)

    ARGUMENTS
    ---------
    ra                    right ascendsion (sexegesimal or decimal degrees)
    dec                   declination (sexegesimal or decimal degrees)
    far                   false alarm rate limit in Hz. Default *1e-5* (~= 1 per day)
    -w <gwid>             the gravitational wave ID (graceDB or human-readable GW forms allowed)
    pathToSettingsFile    path to the yaml settings file
    -c <centerDeg>        the central longitude line (deg)
    pathToMapDirectory    path to a directory containing localisation maps
    ps1ExpId              a panstarrs exposure ID
    mjdStart              start of an MJD range
    mjdEnd                end of the MJD range
    inLastNMins           in the last N number of minutes
    pathToLVMap           path to the LV likelihood map
    sec                   time in seconds
    -t <telescope>        select an individual telescope (default is all telescopes) [ps1|atlas]
    -p <projection>       skymap projection. Default *mercator*. [mercator|gnomonic|mollweide]
    -f <filters>          which exposure filters to show on plots

    FLAGS
    -----
    -h, --help            show this help message
    -s, --settings        the settings file
    -n, --updateNed       update the NED database steam
    -d, --daemon          listen in daemon mode
    -a, --all             plot all timeline plot (including the CPU intensive -21-0 days and all transients/footprints plots)
    -P, --no-pointings    do not update pointings 
    -o, --default-output  output files to the default breaker output location (as set in settings file)
    -e, --exposures       overlay atlas and ps1 exposures

"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
import readline
import glob
import pickle
from docopt import docopt
from fundamentals import tools, times
from breaker.update_ps1_atlas_footprint_tables import update_ps1_atlas_footprint_tables
from breaker.plots.plot_wave_observational_timelines import plot_wave_observational_timelines
from breaker.plots.plot_wave_matched_source_maps import plot_wave_matched_source_maps
from breaker.fakers.generate_faker_catalogue import generate_faker_catalogue
from breaker.stats.survey_footprint import survey_footprint
from breaker.gracedb.listen import listen as mlisten
from astrocalc.times import now as mjdNow
from subprocess import Popen, PIPE, STDOUT
# from ..__init__ import *


def main(arguments=None):
    """
    *The main function used when ``cl_utils.py`` is run as a single script from the cl, or when installed as a cl command*
    """

    # setup the command-line util settings
    su = tools(
        arguments=arguments,
        docString=__doc__,
        logLevel="WARNING",
        options_first=False,
        projectName="breaker"
    )
    arguments, settings, log, dbConn = su.setup()

    # unpack remaining cl arguments using `exec` to setup the variable names
    # automatically
    for arg, val in arguments.iteritems():
        if arg[0] == "-":
            varname = arg.replace("-", "") + "Flag"
        else:
            varname = arg.replace("<", "").replace(">", "")
        if isinstance(val, str) or isinstance(val, unicode):
            exec(varname + " = '%s'" % (val,))
        else:
            exec(varname + " = %s" % (val,))
        if arg == "--dbConn":
            dbConn = val
        log.debug('%s = %s' % (varname, val,))

    ## START LOGGING ##
    startTime = times.get_now_sql_datetime()
    log.info(
        '--- STARTING TO RUN THE cl_utils.py AT %s' %
        (startTime,))

    if init:
        from os.path import expanduser
        home = expanduser("~")
        filepath = home + "/.config/breaker/breaker.yaml"
        try:
            cmd = """open %(filepath)s""" % locals()
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        except:
            pass
        try:
            cmd = """start %(filepath)s""" % locals()
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        except:
            pass

    if not far:
        far = 1e-5

    if gwid and gwid[:2] == "GW":
        for g in settings["gravitational waves"]:
            if settings["gravitational waves"][g]["human-name"] == gwid.strip():
                gwid = g
    if wFlag and wFlag[:2] == "GW":
        for g in settings["gravitational waves"]:
            if settings["gravitational waves"][g]["human-name"] == wFlag.strip():
                wFlag = g

    # CALL FUNCTIONS/OBJECTS
    if update:
        pointingsFlag = True
        if nopointingsFlag:
            pointingsFlag = False
        u = update_ps1_atlas_footprint_tables(
            log=log,
            settings=settings,
            updateNed=updateNedFlag,
            updateAll=allFlag,
            updatePointings=pointingsFlag
        )
        u.get()
    if plot and history:
        p = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            plotType="history"
        )
        p.get()
    if plot and timeline:
        if not pFlag:
            pFlag = "mercator"

        if fFlag:
            filters = list(fFlag)
        else:
            filters = False

        p = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            gwid=wFlag,
            plotType="timeline",
            allPlots=allFlag,
            telescope=tFlag,
            projection=pFlag,
            filters=filters,
            probabilityCut=True
        )
        p.get()
    if plot and sources:
        p = plot_wave_matched_source_maps(
            log=log,
            settings=settings,
            gwid=gwid
        )
        p.get()

    if faker:
        f = generate_faker_catalogue(
            log=log,
            settings=settings,
            ps1ExpId=ps1ExpId,
            gwid=False
        )
        f.get()
    if stats:
        s = survey_footprint(
            log=log,
            settings=settings,
            gwid=gwid,
            telescope=telescope
        )
        s.get()
    if listen and inLastNMins:
        timeNowMjd = mjdNow(
            log=log
        ).get_mjd()
        startMJD = float(timeNowMjd) - float(inLastNMins) / (60 * 60 * 24.)
        this = mlisten(
            log=log,
            settings=settings,
            #label="EM_READY | EM_Selected | ADVOK",
            label="",
            farThreshold=far,
            startMJD=float(startMJD),
            endMJD=float(timeNowMjd) + 1.
        )
        this.get_maps()
    if listen and mjdStart:
        this = mlisten(
            log=log,
            settings=settings,
            # label="EM_READY | EM_Selected | ADVOK",
            label="",
            farThreshold=far,
            startMJD=float(mjdStart),
            endMJD=float(mjdEnd)
        )
        this.get_maps()
    if listen and daemonFlag:
        if sec:
            daemon = float(sec)
        else:
            daemon = True
        this = mlisten(
            log=log,
            settings=settings,
            # label="EM_READY | EM_Selected | ADVOK",
            label="",
            farThreshold=far,
            daemon=daemon
        )
        this.get_maps()

    if skymap:
        if exposuresFlag:
            databaseConnRequired = True
        else:
            databaseConnRequired = False

        plotter = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            databaseConnRequired=databaseConnRequired
        )

        if exposuresFlag:
            plotParameters, ps1Transients, ps1Pointings, atlasPointings, atlasTransients = plotter.get_gw_parameters_from_settings(
                gwid=gwid,
                inFirstDays=(0, 31)
            )
        else:
            ps1Transients = []
            atlasTransients = []
            ps1Pointings = []
            atlasPointings = []

        ps1Transients = []
        atlasTransients = []

        if not cFlag:
            cFlag = 0.
        else:
            cFlag = float(cFlag)

        if defaultoutputFlag:
            outputDirectory = False
        else:
            outputDirectory = "."

        plotter.generate_probability_plot(
            gwid=gwid,
            ps1Transients=ps1Transients,
            atlasTransients=atlasTransients,
            ps1Pointings=ps1Pointings,
            atlasPointings=atlasPointings,
            pathToProbMap=pathToLVMap,
            fileFormats=["pdf", "png"],
            outputDirectory=outputDirectory,
            projection="mollweide",
            plotType="timeline",
            folderName="all_sky_plots",
            fitsImage=False,
            allSky=True,
            center=cFlag
        )

        plotter.generate_probability_plot(
            gwid=gwid,
            pathToProbMap=pathToLVMap,
            fileFormats=["pdf", "png"],
            outputDirectory=outputDirectory,
            projection="cartesian",
            plotType="timeline",
            folderName="all_sky_plots",
            fitsImage=True,
            allSky=True,
            center=cFlag
        )

    if contour:
        from breaker.transients import annotator
        an = annotator(
            log=log,
            settings=settings,
            gwid=gwid
        )

        from astrocalc.coords import unit_conversion
        # ASTROCALC UNIT CONVERTER OBJECT
        converter = unit_conversion(
            log=log
        )
        ra = converter.ra_sexegesimal_to_decimal(
            ra=ra
        )
        dec = converter.dec_sexegesimal_to_decimal(
            dec=dec
        )
        transients = {"cl": (ra, dec)}
        transientNames, probs = an.annotate(transients)
        percentage = probs[0]
        print "The transient lies within the inner %(percentage)s%% likelihood contour of event %(gwid)s" % locals()

    if "dbConn" in locals() and dbConn:
        dbConn.commit()
        dbConn.close()
    ## FINISH LOGGING ##
    endTime = times.get_now_sql_datetime()
    runningTime = times.calculate_time_difference(startTime, endTime)
    log.info('-- FINISHED ATTEMPT TO RUN THE cl_utils.py AT %s (RUNTIME: %s) --' %
             (endTime, runningTime, ))

    return


if __name__ == '__main__':
    main()
