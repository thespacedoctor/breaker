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
    breaker update [-n] [-s <pathToSettingsFile>]
    breaker skymap <gwid> <pathToLVMap>
    breaker plot (timeline|history|sources) [-w <gwid>] [-s <pathToSettingsFile>]
    breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]
    breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
    breaker stats <gwid> [-s <pathToSettingsFile>]
    breaker listen <far> (<mjdStart> <mjdEnd> | <inLastNMins>) [-s <pathToSettingsFile>]
    breaker listen -d <far> [<sec>] [-s <pathToSettingsFile>]

    COMMANDS
    --------
    init                  setup the breaker settings file for the first time
    update                update the PS1 footprint table in breaker database and associate with GW-IDs. Optionally download overlapping NED source and also add to the database
    skymap                generate an all sky FITS & PDF image map given the path to the LV likeihood map (Meractor and Mollweide projections respectively)
    plot                  enter plotting mode
    timeline              plot from the epoch of the wave detection forward in time
    history               plot from now back in time over the last days, weeks and months
    comparison            produce a multi-panel plot to compare wave maps
    stats                 generate some coverage stats for a given wave survey campaign
    sources               overplot map with NED sources found within the wave campaign footprint
    faker                 generate a catalogue of simulated transient sources in PS1 exposure ID footprint
    listen                connect to grace DB and download maps found within the given time range

    ARGUMENTS
    ---------
    far                   false alarm rate limit in Hz (1e-7 Hz ~= 3.2 per year)
    gwid                  the gravitational wave ID
    pathToSettingsFile    path to the yaml settings file
    pathToMapDirectory    path to a directory containing localisation maps
    ps1ExpId              a panstarrs exposure ID
    mjdStart              start of an MJD range
    mjdEnd                end of the MJD range
    inLastNMins           in the last N number of minutes
    pathToLVMap           path to the LV likelihood map
    sec                   time in seconds

    FLAGS
    -----
    -h, --help            show this help message
    -s, --settings        the settings file
    -n, --updateNed       update the NED database steam
    -d, --daemon          listen in daemon mode
    -w, --waveId          a gravitational wave ID

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
from breaker.update_ps1_footprint_table import update_ps1_footprint_table
from breaker.plots.plot_wave_observational_timelines import plot_wave_observational_timelines
from breaker.plots.plot_wave_matched_source_maps import plot_wave_matched_source_maps
from breaker.fakers.generate_faker_catalogue import generate_faker_catalogue
from breaker.stats.survey_footprint import survey_footprint
from breaker.plots.plot_multi_panel_alternate_map_comparison import plot_multi_panel_alternate_map_comparison
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
        logLevel="DEBUG",
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

    # CALL FUNCTIONS/OBJECTS
    if update:
        u = update_ps1_footprint_table(
            log=log,
            settings=settings,
            updateNed=updateNedFlag
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
        p = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            gwid=gwid,
            plotType="timeline"
        )
        p.get()
    if plot and sources:
        p = plot_wave_matched_source_maps(
            log=log,
            settings=settings,
            gwid=gwid
        )
        p.get()
    if plot and comparison:
        p = plot_multi_panel_alternate_map_comparison(
            log=log,
            settings=settings,
            gwid=gwid,
            pathToMapDirectory=pathToMapDirectory
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
            gwid=gwid
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
            label="EM_READY",
            farThreshold=far,
            startMJD=float(startMJD),
            endMJD=float(timeNowMjd) + 1.
        )
        this.get_maps()
    if listen and mjdStart:
        this = mlisten(
            log=log,
            settings=settings,
            label="EM_READY",
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
            label="EM_READY",
            farThreshold=far,
            daemon=daemon
        )
        this.get_maps()

    if skymap:
        plotter = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            databaseConnRequired=False
        )

        plotter.generate_probability_plot(
            gwid=gwid,
            pathToProbMap=pathToLVMap,
            fileFormats=["pdf"],
            outputDirectory=".",
            projection="mollweide",
            plotType="timeline"
        )

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
