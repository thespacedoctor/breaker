#!/usr/local/bin/python
# encoding: utf-8
"""
*The CL tools for breaker*

:Author:
    David Young

:Date Created:
    October 29, 2015

Usage:
    breaker update [-n] [-s <pathToSettingsFile>]
    breaker plot (timeline|history|sources [<gwid>]) [-s <pathToSettingsFile>]
    breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]
    breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
    breaker stats <gwid> [-s <pathToSettingsFile>]
    breaker listen <far> (<mjdStart> <mjdEnd> | <inLastNMins>) [-s <pathToSettingsFile>]
    breaker listen -d <far> [-s <pathToSettingsFile>]

    -h, --help            show this help message
    -s, --settings        the settings file
    -n, --updateNed       update the NED database steam
    -d, --daemon          listen in daemon mode
    far                   false alarm rate limit in Hz (1e-7 Hz ~= 3.2 per year)
    plot                  update the gravitational wave plots
    timeline              observations looking forward from date of GW detection
    history               observations from the past x days
    faker                 generate a catalogue of simulated transient source in ps1 FP
    stats                 output some stats for the GW surveys
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
# from ..__init__ import *


def tab_complete(text, state):
    return (glob.glob(text + '*') + [None])[state]


def main(arguments=None):
    """
    *The main function used when ``cl_utils.py`` is run as a single script from the cl, or when installed as a cl command*

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

    # setup the command-line util settings
    su = tools(
        arguments=arguments,
        docString=__doc__,
        logLevel="DEBUG",
        options_first=False,
        projectName="breaker"
    )
    arguments, settings, log, dbConn = su.setup()

    # tab completion for raw_input
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(tab_complete)

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
        this = mlisten(
            log=log,
            settings=settings,
            label="EM_READY",
            farThreshold=far,
            daemon=True
        )
        this.get_maps()

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
