Command-Line Usage
==================

.. code-block:: bash 
   
    
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
    
    
