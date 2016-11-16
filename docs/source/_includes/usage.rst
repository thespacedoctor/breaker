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
        breaker update [-n] [-s <pathToSettingsFile>]
        breaker skymap <gwid> <pathToLVMap>
        breaker plot (timeline|history|sources) [-w <gwid>] [-s <pathToSettingsFile>]
        breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]
        breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
        breaker stats <gwid> [-s <pathToSettingsFile>]
        breaker listen <far> (<mjdStart> <mjdEnd> | <inLastNMins>) [-s <pathToSettingsFile>]
        breaker listen -d <far> [-s <pathToSettingsFile>]
    
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
    
        FLAGS
        -----
        -h, --help            show this help message
        -s, --settings        the settings file
        -n, --updateNed       update the NED database steam
        -d, --daemon          listen in daemon mode
        -w, --waveId          a gravitational wave ID
    
    
