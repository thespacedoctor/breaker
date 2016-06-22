Usage
======

.. code-block:: bash 
   
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
    