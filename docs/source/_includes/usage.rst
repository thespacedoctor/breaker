Usage
======

.. code-block:: bash 
   
    breaker update [-s <pathToSettingsFile>]
    breaker plot (timeline|history|sources [<gwid>]) [-s <pathToSettingsFile>]
    breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]
    breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
    breaker stats <gwid> [-s <pathToSettingsFile>]

    -h, --help            show this help message
    -s, --settings        the settings file
    plot                  update the gravitational wave plots
    timeline              observations looking forward from date of GW detection
    history               observations from the past x days
    faker                 generate a catalogue of simulated transient source in ps1 FP
    stats                 output some stats for the GW surveys
    