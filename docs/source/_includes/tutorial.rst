Tutorial
========

Example Maps
^^^^^^^^^^^^

To follow along with this tutorial, you may want to grab some of the publically available wave maps. Here are some maps for `GW150914 (aka G184098) <https://losc.ligo.org/events/GW150914/>`_ and `GW151226 (aka G211117) <https://losc.ligo.org/events/GW151226/>`_.

Settings File
^^^^^^^^^^^^^

Most settings for ``breaker`` are found and set in the ``breaker.yaml`` settings file. By default these are placed in ``~/.config/breaker/breaker.yaml`` and running breaker for the first time with generate a default settings file at this location.

As well as some fundamental logging settings, the settings file needs to contain certain database info and ssh tunnel setups, which obvious need to remain private. If you need these settings, ask Dave Young.

**General settings** include where ``breaker`` is to store maps and where to place the output files (plots, faker lists, FITS stamps etc).

.. code-block:: yaml  

    # I/O SETTINGS
    gw maps directory: /Users/Dave/.config/breaker/maps
    output directory: "/Users/Dave/Desktop/breaker-output"

There are also **wave specific settings**, indicting where to find the most accurate map for the wave for plots, timeline ranges and details of skyareas to show in the plots. Here's an exmaple for the first burst:

.. code-block:: yaml

    G184098:
        time:
            mjdStart: 57279.90
            mjdEnd: 57369.90
        plot:
            raRange: 48.  # CENTRAL WIDTH IN DEGREES
            decRange: 45.  # CENTRAL HEIGHT IN DEGREES
            centralCoordinate: [141., 0.0]
        mapPath: "/Users/Dave/.config/breaker/maps/G184098/LALInference_skymap.fits"

Refreshing the Database
^^^^^^^^^^^^^^^^^^^^^^^

To update the PS1 footprint table in breaker database and associate these footprints with the GWs run the following command:

.. code-block:: bash  

    breaker update

and to also download all the overlapping NED sources and add them to the database, use the ``-n`` flag:

.. code-block:: bash  

    breaker update -n

Plots
^^^^^

Once you have the settings file organised and some skymaps maps downloaded from graceDB you can start plotting.

Timeline and History Plots
--------------------------

It's possible to plot a timeline of obsevations over the likelihood map for each wave. By choose the ``breaker plot timeline`` command, the code plots from the epoch of the wave detection (in settings file) forward in time. Alternatively by choosing the ``breaker plot history`` command, the code will plot from now back in time over the last days, weeks and months. 

For example the follow command will produce a set of plots for the wave G184098 = GW150914:

.. code-block:: bash 
     
    breaker plot timeline -w G184098

The plots produced in the output directory (from settings file) are:

.. code-block:: bash 
    
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_in_First_3_Days_of_Wave_Detection_tan.png
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_Between_3-10_Days_of_Wave_Detection_tan.png
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_Between_10-17_Days_of_Wave_Detection_tan.png
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_Between_17-24_Days_of_Wave_Detection_tan.png
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_Between_24-31_Days_of_Wave_Detection_tan.png
    G184098_Probability_Map_PS1_Footprints_and_Transients_Discovered_gt_31_Days_of_Wave_Detection_tan.png
 
and look similar to this:

.. image:: https://i.imgur.com/EC0oyhq.png
        :width: 800px
        :alt: Example Timeline Plot

To run the history command for the same wave:

.. code-block:: bash 
     
    breaker plot history -w G184098

Note running either of these commands without a GWID will generate the timeline/history plots for *all* waves found in your settings file:

.. code-block:: bash 
     
    breaker plot timeline

Alongside the PNG plots, a FITS image is also generated showing the same cutout sky-area as the plots. The signal in the FITS image scales with the probability in the Healpix map.

.. image:: https://i.imgur.com/PXcsfmw.png
        :width: 1000px
        :alt: FITS image of Healpix map

Overplotting NED Sources
------------------------

If the database tables are brought up-to-date using the ``breaker -n update`` command, it is possible to overplot NED sources found within the wave campaign footprint. More fine-grained control of these plots can be gained by scripting solutions by importing ``breaker`` into your own python code. But running the command:

.. code-block:: bash

    breaker plot sources -w G184098

produces this plot:

.. image:: https://i.imgur.com/vn8tTJy.png
        :width: 800px
        :alt: NED source found in wave footprint 

    

Mutli-Panel Comparison Plots
----------------------------

The localisation maps for each wave come in various flavours at different stages of processing and with varying degrees of accuracy. It can be useful to produce a multi-panel plot of these maps to compare them. The following command with generate this plot, with a normalise colour range so the probabilities on each map can be directly compared.

.. code-block:: bash 

    breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]

So for example:

.. code-block:: bash 

    breaker plot comparison G211117 /Users/Dave/git_repos/breaker/breaker/plots/tests/input

produces the following plot in the output directory found in the settings file.

.. image:: https://i.imgur.com/9jubCq2.png
        :width: 1000px
        :alt: GW151226 4 Panel Comparison Plot

Fake Source Catalogues
^^^^^^^^^^^^^^^^^^^^^^

It might be useful at some point to determine the completeness of our campaigns. The ``faker`` command will take a PS1 exposure and extract out all NED galaxy sources with redshift and semi-major axis measurements in the FOV of that exposure. For each of those galaxies a fake transient is placed at a random loction within the galaxy semi-major axes. An extra 17.6% locations are then randomly distributed throughout the area of the exposure to give a overall total of 85% galaxy associations and 15% 'orphans'. Two versions of the fake source catalogue are output, *trimmed* and *complete*, which can then be used to test our pipelines end-to-end.

**Trimmed** example:

.. code-block:: bash 
    
    index,ra,dec,i-mag
    0001,132.76954,4.56831,17.50
    0002,132.70450,4.55963,18.76
    0003,132.81176,4.58280,18.86
    0004,132.74161,4.49493,17.46
    0005,132.82488,4.48862,18.99
    0006,132.71868,4.45854,19.31
    0007,132.60267,4.61480,18.18
    0008,132.59662,4.60154,17.76
    ...

**Complete** example:

.. code-block:: bash 
    
    index,ra,dec,i-mag,redshift,galaxy-id,2mass-k-mag,2mass-k-mag-error
    0001,132.76954,4.56831,17.50,0.073,"SDSS J085105.10+043414.0",15.00,0.14
    0002,132.70450,4.55963,18.76,0.095,"SDSS J085048.39+043335.7",14.45,null
    0003,132.81176,4.58280,18.86,0.071,"SDSS J085114.79+043453.7",14.58,null
    0004,132.74161,4.49493,17.46,0.095,"SDSS J085057.98+042943.8",14.79,0.12
    0005,132.82488,4.48862,18.99,0.071,"SDSS J085118.00+042918.8",null,null
    0006,132.71868,4.45854,19.31,0.077,"SDSS J085052.02+042732.4",null,null
    0007,132.60267,4.61480,18.18,0.097,"SDSS J085024.94+043654.9",15.16,0.17
    0008,132.59662,4.60154,17.76,0.077,"SDSS J085023.19+043602.4",null,null
    ...


Campaign Stats
^^^^^^^^^^^^^^

The ``stats`` command can be run to generate some stats for a given wave survey campaign. For example the command:

.. code-block:: bash 
    
    breaker stats G211117

will rattle throught the ATLAS and PS1 footprints in cronological order and determine some cummulative stats, including the total sky-area covered (squ. deg.) and the total likelihood covered (in 2-dimensions only):

.. code-block:: bash

    0/1449.  MJD: 57382.29419. AREA: 30.67. PROB: 0.00923. SURVEY: atlas
    1/1449.  MJD: 57382.302442. AREA: 59.51. PROB: 0.02116. SURVEY: atlas
    2/1449.  MJD: 57382.313403. AREA: 87.18. PROB: 0.02246. SURVEY: atlas
    3/1449.  MJD: 57384.216272. AREA: 87.18. PROB: 0.02246. SURVEY: ps1
    4/1449.  MJD: 57384.216771. AREA: 87.18. PROB: 0.02246. SURVEY: ps1
    5/1449.  MJD: 57384.221982. AREA: 87.18. PROB: 0.02246. SURVEY: ps1 
    ...
    ...

Download Recentlt Detected Wave Maps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``listen`` command is used to connect to `graceDB <https://gracedb.ligo.org>`_ and download the maps from recently detected waves. You can connect either once and download all maps within a time range, or connect in daemon mode to ping graceDB every 60 secs for new maps.

To connect and download maps between MJDs 57382. and 57384. with a false alarm rate lower limit of 1e-7 Hz:

.. code-block:: bash 
 
    > breaker listen 1e-7 57382. 57384.
    Downloading bayestar.fits.gz for GW event G211117

Or to download maps within the last 15 mins:

.. code-block:: bash 
 
    > breaker listen 1e-7 15
    
To connect in daemon mode:

.. code-block:: bash 

    > breaker listen -d 1e-7
    Downloading bayestar.fits.gz for GW event G211117
    Downloading skymap.fits.gz for GW event G194575
    2 recent events found, will try again in 60 secs
    No recent events, will try again in 60 secs
    
Note the first time ``breaker`` connects to graceDB in daemon mode it downloads all maps from the beginning of time.




