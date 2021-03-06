Command-Line Tutorial
=====================

This is a tutorial for using the CL ``breaker`` tools. To use ``breaker`` in your python scripts, see the full documentation.

Example Maps
^^^^^^^^^^^^

To follow along with this tutorial, you may want to grab some of the publically available wave maps. Here are some maps for `GW150914 (aka G184098) <https://losc.ligo.org/events/GW150914/>`_ and `GW151226 (aka G211117) <https://losc.ligo.org/events/GW151226/>`_.

Settings File
^^^^^^^^^^^^^

Most settings for ``breaker`` are found and set in the ``breaker.yaml`` settings file. By default these are placed in ``~/.config/breaker/breaker.yaml`` and running breaker for the first time with generate a default settings file at this location.

As well as some fundamental logging settings, the settings file needs to contain certain database info and ssh tunnel setups, which obviously need to remain private. If you need these settings, ask Dave Young.

**General settings** include where ``breaker`` is to store maps and where to place the output files (plots, faker lists, FITS stamps etc).

.. code-block:: yaml  

    # I/O SETTINGS
    gw maps directory: /Users/Dave/.config/breaker/maps
    output directory: "/Users/Dave/Desktop/breaker-output"

There are also **wave specific settings**, indicting where to find the most accurate map for the wave for plots, timeline ranges and details of sky-areas to show in the plots. Here's an example for the first burst:

.. code-block:: yaml

    G184098:
        human-name: GW150914    
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

To update the PS1 footprint table in the breaker database and associate these footprints with the GWs run the following command:

.. code-block:: bash  

    breaker update

and to also download all the overlapping NED sources and add them to the database, use the ``-n`` flag:

.. code-block:: bash  

    breaker update -n

Skymaps
^^^^^^^

The Ligo-Virgo likeihood maps are a Healpix rendering of the sky and delivered in a FITS binary table format. To convert these maps to a PDF image (molleweide projection) and an all-sky FITS image (normal mercator projection) use the `skymap` command. For example to convert the *G211117/GW151226* event baystar-map run the following:

.. code-block:: bash 
     
    breaker skymap G211117 ~/.config/breaker/maps/G211117/bayestar.fits.gz
    
This will generate the FITS and PDF images in the CWD:

.. image:: https://i.imgur.com/n5NDDZj.png
        :width: 800px
        :alt: GW151226 FITS image

.. image:: https://i.imgur.com/GynPdBY.png
        :width: 800px
        :alt: GW151226 PDF Mollweide Projection

Plots
^^^^^

Once you have the settings file organised and some sky-maps maps downloaded from graceDB you can start plotting.

Timeline and History Plots
--------------------------

It's possible to plot a timeline of observations over the likelihood map for each wave. By choosing the ``breaker plot timeline`` command, the code plots from the epoch of the wave detection (in settings file) forward in time. Alternatively by choosing the ``breaker plot history`` command, the code will plot from now back in time over the last days, weeks and months. 

For example the following command will produce a set of plots for the wave G184098 = GW150914:

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

Over-plotting NED Sources
------------------------=

If the database tables are brought up-to-date using the ``breaker -n update`` command, it is possible to overplot NED sources found within the wave campaign footprint. More fine-grained control of these plots can be gained by scripting solutions by importing ``breaker`` into your own python code. But running the command:

.. code-block:: bash

    breaker plot sources -w G184098

produces this plot:

.. image:: https://i.imgur.com/vn8tTJy.png
        :width: 800px
        :alt: NED source found in wave footprint 

    

Multi-Panel Comparison Plots
----------------------------

The localisation maps for each wave come in various flavours at different stages of processing and with varying degrees of accuracy. It can be useful to produce a multi-panel plot of these maps to compare them. The following command will generate this plot, with a normalise colour range so the probabilities on each map can be directly compared.

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

It might be useful at some point to determine the completeness of our campaigns. The ``faker`` command will take a PS1 exposure and extract out all NED galaxy sources with redshift and semi-major axis measurements in the FOV of that exposure. For each of those galaxies a fake transient is placed at a random location within the galaxy semi-major axes. An extra 17.6% locations are then randomly distributed throughout the area of the exposure to give a overall total of 85% galaxy associations and 15% 'orphans'. Two versions of the fake source catalogue are output, *trimmed* and *complete*, which can then be used to test our pipelines end-to-end.

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

will rattle through the ATLAS and PS1 footprints in chronological order and determine some cumulative stats, including the total sky-area covered (squ. deg.) and the total likelihood covered (in 2-dimensions only):

.. code-block:: bash

    0/1449.  MJD: 57382.29419. AREA: 30.67. PROB: 0.00923. SURVEY: atlas
    1/1449.  MJD: 57382.302442. AREA: 59.51. PROB: 0.02116. SURVEY: atlas
    2/1449.  MJD: 57382.313403. AREA: 87.18. PROB: 0.02246. SURVEY: atlas
    3/1449.  MJD: 57384.216272. AREA: 87.18. PROB: 0.02246. SURVEY: ps1
    4/1449.  MJD: 57384.216771. AREA: 87.18. PROB: 0.02246. SURVEY: ps1
    5/1449.  MJD: 57384.221982. AREA: 87.18. PROB: 0.02246. SURVEY: ps1 
    ...
    ...

Download Recently Detected Wave Maps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before running the ``listen`` command you need to create a ``.netrc`` file with your GraceDb credentials (with 600 permissions). `See here for a tutorial <https://dcc.ligo.org/public/0118/G1500442/010/ligo-virgo-emfollowup-tutorial.html>`_

Alternatively you can add the GraceDB robot credentials into breaker's settings file. Just take the username and password found in your ``.netrc`` and add them to ``breaker.yaml`` as follows:

.. code-block:: yaml 
    
    graceDB robot credentials: 
        username: <yourLigoUsername>
        password: <yourLigoRobotPassword>

Breaker will first check its own settings file for the GraceDB credentials and then the ``.netrc`` file in your home directory, in that order.
        
The ``listen`` command is used to connect to `graceDB <https://gracedb.ligo.org>`_ and download the maps from recently detected waves. You can connect either once and download all maps within a time range, or connect in daemon mode to ping graceDB every 60 secs for new maps.

To connect and download maps between MJDs 57382. and 57384. with a false alarm rate lower limit of 1e-5 Hz:

.. code-block:: bash 
 
    > breaker listen 1e-5 57382. 57384.
    NEW GRAVITATIONAL WAVE EVENT FOUND ...
        GraceDB ID: G211117
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading LALInference_skymap.fits.gz
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading bayestar.fits.gz
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading LIB_skymap.fits.gz

    METADATA FOR G211117 ...
    Date Added to GraceDB: 2015-12-26 03:40:00 UTC
    Detection Interferometers: H1,L1
    Detection Pipeline: gstlal
    Discovery Group: CBC
    Discovery Search Type: HighMass
    Event Submitter: gstlalcbc
    False Alarm Rate: 3.33262857227e-11 Hz
    GPS Event Time: 1135136350.647758
    GraceDB ID: G211117
    Hanford MJD: 57382.152009812
    Livingston MJD: 57382.1520098019
    MJD Difference Seconds: 0.0008749962
    Maps:
      LALInference3d.fits.gz: false
      LALInference_skymap.fits.gz: true
      LIB_skymap.fits.gz: false
      bayestar.fits.gz: true
      bayestar3d.fits.gz: false
      skymap.fits.gz: false

Or to download maps within the last 15 mins:

.. code-block:: bash 
 
    > breaker listen 1e-5 15
    
To connect in daemon mode:

.. code-block:: bash 

    > breaker listen -d 1e-5
    NEW GRAVITATIONAL WAVE EVENT FOUND ...
        GraceDB ID: G211117
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading LALInference_skymap.fits.gz
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading bayestar.fits.gz
    NEW MAP FOUND FOR GW EVENT G211117 ...
        Downloading LIB_skymap.fits.gz

    METADATA FOR G211117 ...
    Date Added to GraceDB: 2015-12-26 03:40:00 UTC
    Detection Interferometers: H1,L1
    Detection Pipeline: gstlal
    Discovery Group: CBC
    Discovery Search Type: HighMass
    Event Submitter: gstlalcbc
    False Alarm Rate: 3.33262857227e-11 Hz
    GPS Event Time: 1135136350.647758
    GraceDB ID: G211117
    Hanford MJD: 57382.152009812
    Livingston MJD: 57382.1520098019
    MJD Difference Seconds: 0.0008749962
    Maps:
      LALInference3d.fits.gz: false
      LALInference_skymap.fits.gz: true
      LIB_skymap.fits.gz: false
      bayestar.fits.gz: true
      bayestar3d.fits.gz: false
      skymap.fits.gz: false

    NEW GRAVITATIONAL WAVE EVENT FOUND ...
        GraceDB ID: G194575
    NEW MAP FOUND FOR GW EVENT G194575 ...
        Downloading skymap.fits.gz

    METADATA FOR G194575 ...
    Date Added to GraceDB: 2015-10-22 13:35:44 UTC
    Detection Interferometers: H1,L1
    Detection Pipeline: gstlal
    Discovery Group: CBC
    Discovery Search Type: LowMass
    Event Submitter: gstlalcbc
    False Alarm Rate: 9.65424329993e-08 Hz
    GPS Event Time: 1129556016.942353
    GraceDB ID: G194575
    Hanford MJD: 57317.5648143102
    Livingston MJD: 57317.5648141476
    MJD Difference Seconds: 0.0140454769
    Maps:
      LALInference3d.fits.gz: false
      LALInference_skymap.fits.gz: false
      LALInference_skymap.fits.gz: false
      LIB_skymap.fits.gz: false
      bayestar.fits.gz: false
      skymap.fits.gz: true

    0 archived and 2 events found, will try again in 60 secs
    2 archived and 0 events found, will try again in 60 secs
    2 archived and 0 events found, will try again in 60 secs
    ...
    
Note the first time ``breaker`` connects to graceDB in daemon mode it downloads all maps from the beginning of operations (2015-09-01 00:00:00 UTC).

Maps are downloaded to whatever directory you have set as ``gw maps directory`` in the breaker settings file.

.. image:: https://i.imgur.com/kkOlSlp.png
        :width: 800px
        :alt: maps and metadata

Alongside the maps you will find a ``meta.yaml`` file containing some pertinent data about the event as reported in GraceDB.

.. code-block:: yaml 
    
    Date Added to GraceDB: 2015-12-26 03:40:00 UTC
    Detection Interferometers: H1,L1
    Detection Pipeline: gstlal
    Discovery Group: CBC
    Discovery Search Type: HighMass
    Event Submitter: gstlalcbc
    False Alarm Rate: 3.33262857227e-11 Hz
    GPS Event Time: 1135136350.647758
    GraceDB ID: G211117
    Hanford MJD: 57382.152009812
    Livingston MJD: 57382.1520098019
    MJD Difference Seconds: 0.0008749962
    Maps:
      LALInference3d.fits.gz: false
      LALInference_skymap.fits.gz: true
      LIB_skymap.fits.gz: false
      bayestar.fits.gz: true
      bayestar3d.fits.gz: false
      skymap.fits.gz: false


Transient Location Likelihoods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To calculate which likelihood contour a transient/location lies within for a given gravity event run the command:

.. code-block:: bash 
    
    breaker contour <gwid> <ra> <dec> 

So for the event 'G211117' (GW151226) if running:

.. code-block:: bash 
    
    breaker contour G211117 60.264 41.7032 

gives the following output:

.. code-block:: text

    The transient lies within the inner 10% likelihood contour of event G211117 




