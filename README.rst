=======
breaker 
=======

*A python package and command-line tools for the PanSTARRS & ATLAS LIGO-VIRGO (PSAT) group to aid surveys of the likely sky-locations of LIGO-VIRGO discovered Gravitational Waves*.

Here's a summary of what's included in the python package:

.. include:: /classes_and_functions.rst

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
    
    

Documentation
=============

Documentation for breaker is hosted by `Read the Docs <http://breaker.readthedocs.org/en/stable/>`__ (last `stable version <http://breaker.readthedocs.org/en/stable/>`__ and `latest version <http://breaker.readthedocs.org/en/latest/>`__).

Installation
============

The easiest way to install breaker us to use ``pip``:

.. code:: bash

    pip install breaker

Or you can clone the `github repo <https://github.com/thespacedoctor/breaker>`__ and install from a local version of the code:

.. code:: bash

    git clone git@github.com:thespacedoctor/breaker.git
    cd breaker
    python setup.py install

To upgrade to the latest version of breaker use the command:

.. code:: bash

    pip install breaker --upgrade

Troubleshooting
^^^^^^^^^^^^^^^

If you're having trouble with the installation here are a few things to try:

**Astropy and Clang**. On Mac OS you may have to set your C-compiler to clang before astropy will install. So before the breaker installation, try:

.. code:: bash

    setenv CC clang


or, for bash:

.. code:: bash

    export CC=clang


Then try and install breaker again.

**healpy**. If you're having trouble installing healpy try installing the `latest version from github <https://github.com/healpy/healpy/releases>`_. Download and extract the tarball.

Untar, set your ``MACOSX_DEPLOYMENT_TARGET`` environment variable and install:

.. code:: bash

    tar -xvf healpy-1.9.0.tar.gz
    cd healpy-1.9.0
    setenv MACOSX_DEPLOYMENT_TARGET 10.11
    python setup.py install

    




Development
-----------

If you want to tinker with the code, then install in development mode.
This means you can modify the code from your cloned repo:

.. code:: bash

    git clone git@github.com:thespacedoctor/breaker.git
    cd breaker
    python setup.py develop

`Pull requests <https://github.com/thespacedoctor/breaker/pulls>`__
are welcomed!


Issues
------

Please report any issues
`here <https://github.com/thespacedoctor/breaker/issues>`__.

License
=======

Copyright (c) 2016 David Young

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

