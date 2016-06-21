breaker 
=========================

*Tools used by the PanSTARRS & ATLAS teams when surveying the likely sky-locations of LIGO-VIRGO discovered Gravitational Waves*.

Usage
======

.. code-block:: bash 
   
    breaker update [-n] [-s <pathToSettingsFile>]
    breaker plot (timeline|history|sources [<gwid>]) [-s <pathToSettingsFile>]
    breaker plot comparison <gwid> <pathToMapDirectory> [-s <pathToSettingsFile>]
    breaker faker <ps1ExpId> [-s <pathToSettingsFile>]
    breaker stats <gwid> [-s <pathToSettingsFile>]

    -h, --help            show this help message
    -s, --settings        the settings file
    -n, --updateNed       update the NED database steam
    plot                  update the gravitational wave plots
    timeline              observations looking forward from date of GW detection
    history               observations from the past x days
    faker                 generate a catalogue of simulated transient source in ps1 FP
    stats                 output some stats for the GW surveys
    
Documentation
=============

Documentation for breaker is hosted by `Read the Docs <http://breaker.readthedocs.org/en/stable/>`__ (last `stable version <http://breaker.readthedocs.org/en/stable/>`__ and `latest version <http://breaker.readthedocs.org/en/latest/>`__).

Installation
============

Currently you need to manually install ``dryxPython`` and ``HMpTy`` packages before installing breaker.

.. todo::

    - remove dependencies on ``dryxPython`` and ``HMpTy``

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

**healpy**. If you're having trouble installing healpy try installing the `lastest version from github <https://github.com/healpy/healpy/releases>`_. Download and extract the tarball.

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

