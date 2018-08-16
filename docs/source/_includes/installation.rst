Installation
============

The easiest way to install breaker is via Anaconda. For some instructions for installing Anaconda `see here <http://astronotes.co.uk/blog/2017/10/04/An-Astronomer's-Guide-to-dotstar-Conda.html>`__.

Once you have Anaconda installed, then create and activate a new conda environment:

.. code:: bash

    conda create -n breaker python=2.7 pip
    source activate breaker

Now do a conda install of healpy before installing breaker (the pip install of healpy seems very flaky):

.. code:: bash

    conda install -c conda-forge healpy

Finally install breaker:

.. code:: bash

    pip install breaker

Installing a Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    



