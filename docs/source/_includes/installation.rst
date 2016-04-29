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

```bash
setenv CC clang
```

or, for bash:

```bash
export CC=clang
```

Then try and install breaker again.

**healpy**. If you're having trouble installing healpy try installing the `lastest version from github <https://github.com/healpy/healpy/releases>`_. Download and extract the tarball.

Untar, set your ``MACOSX_DEPLOYMENT_TARGET`` environment variable and install:

```
tar -xvf healpy-1.9.0.tar.gz
cd healpy-1.9.0
setenv MACOSX_DEPLOYMENT_TARGET 10.11
python setup.py install
```



