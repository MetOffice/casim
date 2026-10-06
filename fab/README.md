# Fab Build Scripts for Casim

This directory contains the files for building Casim with Fab. It
needs at least Fab version 2.3.0.


## Building
The build script is a Python script that relies on Fab.
Casim can be compiled by providing the required
command line options to the script.

Please read the Fab documentation (esp the 
[introduction to the Fab base class](https://metoffice.github.io/fab/fab_base/index.html)
for details. See also the next section ('setup') if you need site-specific
modifications.

The actual build is defined by command line parameters to the
``fab_Casim.py`` script. Use the ``-h`` option to list all available
options.

In the ``site_specific`` directory are various site-specific
setups. The main one is called ``default``, and it contains setting
for any compiler currently supported by Fab. But each site can
modify the settings. Have a look at the existing configurations
already contained in Casim, and check the
[Fab documentation](https://metoffice.github.io/fab/fab_base/config.html)
for a full explanation of the available options. 

An example build can be done as follows::

    ./fab_casim.py --site nci --platform gadi --suite gnu  --profile debug

The compilers are selected by specifying a suite (Fab supports out of
the box ``gnu``, ``intel-classic``, ``intel-llvm``, ``nvidia`` and
``cray``), and it will use the corresponding Fortran and C compiler.
If MPI is enabled (which is the default), Fab will search for
corresponding ``mpif90`` and ``mpicc`` compiler wrapper, and verify
that they are indeed of the right suite. If you need to use say
a different C compiler (e.g. use ``gcc`` in an otherwise Intel build),
use the ``-cc`` command line option.

Any Fab script also supports the ``--available-compilers`` flag, which
just lists all compilers that Fab knows about that are available on the
system. Example output (of a system that has gfortran and mpif90 as a
wrapper for gfortran)::

    ----- Available compiler and linkers -----
    Gcc - gcc: gcc
    Mpicc - mpicc-gcc: mpicc
    Gfortran - gfortran: gfortran
    Mpif90 - mpif90-gfortran: mpif90
    Linker - linker-gcc: gcc
    Linker - linker-gfortran: gfortran
    Linker - linker-mpif90-gfortran: mpif90
    Linker - linker-mpicc-gcc: mpicc

The Fab workspace defaults to ``./fab-workspace``, but this can be
changed using the ``--fab-workspace`` command line option.

If the build finished successfully, the library will be in a directory
like ``fab-workspace/casim-debug-gfortran/``, it is
called ``libcasim.a``. The actual directory name will depend on the
options specified of course.

## Setting up site-specific options
If you need site-specific options (e.g. you want to change the
default compiler flags used for your compiler), create
a directory with the name of your site and platform under
``site-specific``. Please check the existing
[Fab documentation](https://metoffice.github.io/fab/fab_base/config.html)
for examples on setting up options, or the existing
site-specific setups under ``fab/site-specific``.
