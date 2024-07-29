ReadMe
========

For installation instructions see ``INSTALL.rst``

CMake
-------

The folder ``CMake/`` contains a CMake based alternative build system for FVCOM. The system is contributed as is, and is not officially maintained. 
The system has been setup for FVCOM v3.2. For newer versions the ``src/CMakeLists.txt`` would most likely have to be edited, and new source 
files would have to be added.

Compilation and running the simulations on Betzy
---------------------------------------------------

Loading the required modules
.................................

Load the following modules or add them to your ``bashrc``. Below, we opt to add them to our ``bashrc``::

  $ module load netCDF/4.7.4-iompi-2020a
  $ module load netCDF-Fortran/4.5.2-iompi-2020a
  $ module load CMake/3.16.4-GCCcore-9.3.0
  $ alias mkbetzy=' FC=mpifort CC=mpicc CXX=mpicpc make '

Note the alias ``mkbetzy`` which will be called later (see below), and can be modified as per your preference.

Get the source code
.....................

Clone the ``master`` branch (source code) of this repository into your ``home`` directory::

  $ cd ~
  $ git clone link/to/master/branch/of/this/repository

Cmake 
.......

Copy the ``fvcom-cmake`` directory, i.e., the contents of the ``betzyfv4fisoc`` branch, 
to your ``home`` directory and make the ``src`` point to your source code location::

  $ cd ~
  $ cd ./fvcom-cmake
  $ rm -f -r src
  $ ln -s /cluster/home/user/fvcom4_fisoc   /cluster/home/user/fvcom-cmake/src

Prepare the experiment directory
..................................

Create a new ``experiment`` directory (say ``fv4ice_runs``)::

  $ cd ~
  $ mkdir fv4ice_runs 

Copy the `fv3_new <https://github.com/abhay26992/FVCOM_Petermann_run_utils>`_ script in this directory. 
Check this file to ensure that ``SRC_DEF`` and ``SRC_RT`` are correct. Also check that the cmake directory is correct. 
Note that unless any changes are made to the steps described above, these should be okay. Create the new experiment as::  

  $ cd ~/fv4ice_runs/
  $ ./fv3_new PF_Q_sg

where, (say) ``PF_Q_sg`` is the experiment name. 

Make 
.....

Go to the ``build`` directory within ``experiment directory/experiment name``::

  $ cd ~/fv4ice_runs/PF_Q_sg/build

You can change the default compiler flags and/or the configuration using the ``ccmake`` tool::

  $ ccmake .

This will start a GUI editor which allows you to change any configuration (except the compiler(s)). 
Use the `CMakeCache.txt  <https://github.com/abhay26992/FVCOM_Petermann_run_utils>`_ file to ensure 
that the correct flags have been toggled ``ON`` before you run the numerical experiments. Then make 
using the alias (``mkbetzy``) that we had defined above in our ``bashrc``::

  $ mkbetzy

Running the numerical experiments
...................................

Copy `fv3_run_slurm <https://github.com/abhay26992/FVCOM_Petermann_run_utils>`_ 
to the ``experiment`` directory (``~/fv4ice_runs``). In this script, ensure that 
the ``RESULTDIR`` (and thus, ``JOBDIR``) are correctly specified.

Create a ``run`` directory within ``experiment directory/experiment name``::

  $ cd ~/fv4ice_runs/PF_Q_sg/
  $ mkdir exp_a

This ``run`` directory contains a ``job submission script``, and the ``start file`` (``namelist file``) and the ``input files`` for that run.

A ``job submission script`` template can be found `here  <https://github.com/abhay26992/FVCOM_Petermann_run_utils>`_. 
In this script, specify the ``SBATCH options``. Again, ensure that the ``RESULTDIR`` (and thus, ``JOBDIR``) are 
correctly specified. 

The ``namelist files`` and ``input files`` required to conduct the subglacial discharge experiments are made publicly available via `Zenodo <https://zenodo.org/records/12803094>`_. 
Note that ``INPUT_DIR = './input/',`` under ``&NML_IO`` in the ``namelist files`` implies that all the input files required to run the numerical experiment in question are contained 
within a single ``input`` directory within the ``run`` directory ``~/fv4ice_runs/PF_Q_sg/exp_a``. Specifically, owing to the large file sizes, the input files are stored in 
the ``projects`` directory (``/cluster/projects/your_project_id/your_directory_name``), and the ``run`` directory (in the ``home`` directory) ``~/fv4ice_runs/PF_Q_sg/exp_a`` contains 
an ``input`` symlink. This, however, can be adjusted according to your preference.

Finally, the simulation is run as::

  $ cd ~/fv4ice_runs
  $ ./fv3_run_slurm PF_Q_sg exp_a run

i.e., ``./fv3_run_slurm experiment_name run_directory_name run`` which starts a *virgin run*. Read the ``./fv3_run_slurm`` script 
to learn about the different run types (e.g., *virgin run*, *restart run*, etc.).




