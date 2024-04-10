=======
ph3pywf
=======


Customized Atomate workflow for thermal conductivity calculation using VASP and Phono3py.


Overview
===========

Ph3pyWF is a software developed to extend the capabilities of the ``atomate`` computational workflow framework for high-throughput lattice thermal conductivity analysis in ceramic materials. 
It provides pre-defined density functional theory (DFT) workflow templates that users can use within the atomate framework to streamline their analyses. 
It also contains utility for analyzing computation workflow result and graph (conductivity vs temperature, density of states, phonon dispersion band structures) plotting.

Ph3pyWF uses finite displacement method for lattice thermal conductivity calculation. It dynamically updates the workflow structure by inserting static DFT calculation subtasks generated from the result of unit-cell optimization task. 

Installation
============

Please install ``phonopy`` and ``phono3py`` before installing ph3pywf 
(please see `phonopy's documentation <https://phonopy.github.io/phonopy/install.html>`_ 
and `phono3py's documentation <https://phonopy.github.io/phono3py/install.html>`_). 

Installation from source code::

    $ git clone https://github.com/MatFrontier/ph3pywf.git
    $ cd ph3pywf
    $ pip install -e .


Requirements
============

Python version >=3.8

Ph3pyWF is used as an extension of ``atomate`` computational framework 
and is intended to run on linux-based HPC system. 
For more user-friendly interaction, it is recommended to use this software in Jupyter development environment. 

<Add instruction for configuring launchpad and MongoDB>

Usage
=====

<WIP>




.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.0.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
