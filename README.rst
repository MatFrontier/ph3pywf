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

Installation from source code:

.. code-block::

    $ git clone https://github.com/MatFrontier/ph3pywf.git
    $ cd ph3pywf
    $ pip install -e .


Requirements
============

Python version >= 3.8

Ph3pyWF is used as an extension of ``atomate`` computational framework 
and is intended to run on linux-based HPC system. 
For more user-friendly interaction, it is recommended to use this software in Jupyter development environment. 

Atomate, FirWorks installation and configuration guide: `link <https://atomate.org/installation.html>`_

Usage (Example Jupyter Notebook)
================================

Create workflow for monoclinic ZrO2:

.. code-block:: python

    import numpy as np
    from pymatgen.ext.matproj import MPRester
    from pymatgen.core import Structure
    from fireworks import LaunchPad
    import os
    from ph3pywf.workflows.core import wf_phono3py

    # Required parameters
    material_id = "mp-2858" # Material id to be searched on MP
    supercell_size_fc3 = (2,2,2)
    supercell_size_fc2 = None
    cutoff_pair_distance = 4.0

    # Get working dir 
    working_dir = os.getcwd()

    # Get materials structure
    # If no local structure file, get from MP and save unitcell
    # structure as "struct_unitcell.json"
    if f"{material_id}_unitcell.json" in os.listdir(working_dir):
        # Read saved local structure file
        print("Found local structure file")
        struct_unitcell = Structure.from_file(f"{material_id}_unitcell.json")

    else:
        print("Getting structure from MaterialsProject")
        # api_key = ''
        with MPRester() as mpr:
            struct_unitcell = mpr.get_structure_by_material_id(material_id)

        # Save a local structure file to avoid accessing MP every time
        with open(f"{material_id}_unitcell.json", "w") as fh:
            fh.write(struct_unitcell.to_json())
    
    # Specify more detailed calculation parameters
    # Create the workflow
    c = {
        "supercell_size_fc3": supercell_size_fc3, 
        "supercell_size_fc2": supercell_size_fc2,
        "cutoff_pair_distance": cutoff_pair_distance,
        "is_nac": True,
        "USER_INCAR_SETTINGS": {
            "GGA":"CA",
            "EDIFF":1.0e-09,
            "EDIFFG":-1.0e-05
        },
        "USER_INCAR_SETTINGS_STATIC": {
            "EDIFF":1.0e-09
        },
        "USER_KPOINTS_SETTINGS": {"reciprocal_density": 64},
        "USER_KPOINTS_SETTINGS_STATIC": {"reciprocal_density": 32},
        "metadata": {
            "tags": [
                "validation",
                "primitive unitcell",
                "NAC"
            ],
        }
    }

    print("Creating workflow...")
    workflow = wf_phono3py(structure=struct_unitcell, 
                        c=c,
                        )
    print("Created workflow")

    # Initialize the launchpad and add our workflow
    print("Sending to LaunchPad...")
    launchpad = LaunchPad.auto_load()
    launchpad.add_wf(workflow)
    print("Sent to LaunchPad")
    print("=== Done ===")


Output of above code below. Use tag(task_label) to query submitted workflow and corresponding calculation results.

.. code-block::

    Creating workflow...
    tag = "2022-03-13-02-15-33-975230"
    {task_label: {$regex:"2022-03-13-02-15-33-975230"}}
    Created workflow
    Sending to LaunchPad...


Post analysis:

.. code-block:: python

    from ph3pywf.utils.post_analysis import Ph3py_Result

    task_label = "2022-03-13-02-15-33-975230"
    path_to_db_json = "/home/jovyan/atomate/config/db.json"
    ref_filenames = ["SAMPLE_THERMAL_CONDUCTIVITY_REF.csv"]
    ref_labels = ["SAMPLE_REF_LABEL"]
    plot_initial = True
    plot_dircs = False
    ymax = 500
    fig_size = (12,9)

    result = Ph3py_Result(task_label, path_to_db_json)
    result.plot_thermal_conductivity(
        ref_filenames=ref_filenames, 
        ref_labels=ref_labels, 
        plot_initial=plot_initial, 
        plot_dircs=plot_dircs, 
        ymax=ymax, 
        fig_size=fig_size
    )
    result.plot_bs()
    result.plot_dos()


Project Structure
=================
.. code-block:: text

    Ph3pyWF/
    ├── docs/                       # Documentation files
    │
    ├── src/ph3pywf/
    │   ├── firetasks
    │   │   └── core.py             # Fundamental building blocks for workflows.
    │   │                             Contains task for dynamically mutate workflow, and task for post-DFT LTC calculation
    │   ├── fireworks
    │   │   └── core.py             # Firetasks to execute in sequence
    │   │                             Contains Fireworks modified from Atomate pre-defined FWs
    │   │
    │   ├── utils                   # Helper functions and configs
    │   │   ├── guard_once.py       # CLI script to check FWs progress and rerun if there's transient error
    │   │   ├── guard_ph3pywf.py    # CLI script to continuously check FWs progress and rerun if there's transient error
    │   │   ├── mission_control.py  # Helpers for guard_once and guard_ph3pywf
    │   │   ├── ph3py.py            # Helpers for phono3py, database, and in-workflow post-analysis
    │   │   ├── Ph3pyRelaxSet.yaml  # VASP input parameter set for lattice relaxation DFT job
    │   │   ├── Ph3pyStaticSet.yaml # VASP input parameter set for static DFT job
    │   │   ├── post_analysis.py    # Helpers for post-workflow analysis
    │   │   ├── sets.py             # VASP InputSet classes
    │   │   ├── tmp_fix.py          # Script to update band structure
    │   │   └── VASPIncarBase.yaml  # Base VASP input parameter
    │   │
    │   └── workflows
    │           └── core.py         # Phono3py calculation workflow
    │
    ├── tests/                      # Test scripts
    │
    ├── ...                         # Scaffold files
    │
    ├── README.md                   # Project information and usage instructions
    ├── setup.cfg                   # The actual config used by setup.py
    └── setup.py                    # Setup script for installing the package


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.0.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
