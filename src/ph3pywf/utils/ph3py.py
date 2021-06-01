# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import numpy as np
from monty.dev import requires
from monty.serialization import loadfn

from pymatgen.core import Lattice, Structure
from pymatgen.phonon.bandstructure import (
    PhononBandStructure,
    PhononBandStructureSymmLine,
)
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.symmetry.bandstructure import HighSymmKpath

# try:
#     from phonopy import Phonopy
#     from phonopy.file_IO import write_disp_yaml
#     from phonopy.structure.atoms import PhonopyAtoms
# except ImportError:
#     Phonopy = None

### Additional import is added here

try:
    from phono3py import Phono3py
    from phono3py.file_IO import write_disp_fc3_yaml
except ImportError:
    Phono3py = None
    
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

@requires(Phono3py, "phono3py not installed!")
def get_displaced_structures(structure, atom_disp=0.01, supercell_matrix=None, yaml_fname=None, **kwargs):
    """
    Modified based on pymatgen.io.phonopy.get_displaced_structures()
    Replaced all phonopy features with phono3py
    """

    cutoff_pair_distance = kwargs.get("cutoff_pair_distance", None)
    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)

    ph_structure = get_phonopy_structure(structure)

    if supercell_matrix is None:
        supercell_matrix = np.eye(3) * np.array((2, 2, 2))

    phonon = Phono3py(unitcell=ph_structure, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(
        distance=atom_disp,
        cutoff_pair_distance=cutoff_pair_distance,
        is_plusminus=is_plusminus,
        is_diagonal=is_diagonal,
    )

    if yaml_fname is not None:
        dataset = phonon.dataset
        write_disp_fc3_yaml(
            dataset=dataset,
            supercell=phonon.get_supercell(),
            filename=yaml_fname,
        )

    # Supercell structures with displacement
    disp_supercells = phonon.get_supercells_with_displacements()
    # Perfect supercell structure
    init_supercell = phonon.get_supercell()
    # Structure list to be returned
    structure_list = [get_pmg_structure(init_supercell)]

    for c in disp_supercells:
        if c is not None:
            structure_list.append(get_pmg_structure(c))

    return structure_list


from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_BORN
from phonopy.units import Bohr, Hartree
from phonopy.harmonic.force_constants import show_drift_force_constants
from phono3py.phonon3.fc3 import show_drift_fc3
from phono3py import Phono3py
from phono3py.file_IO import (parse_disp_fc3_yaml,
                              parse_disp_fc2_yaml,
                              parse_FORCES_FC2,
                              parse_FORCES_FC3,
                              read_fc3_from_hdf5,
                              read_fc2_from_hdf5)

def run_thermal_conductivity(phono3py, t_min=0, t_max=1001, t_step=10):
    # Create fc3 and fc2 from disp_fc3.yaml and FORCES_FC3
    disp_dataset = parse_disp_fc3_yaml(filename="phonopy_disp.yaml")
    forces_fc3 = parse_FORCES_FC3(disp_dataset, filename="FORCES_FC3")
    phono3py.produce_fc3(forces_fc3,
                         displacement_dataset=disp_dataset,
                         symmetrize_fc3r=True)
    print("INFO: phono3py.produce_fc3: symmetrize_fc3r=True")
    fc3 = phono3py.get_fc3()
    fc2 = phono3py.get_fc2()

    show_drift_fc3(fc3)
    show_drift_force_constants(fc2, name='fc2')
    
    phono3py.init_phph_interaction()
    
    phono3py.run_thermal_conductivity(
        temperatures=range(t_min, t_max, t_step), 
        boundary_mfp=1e6, # This is to avoid divergence of phonon life time.
        write_kappa=True)

    # Conductivity_RTA object (https://git.io/vVRUW)
    cond_rta = phono3py.get_thermal_conductivity()
    
    return cond_rta