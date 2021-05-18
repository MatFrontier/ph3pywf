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
def get_displaced_structures(pmg_structure, atom_disp=0.01, supercell_matrix=None, yaml_fname=None, **kwargs):
    r"""
    Modified based on pymatgen.io.phonopy.get_displaced_structures()
    Replaced all phonopy features with phono3py
    """

    cutoff_pair_distance = kwargs.get("cutoff_pair_distance", None)
    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)

    ph_structure = get_phonopy_structure(pmg_structure)

    if supercell_matrix is None:
        supercell_matrix = np.eye(3) * np.array((1, 1, 1))

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
