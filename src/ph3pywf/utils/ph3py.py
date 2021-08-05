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
    from phono3py.file_IO import write_disp_fc3_yaml, write_disp_fc2_yaml
except ImportError:
    Phono3py = None
    
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

@requires(Phono3py, "phono3py not installed!")
def get_displaced_structures(structure, 
                             atom_disp=0.03, 
                             supercell_matrix_fc3=None, 
                             supercell_matrix_fc2=None,
                             yaml_fname_fc3=None, 
                             yaml_fname_fc2=None, 
                             **kwargs):
    """
    Modified based on pymatgen.io.phonopy.get_displaced_structures()
    Replaced all phonopy features with phono3py
    supercell_matrix_fc2: for fc2 calculation
    """

    cutoff_pair_distance = kwargs.get("cutoff_pair_distance", None)
    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)
    primitive_matrix = kwargs.get("primitive_matrix", None),

    ph_structure = get_phonopy_structure(structure)

    if supercell_matrix_fc3 is None:
        supercell_matrix_fc3 = np.eye(3) * np.array((2, 2, 2))

    phonon = Phono3py(
        unitcell=ph_structure, 
        supercell_matrix=supercell_matrix_fc3,
        phonon_supercell_matrix=supercell_matrix_fc2,
    )
    
    phonon.generate_displacements(
        distance=atom_disp,
        cutoff_pair_distance=cutoff_pair_distance,
        is_plusminus=is_plusminus,
        is_diagonal=is_diagonal,
    )
    
    if yaml_fname_fc3 is not None:
        write_disp_fc3_yaml(
            dataset=phonon.dataset,
            supercell=phonon.supercell,
            filename=yaml_fname_fc3,
        )

    # Supercell structures with displacement
    disp_supercells_fc3 = phonon.get_supercells_with_displacements()
    # Original supercell structure
    init_supercell_fc3 = phonon.supercell
    # Structure list to be returned
    structure_list_fc3 = [get_pmg_structure(init_supercell_fc3)]

    for c in disp_supercells_fc3:
        if c is not None:
            structure_list_fc3.append(get_pmg_structure(c))
            
    # For 2nd order (fc2)
    # Original supercell structure
    init_supercell_fc2 = phonon.phonon_supercell
    # Structure list to be returned
    structure_list_fc2 = [get_pmg_structure(init_supercell_fc2)]
    
    if supercell_matrix_fc2 is not None:
        
        phonon.generate_fc2_displacements(
            distance=atom_disp,
            is_plusminus=is_plusminus,
            is_diagonal=is_diagonal,
        )
        
        if yaml_fname_fc2 is not None:
            write_disp_fc2_yaml(
                dataset=phonon.phonon_dataset,
                supercell=phonon.phonon_supercell,
                filename=yaml_fname_fc2,
            )
            
        disp_supercells_fc2 = phonon.get_phonon_supercells_with_displacements()
        for c in disp_supercells_fc2:
            if c is not None:
                structure_list_fc2.append(get_pmg_structure(c))

    return structure_list_fc3, structure_list_fc2


from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_BORN
from phonopy.units import Bohr, Hartree
from phonopy.harmonic.force_constants import show_drift_force_constants
from phono3py.phonon3.fc3 import show_drift_fc3
from phono3py import Phono3py
from phono3py.cui.create_force_constants import parse_forces
from phono3py.file_IO import (parse_disp_fc3_yaml,
                              parse_disp_fc2_yaml,
                              parse_FORCES_FC2,
                              parse_FORCES_FC3,
                              write_fc3_to_hdf5,
                              write_fc2_to_hdf5,
                              read_fc3_from_hdf5,
                              read_fc2_from_hdf5)

def run_thermal_conductivity(phono3py, t_min=0, t_max=1001, t_step=10):
    """
    Evaluate lattice thermal conductivity and write kappa-*.hdf5.
    Phono3py object must be instantiated.
    This function reads disp_fcx.yaml and FORCES_FCx.
    """
    
    # create fc3 and fc2
    if phono3py.phonon_supercell_matrix is not None:
        # if fc2 supercell size is specified, 
        # prepare dataset of fc3 from disp_fc3.yaml and FORCES_FC3
        # and fc2 from disp_fc2.yaml and FORCES_FC2
        print("INFO: preparing fc3 from disp_fc3.yaml and FORCES_FC3")
        phono3py.dataset = parse_forces(phono3py,
                                        # cutoff_pair_distance=None, # optional
                                        force_filename="FORCES_FC3",
                                        disp_filename="disp_fc3.yaml",
                                        fc_type="fc3")
        print("INFO: preparing fc2 from disp_fc2.yaml and FORCES_FC2")
        phono3py.phonon_dataset = parse_forces(phono3py,
                                               force_filename="FORCES_FC2",
                                               disp_filename="disp_fc2.yaml",
                                               fc_type="phonon_fc2")
    else:
        # if fc2 supercell size is NOT specified, 
        # prepare dataset of fc3 and fc2 from disp_fc3.yaml and FORCES_FC3
        print("INFO: preparing fc3 and fc2 from disp_fc3.yaml and FORCES_FC3")
        phono3py.dataset = parse_forces(phono3py,
                                        # cutoff_pair_distance=None, # optional
                                        force_filename="FORCES_FC3",
                                        disp_filename="disp_fc3.yaml",
                                        fc_type="fc3")
        phono3py.phonon_dataset = parse_forces(phono3py,
                                               force_filename="FORCES_FC3",
                                               disp_filename="disp_fc3.yaml",
                                               fc_type="fc2")
    print("INFO: fc3 and fc2 dataset prepared")
    
    # produce fc3 and fc2
    print("INFO: phono3py.produce_fc3: using default paramters")
    phono3py.produce_fc3()
    print("INFO: phono3py.produce_fc2: using default paramters")
    phono3py.produce_fc2()
    
    fc3 = phono3py.fc3
    fc2 = phono3py.fc2
    write_fc3_to_hdf5(fc3)
    write_fc2_to_hdf5(fc2)
    
    show_drift_fc3(fc3)
    show_drift_force_constants(fc2, name="fc2")
    
    phono3py.init_phph_interaction()
    
    phono3py.run_thermal_conductivity(
        temperatures=range(t_min, t_max, t_step), 
        is_isotope=True, 
        boundary_mfp=1e6, # This is to avoid divergence of phonon life time.
        write_kappa=True)

    # Conductivity_RTA object (https://git.io/vVRUW)
    cond_rta = phono3py.get_thermal_conductivity()
    
    return cond_rta


from phonopy.cui.phonopy_script import file_exists
from phono3py.file_IO import parse_disp_fc2_yaml, parse_FORCES_FC2, get_length_of_first_line
from phonopy.file_IO import write_FORCE_SETS

def create_FORCE_SETS_from_FORCES_FCx(forces_filename="FORCES_FC3", disp_filename="disp_fc3.yaml", log_level=1):
    with open(forces_filename, 'r') as f:
        len_first_line = get_length_of_first_line(f)

    if len_first_line == 3:
        file_exists(disp_filename, log_level)
        disp_dataset = parse_disp_fc2_yaml(filename=disp_filename)
        file_exists(forces_filename, log_level)
        parse_FORCES_FC2(disp_dataset, filename=forces_filename)
        if log_level:
            print("Displacement dataset was read from \"%s\"."
                  % disp_filename)
        write_FORCE_SETS(disp_dataset)

        if log_level:
            print("FORCE_SETS has been created.")


from phonopy.file_IO import write_BORN
from atomate.vasp.database import VaspCalcDb

def create_BORN_file_from_tag(tag, db_file):
    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
    
    # get relaxation task document 
    print(f"Extracting doc with keyword: {tag} structure optimization")
    doc_relaxation = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} structure optimization"},
        }
    )
    
    # extract unitcell structure, borns, and epsilon
    unitcell = Structure.from_dict(doc_relaxation["calcs_reversed"][0]["output"]["structure"])
    ph_unitcell = get_phonopy_structure(unitcell)
    borns = np.array(doc_relaxation["calcs_reversed"][0]["output"]["outcar"]["born"])
    epsilon = np.array(doc_relaxation["calcs_reversed"][0]["output"]["epsilon_static"])
    
    # write BORN file
    print("Writing BORN file")
    write_BORN(ph_unitcell, borns, epsilon)