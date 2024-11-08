# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

# TODO: keep only the phono3py auxilary functions here, and rename as "ph3py_aux" perhaps

import numpy as np
from monty.dev import requires

from pymatgen.core import Structure
from pymatgen.phonon.bandstructure import (
    PhononBandStructureSymmLine,
)
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.symmetry.bandstructure import HighSymmKpath

import gridfs
import zlib
import json
import os
from monty.json import MontyEncoder

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
def get_displaced_structures(
    structure,
    atom_disp=0.03,
    supercell_matrix_fc3=None,
    supercell_matrix_fc2=None,
    yaml_fname_fc3=None,
    yaml_fname_fc2=None,
    **kwargs,
):
    """
    Modified based on pymatgen.io.phonopy.get_displaced_structures()
    Replaced all phonopy features with phono3py
    supercell_matrix_fc2: for fc2 calculation
    """

    cutoff_pair_distance = kwargs.get("cutoff_pair_distance", None)
    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)
    primitive_matrix = kwargs.get("primitive_matrix", None)
    is_symmetry = kwargs.get("is_symmetry", True)
    symprec = kwargs.get("symprec", 1e-5)

    ph_structure = get_phonopy_structure(structure)

    if supercell_matrix_fc3 is None:
        supercell_matrix_fc3 = np.eye(3) * np.array((2, 2, 2))

    phonon = Phono3py(
        unitcell=ph_structure,
        supercell_matrix=supercell_matrix_fc3,
        phonon_supercell_matrix=supercell_matrix_fc2,
        primitive_matrix=primitive_matrix,
        is_symmetry=is_symmetry,
        symprec=symprec,
    )

    print("get_displaced_structures: symprec = {}".format(phonon.symmetry.tolerance))
    print("get_displaced_structures: is_symmetry = {}".format(phonon._is_symmetry))

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
    disp_supercells_fc3 = phonon.supercells_with_displacements
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

        disp_supercells_fc2 = phonon.phonon_supercells_with_displacements
        for c in disp_supercells_fc2:
            if c is not None:
                structure_list_fc2.append(get_pmg_structure(c))

    return structure_list_fc3, structure_list_fc2


from phonopy.harmonic.force_constants import show_drift_force_constants
from phono3py.phonon3.fc3 import show_drift_fc3
from phono3py import Phono3py
from phono3py.cui.create_force_constants import parse_forces
from phono3py.file_IO import (
    parse_disp_fc2_yaml,
    parse_FORCES_FC2,
    write_fc3_to_hdf5,
    write_fc2_to_hdf5,
)


def run_thermal_conductivity(phono3py, is_LBTE=False, t_min=0, t_max=1001, t_step=10):
    """
    Evaluate lattice thermal conductivity and write kappa-*.hdf5.
    Phono3py object must be instantiated.
    This function reads disp_fcx.yaml and FORCES_FCx.
    """

    # create fc3 and fc2
    if not (phono3py.phonon_supercell_matrix == phono3py.supercell_matrix).all():
        # if fc2 supercell size is specified,
        # prepare dataset of fc3 from disp_fc3.yaml and FORCES_FC3
        # and fc2 from disp_fc2.yaml and FORCES_FC2
        print("INFO: preparing fc3 from disp_fc3.yaml and FORCES_FC3")
        phono3py.dataset = parse_forces(
            phono3py,
            # cutoff_pair_distance=None, # optional
            force_filename="FORCES_FC3",
            disp_filename="disp_fc3.yaml",
            fc_type="fc3",
        )
        print("INFO: preparing fc2 from disp_fc2.yaml and FORCES_FC2")
        phono3py.phonon_dataset = parse_forces(
            phono3py,
            force_filename="FORCES_FC2",
            disp_filename="disp_fc2.yaml",
            fc_type="phonon_fc2",
        )
    else:
        # if fc2 supercell size is NOT specified,
        # prepare dataset of fc3 and fc2 from disp_fc3.yaml and FORCES_FC3
        print("INFO: preparing fc3 and fc2 from disp_fc3.yaml and FORCES_FC3")
        phono3py.dataset = parse_forces(
            phono3py,
            # cutoff_pair_distance=None, # optional
            force_filename="FORCES_FC3",
            disp_filename="disp_fc3.yaml",
            fc_type="fc3",
        )
        phono3py.phonon_dataset = parse_forces(
            phono3py,
            force_filename="FORCES_FC3",
            disp_filename="disp_fc3.yaml",
            fc_type="fc2",
        )
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
        is_LBTE=is_LBTE,
        temperatures=range(t_min, t_max, t_step),
        is_isotope=True,
        boundary_mfp=1e6,  # This is to avoid divergence of phonon life time.
        write_kappa=True,
    )

    # Conductivity_RTA object (https://git.io/vVRUW)
    cond_rta = phono3py.thermal_conductivity

    return cond_rta


from phonopy.cui.phonopy_script import file_exists
from phono3py.file_IO import (
    parse_disp_fc2_yaml,
    parse_FORCES_FC2,
    get_length_of_first_line,
)
from phonopy.file_IO import write_FORCE_SETS


def create_FORCE_SETS_from_FORCES_FCx(
    forces_filename="FORCES_FC3", disp_filename="disp_fc3.yaml", log_level=1
):
    with open(forces_filename, "r") as f:
        len_first_line = get_length_of_first_line(f)

    if len_first_line == 3:
        file_exists(disp_filename, log_level)
        disp_dataset = parse_disp_fc2_yaml(filename=disp_filename)
        file_exists(forces_filename, log_level)
        parse_FORCES_FC2(disp_dataset, filename=forces_filename)
        if log_level:
            print('Displacement dataset was read from "%s".' % disp_filename)
        write_FORCE_SETS(disp_dataset)

        if log_level:
            print("FORCE_SETS has been created.")


from phonopy.file_IO import write_BORN
from atomate.vasp.database import VaspCalcDb


def create_BORN_file_from_tag(tag, db_file):
    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

    # get relaxation task document
    print(f"Extracting doc with keyword: {tag} BORN")
    doc_relaxation = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} BORN"},
        }
    )

    # extract unitcell structure, borns, and epsilon
    unitcell = Structure.from_dict(
        doc_relaxation["calcs_reversed"][0]["output"]["structure"]
    )
    ph_unitcell = get_phonopy_structure(unitcell)
    borns = np.array(doc_relaxation["calcs_reversed"][0]["output"]["outcar"]["born"])
    epsilon = np.array(
        doc_relaxation["calcs_reversed"][0]["output"]["outcar"]["dielectric_tensor"]
    )

    # write BORN file
    print("Writing BORN file")
    write_BORN(ph_unitcell, borns, epsilon)


from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from phonopy import Phonopy
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.phonon.dos import CompletePhononDos, PhononDos


def get_phonon_band_structure_symm_line_ph3pywf(
    phonon: Phonopy,
    has_nac: bool = False,
    filename: str = None,
    line_density: float = 20.0,
    symprec: float = 0.01,
) -> PhononBandStructureSymmLine:
    """
    Get a phonon band structure along a high symmetry path from an initialized phonopy
    instance.
    Args:
        phonon: A Phonopy instance.
        has_nac: specify if the band structure has been produced taking into account
                non-analytical corrections at Gamma. If True frequenciens at Gamma from
                diffent directions will be stored in naf. Default False.
        filename : str, optional
            File name used to write ``band.yaml`` like file.
        line_density: The density along the high symmetry path.
        symprec: Symmetry precision passed to phonopy and used for determining
            the band structure path.
    Returns:
        The line mode band structure.
    """
    structure_phonopy = phonon.unitcell
    structure = get_pmg_structure(structure_phonopy)
    kpath = HighSymmKpath(structure, symprec=symprec)

    kpoints, labels = kpath.get_kpoints(
        line_density=line_density, coords_are_cartesian=False
    )

    phonon.run_band_structure([kpoints])

    if filename is not None:
        phonon.write_yaml_band_structure(filename=filename)

    frequencies = np.array(phonon.band_structure.frequencies)
    frequencies = frequencies[0, :, :].T

    labels_dict = {a: k for a, k in zip(labels, kpoints) if a != ""}

    return PhononBandStructureSymmLine(
        kpoints,
        frequencies,
        structure.lattice,
        has_nac=has_nac,
        labels_dict=labels_dict,
        coords_are_cartesian=True,
    )


def get_phonon_dos_ph3pywf(
    phonon: Phonopy,
    mesh_density: float = 100.0,
    num_dos_steps: int = 200,
) -> CompletePhononDos:
    """
    Get a projected phonon density of states from an initialized phonopy instance.
    Args:
        phonon: A Phonopy instance.
        mesh_density: The density of the q-point mesh. See the docstring
            for the ``mesh`` argument in Phonopy.init_mesh() for more details.
        num_dos_steps: Number of frequency steps in the energy grid.
    Returns:
        The density of states.
    """
    structure_phonopy = phonon.unitcell
    structure = get_pmg_structure(structure_phonopy)
    phonon.run_mesh(
        mesh_density,
        is_mesh_symmetry=False,
        with_eigenvectors=True,
        is_gamma_center=True,
    )

    # get min, max, step frequency
    frequencies = phonon.get_mesh_dict()["frequencies"]
    freq_min = frequencies.min()
    freq_max = frequencies.max()
    freq_pitch = (freq_max - freq_min) / num_dos_steps

    phonon.run_projected_dos(
        freq_min=freq_min, freq_max=freq_max, freq_pitch=freq_pitch
    )

    dos_raw = phonon.projected_dos.get_partial_dos()
    pdoss = dict(zip(structure, dos_raw[1]))

    total_dos = PhononDos(dos_raw[0], dos_raw[1].sum(axis=0))
    return CompletePhononDos(structure, total_dos, pdoss)


import yaml


def write_yaml_from_dict(d, filename):
    """
    Write dict to .yaml file.
    d (dict)
    filename (str)
    """
    with open(filename, "w") as outfile:
        yaml.dump(d, outfile, default_flow_style=False)


# DATABASE HELPERS #

# TODO: move these to "db"
def insert_gridfs_file(filename, db, collection):
    """
    Insert the file with given filename into GridFS.
    Args:
        filename (str): the path to file
        db (pymongo.database.Database): usually mmdb.db
        collection (str): the GridFS collection name
    Returns:
        dictionary of gridfs info for retrieval
    """
    # standard gridfs preparation
    fs = gridfs.GridFS(db, collection=collection)

    # read file
    with open(filename, "rb") as fh:
        contents = fh.read()

    # upload
    fs_id = fs.put(contents)

    # prepare dict to return
    fs_info = {
        "fs_id": fs_id,
        "type": "file",
        "filename": os.path.basename(filename),
        "collection": collection,
    }

    return fs_info


def insert_gridfs_dict(d, db, collection):
    """
    Insert the given document into GridFS.
    Args:
        d (dict): the document
        db (pymongo.database.Database): usually mmdb.db
        collection (str): the GridFS collection name
    Returns:
        dictionary of gridfs info for retrieval
    """
    # standard gridfs preparation
    fs = gridfs.GridFS(db, collection=collection)

    # always perform the string conversion when inserting directly to gridfs
    d = json.dumps(d, cls=MontyEncoder)
    d = zlib.compress(d.encode(), True)

    # upload
    fs_id = fs.put(d)

    # prepare dict to return
    fs_info = {
        "fs_id": fs_id,
        "type": "dict",
        "collection": collection,
    }

    return fs_info


def get_file_from_gridfs(fs_id, db, collection, filename):
    # standard gridfs preparation
    fs = gridfs.GridFS(db, collection=collection)

    # retrieve file contents from gridfs
    contents = fs.get(fs_id).read()

    # write file
    with open(filename, "wb") as fh:
        fh.write(contents)
    return


def get_dict_from_gridfs(fs_id, db, collection):
    # standard gridfs preparation
    fs = gridfs.GridFS(db, collection=collection)

    # retrieve dict, and decompress
    bs_json = zlib.decompress(fs.get(fs_id).read())
    obj_dict = json.loads(bs_json.decode())

    return obj_dict


def get_from_gridfs(db, task_doc, key, filename=None):
    """
    Retrieve the object (or write the file) of type key associated with given task document.
    Args:
        db (pymongo.database.Database): usually mmdb.db
        task_doc (dict): document in db collection
        key (str): key to find gridfs info from task_doc
        filename (str): if specified, override the original filename when write,
            otherwise the original filename is used
    Returns:
        The dictionary stored or nothing for file object
    """
    # get gridfs info for retrieval
    fs_info = task_doc[key]
    collection = fs_info["collection"]
    fs_id = fs_info["fs_id"]

    # retrieve file
    if fs_info["type"] == "file":
        filename = filename or fs_info["filename"]
        return get_file_from_gridfs(fs_id, db, collection, filename)

    # retrieve dict
    if fs_info["type"] == "dict":
        return get_dict_from_gridfs(fs_id, db, collection)


# POST ANALYSIS HELPER #

from pymatgen.io.vasp import Kpoints

# TODO: move these to "post_analysis"
def convtest_preview_mesh(db_file_local, tag, mesh_densities=None):
    mesh_densities = mesh_densities or [256 * k for k in range(1, 150)]

    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file_local, admin=True)

    # read fc3 doc from DB for structure
    doc_fc3 = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} disp_fc3"},
        }
    )

    # get supercell structure for fc3
    structure = Structure.from_dict(doc_fc3["input"]["structure"])

    # if mesh_list is not provided
    # set q-point mesh using Kpoints class method
    mesh_list = []
    _tmp = []
    for mesh_density in mesh_densities:
        qpoints = Kpoints.automatic_density_by_vol(structure, mesh_density)
        #     print(f"{mesh_density = }") # FOR TESTING
        #     print(f"{qpoints.kpts[0] = }") # FOR TESTING
        if any(qpt % 2 == 0 for qpt in qpoints.kpts[0]):
            #         print("\tignore even") # FOR TESTING
            _tmp.append(mesh_density)
            continue
        if len(mesh_list) and qpoints.kpts[0] == mesh_list[-1]:
            #         print("\tignore repeated") # FOR TESTING
            _tmp.append(mesh_density)
            continue

        mesh_list.append(qpoints.kpts[0])

    print("Generated mesh list:")  # FOR TESTING
    for mesh in mesh_list:  # FOR TESTING
        print(f"\t{mesh}")  # FOR TESTING

    for m in _tmp:
        mesh_densities.remove(m)

    print("Reduced mesh densities:")  # FOR TESTING
    for mesh_d in mesh_densities:  # FOR TESTING
        print(f"\t{mesh_d}")  # FOR TESTING


# OTHER HELPER FUNCTIONS #

# TODO: move this to "mission_control"?
def check_time_of_tasks(db_file_local, tag):

    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file_local, admin=True)

    docs_p = mmdb.collection.find(
        {
            "task_label": {"$regex": f"{tag}"},
        }
    )

    docs = []
    for p in docs_p:
        docs.append(p)
    # print(type(docs[0]))
    for d in docs:
        print('{{task_label:"{}"}}'.format(d["task_label"]))
        if "nsites" in d:
            print("nsites = {}".format(d["nsites"]))
        if "input" in d:
            print("EDIFF = {}".format(d["input"]["incar"]["EDIFF"]))
        if "run_stats" in d:
            if "standard" in d["run_stats"]:
                print(
                    "Maximum memory used (kb) = {}".format(
                        d["run_stats"]["standard"]["Maximum memory used (kb)"]
                    )
                )
            print(
                "Elapsed time (sec) = {}".format(
                    d["run_stats"]["overall"]["Elapsed time (sec)"]
                )
            )
        print("")