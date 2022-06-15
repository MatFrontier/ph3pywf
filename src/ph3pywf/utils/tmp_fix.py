from ph3pywf.utils.ph3py import (
    create_BORN_file_from_tag,
    get_phonon_band_structure_symm_line_ph3pywf,
)
import numpy as np
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from pymatgen.io.phonopy import get_phonopy_structure
from phonopy import Phonopy
from datetime import datetime
from ph3pywf.utils.ph3py import get_from_gridfs
from phonopy.file_IO import parse_FORCE_CONSTANTS
import phonopy.cui.load_helper as load_helper
from phonopy.interface.calculator import get_default_physical_units


def update_bs(tag, db_file):
    born_filename = None

    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

    # read addertask_dict from DB
    addertask_dict = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} DisplacedStructuresAdderTask"},
        }
    )

    # get user settings from the adder task
    supercell_size_fc3 = addertask_dict["user_settings"].get("supercell_size_fc3", None)
    supercell_size_fc2 = addertask_dict["user_settings"].get("supercell_size_fc2", None)
    primitive_matrix = addertask_dict["user_settings"].get("primitive_matrix", None)
    is_nac = addertask_dict["user_settings"].get("is_nac", False)
    is_symmetry = addertask_dict["user_settings"].get("is_symmetry", True)
    symprec = addertask_dict["user_settings"].get("symprec", 1e-5)

    # if non-analytical term correction is on, generate BORN file
    if is_nac and born_filename is None:
        create_BORN_file_from_tag(tag, db_file)
        born_filename = "BORN"

    # get supercell_matrices
    supercell_matrix_fc3 = np.eye(3) * np.array(supercell_size_fc3)
    if supercell_size_fc2 is not None:
        supercell_matrix_fc2 = np.eye(3) * np.array(supercell_size_fc2)
    else:
        supercell_matrix_fc2 = supercell_matrix_fc3

    # get relaxation task document
    doc_relaxation = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} structure optimization"},
        }
    )

    # get optimized unitcell structure from relaxation output
    # if relaxation doc not found, get structure from Addertask input
    if doc_relaxation is None:
        doc_adderfw = mmdb.db["fireworks"].find_one(
            {
                "name": {"$regex": f"{tag} DisplacedStructuresAdderTask"},
            }
        )
        unitcell = Structure.from_dict(
            doc_adderfw["spec"]["_tasks"][0]["struct_unitcell"]
        )
    else:
        unitcell = Structure.from_dict(
            doc_relaxation["calcs_reversed"][0]["output"]["structure"]
        )

    ph_unitcell = get_phonopy_structure(unitcell)

    # get FORCE_CONSTANTS
    ph3py_entry = mmdb.db["ph3py_tasks"].find_one({"task_label": tag})
    get_from_gridfs(mmdb.db, ph3py_entry, "FORCE_CONSTANTS")

    # prepare Phonopy instance
    phonon = Phonopy(
        ph_unitcell,
        supercell_matrix=supercell_matrix_fc2,
        primitive_matrix=primitive_matrix,
        is_symmetry=is_symmetry,
        symprec=symprec,
    )
    if is_nac:
        units = get_default_physical_units(None)
        phonon.nac_params = load_helper.get_nac_params(
            primitive=phonon.primitive,
            born_filename=born_filename,
            is_nac=is_nac,
            nac_factor=units["nac_factor"],
        )
    phonon.set_force_constants(parse_FORCE_CONSTANTS())

    bs = get_phonon_band_structure_symm_line_ph3pywf(
        phonon,
        has_nac=is_nac,
        filename="band.yaml",
        line_density=80,
    )

    mmdb.db["ph3py_tasks"].update_one(
        {"task_label": tag},
        {
            "$set": {
                "band_structure": bs.as_dict(),
                "last_updated": datetime.utcnow(),
            }
        },
    )
