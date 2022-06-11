from ph3pywf.utils.ph3py import (
    create_FORCE_SETS_from_FORCES_FCx,
    create_BORN_file_from_tag,
    get_phonon_band_structure_symm_line_ph3pywf,
    write_yaml_from_dict,
)
from phono3py.file_IO import (
    parse_disp_fc2_yaml,
    parse_disp_fc3_yaml,
    write_FORCES_FC2,
    write_FORCES_FC3,
)
import numpy as np
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from pymatgen.io.phonopy import get_phonopy_structure
import phonopy
from datetime import datetime


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

    # get force_sets from the disp_fc3-* runs in DB
    force_sets_fc3 = []
    docs_p_fc3 = mmdb.collection.find(
        {
            "task_label": {"$regex": f"{tag} disp_fc3-*"},
        }
    )
    docs_disp_fc3 = []
    for p in docs_p_fc3:
        docs_disp_fc3.append(p)

    docs_disp_fc3 = sorted(docs_disp_fc3, key=lambda i: i["task_label"])
    for d in docs_disp_fc3:
        forces = np.array(d["output"]["forces"])
        force_sets_fc3.append(forces)

    # get force_sets_fc2 from the disp_fc2-* runs in DB
    if supercell_size_fc2 is not None:
        force_sets_fc2 = []
        docs_p_fc2 = mmdb.collection.find(
            {
                "task_label": {"$regex": f"{tag} disp_fc2-*"},
            }
        )
        docs_disp_fc2 = []
        for p in docs_p_fc2:
            docs_disp_fc2.append(p)

        docs_disp_fc2 = sorted(docs_disp_fc2, key=lambda i: i["task_label"])
        for d in docs_disp_fc2:
            forces = np.array(d["output"]["forces"])
            force_sets_fc2.append(forces)

    # get disp_fc3.yaml from DB
    # and get disp_dataset_fc3 from disp_fc3.yaml
    write_yaml_from_dict(addertask_dict["yaml_fc3"], "disp_fc3.yaml")

    disp_dataset_fc3 = parse_disp_fc3_yaml(filename="disp_fc3.yaml")

    # generate FORCES_FC3
    write_FORCES_FC3(disp_dataset_fc3, force_sets_fc3, filename="FORCES_FC3")

    # get disp_fc2.yaml from DB
    # and get disp_dataset_fc2 from disp_fc2.yaml
    if supercell_size_fc2 is not None:
        write_yaml_from_dict(addertask_dict["yaml_fc2"], "disp_fc2.yaml")

        disp_dataset_fc2 = parse_disp_fc2_yaml(filename="disp_fc2.yaml")

        # generate FORCES_FC2
        write_FORCES_FC2(disp_dataset_fc2, force_sets_fc2, filename="FORCES_FC2")

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

    # create phonopy FORCE_SETS
    if supercell_size_fc2 is not None:
        # create FORCE_SETS from fc2
        create_FORCE_SETS_from_FORCES_FCx(
            forces_filename="FORCES_FC2",
            disp_filename="disp_fc2.yaml",
        )
    else:
        # create FORCE_SETS from fc3
        create_FORCE_SETS_from_FORCES_FCx(
            forces_filename="FORCES_FC3",
            disp_filename="disp_fc3.yaml",
        )

    # prepare Phonopy instance
    phonon = phonopy.load(
        supercell_matrix=supercell_matrix_fc2,
        primitive_matrix=primitive_matrix,
        is_nac=is_nac,
        unitcell=ph_unitcell,
        is_symmetry=is_symmetry,
        symprec=symprec,
        force_sets_filename="FORCE_SETS",
        born_filename=born_filename if is_nac else None,
    )

    bs = get_phonon_band_structure_symm_line_ph3pywf(
        phonon, has_nac=is_nac, filename="band.yaml"
    )

    mmdb.db["ph3py_tasks"].update_one(
        {"task_label": "2022-02-23-17-49-45-226766"},
        {
            "$set": {
                "band_structure": bs.as_dict(),
                "last_updated": datetime.utcnow(),
            }
        },
    )

    return bs
