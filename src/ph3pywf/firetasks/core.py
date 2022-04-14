# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, FiretaskBase
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.utils.ph3py import (
    get_displaced_structures, 
    run_thermal_conductivity, 
    create_FORCE_SETS_from_FORCES_FCx, 
    create_BORN_file_from_tag, 
    get_phonon_band_structure_symm_line_ph3pywf,
    get_phonon_dos_ph3pywf,
    write_yaml_from_dict,
    insert_gridfs_file,
    insert_gridfs_dict,
)
from fireworks import explicit_serialize
from atomate.utils.utils import env_chk
import numpy as np
from atomate.vasp.fireworks.core import StaticFW
from pymatgen.io.vasp import Kpoints
import yaml
import os
import sys
from atomate.utils.utils import get_logger
from phono3py import Phono3py
from pymatgen.io.phonopy import get_phonopy_structure
import h5py
from phono3py.file_IO import parse_disp_fc2_yaml, parse_disp_fc3_yaml, write_FORCES_FC2, write_FORCES_FC3
from phonopy.file_IO import write_FORCE_CONSTANTS
import phonopy

logger = get_logger(__name__)

@explicit_serialize
class DisplacedStructuresAdderTask(FiretaskBase):
    """
    Intermediate Firetask that reads optimized structure, 
    generate displaced structures, 
    and dynamically add StaticFWs to WF.
    
    Required params: 
        tag (str): unique label to identify contents related to this WF.
        db_file (str): path to file containing the database credentials. Supports env_chk.
        
    Optional params: 
        supercell_size_fc3 (tuple): Supercell dimension for 3rd order force constants. (2, 2, 2) by default. 
        supercell_size_fc2 (tuple): Optional supercell dimension for 2nd order force constants. 
        cutoff_pair_distance (float): set to reduce the number of supercells with displacements to be calculated.
        atom_disp (float): atomic displacement. Default is 0.03 Angstrom.
        struct_unitcell (Structure): optimized unitcell structure. 
        vis_static (VaspInputSet): input set for static VASP jobs.
        user_settings (dict): c
        is_nac (bool): If True, non-analytical term correction is on.
    """
    required_params = ["tag", "db_file"]
    optional_params = ["supercell_size_fc3",
                       "supercell_size_fc2",
                       "cutoff_pair_distance", 
                       "atom_disp", 
                       "struct_unitcell", 
                       "vis_static", 
                       "primitive_matrix",
                       "user_settings",
                       "is_nac",
                       "is_symmetry",
                       "symprec"]
        
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        supercell_size_fc3 = self.get("supercell_size_fc3", (2,2,2))
        supercell_size_fc2 = self.get("supercell_size_fc2", None)
        cutoff_pair_distance = self.get("cutoff_pair_distance", None)
        atom_disp = self.get("atom_disp", 0.03)
        struct_unitcell = self.get("struct_unitcell", None)
        vis_static = self.get("vis_static", None)
        primitive_matrix = self.get("primitive_matrix", None)
        is_nac = self.get("is_nac", False)
        is_symmetry = self.get("is_symmetry", True)
        symprec = self.get("symprec", 1e-5)
        
        logger.info("Adder: DEBUG VER 05/27 10:52")
        
        # initialize document to be saved in DB
        addertask_dict = {}
        addertask_dict["user_settings"] = self.get("user_settings", {})
        
        # read optimized structures
        mmdb = VaspCalcDb.from_db_file(db_file)
        if struct_unitcell is None:
            doc = mmdb.collection.find_one(
                {
                    "task_label": {"$regex": f"{tag} structure optimization"},
                },
#                 {"calcs_reversed": 1}, # appears to be "projection", not sure if necessary
            )
            logger.info("Adder: Found doc task_id = {}".format(doc["task_id"]))
            struct_unitcell = Structure.from_dict(doc["calcs_reversed"][0]["output"]["structure"])
        
        # generate displaced structures from optimized structure
        supercell_matrix_fc3 = np.eye(3) * np.array(supercell_size_fc3)
        if supercell_size_fc2 is not None:
            supercell_matrix_fc2 = np.eye(3) * np.array(supercell_size_fc2)
        else:
            supercell_matrix_fc2 = supercell_matrix_fc3
        
        struct_displaced_fc3, struct_displaced_fc2 = get_displaced_structures(
            structure=struct_unitcell,
            atom_disp=atom_disp,
            supercell_matrix_fc3=supercell_matrix_fc3,
            supercell_matrix_fc2=supercell_matrix_fc2,
            yaml_fname_fc3="disp_fc3.yaml",
            yaml_fname_fc2="disp_fc2.yaml",
            cutoff_pair_distance=cutoff_pair_distance,
            primitive_matrix=primitive_matrix,
            is_symmetry=is_symmetry,
            symprec=symprec,
        )
        
        # save disp_fc3.yaml in DB collection
        with open("disp_fc3.yaml", "r") as fh:
            addertask_dict["yaml_fc3"] = yaml.load(fh, Loader=yaml.SafeLoader)
        
        # save disp_fc2.yaml in DB collection
        if supercell_size_fc2 is not None:
            with open("disp_fc2.yaml", "r") as fh:
                addertask_dict["yaml_fc2"] = yaml.load(fh, Loader=yaml.SafeLoader)
        
        calc_dir = os.getcwd()
        fullpath = os.path.abspath(calc_dir)
        addertask_dict["dir_name"] = fullpath
        addertask_dict["last_updated"] = datetime.utcnow()
        addertask_dict["task_label"] = f"{tag} DisplacedStructuresAdderTask"
        mmdb.insert_task(addertask_dict)
        
        # initialize new fws
        new_fws = []
        
        # update vis_static preparation
        vis_dict = vis_static.as_dict()
        
        # append StaticFW for BORN file generation
        if is_nac:
            vis_dict_born = vis_dict.copy()
            vis_dict_born["structure"] = struct_unitcell.as_dict() # update vis_static_born
            vis_dict_born["lepsilon"] = True # update vis_static_born
            vis_static_born = vis_static.__class__.from_dict(vis_dict_born) # update vis_static_born
            fw = StaticFW(structure=struct_unitcell,
                          vasp_input_set=vis_static_born, 
                          name=f"{tag} BORN",
                          prev_calc_loc=None,
                          prev_calc_dir=None,
                         )
            logger.info("Adder: adding BORN FW")
            new_fws.append(fw)
        
        # for each structure in the displaced structures append a StaticFW
        for i, structure in enumerate(struct_displaced_fc3):
            if i==0: 
                continue # Skip undeformed supercell
            disp_id = f"{i:05d}"
#             logger.info(f"Adder: Before update: vis.structure.num_sites={vis_static.structure.num_sites}")
            vis_dict["structure"] = structure.as_dict() # update vis_static
            vis_static = vis_static.__class__.from_dict(vis_dict) # update vis_static
#             logger.info(f"Adder: After update: vis.structure.num_sites={vis_static.structure.num_sites}")
            fw = StaticFW(structure=structure,
                          vasp_input_set=vis_static, 
                          name=f"{tag} disp_fc3-{disp_id}",
                          prev_calc_loc=None,
                          prev_calc_dir=None,
                         )
            logger.info(f"Adder: adding fc3 FW {disp_id}")
            new_fws.append(fw)
        
        # append StaticFWs for fc2
        if supercell_size_fc2 is not None:
            for i, structure in enumerate(struct_displaced_fc2):
                if i==0: 
                    continue # Skip undeformed supercell
                disp_id = f"{i:05d}"
                vis_dict["structure"] = structure.as_dict() # update vis_static
                vis_static = vis_static.__class__.from_dict(vis_dict) # update vis_static
#                 logger.info(f"Adder: After update: vis.structure.num_sites={vis_static.structure.num_sites}")
                fw = StaticFW(structure=structure,
                              vasp_input_set=vis_static, 
                              name=f"{tag} disp_fc2-{disp_id}",
                              prev_calc_loc=None,
                              prev_calc_dir=None,
                             )
                logger.info(f"Adder: adding fc2 FW {disp_id}")
                new_fws.append(fw)
        
        # return additions of combined FWs
        if len(new_fws) != 0:
            logger.info("Adder: returning FWAction")
            return FWAction(detours=new_fws)
        

@explicit_serialize
class Phono3pyAnalysisToDb(FiretaskBase):
    """
    Compute thermal conductivity and store in db.
    
    Required params: 
        tag (str): unique label to identify contents related to this WF.
        db_file (str): path to file containing the database credentials. Supports env_chk.
        
    Optional params: 
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        t_step (float): temperature step (in K)
        mesh (list): sampling mesh numbers in reciprocal space.
        born_filename (str): filename corresponding to "BORN", a file contains non-analytical term correction parameters.
            Specify to use user provided BORN file.
        metadata (dict): meta data.
    
    """
    required_params = ["tag", "db_file"]
    optional_params = ["t_min",
                       "t_max",
                       "t_step",
                       "mesh",
                       "born_filename"]
        
    def run_task(self, fw_spec):
        # initialize doc
        ph3py_dict = {}
        
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        t_min = self.get("t_min", 100)
        t_max = self.get("t_max", 1301)
        t_step = self.get("t_step", 50)
        mesh = self.get("mesh", [11, 11, 11])
        born_filename = self.get("born_filename", None)
        
        ph3py_dict["task_label"] = tag
        
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
        # copy user settings
        ph3py_dict["user_settings"] = addertask_dict["user_settings"]
        # update user settings
        ph3py_dict["user_settings"]["t_min"] = t_min
        ph3py_dict["user_settings"]["t_max"] = t_max
        ph3py_dict["user_settings"]["t_step"] = t_step
        ph3py_dict["user_settings"]["mesh"] = mesh
        ph3py_dict["user_settings"]["born_filename"] = born_filename
        
        # get force_sets from the disp_fc3-* runs in DB
        force_sets_fc3 = []
        logger.info("PostAnalysis: Extracting fc3 docs from DB")
        docs_p_fc3 = mmdb.collection.find(
            {
                "task_label": {"$regex": f"{tag} disp_fc3-*"},
            }
        )
        docs_disp_fc3 = []
        for p in docs_p_fc3:
            docs_disp_fc3.append(p)
        
        docs_disp_fc3 = sorted(docs_disp_fc3, key = lambda i: i["task_label"])
        logger.info("PostAnalysis: Read {} docs".format(len(docs_disp_fc3)))
        logger.info("PostAnalysis: Generating fc3 force sets")
        for d in docs_disp_fc3:
            logger.info("PostAnalysis: Reading doc: {}".format(d["task_label"]))
            forces = np.array(d["output"]["forces"])
            force_sets_fc3.append(forces)
        
        # get force_sets_fc2 from the disp_fc2-* runs in DB
        if supercell_size_fc2 is not None:
            force_sets_fc2 = []
            logger.info("PostAnalysis: Extracting fc2 docs from DB")
            docs_p_fc2 = mmdb.collection.find(
                {
                    "task_label": {"$regex": f"{tag} disp_fc2-*"},
                }
            )
            docs_disp_fc2 = []
            for p in docs_p_fc2:
                docs_disp_fc2.append(p)

            docs_disp_fc2 = sorted(docs_disp_fc2, key = lambda i: i["task_label"])
            logger.info("PostAnalysis: Read {} docs".format(len(docs_disp_fc2)))
            logger.info("PostAnalysis: Generating fc2 force sets")
            for d in docs_disp_fc2:
                logger.info("PostAnalysis: Reading doc: {}".format(d["task_label"]))
                forces = np.array(d["output"]["forces"])
                force_sets_fc2.append(forces)
            
        # get disp_fc3.yaml from DB
        # and get disp_dataset_fc3 from disp_fc3.yaml
        logger.info("PostAnalysis: Writing disp_fc3.yaml")
        write_yaml_from_dict(addertask_dict["yaml_fc3"], "disp_fc3.yaml")

        disp_dataset_fc3 = parse_disp_fc3_yaml(filename="disp_fc3.yaml")

        # generate FORCES_FC3
        logger.info("PostAnalysis: Writing FORCES_FC3")
        write_FORCES_FC3(disp_dataset_fc3, force_sets_fc3, filename="FORCES_FC3")
        
        # get disp_fc2.yaml from DB
        # and get disp_dataset_fc2 from disp_fc2.yaml
        if supercell_size_fc2 is not None:
            logger.info("PostAnalysis: Writing disp_fc2.yaml")
            write_yaml_from_dict(addertask_dict["yaml_fc2"], "disp_fc2.yaml")
            
            disp_dataset_fc2 = parse_disp_fc2_yaml(filename="disp_fc2.yaml")
            
            # generate FORCES_FC2
            logger.info("PostAnalysis: Writing FORCES_FC2")
            write_FORCES_FC2(disp_dataset_fc2, force_sets_fc2, filename="FORCES_FC2")
            
        # get relaxation task document 
        doc_relaxation = mmdb.collection.find_one(
            {
                "task_label": {"$regex": f"{tag} structure optimization"},
            }
        )
        
        # get optimized unitcell structure
        unitcell = Structure.from_dict(doc_relaxation["calcs_reversed"][0]["output"]["structure"])
        ph_unitcell = get_phonopy_structure(unitcell)
        unitcell.to(fmt="poscar", filename="POSCAR-unitcell") # FOR TESTING
        
        # if non-analytical term correction is on, generate BORN file
        if is_nac and born_filename is None:
            logger.info("PostAnalysis: generating BORN file")
            create_BORN_file_from_tag(tag, db_file)
            born_filename = "BORN"
        
        # get supercell_matrices
        supercell_matrix_fc3 = np.eye(3) * np.array(supercell_size_fc3)
        if supercell_size_fc2 is not None:
            supercell_matrix_fc2 = np.eye(3) * np.array(supercell_size_fc2)
        else:
            supercell_matrix_fc2 = supercell_matrix_fc3
        
        # prepare Phono3py instance
        logger.info("PostAnalysis: Preparing Phono3py instance")
        phono3py = Phono3py(unitcell=ph_unitcell,
                            supercell_matrix=supercell_matrix_fc3,
                            phonon_supercell_matrix=supercell_matrix_fc2,
                            primitive_matrix=primitive_matrix,
                            is_symmetry=is_symmetry,
                            symprec=symprec,
                            mesh=mesh,
                            log_level=1, # log_level=0 make phono3py quiet
                           )

        # use run_thermal_conductivity()
        # which will read disp_fc3.yaml and FORCES_FC3
        # this operation will generate file kappa-*.hdf5
        logger.info("PostAnalysis: Evaluating thermal conductivity")
        if "kappa-m{}{}{}.hdf5".format(*mesh) in os.listdir(os.getcwd()):
            logger.info("PostAnalysis: Deleting previous kappa-m{}{}{}.hdf5 file".format(*mesh))
            os.remove("kappa-m{}{}{}.hdf5".format(*mesh))
        run_thermal_conductivity(phono3py, t_min, t_max, t_step)
        
#         # save fc2 and fc3
#         logger.info("PostAnalysis: Saving fc2 and fc3")
#         ph3py_dict["fc3"] = phono3py.fc3.tolist()
#         ph3py_dict["fc2"] = phono3py.fc2.tolist()

        # create phonopy FORCE_SETS
        if supercell_size_fc2 is not None: 
            # create FORCE_SETS from fc2
            logger.info("PostAnalysis: Creating FORCE_SETS from FORCES_FC2")
            create_FORCE_SETS_from_FORCES_FCx(forces_filename="FORCES_FC2", 
                                              disp_filename="disp_fc2.yaml",
                                             )
        else: 
            # create FORCE_SETS from fc3
            logger.info("PostAnalysis: Creating FORCE_SETS from FORCES_FC3")
            create_FORCE_SETS_from_FORCES_FCx(forces_filename="FORCES_FC3", 
                                              disp_filename="disp_fc3.yaml",
                                             )
        
        # prepare Phonopy instance
        logger.info("PostAnalysis: Preparing Phonopy instance")
        logger.info(f"PostAnalysis: is_nac = {is_nac}")
        if is_nac:
            logger.info(f"PostAnalysis: Reading nac_params from file: \"{born_filename}\"")
        phonon = phonopy.load(supercell_matrix=supercell_matrix_fc2,
                              primitive_matrix=primitive_matrix,
                              is_nac=is_nac,
                              unitcell=ph_unitcell,
                              is_symmetry=is_symmetry,
                              symprec=symprec,
                              force_sets_filename="FORCE_SETS",
                              born_filename=born_filename if is_nac else None)
        
        # write FORCE_CONSTANTS
        logger.info("PostAnalysis: Creating FORCE_CONSTANTS")
        write_FORCE_CONSTANTS(phonon.force_constants,
                              filename="FORCE_CONSTANTS")
        
        # evaluate and save phonon dispersion band structure
        logger.info("PostAnalysis: Evaluating phonon dispersion band structure")
        bs = get_phonon_band_structure_symm_line_ph3pywf(phonon,
                                                         has_nac=is_nac,
                                                         filename="band.yaml")
        
        ph3py_dict["band_structure"] = bs.as_dict()
        write_yaml_from_dict(ph3py_dict["band_structure"], "band_structure.yaml") # FOR TESTING
        
        # parse phonon DOS
        logger.info("PostAnalysis: Parsing phonon DOS")
        dos = get_phonon_dos_ph3pywf(phonon)
        
        ph3py_dict["dos"] = dos.as_dict()
        write_yaml_from_dict(ph3py_dict["dos"], "dos.yaml") # FOR TESTING

        # parse kappa-*.hdf5
        logger.info("PostAnalysis: Parsing kappa-*.hdf5")
        f = h5py.File("kappa-m{}{}{}.hdf5".format(*mesh))
        kappa_dict = {}
        for item in list(f):
            logger.info("PostAnalysis: Reading property: {}".format(item))
            if item == "kappa_unit_conversion":
                ph3py_dict[item] = f[item][()]
                kappa_dict[item] = f[item][()]
                continue
            if item in ["temperature", "kappa", "mesh"]:
                ph3py_dict[item] = f[item][:].tolist()
            kappa_dict[item] = f[item][:].tolist()
#             write_yaml_from_dict(ph3py_dict[item], "kappa."+item+".yaml") # FOR TESTING
        
        # insert kappa information to gridfs
        ph3py_dict["kappa_hdf5"] = insert_gridfs_dict(kappa_dict, mmdb.db, "ph3py_tasks_fs")
        
        # insert designated files to gridfs
        filenames = ["fc2.hdf5",
                     "fc3.hdf5",
                     "FORCE_CONSTANTS", 
                    ]
        for filename in filenames:
            ph3py_dict[filename] = insert_gridfs_file(filename, mmdb.db, "ph3py_tasks_fs")
        
        # add more informations in ph3py_dict
        ph3py_dict["structure"] = unitcell.as_dict()
        ph3py_dict["formula_pretty"] = unitcell.composition.reduced_formula
        ph3py_dict["success"] = True
        ph3py_dict["supercell_size_fc3"] = supercell_size_fc3
        ph3py_dict["supercell_size_fc2"] = supercell_size_fc2
        
        calc_dir = os.getcwd()
        fullpath = os.path.abspath(calc_dir)
        ph3py_dict["dir_name"] = fullpath

        # store results in ph3py_tasks collection
        coll = mmdb.db["ph3py_tasks"]
        # if coll.find_one({"task_label": {"$regex": f"{tag}"}}) is not None:
        ph3py_dict["last_updated"] = datetime.utcnow()
        
        write_yaml_from_dict(ph3py_dict, "ph3py_dict.yaml") # FOR TESTING
            
        coll.update_one(
            {"task_label": {"$regex": tag}}, 
            {"$set": ph3py_dict}, 
            upsert=True,
        )
        
        return FWAction()
        
        
    
#########################
# TESTING MODULES BELOW #
#########################

@explicit_serialize
class StoreStructureTask(FiretaskBase):
    required_params = ["tag", "db_file", "structure", "name", "terminate"]
    
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        structure = self.get("structure", None)
        name = self.get("name", "NoName")
        terminate = bool(self.get("terminate", False))
        
        # connect to DB
        mmdb = VaspCalcDb.from_db_file(db_file)
        coll = mmdb.collection
#         coll = mmdb.db.get_collection("_IOtest")
        
        # create structure dict
        struct_dict = {"calcs_reversed":[{"output":{"structure":None}}],"task_label":None}
        struct_dict["task_label"] = f"{tag} {name}"
        struct_dict["calcs_reversed"][0]["output"]["structure"] = structure.as_dict()
        
        # store structure in DB collection
        coll.insert_one(struct_dict)
        
        # terminate FW dependently
        if terminate:
            return FWAction()

@explicit_serialize
class Phono3pyMeshConvergenceToDb(FiretaskBase):
    """
    Compute thermal conductivity using different mesh density and store in db.
    
    Required params: 
        tag (str): unique label to identify contents related to this WF.
        db_file (str): path to file containing the database credentials. Supports env_chk.

    Optional params: 
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        t_step (float): temperature step (in K)
        mesh_densities (list): list of reciprocal densities for q-points.
        mesh_list (list): list of q-points meshes. 
            if provided, will skip automatic generation of meshes using mesh_densities.
        tag_for_copy (str): different label for a copy, if provided, will create a copy in 
            ph3py_tasks_convergence_test instead of overwriting the existing document.
    
    """
    required_params = ["tag", "db_file"]
    optional_params = ["t_min",
                       "t_max",
                       "t_step",
                       "mesh_densities",
                       "mesh_list",
                       "tag_for_copy"]
        
    def run_task(self, fw_spec):
        # initialize doc
        ph3py_dict = {}
        
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        tag_for_copy = self.get("tag_for_copy", None)
        t_min = self.get("t_min", 100)
        t_max = self.get("t_max", 1301)
        t_step = self.get("t_step", 50)
        mesh_densities = self.get("mesh_densities", [256 * k for k in range(1,40)])
        mesh_list = self.get("mesh_list", None)
        
        ph3py_dict["task_label"] = tag_for_copy if tag_for_copy else tag
                
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
        is_symmetry = addertask_dict["user_settings"].get("is_symmetry", True)
        symprec = addertask_dict["user_settings"].get("symprec", 1e-5)
        # copy user settings
        ph3py_dict["user_settings"] = addertask_dict["user_settings"]
        # update user settings
        ph3py_dict["user_settings"]["t_min"] = t_min
        ph3py_dict["user_settings"]["t_max"] = t_max
        ph3py_dict["user_settings"]["t_step"] = t_step
        
        # add more informations in ph3py_dict
        calc_dir = os.getcwd()
        fullpath = os.path.abspath(calc_dir)
        ph3py_dict["dir_name"] = fullpath
        ph3py_dict["success"] = False
        
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
        if not mesh_list:
            mesh_list = []
            _tmp = []
            for mesh_density in mesh_densities:
                qpoints = Kpoints.automatic_density_by_vol(
                    structure, mesh_density
                )
                print(f"{mesh_density = }") # FOR TESTING
                print(f"{qpoints.kpts[0] = }") # FOR TESTING
                if any(qpt%2==0 for qpt in qpoints.kpts[0]):
                    print("\tignore even") # FOR TESTING
                    _tmp.append(mesh_density)
                    continue
                if len(mesh_list) and qpoints.kpts[0] == mesh_list[-1]:
                    print("\tignore repeated") # FOR TESTING
                    _tmp.append(mesh_density)
                    continue

                mesh_list.append(qpoints.kpts[0])
        
            logger.info("Generated mesh list:")
            for mesh in mesh_list:
                print(f"\t{mesh}")

            for m in _tmp:
                mesh_densities.remove(m)
                

        ph3py_dict["user_settings"]["mesh_densities"] = mesh_densities
        ph3py_dict["mesh_list"] = mesh_list
 
        # get force_sets from the disp_fc3-* runs in DB
        force_sets_fc3 = []
        logger.info("PostAnalysis: Extracting fc3 docs from DB")
        docs_p_fc3 = mmdb.collection.find(
            {
                "task_label": {"$regex": f"{tag} disp_fc3-*"},
            }
        )
        docs_disp_fc3 = []
        for p in docs_p_fc3:
            docs_disp_fc3.append(p)
        
        docs_disp_fc3 = sorted(docs_disp_fc3, key = lambda i: i["task_label"])
        logger.info("PostAnalysis: Read {} docs".format(len(docs_disp_fc3)))
        logger.info("PostAnalysis: Generating fc3 force sets")
        for d in docs_disp_fc3:
            logger.info("PostAnalysis: Reading doc: {}".format(d["task_label"]))
            forces = np.array(d["output"]["forces"])
            force_sets_fc3.append(forces)
        
        # get force_sets_fc2 from the disp_fc2-* runs in DB
        if supercell_size_fc2 is not None:
            force_sets_fc2 = []
            logger.info("PostAnalysis: Extracting fc2 docs from DB")
            docs_p_fc2 = mmdb.collection.find(
                {
                    "task_label": {"$regex": f"{tag} disp_fc2-*"},
                }
            )
            docs_disp_fc2 = []
            for p in docs_p_fc2:
                docs_disp_fc2.append(p)

            docs_disp_fc2 = sorted(docs_disp_fc2, key = lambda i: i["task_label"])
            logger.info("PostAnalysis: Read {} docs".format(len(docs_disp_fc2)))
            logger.info("PostAnalysis: Generating fc2 force sets")
            for d in docs_disp_fc2:
                logger.info("PostAnalysis: Reading doc: {}".format(d["task_label"]))
                forces = np.array(d["output"]["forces"])
                force_sets_fc2.append(forces)
            
        # get disp_fc3.yaml from DB
        # and get disp_dataset_fc3 from disp_fc3.yaml
        logger.info("PostAnalysis: Writing disp_fc3.yaml")
        write_yaml_from_dict(addertask_dict["yaml_fc3"], "disp_fc3.yaml")

        disp_dataset_fc3 = parse_disp_fc3_yaml(filename="disp_fc3.yaml")

        # generate FORCES_FC3
        logger.info("PostAnalysis: Writing FORCES_FC3")
        write_FORCES_FC3(disp_dataset_fc3, force_sets_fc3, filename="FORCES_FC3")
        
        # get disp_fc2.yaml from DB
        # and get disp_dataset_fc2 from disp_fc2.yaml
        if supercell_size_fc2 is not None:
            logger.info("PostAnalysis: Writing disp_fc2.yaml")
            write_yaml_from_dict(addertask_dict["yaml_fc2"], "disp_fc2.yaml")
            
            disp_dataset_fc2 = parse_disp_fc2_yaml(filename="disp_fc2.yaml")
            
            # generate FORCES_FC2
            logger.info("PostAnalysis: Writing FORCES_FC2")
            write_FORCES_FC2(disp_dataset_fc2, force_sets_fc2, filename="FORCES_FC2")

        # get relaxation task document 
        doc_relaxation = mmdb.collection.find_one(
            {
                "task_label": {"$regex": f"{tag} structure optimization"},
            }
        )
        
        # get optimized unitcell structure
        unitcell = Structure.from_dict(doc_relaxation["calcs_reversed"][0]["output"]["structure"])
        ph_unitcell = get_phonopy_structure(unitcell)
        unitcell.to(fmt="poscar", filename="POSCAR-unitcell") # FOR TESTING
        
        # get supercell_matrices
        supercell_matrix_fc3 = np.eye(3) * np.array(supercell_size_fc3)
        if supercell_size_fc2 is not None:
            supercell_matrix_fc2 = np.eye(3) * np.array(supercell_size_fc2)
        else:
            supercell_matrix_fc2 = supercell_matrix_fc3
        
        # prepare Phono3py instance
        logger.info("PostAnalysis: Preparing Phono3py instance")
        phono3py = Phono3py(unitcell=ph_unitcell,
                            supercell_matrix=supercell_matrix_fc3,
                            phonon_supercell_matrix=supercell_matrix_fc2,
                            primitive_matrix=primitive_matrix,
                            is_symmetry=is_symmetry,
                            symprec=symprec,
                            mesh=[3,3,3],
                            log_level=1, # log_level=0 make phono3py quiet
                           )
        
        # initialize list for kappa storage
        ph3py_dict["convergence_test"] = []
        ph3py_dict["temperature"] = list(range(t_min, t_max, t_step))
        
        # change to ph3py_tasks_convergence_test collection
        coll = mmdb.db["ph3py_tasks_convergence_test"]
        
        # run thermal conductivity using different mesh numbers
        for mesh in mesh_list:
            # update mesh number
            logger.info(f"PostAnalysis: {mesh = }")
            phono3py.mesh_numbers = mesh
            
            # use run_thermal_conductivity()
            # which will read disp_fc3.yaml and FORCES_FC3
            # this operation will generate file kappa-*.hdf5
            logger.info("PostAnalysis: Evaluating thermal conductivity")
            if "kappa-m{}{}{}.hdf5".format(*mesh) in os.listdir(os.getcwd()):
                logger.info("PostAnalysis: Deleting previous kappa-m{}{}{}.hdf5 file".format(*mesh))
                os.remove("kappa-m{}{}{}.hdf5".format(*mesh))
            run_thermal_conductivity(phono3py, t_min, t_max, t_step)
            
            # parse kappa-*.hdf5
            logger.info("PostAnalysis: Parsing kappa-m{}{}{}.hdf5".format(*mesh))
            f = h5py.File("kappa-m{}{}{}.hdf5".format(*mesh))
            kappa_dict = {"mesh": mesh, "kappa": f["kappa"][:].tolist()}
            ph3py_dict["convergence_test"].append(kappa_dict)
            
            # store result of progress so far (possible to optimize but not critical)
            ph3py_dict["last_updated"] = datetime.utcnow()
            if tag_for_copy:
                coll.update_one(
                    {"task_label": tag_for_copy},
                    {"$set": ph3py_dict},
                    upsert=True,
                )
            else:
                coll.update_one(
                    {"task_label": tag},
                    {"$set": ph3py_dict},
                    upsert=True,
                )
        

        # store results in collection
        ph3py_dict["last_updated"] = datetime.utcnow()
        ph3py_dict["success"] = True
        
        if tag_for_copy:
            coll.update_one(
                {"task_label": tag_for_copy},
                {"$set": ph3py_dict},
                upsert=True,
            )
        else:
            coll.update_one(
                {"task_label": tag},
                {"$set": ph3py_dict},
                upsert=True,
            )
        
        return FWAction()

@explicit_serialize
class Phono3pyEvaluateKappaFromConvTest(FiretaskBase):
    """
    Compute thermal conductivity using different mesh density and store in db.
    
    Required params: 
        tag (str): unique label to identify contents related to this WF.
        db_file (str): path to file containing the database credentials. Supports env_chk.
    
    Optional params: 
        tag_for_copy (str): different label for a copy, if provided, will create a copy in 
            ph3py_tasks_convergence_test instead of overwriting the existing document.
    
    """
    required_params = ["tag", "db_file"]
    optional_params = ["tag_for_copy"]
    
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        tag_for_copy = self.get("tag_for_copy", None)
        
        # define fitting exp function
        def _exp_func(x, kappa_inf, epsilon):
            # kappa_inf is the key, as the value at inf x
            y = kappa_inf * (1 - np.exp(-x / epsilon))
            return y
        
        # connect to DB
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        # read conv_test_doc from DB
        if tag_for_copy:
            conv_test_doc = mmdb.db["ph3py_tasks_convergence_test"].find_one(
                {"task_label": tag_for_copy}
            )
        else:
            conv_test_doc = mmdb.db["ph3py_tasks_convergence_test"].find_one(
                {"task_label": tag}
            )
        
        # get calculation result
        temperature = np.array(conv_test_doc["temperature"])
        mesh_densities = conv_test_doc["user_settings"]["mesh_densities"]
        mesh_list = []
        kappa_list = []
        for each_mesh in conv_test_doc["convergence_test"]:
            mesh_list.append(each_mesh["mesh"])
            kappa_list.append(np.array(each_mesh["kappa"])[:,:])

        # curve fitting
        from scipy.optimize import curve_fit
        bounds = (0, np.inf)
        kappa_fitted = np.zeros((len(temperature), 6))
        for direction in range(3):
            for t in range(len(temperature)):
                logger.info(f"Now fitting direction = {direction}, temperature = {temperature[t]}")
                xdata = np.array(
                    [conv_test_doc["convergence_test"][i]["mesh"][direction] for i in range(len(mesh_list))]
                )
                ydata = np.array(
                    [kappa_list[mesh][t, direction] for mesh in range(len(mesh_list))]
                )
                popt, _ = curve_fit(_exp_func, xdata, ydata, bounds=bounds)
                print(f"kappa_inf = {popt[0]}") # FOR TESTING
                kappa_fitted[t, direction] = popt[0]
        
        for direction in range(3,6):
            for t in range(len(temperature)):
                # directly use the value from max mesh
                kappa_fitted[t, direction] = kappa_list[-1][t, direction]

        # initialize doc
        ph3py_dict = {}
        
        # update evaluated kappa
        ph3py_dict["temperature_fitted"] = temperature.tolist()
        ph3py_dict["kappa_fitted"] = kappa_fitted.tolist()
        ph3py_dict["mesh_list"] = mesh_list
        
        # store results in ph3py_tasks collection
        coll = mmdb.db["ph3py_tasks"]
        ph3py_dict["last_updated"] = datetime.utcnow()
        
        # write_yaml_from_dict(ph3py_dict, "ph3py_dict.yaml") # FOR TESTING
        
        if tag_for_copy:
            coll.update_one(
                {"task_label": tag_for_copy},
                {"$set": ph3py_dict},
                upsert=True,
            )
        else:
            coll.update_one(
                {"task_label": tag},
                {"$set": ph3py_dict},
                upsert=True,
            )
        
        return FWAction()