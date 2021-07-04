# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.utils.ph3py import get_displaced_structures, run_thermal_conductivity, create_FORCE_SETS_from_FORCES_FC3
from fireworks import explicit_serialize
from atomate.utils.utils import env_chk
import numpy as np
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
import yaml
import os
import sys
from atomate.utils.utils import get_logger
from phono3py import Phono3py
from pymatgen.io.phonopy import get_phonopy_structure, get_phonon_band_structure_symm_line_from_fc, get_phonon_dos_from_fc
import h5py
from phono3py.file_IO import parse_disp_fc3_yaml, write_FORCES_FC3
from phonopy.file_IO import parse_FORCE_CONSTANTS, write_FORCE_CONSTANTS
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
        supercell_size (tuple): (2, 2, 2) by default. 
        cutoff_pair_distance (float): set to reduce the number of supercells with displacements to be calculated.
        atom_disp (float): atomic displacement. Default is 0.01 Angstrom.
        struct_unitcell (Structure): optimized unitcell structure. 
        vis_static (VaspInputSet): input set for static VASP jobs.
        is_reduced_test (bool): if set to True, there'll be only a few staticFW generated.
    """
    required_params = ["tag", "db_file"]
    optional_params = ["supercell_size", 
                       "cutoff_pair_distance", 
                       "atom_disp", 
                       "struct_unitcell", 
                       "vis_static", 
                       "is_reduced_test"]
        
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        supercell_size = self.get("supercell_size", (2,2,2))
        cutoff_pair_distance = self.get("cutoff_pair_distance", None)
        atom_disp = self.get("atom_disp", 0.03)
        struct_unitcell = self.get("struct_unitcell", None)
        vis_static = self.get("vis_static", None)
#         name = self.get("name", "DisplacedStructuresAdderTask")
        is_reduced_test = self.get("is_reduced_test", False)
        
        logger.info("Adder: DEBUG VER 05/27 10:52")
        
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
        supercell_matrix = np.eye(3) * np.array(supercell_size) 
        struct_displaced = get_displaced_structures(
            structure=struct_unitcell,
            atom_disp=atom_disp,
            supercell_matrix=supercell_matrix,
            yaml_fname="disp_fc3.yaml",
            cutoff_pair_distance=cutoff_pair_distance,
        )
        
        # save disp_fc3.yaml to DB collection
        phonopy_disp_dict = {}
        with open("disp_fc3.yaml", "r") as fh:
            phonopy_disp_dict["yaml"] = yaml.load(fh, Loader=yaml.SafeLoader)
        
        calc_dir = os.getcwd()
        fullpath = os.path.abspath(calc_dir)
        phonopy_disp_dict["dir_name"] = fullpath
        phonopy_disp_dict["last_updated"] = datetime.utcnow()
        phonopy_disp_dict["task_label"] = f"{tag} DisplacedStructuresAdderTask"
        mmdb.insert_task(phonopy_disp_dict)
        
        # initialize new fws
        new_fws = []
        
        # update vis_static preparation
        vis_dict = vis_static.as_dict()
        
        # for each structure in the displaced structures append a StaticFW
        for i, structure in enumerate(struct_displaced):
            if i==0: 
                continue # Skip undeformed supercell
            if is_reduced_test and i==5: # For dynamic wf testing
                logger.info("Adder: Stop FW generation for dynamic WF testing")
                break # For dynamic wf testing
            disp_id = f"{i:05d}"
#             logger.info(f"Adder: Before update: vis.structure.num_sites={vis_static.structure.num_sites}")
            vis_dict["structure"] = structure.as_dict() # update vis_static
            vis_static = vis_static.from_dict(vis_dict) # update vis_static
            logger.info(f"Adder: After update: vis.structure.num_sites={vis_static.structure.num_sites}")
            fw = StaticFW(structure=structure,
                          vasp_input_set=vis_static, 
                          name=f"{tag} disp-{disp_id}",
                          prev_calc_loc=None,
                          prev_calc_dir=None,
                         )
            logger.info(f"Adder: adding FW {disp_id}")
            new_fws.append(fw)
        
        # return additions of combined FWs
        if len(new_fws) != 0:
            logger.info("Adder: returning FWAction")
            return FWAction(detours=new_fws)
        

@explicit_serialize
class Phono3pyAnalysisToDb(FiretaskBase):
    """
    Compute thermal conductivity?
    
    Required params: 
        tag (str): unique label to identify contents related to this WF.
        db_file (str): path to file containing the database credentials. Supports env_chk.
        
    Optional params: 
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        t_step (float): temperature step (in K)
        supercell_size (tuple)
        primitive_matrix (ndarray): transformation matrix to primitive cell from unit cell.
            Primitive matrix with respect to unit cell.
            shape=(3, 3), dtype='double', order='C'
        mesh (list): sampling mesh numbers in reciprocal space.
        metadata (dict): meta data.
        user_settings (dict): c
    
    """
    required_params = ["tag", "db_file"]
    optional_params = ["t_min",
                       "t_max",
                       "t_step",
                       "supercell_size", 
                       "primitive_matrix",
                       "mesh",
                       "metadata",
                       "user_settings"]
        
    def run_task(self, fw_spec):
        # initialize doc
        ph3py_dict = {}
        
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        t_min = self.get("t_min", 0)
        t_max = self.get("t_max", 1001)
        t_step = self.get("t_step", 10)
        supercell_size = self.get("supercell_size", (2,2,2))
        primitive_matrix = self.get("primitive_matrix", None)
        mesh = self.get("mesh", [11, 11, 11])
        ph3py_dict["metadata"] = self.get("metadata", {})
        ph3py_dict["user_settings"] = self.get("user_settings", {})
        ph3py_dict["task_label"] = tag
        
        # get force_sets from the disp-* runs in DB
        force_sets = []
        logger.info("PostAnalysis: Extracting docs from DB")
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        docs_p = mmdb.collection.find(
            {
                "task_label": {"$regex": f"{tag} disp*"},
            }
        )
        docs_disp = []
        for p in docs_p:
            docs_disp.append(p)

        docs_disp = sorted(docs_disp, key = lambda i: i["task_label"])
        logger.info("PostAnalysis: Read {} docs".format(len(docs_disp)))
        logger.info("PostAnalysis: Generating force sets")
        for d in docs_disp:
            forces = np.array(d["output"]["forces"])
            force_sets.append(forces)

        # get disp_fc3.yaml from DB
        # and get disp_dataset from disp_fc3.yaml
        phonopy_disp_dict = mmdb.collection.find_one(
            {
                "task_label": {"$regex": f"{tag} DisplacedStructuresAdderTask"},
            }
        )
        logger.info("PostAnalysis: Writing disp_fc3.yaml")
        with open("disp_fc3.yaml", "w") as outfile:
            yaml.dump(phonopy_disp_dict["yaml"], outfile, default_flow_style=False)

        disp_dataset = parse_disp_fc3_yaml(filename="disp_fc3.yaml")

        # generate FORCES_FC3
        logger.info("PostAnalysis: Writing FORCES_FC3")
        write_FORCES_FC3(disp_dataset, force_sets, filename="FORCES_FC3")

        # prepare Phono3py object
        doc_opt = mmdb.collection.find_one(
            {
                "task_label": {"$regex": f"{tag} structure optimization"},
            }
        )

        unitcell = Structure.from_dict(doc_opt["calcs_reversed"][0]["output"]["structure"])
        ph_unitcell = get_phonopy_structure(unitcell)
        supercell_matrix = np.eye(3) * np.array(supercell_size) 
        phono3py = Phono3py(unitcell=ph_unitcell,
                            supercell_matrix=supercell_matrix,
                            primitive_matrix=primitive_matrix,
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
        logger.info("PostAnalysis: Creating FORCE_SETS")
        create_FORCE_SETS_from_FORCES_FC3(forces_filename="FORCES_FC3", 
                                          disp_filename="disp_fc3.yaml",
                                         )

        # create FORCE_CONSTANTS
        logger.info("PostAnalysis: Creating FORCE_CONSTANTS")
        phonon = phonopy.load(supercell_matrix=supercell_matrix,
                              primitive_matrix=primitive_matrix,
                              unitcell=ph_unitcell,
                              force_sets_filename="FORCE_SETS")
        write_FORCE_CONSTANTS(phonon.get_force_constants(),
                              filename="FORCE_CONSTANTS")

        # save phonon dispersion band structure object
        logger.info("PostAnalysis: Evaluating phonon dispersion band structure")
        force_constants = parse_FORCE_CONSTANTS()
        bs = get_phonon_band_structure_symm_line_from_fc(unitcell, 
                                                         supercell_matrix, 
                                                         force_constants, 
                                                         primitive_matrix=primitive_matrix)
        ph3py_dict["band_structure"] = bs.as_dict()
        with open("band_structure.yaml", "w") as outfile: # FOR TESTING
            yaml.dump(ph3py_dict["band_structure"], outfile, default_flow_style=False) # FOR TESTING
#         plotter = PhononBSPlotter(bs)
#         plotter.save_plot("plot.png","png")
        
        # parse phonon DOS
        logger.info("PostAnalysis: Parsing phonon DOS")
        dos = get_phonon_dos_from_fc(unitcell, 
                                     supercell_matrix, 
                                     force_constants, 
                                     primitive_matrix=primitive_matrix)
        
        ph3py_dict["dos"] = dos.as_dict()
        with open("dos.yaml", "w") as outfile: # FOR TESTING
            yaml.dump(ph3py_dict["dos"], outfile, default_flow_style=False) # FOR TESTING

        # parse kappa-*.hdf5
        logger.info("PostAnalysis: Parsing kappa-*.hdf5")
        f = h5py.File("kappa-m{}{}{}.hdf5".format(*mesh))
        for item in list(f):
            logger.info("PostAnalysis: Reading property: {}".format(item))
            if item == "kappa_unit_conversion":
                ph3py_dict[item] = f[item][()]
                continue
            if item == "mode_kappa": # skip mode_kappa to reduce size of document
                logger.info("PostAnalysis: Skipping property: mode_kappa")
                continue
            ph3py_dict[item] = f[item][:].tolist()
        
        
        # add more informations in ph3py_dict
        ph3py_dict["structure"] = unitcell.as_dict()
        ph3py_dict["formula_pretty"] = unitcell.composition.reduced_formula
        ph3py_dict["success"] = True
        ph3py_dict["supercell_size"] = supercell_size
        
        calc_dir = os.getcwd()
        fullpath = os.path.abspath(calc_dir)
        ph3py_dict["dir_name"] = fullpath

        # store results in ph3py_tasks collection
        coll = mmdb.db["ph3py_tasks"]
        # if coll.find_one({"task_label": {"$regex": f"{tag}"}}) is not None:
        ph3py_dict["last_updated"] = datetime.utcnow()
        
        with open("ph3py_dict.yaml", "w") as outfile: # FOR TESTING
            yaml.dump(ph3py_dict, outfile, default_flow_style=False) # FOR TESTING
            
        coll.update_one(
            {"task_label": {"$regex": f"{tag}"}}, {"$set": ph3py_dict}, upsert=True
        )
        
        
        
        
    
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
        