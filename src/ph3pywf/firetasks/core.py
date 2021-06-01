# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.utils.ph3py import get_displaced_structures, run_thermal_conductivity
from fireworks import explicit_serialize
from atomate.utils.utils import env_chk
import numpy as np
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
import yaml
import os
from atomate.utils.utils import get_logger
from phono3py import Phono3py
from pymatgen.io.phonopy import get_phonopy_structure
import h5py

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
        cutoff_pair_distance (float)
        struct_unitcell (Structure): optimized unitcell structure
        vis_static (VaspInputSet)
        name (str)
        is_reduced_test (bool): if switch on, there'll be only a few staticFW generated
    """
    required_params = ["tag", "db_file"]
    optional_params = ["supercell_size", 
                       "cutoff_pair_distance", 
                       "struct_unitcell", 
                       "vis_static", 
#                        "name", 
                       "is_reduced_test"]
        
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        supercell_size = self.get("supercell_size", (2,2,2))
        cutoff_pair_distance = self.get("cutoff_pair_distance", None)
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
            supercell_matrix=supercell_matrix,
            yaml_fname="phonopy_disp.yaml",
            cutoff_pair_distance=cutoff_pair_distance,
        )
        
        # save phonopy_disp.yaml to DB collection
        phonopy_disp_dict = {}
        with open("phonopy_disp.yaml", "r") as fh:
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
        
    
    """
    required_params = ["tag", "db_file"]
    optional_params = ["t_min",
                       "t_max",
                       "t_step",
                       "supercell_size", 
                       "mesh",
                       "metadata"]
        
    def run_task(self, fw_spec):
        # initialize doc
        ph3py_dict = {}
        
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        t_min = self.get("t_min", 0)
        t_max = self.get("t_max", 1001)
        t_step = self.get("t_step", 10)
        supercell_size = self.get("supercell_size", (2,2,2))
        mesh = self.get("mesh", [11, 11, 11])
        ph3py_dict["metadata"] = self.get("metadata", {})
        ph3py_dict["metadata"].update({"task_label_tag": tag})
        
        # get force_sets from the disp-* runs in DB
        force_sets = []
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        docs_p = mmdb.collection.find(
            {
                "task_label": {"$regex": f"{tag} disp*"},
            }
        )
        docs_disp = []
        for p in docs_p:
            docs_disp.append(p)
        
        for d in docs_disp:
            forces = np.array(d["output"]["forces"])
            force_sets.append(forces)
        
        # get phonopy_disp.yaml from DB
        # and get disp_dataset from phonopy_disp.yaml
        phonopy_disp_dict = mmdb.collection.find_one(
            {
                "task_label": {"$regex": f"{tag} DisplacedStructuresAdderTask"},
            }
        )
        
        with open("phonopy_disp.yaml", "w") as outfile:
            yaml.dump(phonopy_disp_dict["yaml"], outfile, default_flow_style=False)
        
        disp_dataset = parse_disp_fc3_yaml(filename="phonopy_disp.yaml")
        
        # generate FORCES_FC3
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
                            mesh=mesh,
                            log_level=1, # log_level=0 make phono3py quiet
                           )
        
        # use run_thermal_conductivity()
        # which will read phonopy_disp.yaml and FORCES_FC3
        # this operation will generate file kappa-*.hdf5
        run_thermal_conductivity(phono3py, t_min, t_max, t_step)
        
        # parse kappa-*.hdf5
        f = h5py.File("kappa-m{}{}{}.hdf5".format(*mesh))
        ph3py_dict["temperature"] = f["temperature"][:].tolist()
        ph3py_dict["kappa"] = f["kappa"][:].tolist()
        ph3py_dict["mesh"] = f["mesh"][:].tolist()
        
        # add more informations in ph3py_dict
        ph3py_dict["structure"] = unitcell.as_dict()
        ph3py_dict["formula_pretty"] = unitcell.composition.reduced_formula
        ph3py_dict["success"] = True
        ph3py_dict["supercell_size"] = supercell_size
        
        # store results in ph3py_tasks collection
        coll = mmdb.db["ph3py_tasks"]
        coll.insert_one(ph3py_dict)
        
        
        
        
    
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
        