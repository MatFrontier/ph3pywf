# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.utils.ph3py import get_displaced_structures
from fireworks import explicit_serialize
from atomate.utils.utils import env_chk
import numpy as np
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

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
    """
    required_params = ["tag", "db_file"]
    optional_params = ["supercell_size", "cutoff_pair_distance", "struct_unitcell", "vis_static", "name"]
        
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        supercell_size = self.get("supercell_size", (2,2,2))
        cutoff_pair_distance = self.get("cutoff_pair_distance", None)
        struct_unitcell = self.get("struct_unitcell", None)
        vis_static = self.get("vis_static", MPStaticSet)
        name = self.get("name", "DisplacedStructuresAdderTask")
        
        # read optimized structures
        if struct_unitcell is None:
            atomate_db = VaspCalcDb.from_db_file(db_file)
            coll = atomate_db.db.get_collection("tasks")
            doc = coll.find_one(
                {
                    "task_label": {"$regex": f"{tag} structure optimization"},
                },
#                 {"calcs_reversed": 1}, # appears to be "projection", not sure if necessary
            )
            
            struct_unitcell = Structure.from_dict(doc["calcs_reversed"][0]["output"]["structure"])
        
        # generate displaced structures from optimized structure
        supercell_matrix = np.eye(3) * np.array(supercell_size) 
        struct_displaced = get_displaced_structures(
            structure=struct_unitcell,
            supercell_matrix=supercell_matrix,
            yaml_fname="phonopy_disp.yaml",
            cutoff_pair_distance=cutoff_pair_distance,
        )
        # TODO: save phonopy_disp.yaml to DB collection
        # Should it be stored in "fireworks" collection
        # or "tasks" collection?
        # Perhaps this should be another firetask?
        
        # initialize new fws
        new_fws = []
        
        # for each structure in the displaced structures append a StaticFW
        for i, structure in enumerate(struct_displaced):
            if i==0: 
                continue # Skip undeformed supercell
            if i>10: # JUST FOR TESTING !!!!!!!!!!!!!!!!
                break # JUST FOR TESTING !!!!!!!!!!!!!!!!
            disp_id = f"{i:05d}"
            fw = StaticFW(structure=structure,
                          vasp_input_set=vis_static, 
                          name=f"{tag} disp-{disp_id}",
                         )
            print(f"adding FW {disp_id}")
            print(structure.composition.reduced_formula)
            new_fws.append(fw)
        
        # return additions of combined FWs
        if len(new_fws) != 0:
            print("returning FWAction")
            return FWAction(detours=new_fws)
        
        
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
        atomate_db = VaspCalcDb.from_db_file(db_file)
#         coll = atomate_db.db.get_collection("_IOtest")
        coll = atomate_db.db.get_collection("tasks")
        
        # create structure dict
        struct_dict = {"calcs_reversed":[{"output":{"structure":None}}],"task_label":None}
        struct_dict["task_label"] = f"{tag} {name}"
        struct_dict["calcs_reversed"][0]["output"]["structure"] = structure.as_dict()
        
        # store structure in DB collection
        coll.insert_one(struct_dict)
        
        # terminate FW dependently
        if terminate:
            return FWAction()