from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.utils.ph3py import get_displaced_structures
from fireworks import explicit_serialize
from atomate.utils.utils import env_chk
import numpy as np

### Intermediate Firetask that reads optimized structure, generate displaced structures, and make addition of StaticFWs.
@explicit_serialize
class DisplacedStructuresAdderTask(FiretaskBase):
    """
    Intermediate Firetask that reads optimized structure, 
    generate displaced structures, 
    and make addition of StaticFWs.
    """
    required_params = ["tag", "db_file"]
    optional_params = ["supercell_size", "cutoff_pair_distance", "name", "struct_unitcell", "vis_static"]
        
    def run_task(self, fw_spec):
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        supercell_size = self.get("supercell_size", (2,2,2))
        cutoff_pair_distance = self.get("cutoff_pair_distance", None)
        name = self.get("name", "DisplacedStructuresAdderTask")
        struct_unitcell = self.get("struct_unitcell", None)
        vis_static = self.get("vis_static", MPStaticSet)
        
        # read optimized structures
        if struct_unitcell is None:
            atomate_db = VaspCalcDb.from_db_file(db_file)
            coll = atomate_db.db.get_collection("tasks")
            doc = coll.find_one(
                {
                    "task_label": {"$regex": f"{tag} structure optimization"},
                    "formula_pretty": struct_unitcell.composition.reduced_formula,
                },
#                 {"calcs_reversed": 1}, # not sure what it does?
            )
            
            struct_unitcell = Structure.from_dict(doc["calcs_reversed"][0]["output"]["structure"])
        
        # generate displaced structures from optimized structure
        supercell_matrix = np.eye(3) * np.array(supercell_size) # supercell size = 2*2*2
        struct_displaced = get_displaced_structures(
            structure=struct_unitcell,
            supercell_matrix=supercell_matrix,
            yaml_fname="phonopy_disp.yaml",
            cutoff_pair_distance=cutoff_pair_distance,
        )
        # TODO: save phonopy_disp.yaml to DB collection?
        
        # initialize new fws
        new_fws = []
        
        # for each structure in the displaced structures create a StaticFW
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
            new_fws.append(fw)
        
        # return WF of combined FW
        wf = Workflow(new_fws)
        if len(new_fws) != 0:
            return FWAction(addition=wf)
        
        
        