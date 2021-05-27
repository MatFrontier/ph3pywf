# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.firetasks.core import DisplacedStructuresAdderTask
# from ph3pywf.firetasks.core import Phono3pyAnalysisToDb
from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

def wf_phono3py(structure, 
                name="phono3py wf",
                c=None,
               ):
    """
    Returns Phono3py calculation workflow.
    Args:
        structure (Structure): input structure.
        supercell_size (tuple): (2, 2, 2) by default. 
        cutoff_pair_distance (float)
        vasp_input_set_relax (VaspInputSet)
        vasp_input_set_static (VaspInputSet)
        tag (str): unique label to identify contents related to this WF.
    Returns:
        Workflow
    """
    
    c = c or {}
    tag = c.get("tag", datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f'))
    supercell_size = c.get("supercell_size", None)
    cutoff_pair_distance = c.get("cutoff_pair_distance", None)
    vasp_input_set_relax = c.get("vasp_input_set_relax", None)
    vasp_input_set_static = c.get("vasp_input_set_static", None)
    db_file = c.get("db_file", DB_FILE)
    metadata = c.get("metadata", None)
    spec = c.get("spec", None)
    is_reduced_test = c.get("is_reduced_test", False)
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    user_potcar_settings = c.get("USER_POTCAR_SETTINGS", {})
    user_potcar_functional = c.get("USER_POTCAR_FUNCTIONAL", None)
    
    # update vasp_input_set_relax
    vasp_input_set_relax = vasp_input_set_relax or MPRelaxSet(structure,
                                                              user_incar_settings=user_incar_settings,
                                                              user_potcar_settings=user_potcar_settings,
                                                              user_potcar_functional=user_potcar_functional,
                                                             )
    
    # structure optimization firework
    fws = [OptimizeFW(structure=structure, 
                      vasp_input_set=vasp_input_set_relax, 
                      name=f"{tag} structure optimization", 
                     )]
    
    # update vasp_input_set_static
    vasp_input_set_static = vasp_input_set_static or MPStaticSet(structure,
                                                                 user_incar_settings=user_incar_settings,
                                                                 user_potcar_settings=user_potcar_settings,
                                                                 user_potcar_functional=user_potcar_functional,
                                                                )
    
    # convert GetDisplacedStructuresFWs to FW and add to FW list
    parents = fws[0]
    
    fw_name = "{}:{} DisplacedStructuresAdderTask".format(
        structure.composition.reduced_formula if structure else "unknown", 
        tag, 
    )
    
    fw = Firework(
        DisplacedStructuresAdderTask(
            tag=tag, 
            db_file=db_file, 
            supercell_size=supercell_size, 
            cutoff_pair_distance=cutoff_pair_distance, 
            vis_static=vasp_input_set_static, 
            is_reduced_test=is_reduced_test,
        ), 
        name=fw_name, 
        parents=parents, 
    )
    
    fws.append(fw)
    
#     # post analysis 
#     parents = fws[-1]
#     fw_name = "{}:{} Phono3pyAnalysisToDb".format(
#         structure.composition.reduced_formula if structure else "unknown", 
#         tag, 
#     )
    
#     fw = Firework(
#         Phono3pyAnalysisToDb(
#             tag=tag, 
#             db_file=db_file,
#         ),
#         name=fw_name, 
#         parents=parents,
#     )
    
#     fws.append(fw)
    
    # create the workflow
    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    
    return Workflow(fws, name=wfname, metadata=metadata)



#########################
# TESTING MODULES BELOW #
#########################

def wf_disp_from_optimized(structure, 
                           name="phono3py wf from optimized struct", 
                           c=None, 
                          ):
    c = c or {}
    tag = c.get("tag", datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f'))
    supercell_size = c.get("supercell_size", None)
    cutoff_pair_distance = c.get("cutoff_pair_distance", None)
    vasp_input_set_relax = c.get("vasp_input_set_relax", None)
    vasp_input_set_static = c.get("vasp_input_set_static", None)
    db_file = c.get("db_file", DB_FILE)
    metadata = c.get("metadata", None)
    spec = c.get("spec", None)
    is_reduced_test = c.get("is_reduced_test", False)
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    user_potcar_settings = c.get("USER_POTCAR_SETTINGS", {})
    user_potcar_functional = c.get("USER_POTCAR_FUNCTIONAL", None)
    
    from ph3pywf.firetasks.core import Phono3pyAnalysisToDb # just to check if detours work
    
    # update vasp_input_set_static
    vasp_input_set_static = vasp_input_set_static or MPStaticSet(structure,
                                                                 user_incar_settings=user_incar_settings,
                                                                 user_potcar_settings=user_potcar_settings,
                                                                 user_potcar_functional=user_potcar_functional,
                                                                )
    
    # call adder FW
    fw_name = "{}:{} DisplacedStructuresAdderTask".format(
        structure.composition.reduced_formula if structure else "unknown", 
        tag, 
    )
    
    fw = Firework(
        DisplacedStructuresAdderTask(
            tag=tag, 
            db_file=db_file, 
            supercell_size=supercell_size, 
            cutoff_pair_distance=cutoff_pair_distance, 
            struct_unitcell=structure, 
            vis_static=vasp_input_set_static, 
            is_reduced_test=is_reduced_test,
        ), 
        name=fw_name, 
#         parents=parents, 
    )
    
    fws = [fw]
    
    # post analysis 
    parents = fws[-1]
    fw_name = "{}:{} Phono3pyAnalysisToDb".format(
        structure.composition.reduced_formula if structure else "unknown", 
        tag, 
    )
    
    fw = Firework(
        Phono3pyAnalysisToDb(
            tag=tag, 
            db_file=db_file,
        ),
        name=fw_name, 
        parents=parents,
    )
    
    fws.append(fw)

    # create the workflow
    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    
    return Workflow(fws, name=wfname, metadata=metadata)



def wf_disp_from_dynatest(structure, 
                          supercell_size=None, 
                          cutoff_pair_distance=None, 
                          vasp_input_set_static=None, 
                          tag=None, 
                          db_file=DB_FILE, 
                         ):
    from ph3pywf.firetasks.core import StoreStructureTask
    
    tag = tag or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    
    fws = [Firework(
        StoreStructureTask(
            tag=tag, 
            db_file=db_file, 
            structure=structure, 
            name="structure optimization",
            terminate=False,
        ),
        name="pseudo structure optimization",
    )]
    
    # get the input set for the static calculations and update it if we passed custom settings
    vis_static = vasp_input_set_static or MPStaticSet(structure, force_gamma=True)
    
    # call adder FW
    fw_name = "{}:{} DisplacedStructuresAdderTask".format(
        structure.composition.reduced_formula if structure else "unknown", 
        tag, 
    )
    
    parents = fws[0]
    
    fw = Firework(
        DisplacedStructuresAdderTask(
            tag=tag, 
            db_file=db_file, 
            supercell_size=supercell_size, 
            cutoff_pair_distance=cutoff_pair_distance, 
            vis_static=vasp_input_set_static, 
        ), 
        name=fw_name, 
        parents=parents,
    )
    fws.append(fw)
    
    # create the workflow
    wf = Workflow(fws)
    wf.name = "{}:{}".format(structure.composition.reduced_formula, 
                             "phono3py calculation using pseudo struct opt")
    
    return wf
