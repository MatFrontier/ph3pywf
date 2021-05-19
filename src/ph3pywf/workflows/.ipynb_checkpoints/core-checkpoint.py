from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.firetasks.core import DisplacedStructuresAdderTask
from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

def wf_phono3py(structure, 
                supercell_size=None, 
                cutoff_pair_distance=None, 
                vasp_input_set_relax=None, 
                vasp_input_set_static=None, 
                tag=None, 
                db_file=DB_FILE, 
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
    tag = tag or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    
    # get the input set for the optimization and update it if we passed custom settings
    vis_relax = vasp_input_set_relax or MPRelaxSet(structure, force_gamma=True)
    
    # Structure optimization firework
    fws = [OptimizeFW(structure=structure, 
                      vasp_input_set=vis_relax, 
                      name=f"{tag} structure optimization", 
                     )]
    
    # get the input set for the static calculations and update it if we passed custom settings
    vis_static = vasp_input_set_static or MPStaticSet(structure, force_gamma=True)
    
    # convert GetDisplacedStructuresFWs to FW and add to FW list
    parents = fws[0]
    
    fw_name = "{}-{} DisplacedStructuresAdderTask".format(
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
        ), 
        name=fw_name, 
        parents=parents, 
    )
    
    fws.append(fw)
    
    # create the workflow
    wf_ph3py = Workflow(fws)
    wf_ph3py.name = "{}:{}".format(structure.composition.reduced_formula, "phono3py calculation")
        
    # post analysis 
    # TODO: write post analysis FW
    
    return wf_ph3py


def wf_disp_from_optimized(structure, 
                           supercell_size=None, 
                           cutoff_pair_distance=None, 
                           vasp_input_set_static=None, 
                           tag=None, 
                           db_file=DB_FILE, 
                          ):
    tag = tag or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    
    # get the input set for the static calculations and update it if we passed custom settings
    vis_static = vasp_input_set_static or MPStaticSet(structure, force_gamma=True)
    
    # call adder FW
    fw_name = "{}-{} DisplacedStructuresAdderTask".format(
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
        ), 
        name=fw_name, 
    )
    
    # create the workflow
    wf = Workflow([fw])
    wf.name = "{}:{}".format(structure.composition.reduced_formula, 
                             "phono3py calculation from optimized struct")
    
    return wf