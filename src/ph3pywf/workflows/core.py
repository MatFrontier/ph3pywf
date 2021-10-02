# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import FWAction, Firework, FiretaskBase, Workflow
from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
from ph3pywf.firetasks.core import DisplacedStructuresAdderTask, Phono3pyAnalysisToDb
from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from ph3pywf.utils.sets import Ph3pyRelaxSet, Ph3pyStaticSet

def wf_phono3py(structure, 
                name="phono3py wf",
                c=None,
               ):
    """
    Returns Phono3py calculation workflow.
    Args:
        structure (Structure): input structure.
        c (dict): workflow config dict:
            tag (str): unique label to identify contents related to this WF.
            supercell_size_fc3 (tuple): Supercell dimension for 3rd order force constants. (2, 2, 2) by default. 
            supercell_size_fc2 (tuple): Optional supercell dimension for 2nd order force constants. 
            cutoff_pair_distance (float): set to reduce the number of supercells with displacements to be calculated.
            atom_disp (float): atomic displacement. Default is 0.03 Angstrom.
            vasp_input_set_relax (VaspInputSet): input set for optimization VASP job. 
            vasp_input_set_static (VaspInputSet): input set for static VASP jobs.
            db_file (str): path to the db file.
            metadata (dict): meta data.
            is_reduced_test (bool): if set to True, there'll be only a few staticFW generated.
            USER_*_SETTINGS (dict): VASP input parameters to override default VaspInputSet settings.
            t_min (float): min temperature (in K)
            t_max (float): max temperature (in K)
            t_step (float): temperature step (in K)
            primitive_matrix (ndarray): transformation matrix to primitive cell from unit cell.
                Primitive matrix with respect to unit cell.
                shape=(3, 3), dtype='double', order='C'
            is_nac (bool): If True, non-analytical term correction is on. 
    Returns:
        Workflow
    """
    
    c = c or {}
    tag = c.get("tag", datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f'))
    supercell_size_fc3 = c.get("supercell_size_fc3", (2,2,2))
    supercell_size_fc2 = c.get("supercell_size_fc2", None)
    cutoff_pair_distance = c.get("cutoff_pair_distance", None)
    atom_disp = c.get("atom_disp", 0.03)
    vasp_input_set_relax = c.get("vasp_input_set_relax", None)
    vasp_input_set_static = c.get("vasp_input_set_static", None)
    db_file = c.get("db_file", DB_FILE)
    metadata = c.get("metadata", {})
    is_reduced_test = c.get("is_reduced_test", False)
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    user_incar_settings_static = c.get("USER_INCAR_SETTINGS_STATIC", {})
    user_potcar_settings = c.get("USER_POTCAR_SETTINGS", {})
    user_potcar_functional = c.get("USER_POTCAR_FUNCTIONAL", None)
    user_kpoints_settings = c.get("USER_KPOINTS_SETTINGS", None)
    t_min = c.get("t_min", 0)
    t_max = c.get("t_max", 1001)
    t_step = c.get("t_step", 10)
    primitive_matrix = c.get("primitive_matrix", None)
    mesh = c.get("mesh", [20,20,20])
    is_nac = c.get("is_nac", False)
    
    # store tag in metadata
    metadata["label"] = tag
    print(f"tag: \"{tag}\"")
    print(f"{{task_label: {{$regex:\"{tag}\"}}}}")
    
    # update vasp_input_set_relax
    vasp_input_set_relax = vasp_input_set_relax or Ph3pyRelaxSet(structure,
                                                                user_incar_settings=user_incar_settings,
                                                                user_potcar_settings=user_potcar_settings,
                                                                user_potcar_functional=user_potcar_functional,
                                                                user_kpoints_settings=user_kpoints_settings,
                                                               )
    

    # structure optimization firework
    fws = [OptimizeFW(structure=structure, 
                      vasp_input_set=vasp_input_set_relax, 
                      name=f"{tag} structure optimization", 
                     )]
    
    # update vasp_input_set_static
    vasp_input_set_static = vasp_input_set_static or Ph3pyStaticSet(structure,
                                                                    user_incar_settings=user_incar_settings_static,
                                                                    user_potcar_settings=user_potcar_settings,
                                                                    user_potcar_functional=user_potcar_functional,
                                                                    user_kpoints_settings=user_kpoints_settings,
                                                                   )
    
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
            supercell_size_fc3=supercell_size_fc3, 
            supercell_size_fc2=supercell_size_fc2, 
            cutoff_pair_distance=cutoff_pair_distance, 
            atom_disp=atom_disp, 
            vis_static=vasp_input_set_static,
            user_settings=c,
            is_nac=is_nac,
        ), 
        name=fw_name, 
        parents=parents, 
    )
    
    fws.append(fw)
    
    # post analysis 
    parents = fws[-1]
    fw_name = "{}-{} Phono3pyAnalysisToDb".format(
        structure.composition.reduced_formula if structure else "unknown", 
        tag, 
    )
    
    fw = Firework(
        Phono3pyAnalysisToDb(
            tag=tag, 
            db_file=db_file,
            t_min=t_min,
            t_max=t_max,
            t_step=t_step,
            mesh=mesh,
        ),
        name=fw_name, 
        parents=parents,
    )
    
    fws.append(fw)
    
    # create the workflow
    wfname = "{}-{}".format(structure.composition.reduced_formula, name)
    
    return Workflow(fws, name=wfname, metadata=metadata)

def wf_ph3py_post_analysis(tag,
                           db_file_local,
                           name="phono3py post analysis only wf",
                           c=None,
                          ):
    """
    Rerun Phono3pyAnalysisToDb task
    and overwrite the previous post ph3py_task doc
    """
    c = c or {}
    db_file = c.get("db_file", DB_FILE)
    t_min = c.get("t_min", 10)
    t_max = c.get("t_max", 1001)
    t_step = c.get("t_step", 10)
    mesh = c.get("mesh", None)
    born_filename = c.get("born_filename", None)
    
    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file_local, admin=True)
    
    # read addertask_dict from DB
    opt_dict = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} structure optimization"},
        }
    )
    formula_pretty = opt_dict["formula_pretty"]
    
    # prepare FW
    fw_name = "{}-{} Phono3pyAnalysisToDb".format(
        formula_pretty, 
        tag, 
    )
    
    fw = Firework(
        Phono3pyAnalysisToDb(
            tag=tag, 
            db_file=db_file,
            t_min=t_min,
            t_max=t_max,
            t_step=t_step,
            mesh=mesh,
            born_filename=born_filename,
        ),
        name=fw_name, 
    )
    
    wfname = "{}-{}".format(formula_pretty, name)

    return Workflow([fw], name=wfname)

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
