# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from datetime import datetime
from fireworks.core.firework import Firework, Workflow
from atomate.vasp.database import VaspCalcDb
from ph3pywf.firetasks.core import (
    DisplacedStructuresAdderTask,
    Phono3pyAnalysisToDb,
    Phono3pyMeshConvergenceToDb,
    Phono3pyEvaluateKappaFromConvTest,
)
from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from ph3pywf.utils.sets import Ph3pyRelaxSet, Ph3pyStaticSet
from ph3pywf.fireworks.core import ForceSymmOptimizeFW, ForceSymmStaticFW


def wf_phono3py(
    structure,
    skip_relax=False,
    name="phono3py wf",
    c=None,
):
    """
    Returns Phono3py calculation workflow.
    Args:
        structure (Structure): input structure.
        skip_relax (bool): if True, skip relaxation and directly feed structure to Adder task.
        c (dict): workflow config dict:
            tag (str): unique label to identify contents related to this WF.
            supercell_size_fc3 (tuple): Supercell dimension for 3rd order force constants.
                (2, 2, 2) by default.
            supercell_size_fc2 (tuple): Optional supercell dimension for 2nd order force constants.
            cutoff_pair_distance (float): set to reduce the number of supercells
                with displacements to be calculated.
            atom_disp (float): atomic displacement. Default is 0.03 Angstrom.
            vasp_input_set_relax (VaspInputSet): input set for optimization VASP job.
            vasp_input_set_static (VaspInputSet): input set for static VASP jobs.
            db_file (str): path to the db file.
            metadata (dict): meta data.
            USER_*_SETTINGS (dict): VASP input parameters to override default VaspInputSet settings.
            t_min (float): min temperature (in K)
            t_max (float): max temperature (in K)
            t_step (float): temperature step (in K)
            primitive_matrix (ndarray): transformation matrix to primitive cell from unit cell.
                Primitive matrix with respect to unit cell.
                shape=(3, 3), dtype='double', order='C'
            mesh (list): sampling mesh numbers in reciprocal space.
            is_nac (bool): if True, non-analytical term correction is on.
            is_symmetry (bool): parameter for phono3py instance.
            symprec (float): structural symmetry tolerance.
                If specified, will override symprec parameter for phono3py instance and VASP setting.
    Returns:
        Workflow
    """

    c = c or {}
    tag = c.get("tag", datetime.utcnow().strftime("%Y-%m-%d-%H-%M-%S-%f"))
    supercell_size_fc3 = c.get("supercell_size_fc3", (2, 2, 2))
    supercell_size_fc2 = c.get("supercell_size_fc2", None)
    cutoff_pair_distance = c.get("cutoff_pair_distance", None)
    atom_disp = c.get("atom_disp", 0.03)
    vasp_input_set_relax = c.get("vasp_input_set_relax", None)
    vasp_input_set_static = c.get("vasp_input_set_static", None)
    db_file = c.get("db_file", DB_FILE)
    metadata = c.get("metadata", {})
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    user_incar_settings_static = c.get("USER_INCAR_SETTINGS_STATIC", {})
    user_potcar_settings = c.get("USER_POTCAR_SETTINGS", {})
    user_potcar_functional = c.get("USER_POTCAR_FUNCTIONAL", None)
    user_kpoints_settings = c.get("USER_KPOINTS_SETTINGS", None)
    user_kpoints_settings_static = c.get("USER_KPOINTS_SETTINGS_STATIC", None)
    t_min = c.get("t_min", 200)
    t_max = c.get("t_max", 1401)
    t_step = c.get("t_step", 50)
    primitive_matrix = c.get("primitive_matrix", None)
    mesh = c.get("mesh", [11, 11, 11])
    is_nac = c.get("is_nac", False)
    is_symmetry = c.get("is_symmetry", True)
    symprec = c.get("symprec", 1e-5)

    # store tag in metadata
    metadata["label"] = tag
    print(f'tag = "{tag}"')
    print(f'{{task_label: {{$regex:"{tag}"}}}}')

    # update symprec in vis user settings
    if not "SYMPREC" in user_incar_settings:
        user_incar_settings["SYMPREC"] = symprec
    if not "SYMPREC" in user_incar_settings_static:
        user_incar_settings_static["SYMPREC"] = symprec

    # update vasp_input_set_relax
    vasp_input_set_relax = vasp_input_set_relax or Ph3pyRelaxSet(
        structure,
        user_incar_settings=user_incar_settings,
        user_potcar_settings=user_potcar_settings,
        user_potcar_functional=user_potcar_functional,
        user_kpoints_settings=user_kpoints_settings,
    )

    fws = []

    # structure optimization firework
    if not skip_relax:
        fws.append(
            ForceSymmOptimizeFW(
                structure=structure,
                vasp_input_set=vasp_input_set_relax,
                name=f"{tag} structure optimization",
                job_type="normal",
            )
        )

    # update vasp_input_set_static
    vasp_input_set_static = vasp_input_set_static or Ph3pyStaticSet(
        structure,
        user_incar_settings=user_incar_settings_static,
        user_potcar_settings=user_potcar_settings,
        user_potcar_functional=user_potcar_functional,
        user_kpoints_settings=user_kpoints_settings_static,
    )

    # convert GetDisplacedStructuresFWs to FW and add to FW list
    parents = fws[0] if not skip_relax else None

    fw_name = "{}-{} DisplacedStructuresAdderTask".format(
        structure.composition.reduced_formula if structure else "unknown",
        tag,
    )

    fw = Firework(
        DisplacedStructuresAdderTask(
            tag=tag,
            db_file=db_file,
            struct_unitcell=structure if skip_relax else None,
            supercell_size_fc3=supercell_size_fc3,
            supercell_size_fc2=supercell_size_fc2,
            cutoff_pair_distance=cutoff_pair_distance,
            atom_disp=atom_disp,
            vis_static=vasp_input_set_static,
            primitive_matrix=primitive_matrix,
            user_settings=c,
            is_nac=is_nac,
            is_symmetry=is_symmetry,
            symprec=symprec,
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


def wf_ph3py_post_analysis(
    tag,
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
    t_min = c.get("t_min", 200)
    t_max = c.get("t_max", 1401)
    t_step = c.get("t_step", 50)
    mesh = c.get("mesh", None)
    born_filename = c.get("born_filename", None)

    metadata = c.get("metadata", {})
    metadata["label"] = tag

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

    return Workflow([fw], name=wfname, metadata=metadata)


def wf_ph3py_get_kappa_convergence(
    tag,
    db_file_local,
    name="phono3py get kappa from convergence test wf",
    tag_for_copy=None,
    c=None,
):
    """
    Rerun Phono3pyMeshConvergenceToDb task
    and Phono3pyEvaluateKappaFromConvTest task
    update kappa in ph3py_task
    """
    c = c or {}
    db_file = c.get("db_file", DB_FILE)
    t_min = c.get("t_min", 200)
    t_max = c.get("t_max", 1401)
    t_step = c.get("t_step", 50)
    mesh_densities = c.get("mesh_densities", [256 * k for k in range(1, 100)])

    metadata = c.get("metadata", {})
    metadata["label"] = tag

    # connect to DB
    mmdb = VaspCalcDb.from_db_file(db_file_local, admin=True)

    # read addertask_dict from DB
    opt_dict = mmdb.collection.find_one(
        {
            "task_label": {"$regex": f"{tag} structure optimization"},
        }
    )
    formula_pretty = opt_dict["formula_pretty"]

    # prepare convergence run FW
    fw_name = "{}-{} Phono3pyMeshConvergenceToDb".format(
        formula_pretty,
        tag,
    )

    fw = Firework(
        Phono3pyMeshConvergenceToDb(
            tag=tag,
            db_file=db_file,
            t_min=t_min,
            t_max=t_max,
            t_step=t_step,
            mesh_densities=mesh_densities,
            tag_for_copy=tag_for_copy,
        ),
        name=fw_name,
    )

    fws = [fw]

    # prepare evaluate kappa FW
    fw_name = "{}-{} Phono3pyEvaluateKappaFromConvTest".format(
        formula_pretty,
        tag,
    )

    parents = fws[0]

    fw = Firework(
        Phono3pyEvaluateKappaFromConvTest(
            tag=tag,
            db_file=db_file,
            tag_for_copy=tag_for_copy,
        ),
        name=fw_name,
        parents=parents,
    )

    fws.append(fw)

    wfname = "{}-{}".format(formula_pretty, name)

    return Workflow(fws, name=wfname, metadata=metadata)


#########################
# TESTING MODULES BELOW #
#########################


def wf_ph3py_convergence_test(
    tag,
    db_file_local,
    name="phono3py convergence test wf",
    c=None,
):
    """
    Rerun Phono3pyMeshConvergenceToDb task
    store doc to ph3py_tasks_convergence_test collection
    """
    c = c or {}
    db_file = c.get("db_file", DB_FILE)
    t_min = c.get("t_min", 200)
    t_max = c.get("t_max", 1401)
    t_step = c.get("t_step", 50)
    mesh_densities = c.get("mesh_densities", [256 * k for k in range(1, 100)])
    mesh_list = c.get("mesh_list", None)
    #     print(mesh_densities) # FOR TESTING

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
    fw_name = "{}-{} Phono3pyMeshConvergenceToDb".format(
        formula_pretty,
        tag,
    )

    fw = Firework(
        Phono3pyMeshConvergenceToDb(
            tag=tag,
            db_file=db_file,
            t_min=t_min,
            t_max=t_max,
            t_step=t_step,
            mesh_densities=mesh_densities,
            mesh_list=mesh_list,
        ),
        name=fw_name,
    )

    wfname = "{}-{}".format(formula_pretty, name)

    return Workflow([fw], name=wfname)
