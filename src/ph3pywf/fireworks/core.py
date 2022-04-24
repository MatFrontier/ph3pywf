# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import warnings

from fireworks import Firework
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import (
    MITMDSet,
    MITRelaxSet,
    MPRelaxSet,
    MPScanRelaxSet,
    MPSOCSet,
    MPStaticSet,
)

from atomate.common.firetasks.glue_tasks import (
    CopyFiles,
    DeleteFiles,
    GzipDir,
    PassCalcLocs,
)
from atomate.vasp.config import (
    DB_FILE,
    HALF_KPOINTS_FIRST_RELAX,
    RELAX_MAX_FORCE,
    VASP_CMD,
    VDW_KERNEL_DIR,
)
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs, pass_vasp_result
from atomate.vasp.firetasks.neb_tasks import (
    TransferNEBTask,
    WriteNEBFromEndpoints,
    WriteNEBFromImages,
)
from atomate.vasp.firetasks.parse_outputs import BoltztrapToDb, VaspToDb
from atomate.vasp.firetasks.run_calc import RunBoltztrap, RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import (
    ModifyIncar,
    WriteNormalmodeDisplacedPoscar,
    WriteScanRelaxFromPrev,
    WriteTransmutedStructureIOSet,
    WriteVaspFromIOSet,
    WriteVaspFromIOSetFromInterpolatedPOSCAR,
    WriteVaspHSEBSFromPrev,
    WriteVaspNSCFFromPrev,
    WriteVaspSOCFromPrev,
    WriteVaspStaticFromPrev,
)

class ForceSymmOptimizeFW(Firework):
    def __init__(
        self,
        structure,
        name="structure optimization",
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        ediffg=None,
        db_file=DB_FILE,
        force_gamma=True,
        job_type="double_relaxation_run",
        max_force_threshold=RELAX_MAX_FORCE,
        auto_npar=">>auto_npar<<",
        half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
        parents=None,
        **kwargs,
    ):
        """
        Modified based on atomate.vasp.fireworks.core.OptimizeFW
        Remove "symprec_noise" from errors_subset_to_detect to force symmetry.
        Optimize the given structure.
        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            job_type (str): custodian job type (default "double_relaxation_run")
            max_force_threshold (float): max force on a site allowed at end; otherwise, reject job
            auto_npar (bool or str): whether to set auto_npar. defaults to env_chk: ">>auto_npar<<"
            half_kpts_first_relax (bool): whether to use half the kpoints for the first relaxation
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, force_gamma=force_gamma, **override_default_vasp_params
        )

        if (
            vasp_input_set.incar["ISIF"] in (0, 1, 2, 7)
            and job_type == "double_relaxation"
        ):
            warnings.warn(
                f"A double relaxation run might not be appropriate with ISIF {vasp_input_set.incar['ISIF']}"
            )

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        # modify handler_group
        # remove "symprec_noise" from errors_subset_to_catch
        subset = list(VaspErrorHandler.error_msgs.keys())
        subset.pop("symprec_noise")
        handler_group = [
            VaspErrorHandler(errors_subset_to_catch=subset),
            MeshSymmetryErrorHandler(),
            UnconvergedErrorHandler(),
            NonConvergingErrorHandler(),
            PotimErrorHandler(),
            PositiveEnergyErrorHandler(),
            FrozenJobErrorHandler(),
            StdErrHandler(),
            LargeSigmaHandler(),
            IncorrectSmearingHandler(),
        ]

        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                handler_group=handler_group,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super().__init__(
            t,
            parents=parents,
            name=f"{structure.composition.reduced_formula}-{name}",
            **kwargs,
        )

class ForceSymmStaticFW(Firework):
    def __init__(
        self,
        structure=None,
        name="static",
        vasp_input_set=None,
        vasp_input_set_params=None,
        vasp_cmd=VASP_CMD,
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=DB_FILE,
        vasptodb_kwargs=None,
        parents=None,
        spec_structure_key=None,
        **kwargs,
    ):
        """
        Modified based on atomate.vasp.fireworks.core.StaticFW
        Standard static calculation Firework - either from a previous location or from a structure.
        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        formula = structure.composition.reduced_formula if structure else "unknown"
        fw_name = f"{formula}-{name}"

        if spec_structure_key is not None:
            vasp_input_set = vasp_input_set or MPStaticSet(
                structure, **vasp_input_set_params
            )
            t.append(
                WriteVaspFromIOSet(
                    structure=structure,
                    vasp_input_set=vasp_input_set,
                    spec_structure_key="prev_calc_structure",
                )
            )
        elif prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif parents:
            if prev_calc_loc:
                t.append(
                    CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True)
                )
            t.append(WriteVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif structure:
            vasp_input_set = vasp_input_set or MPStaticSet(
                structure, **vasp_input_set_params
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        # modify handler_group
        # remove "symprec_noise" from errors_subset_to_catch
        subset = list(VaspErrorHandler.error_msgs.keys())
        subset.pop("symprec_noise")
        handler_group = [
            VaspErrorHandler(errors_subset_to_catch=subset),
            MeshSymmetryErrorHandler(),
            UnconvergedErrorHandler(),
            NonConvergingErrorHandler(),
            PotimErrorHandler(),
            PositiveEnergyErrorHandler(),
            FrozenJobErrorHandler(),
            StdErrHandler(),
            LargeSigmaHandler(),
            IncorrectSmearingHandler(),
        ]

        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                handler_group=handler_group,
                auto_npar=">>auto_npar<<",
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super().__init__(t, parents=parents, name=fw_name, **kwargs)