from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Incar
from pathlib import Path
from monty.serialization import loadfn

MODULE_DIR = Path(__file__).resolve().parent

def _load_yaml_config(fname):
    config = loadfn(str(MODULE_DIR / ("%s.yaml" % fname)))
    if "PARENT" in config:
        parent_config = _load_yaml_config(config["PARENT"])
        for k, v in parent_config.items():
            if k not in config:
                config[k] = v
            elif isinstance(v, dict):
                v_new = config.get(k, {})
                v_new.update(v)
                config[k] = v_new
    return config


class Ph3pyRelaxSet(DictSet):
    """
    Implementation of VaspInputSet for relaxation jobs in Ph3py workflow. 
    """

    CONFIG = _load_yaml_config("Ph3pyRelaxSet")

    def __init__(self, structure, **kwargs):
        """
        :param structure: Structure
        :param kwargs: Same as those supported by DictSet.
        """
        super().__init__(structure, Ph3pyRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class Ph3pyStaticSet(DictSet):
    """
    Implementation of VaspInputSet for static jobs in Ph3py workflow. 
    """

    CONFIG = _load_yaml_config("Ph3pyStaticSet")

    def __init__(
        self, 
        structure,
        lepsilon=False,
        lcalcpol=False,
        **kwargs):
        """
        Args:
            structure (Structure): Structure from previous run.
            lepsilon (bool): Whether to add static dielectric calculation.
            lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
                for electronic polarization
            **kwargs: kwargs supported by MPRelaxSet.
        """
        super().__init__(structure, Ph3pyStaticSet.CONFIG, **kwargs)
        self.kwargs = kwargs
        self.lepsilon = lepsilon
        self.lcalcpol = lcalcpol
        
    @property
    def incar(self):
        """
        :return: Incar
        """
        parent_incar = super().incar
        incar = Incar(parent_incar)

        if self.lepsilon:
            incar["IBRION"] = 8
            incar["LEPSILON"] = True

            # LPEAD=T: numerical evaluation of overlap integral prevents
            # LRF_COMMUTATOR errors and can lead to better expt. agreement
            # but produces slightly different results
            incar["LPEAD"] = True

            # Note that DFPT calculations MUST unset NSW. NSW = 0 will fail
            # to output ionic.
            incar.pop("NSW", None)
            incar.pop("NPAR", None)

            # tighter ediff for DFPT
            incar["EDIFF"] = 1e-5

        if self.lcalcpol:
            incar["LCALCPOL"] = True

        for k in ["MAGMOM", "NUPDOWN"] + list(self.user_incar_settings.keys()):
            # For these parameters as well as user specified settings, override
            # the incar settings.
            if parent_incar.get(k, None) is not None:
                incar[k] = parent_incar[k]
            else:
                incar.pop(k, None)

        # use new LDAUU when possible b/c the Poscar might have changed
        # representation
        if incar.get("LDAU"):
            u = incar.get("LDAUU", [])
            j = incar.get("LDAUJ", [])
            if sum([u[x] - j[x] for x, y in enumerate(u)]) > 0:
                for tag in ("LDAUU", "LDAUL", "LDAUJ"):
                    incar.update({tag: parent_incar[tag]})
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in incar:
                incar.update({"LMAXMIX": parent_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        incar["EDIFF"] = min(incar.get("EDIFF", 1), parent_incar["EDIFF"])
        return incar
