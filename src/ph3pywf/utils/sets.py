from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
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

    def __init__(self, structure, **kwargs):
        """
        :param structure: Structure
        :param kwargs: Same as those supported by DictSet.
        """
        super().__init__(structure, Ph3pyRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs
        
