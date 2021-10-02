from pymatgen.io.vasp.sets import _load_yaml_config, DictSet
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent

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
        
