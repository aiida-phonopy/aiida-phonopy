from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.work.run import run, submit, async

from aiida.orm import load_node, DataFactory, WorkflowFactory
from aiida.orm.data.base import Str, Float, Bool, Int
from aiida.work.workchain import _If, _While

import numpy as np

__testing__ = False

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

BandStructureData = DataFactory('phonopy.band_structure')
GruneisenPhonopy = WorkflowFactory('phonopy.gruneisen')
PhononPhonopy = WorkflowFactory('phonopy.phonon')

OptimizeStructure = WorkflowFactory('phonopy.optimize')

class PhononPhono3py(WorkChain):


    @classmethod
    def define(cls, spec):
        super(PhononPhono3py, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("optimize", valid_type=Bool, required=False, default=Bool(True))
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))  # false by default
        spec.input("chunks", valid_type=Int, required=False, default=Int(100))

        spec.outline(cls.harmonic_calculation,
                     _While(cls.not_converged)(cls.get_data_sets,
                                               cls.get_thermal_conductivity),
                     cls.collect_data)

    def not_converged(self):

        pass

    def harmonic_calculation(self):
        pass

    def get_data_sets(self):
        pass

    def get_thermal_conductivity(self):

        pass

    def collect_data(self):

        pass
