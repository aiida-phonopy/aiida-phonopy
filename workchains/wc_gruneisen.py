# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.work.run import run, submit, async

from aiida.orm import Code, CalculationFactory, load_node, DataFactory, DataFactory, WorkflowFactory

from aiida.orm.data.base import Str, Float, Bool

# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows
from aiida.workflows.wc_phonon import PhononPhonopy, get_path_using_seekpath, get_born_parameters

ForceConstantsData = DataFactory('force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

import numpy as np
from generate_inputs import *

def get_phonon(structure, force_constants, ph_settings):
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy import Phonopy

    # Generate phonopy phonon object
    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)

    phonon = Phonopy(bulk,
                     ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_precision)

    phonon.set_force_constants(force_constants.get_array())

    if force_constants.epsilon_and_born_exist():
        phonon.set_nac_params(get_born_parameters(phonon,
                                                  force_constants.get_born_charges(),
                                                  force_constants.get_epsilon(),
                                                  ph_settings.dict.symmetry_precision))
    return phonon

@workfunction
def phonopy_gruneisen(phonon_plus_structure,
                      phonon_plus_fc,
                      phonon_minus_structure,
                      phonon_minus_fc,
                      phonon_origin_structure,
                      phonon_origin_fc,
                      ph_settings):
    from phonopy import PhonopyGruneisen

    phonon_plus = get_phonon(phonon_plus_structure,
                             phonon_plus_fc,
                             ph_settings)

    phonon_minus = get_phonon(phonon_minus_structure,
                              phonon_minus_fc,
                              ph_settings)

    phonon_origin = get_phonon(phonon_origin_structure,
                               phonon_origin_fc,
                               ph_settings)

    gruneisen = PhonopyGruneisen(phonon_origin,  # equilibrium
                                 phonon_plus,  # plus
                                 phonon_minus)  # minus

    gruneisen.set_mesh(ph_settings.dict.mesh, is_gamma_center=False, is_mesh_symmetry=True)

    # BAND STRUCTURE
    band_structure = get_path_using_seekpath(phonon_origin.get_primitive())
    gruneisen.set_band_structure(band_structure.get_band_ranges(),
                                 band_structure.get_number_of_points())

    band_structure.set_band_structure_gruneisen(gruneisen.get_band_structure())

    # mesh
    mesh = gruneisen.get_mesh()
    frequencies_mesh = np.array(mesh.get_frequencies())
    gruneisen_mesh = np.array(mesh.get_gruneisen())
    q_points_mesh = np.array(mesh.get_qpoints())
    weights_mesh = np.array(mesh.get_weights())
    eigenvalues_mesh = np.array(mesh.get_eigenvalues())

    # build mesh
    mesh_array = ArrayData()
    mesh_array.set_array('frequencies', frequencies_mesh)
    mesh_array.set_array('gruneisen', gruneisen_mesh)
    mesh_array.set_array('q_points', q_points_mesh)
    mesh_array.set_array('eigenvalues', eigenvalues_mesh)
    mesh_array.set_array('weights', weights_mesh)


    return {'band_structure': band_structure, 'mesh': mesh_array}


class GruneisenPhonopy(WorkChain):
    """
    Workchain to calculate the mode Gruneisen parameters
    """
    @classmethod
    def define(cls, spec):
        super(GruneisenPhonopy, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("machine", valid_type=ParameterData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        # spec.input("optimize", valid_type=Bool, required=False, default=Bool(True))
        #spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("stress_displacement", valid_type=Float, required=False, default=Float(1e-2))

        spec.outline(cls.create_unit_cell_expansions, cls.calculate_gruneisen)

    def create_unit_cell_expansions(self):

        print('start Gruneisen {}'.format(self.pid))
        print ('start create cell expansions')

        # For testing
        testing = False
        if testing:
            self.ctx._content['plus'] = load_node(9842)
            self.ctx._content['origin'] = load_node(9839)
            self.ctx._content['minus'] = load_node(9845)
            return

        calcs = {}
        for expansions in {'plus': float(self.inputs.stress_displacement),
                           'origin': 0.0,
                           'minus': -float(self.inputs.stress_displacement)}.items():

            future = submit(PhononPhonopy,
                            structure=self.inputs.structure,
                            machine=self.inputs.machine,
                            ph_settings=self.inputs.ph_settings,
                            es_settings=self.inputs.es_settings,
                            pressure=Float(expansions[1]),
                            optimize=Bool(True)
                            )

            calcs[expansions[0]] = future
            print ('phonon workchain: {} {}'.format(expansions[0], future.pid))

        return ToContext(**calcs)

    def calculate_gruneisen(self):

        self.report('calculate gruneisen')
        print ('calculate gruneisen')
        print self.ctx.plus, self.ctx.minus, self.ctx.origin

        gruneisen_results = phonopy_gruneisen(phonon_plus_structure=self.ctx.plus.out.final_structure,
                                              phonon_plus_fc=self.ctx.plus.out.force_constants,
                                              phonon_minus_structure=self.ctx.minus.out.final_structure,
                                              phonon_minus_fc=self.ctx.minus.out.force_constants,
                                              phonon_origin_structure=self.ctx.origin.out.final_structure,
                                              phonon_origin_fc=self.ctx.origin.out.force_constants,
                                              ph_settings=self.inputs.ph_settings)

        self.out('band_structure', gruneisen_results['band_structure'])
        self.out('mesh', gruneisen_results['mesh'])