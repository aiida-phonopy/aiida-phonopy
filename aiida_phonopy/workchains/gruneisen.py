from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.engine import WorkChain, ToContext
from aiida.engine import workfunction
from aiida.engine import run, submit

from aiida.plugins import load_node, DataFactory, WorkflowFactory

from aiida.orm import Str, Float, Bool

# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows

ForceConstantsData = DataFactory('phonopy.force_constants')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

BandStructureData = DataFactory('phonopy.band_structure')
PhononPhonopy = WorkflowFactory('phonopy.phonon')

import numpy as np

__testing__ = False


def get_phonon(structure, force_constants, ph_settings, nac_data=None):
    from phonopy import Phonopy
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    phonon = Phonopy(phonopy_atoms_from_structure(structure),
                     ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_tolerance)

    if force_constants is not None:
        phonon.set_force_constants(force_constants.get_data())

    if nac_data is not None:
            primitive = phonon.get_primitive()
            nac_parameters = nac_data.get_born_parameters_phonopy(primitive_cell=primitive.get_cell())
            phonon.set_nac_params(nac_parameters)

    return phonon


def get_commensurate(structure, ph_settings):
    from phonopy import Phonopy
    from phonopy.harmonic.dynmat_to_fc import DynmatToForceConstants
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    phonon = Phonopy(phonopy_atoms_from_structure(structure),
                     ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_tolerance)

    primitive = phonon.get_primitive()
    supercell = phonon.get_supercell()

    dynmat2fc = DynmatToForceConstants(primitive, supercell)
    commensurate = dynmat2fc.get_commensurate_points()

    return dynmat2fc, commensurate


def get_force_constants(phonon_origin, gruneisen, commensurate, volumes):
    from phonopy.harmonic.dynmat_to_fc import DynmatToForceConstants
    from phonopy.units import VaspToTHz

    from copy import deepcopy
    phonon = deepcopy(phonon_origin)

    phonon.set_qpoints_phonon(commensurate,
                              is_eigenvectors=True)
    frequencies, eigenvectors = phonon.get_qpoints_phonon()

    primitive = phonon.get_primitive()
    supercell = phonon.get_supercell()

    dynmat2fc = DynmatToForceConstants(primitive, supercell)
    volume_ref = phonon.get_unitcell().get_volume()

    force_constants_list = []
    for volume in volumes:
        renormalized_frequencies = []
        for freq, g in zip(frequencies, gruneisen):
            renormalized_frequencies.append(freq + (freq * (np.exp(-g*np.log(volume/volume_ref))-1)))
        renormalized_frequencies = np.array(renormalized_frequencies)

        # Fixing Gamma point data
        renormalized_frequencies[0][0:3] = [0.0, 0.0, 0.0]

        dynmat2fc.set_dynamical_matrices(renormalized_frequencies / VaspToTHz, eigenvectors)
        dynmat2fc.run()
        force_constants_list.append(np.array(dynmat2fc.get_force_constants()))

    return force_constants_list


def get_thermal_properties(structure, ph_settings, force_constants_list):
    from phonopy import Phonopy
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    free_energy_list = []
    entropy_list = []
    cv_list = []
    temperature = None
    for fc in force_constants_list:
        phonon = Phonopy(phonopy_atoms_from_structure(structure),
                         ph_settings.dict.supercell,
                         primitive_matrix=ph_settings.dict.primitive,
                         symprec=ph_settings.dict.symmetry_tolerance)

        # Normalization factor primitive to unit cell
        normalization_factor = phonon.unitcell.get_number_of_atoms()/phonon.primitive.get_number_of_atoms()

        phonon.set_force_constants(fc)
        phonon.set_mesh(ph_settings.dict.mesh, is_eigenvectors=True, is_mesh_symmetry=False)
        phonon.set_thermal_properties()
        temperature, free_energy, entropy, cv = phonon.get_thermal_properties()
        free_energy_list.append(np.array(free_energy) * normalization_factor)
        entropy_list.append(np.array(entropy) * normalization_factor)
        cv_list.append(np.array(cv) * normalization_factor)

    return np.array(free_energy_list), np.array(entropy_list).T, np.array(cv_list).T, temperature


def get_gruneisen_at_list(phonon_origin, phonon_plus, phonon_minus, list_qpoints):

    from phonopy.gruneisen.core import GruneisenBase
    from copy import deepcopy

    gruneisen = GruneisenBase(phonon_origin.get_dynamical_matrix(),
                              phonon_plus.get_dynamical_matrix(),
                              phonon_minus.get_dynamical_matrix())

    gruneisen.set_qpoints(list_qpoints)
    gamma = gruneisen.get_gruneisen()

    factor = phonon_origin.get_unit_conversion_factor(),
    eigenvalues = gruneisen.get_eigenvalues()
    frequencies = np.sqrt(abs(eigenvalues)) * np.sign(eigenvalues) * factor

    phonon = deepcopy(phonon_origin)

    phonon.set_qpoints_phonon(list_qpoints,
                              is_eigenvectors=True)
    frequencies_check, eigenvectors = phonon.get_qpoints_phonon()

    # Make sure that both sets of frequencies are the same (should be!)
    np.testing.assert_almost_equal(frequencies, frequencies_check)

    return gamma, frequencies, eigenvectors


@workfunction
def phonopy_gruneisen(**kwargs):

    from phonopy import PhonopyGruneisen

    phonon_plus_structure = kwargs.pop('phonon_plus_structure')
    phonon_plus_fc = kwargs.pop('phonon_plus_fc')
    phonon_minus_structure = kwargs.pop('phonon_minus_structure')
    phonon_minus_fc = kwargs.pop('phonon_minus_fc')
    phonon_origin_structure = kwargs.pop('phonon_origin_structure')
    phonon_origin_fc = kwargs.pop('phonon_origin_fc')
    ph_settings = kwargs.pop('ph_settings')
    bands = kwargs.pop('bands')

    if 'phonon_plus_nac' in kwargs:
        print ('using nac in guneisen')
        phonon_plus_nac = kwargs.pop('phonon_plus_nac')
        phonon_minus_nac = kwargs.pop('phonon_minus_nac')
        phonon_origin_nac = kwargs.pop('phonon_origin_nac')
    else:
        phonon_plus_nac = None
        phonon_minus_nac = None
        phonon_origin_nac = None

    phonon_plus = get_phonon(phonon_plus_structure,
                             phonon_plus_fc,
                             ph_settings,
                             nac_data=phonon_plus_nac)

    phonon_minus = get_phonon(phonon_minus_structure,
                              phonon_minus_fc,
                              ph_settings,
                              nac_data=phonon_minus_nac)

    phonon_origin = get_phonon(phonon_origin_structure,
                               phonon_origin_fc,
                               ph_settings,
                               nac_data=phonon_origin_nac)

    gruneisen = PhonopyGruneisen(phonon_origin,  # equilibrium
                                 phonon_plus,  # plus
                                 phonon_minus)  # minus

    gruneisen.set_mesh(ph_settings.dict.mesh, is_gamma_center=False, is_mesh_symmetry=True)

    # band structure
    gruneisen.set_band_structure(bands.get_bands())

    band_structure = BandStructureData(bands=bands.get_bands(),
                                       labels=bands.get_labels(),
                                       unitcell=bands.get_unitcell())

    band_structure.set_band_structure_gruneisen(gruneisen.get_band_structure())

    # mesh
    mesh_data = gruneisen.get_mesh()

    mesh_array = ArrayData()
    mesh_array.set_array('frequencies', np.array(mesh_data[2]))
    mesh_array.set_array('gruneisen', np.array(mesh_data[4]))
    mesh_array.set_array('q_points', np.array(mesh_data[0]))
    mesh_array.set_array('eigenvectors', np.array(mesh_data[3]))
    mesh_array.set_array('weights', np.array(mesh_data[1]))

    # commensurate
    dynmat2fc, commensurate_q_points = get_commensurate(phonon_origin_structure, ph_settings)
    commensurate_gruneisen, commensurate_frequencies, eigenvectors = get_gruneisen_at_list(phonon_origin,
                                                                                           phonon_plus,
                                                                                           phonon_minus,
                                                                                           commensurate_q_points)
    commensurate_array = ArrayData()
    commensurate_array.set_array('q_points', commensurate_q_points)
    commensurate_array.set_array('gruneisen', commensurate_gruneisen)
    commensurate_array.set_array('frequencies', commensurate_frequencies)
    commensurate_array.set_array('eigenvectors', eigenvectors)

    return {'band_structure': band_structure, 'mesh': mesh_array, 'commensurate': commensurate_array}


@workfunction
def get_eos(phonon_plus_structure, phonon_origin_structure, phonon_minus_structure,
            phonon_plus_data, phonon_origin_data, phonon_minus_data):

    # Wayaround to VASP  (energy without entropy not correct in vasprun.xml)

    # e_plus = phonon_plus_data.get_array('electronic_steps')[None][-1][-1][-1]['e_wo_entrp']
    # e_origin = phonon_origin_data.get_array('electronic_steps')[None][-1][-1][-1]['e_wo_entrp']
    # e_minus = phonon_minus_data.get_array('electronic_steps')[None][-1][-1][-1]['e_wo_entrp']

    # print 'plus new:', e_plus, phonon_plus_structure.get_cell_volume()
    # print 'origin new:', e_origin, phonon_origin_structure.get_cell_volume()
    # print 'minus new:', e_minus, phonon_minus_structure.get_cell_volume()

    e_plus = phonon_plus_data.dict.energy
    e_origin = phonon_origin_data.dict.energy
    e_minus = phonon_minus_data.dict.energy

    vol_fit = np.polyfit([phonon_plus_structure.get_cell_volume(),
                          phonon_origin_structure.get_cell_volume(),
                          phonon_minus_structure.get_cell_volume()],
                         [e_plus, e_origin, e_minus],
                         2)
    p_energy = np.poly1d(vol_fit)

    stress_fit = np.polyfit([phonon_plus_structure.get_cell_volume(),
                             phonon_origin_structure.get_cell_volume(),
                             phonon_minus_structure.get_cell_volume()],
                            [np.average(np.diag(phonon_plus_data.dict.stress)),
                             np.average(np.diag(phonon_origin_data.dict.stress)),
                             np.average(np.diag(phonon_minus_data.dict.stress))],
                            2)
    p_stress = np.poly1d(stress_fit)

    volumes = np.linspace(phonon_origin_structure.get_cell_volume() * 0.95,
                          phonon_origin_structure.get_cell_volume() * 1.05,
                          num=10)

    eos = ArrayData()
    eos.set_array('volumes', volumes)
    eos.set_array('energies', p_energy(volumes))
    eos.set_array('stresses', p_stress(volumes))

    return {'eos': eos}


@workfunction
def phonopy_qha_prediction(phonon_structure,
                           force_constants,
                           eos,
                           ph_settings,
                           commensurate):

    phonon = get_phonon(phonon_structure,
                        force_constants,
                        ph_settings)

    force_constants_list = get_force_constants(phonon,
                                               commensurate.get_array('gruneisen'),
                                               commensurate.get_array('q_points'),
                                               eos.get_array('volumes'))

    free_energy_list, entropy_list, cv_list, temperature = get_thermal_properties(phonon_structure,
                                                                                  ph_settings,
                                                                                  force_constants_list)

    # testing
    if __testing__:

        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.plot(free_energy_list)

        plt.figure(2)
        plt.plot(eos.get_array('volumes'), eos.get_array('energies'))

        plt.figure(3)
        plt.plot(eos.get_array('stresses'), eos.get_array('energies'))

        plt.figure(4)
        plt.plot(temperature, cv_list)

        plt.figure(5)
        plt.plot(temperature, entropy_list)
        plt.show()
        # end testing

    qha_output = get_qha(eos,
                         temperature,
                         free_energy_list,
                         cv_list,
                         entropy_list)

    stress_fit = np.polyfit(eos.get_array('volumes'),
                            eos.get_array('stresses'), 2)
    p_stress = np.poly1d(stress_fit)

    volumes_test = np.linspace(eos.get_array('volumes')[0]-20, eos.get_array('volumes')[-1]+20, 50)

    # testing
    if __testing__:
        import matplotlib.pyplot as plt
        plt.plot(eos.get_array('volumes'), eos.get_array('stresses'), 'ro')
        plt.plot(volumes_test, p_stress(volumes_test), label='poly')
        plt.plot(volumes_test,
                 fit_pressure_volume(eos.get_array('volumes'), eos.get_array('stresses'), volumes_test),
                 label='inverse')

        plt.legend()
        plt.show()

    # stresses = p_stress(qha_output['volume_temperature'])
    stresses = fit_pressure_volume(eos.get_array('volumes'), eos.get_array('stresses'), qha_output['volume_temperature'])

    prediction = {'temperature': [0.0, qha_output['qha_temperatures'][-1]],
                  'volume_range': [min(qha_output['volume_temperature']), max(qha_output['volume_temperature'])],
                  'stress_range': [max(stresses), min(stresses)]}

    return {'qha_prediction': Dict(dict=prediction)
}


def fit_pressure_volume(volume_data, stress_data, volumes_prediction):
    from scipy.optimize import curve_fit

    def func(x, a, b):
        return a / x + b

    popt, pcov = curve_fit(func, volume_data, stress_data)

    return func(volumes_prediction, *popt)

def get_qha(eos, temperatures, fe_phonon, cv, entropy, t_max=1000):

    from phonopy import PhonopyQHA
    import numpy as np

    phonopy_qha = PhonopyQHA(eos.get_array('volumes'),
                             eos.get_array('energies'),
                             eos="vinet",
                             temperatures=np.array(temperatures),
                             free_energy=np.array(fe_phonon).T,
                             cv=np.array(cv),
                             entropy=np.array(entropy),
                             t_max=t_max,
                             verbose=False)

    # testing
    if __testing__:
        import matplotlib.pyplot as plt
        plt = phonopy_qha.plot_qha()
        plt.show()

    qha_results = {'qha_temperatures': phonopy_qha._qha._temperatures[:phonopy_qha._qha._max_t_index],
                   'helmholtz_volume': phonopy_qha.get_helmholtz_volume(),
                   'thermal_expansion': phonopy_qha.get_thermal_expansion(),
                   'volume_temperature':  phonopy_qha.get_volume_temperature(),
                   'heat_capacity_P_numerical': phonopy_qha.get_heat_capacity_P_numerical(),
                   'volume_expansion': phonopy_qha.get_volume_expansion(),
                   'gibbs_temperature': phonopy_qha.get_gibbs_temperature()}

    return qha_results


class GruneisenPhonopy(WorkChain):
    """
    Workchain to calculate the mode Gruneisen parameters
    """
    @classmethod
    def define(cls, spec):
        super(GruneisenPhonopy, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=Dict)
        spec.input("es_settings", valid_type=Dict)
        # Optional arguments
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))  # in kB
        spec.input("stress_displacement", valid_type=Float, required=False, default=Float(2.0))  # in kB
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))

        spec.outline(cls.create_unit_cell_expansions, cls.calculate_gruneisen)

    def create_unit_cell_expansions(self):

        print('start Gruneisen (pk={})'.format(self.pid))
        print('start create cell expansions')

        # For testing
        testing = False
        if testing:
            self.ctx._content['plus'] = load_node(13603)
            self.ctx._content['origin'] = load_node(13600)
            self.ctx._content['minus'] = load_node(13606)
            return

        calcs = {}
        for expansions in {'plus': float(self.inputs.pressure) + float(self.inputs.stress_displacement),
                           'origin': float(self.inputs.pressure),
                           'minus': float(self.inputs.pressure) - float(self.inputs.stress_displacement)}.items():

            future = submit(PhononPhonopy,
                            structure=self.inputs.structure,
                            ph_settings=self.inputs.ph_settings,
                            es_settings=self.inputs.es_settings,
                            pressure=Float(expansions[1]),
                            optimize=Bool(True),
                            use_nac=self.inputs.use_nac
                            )

            calcs[expansions[0]] = future
            print ('phonon workchain: {} {}'.format(expansions[0], future.pid))

        return ToContext(**calcs)

    def calculate_gruneisen(self):

        self.report('calculate gruneisen')
        print('calculate gruneisen')
        print(self.ctx.plus, self.ctx.minus, self.ctx.origin)

        input_gruneisen = {'phonon_plus_structure' : self.ctx.plus.out.final_structure,
                           'phonon_plus_fc' : self.ctx.plus.out.force_constants,
                           'phonon_minus_structure' : self.ctx.minus.out.final_structure,
                           'phonon_minus_fc': self.ctx.minus.out.force_constants,
                           'phonon_origin_structure' : self.ctx.origin.out.final_structure,
                           'phonon_origin_fc' : self.ctx.origin.out.force_constants,
                           'ph_settings': self.inputs.ph_settings,
                           'bands': self.ctx.origin.out.band_structure}

        if 'nac_data' in self.ctx.origin.get_outputs_dict():
            print ('loading nac (grune)')
            input_gruneisen.update({'phonon_plus_nac': self.ctx.plus.out.nac_data,
                                    'phonon_minus_nac': self.ctx.minus.out.nac_data,
                                    'phonon_origin_nac': self.ctx.origin.out.nac_data})

        gruneisen_results = phonopy_gruneisen(**input_gruneisen)

        self.out('band_structure', gruneisen_results['band_structure'])
        self.out('mesh', gruneisen_results['mesh'])
        self.out('commensurate', gruneisen_results['commensurate'])

        eos = get_eos(phonon_plus_structure=self.ctx.plus.out.final_structure,
                      phonon_origin_structure=self.ctx.origin.out.final_structure,
                      phonon_minus_structure=self.ctx.minus.out.final_structure,
                      phonon_plus_data=self.ctx.plus.out.optimized_data,
                      phonon_origin_data=self.ctx.origin.out.optimized_data,
                      phonon_minus_data=self.ctx.minus.out.optimized_data)

        prediction_results = phonopy_qha_prediction(phonon_structure=self.ctx.origin.out.final_structure,
                                                    force_constants=self.ctx.origin.out.force_constants,
                                                    eos=eos['eos'],
                                                    ph_settings=self.inputs.ph_settings,
                                                    commensurate=gruneisen_results['commensurate'])

        self.out('prediction', prediction_results['qha_prediction'])

        print(prediction_results)
