from aiida.work.workchain import WorkChain, ToContext
from aiida.orm import Code, DataFactory, WorkflowFactory
from aiida.orm.data.base import Float, Bool, Int
from aiida.work.workchain import if_
import numpy as np
from aiida_phonopy.common.generate_inputs import generate_inputs
from aiida_phonopy.common.utils import (phonopy_bulk_from_structure,
                                        phonopy_atoms_to_structure,
                                        get_nac_from_data,
                                        get_path_using_seekpath)

# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows
# from aiida.workflows.wc_optimize import OptimizeStructure

ForceConstantsData = DataFactory('phonopy.force_constants')
ForceSetsData = DataFactory('phonopy.force_sets')
BandStructureData = DataFactory('phonopy.band_structure')
PhononDosData = DataFactory('phonopy.phonon_dos')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')
OptimizeStructure = WorkflowFactory('phonopy.optimize')


class PhononPhonopy(WorkChain):
    """ Workchain to do a phonon calculation using phonopy

    :param structure: StructureData object that contains the crystal structure
        unit cell
    :param ph_settings: ParametersData object that contains a dictionary with
        the data needed to run phonopy:
            'supercell_matrix': [[2,0,0],
                                 [0,2,0],
                                 [0,0,2]],
            'primitive_matrix': [[1.0, 0.0, 0.0],
                                 [0.0, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]],
            'distance': 0.01,
            'mesh': [40, 40, 40],
            # 'code': 'phonopy@boston'
        'code' is include to run phonopy remotely otherwise run phonopy localy.
    :param es_settings: ParametersData object that contains a dictionary with
        the setting needed to calculate the electronic structure. The
        structure of this dictionary strongly depends on the software
        (VASP, QE, LAMMPS, ...)
    :param optimize: Set true to perform a crystal structure optimization
        before the phonon calculation (default: True)
    :param pressure: Set the external pressure (stress tensor) at which the
        optimization is performed in KBar (default: 0)
    """
    @classmethod
    def define(cls, spec):
        super(PhononPhonopy, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("optimize",
                   valid_type=Bool, required=False, default=Bool(True))
        spec.input("pressure",
                   valid_type=Float, required=False, default=Float(0.0))
        spec.input("is_nac",
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input("run_phonopy",
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input("remote_phonopy",
                   valid_type=Bool, required=False, default=Bool(False))

        spec.outline(cls.check_ph_settings,
                     if_(cls.use_optimize)(cls.optimize),
                     cls.set_unitcell,
                     cls.create_displacements,
                     cls.run_force_and_nac_calculations,
                     cls.create_force_sets,
                     if_(cls.is_nac)(cls.create_nac_params),
                     cls.set_default_outputs,
                     cls.prepare_phonopy_inputs,
                     if_(cls.run_phonopy)(
                         if_(cls.remote_phonopy)(
                             cls.create_phonopy_builder,
                             cls.run_phonopy_remote,
                             cls.collect_data,
                         ).else_(
                             cls.run_phonopy_in_workchain,
                         )))

    def use_optimize(self):
        return self.inputs.optimize

    def remote_phonopy(self):
        return self.inputs.remote_phonopy

    def run_phonopy(self):
        return self.inputs.run_phonopy

    def is_nac(self):
        return self.inputs.is_nac

    def check_ph_settings(self):
        self.report('check ph_settings')

        ph_settings_dict = self.inputs.ph_settings.get_dict()

        if 'supercell_matrix' not in ph_settings_dict:
            raise RuntimeError("supercell_matrix is not set.")
        if 'primitive_matrix' not in ph_settings_dict:
            ph_settings_dict['primitive_matrix'] = np.eye(3, dtype='double')
        if 'distance' not in ph_settings_dict:
            ph_settings_dict['distance'] = 0.01
        if 'mesh' not in ph_settings_dict:
            ph_settings_dict['mesh'] = [20, 20, 20]
        if 'symmetry_precision' not in ph_settings_dict:
            ph_settings_dict['symmetry_precision'] = 1e-5
        if 'code' not in ph_settings_dict and self.inputs.remote_phonopy:
            raise RuntimeError("code has to be specified in ph_settings.")
        if 'code' in ph_settings_dict and 'options' not in ph_settings_dict:
            raise RuntimeError(
                "options of code has to be specified in ph_settings.")
        self.ctx.ph_settings = ParameterData(dict=ph_settings_dict)

    def optimize(self):
        self.report('start optimize')
        future = self.submit(OptimizeStructure,
                             structure=self.inputs.structure,
                             es_settings=self.inputs.es_settings,
                             pressure=self.inputs.pressure)
        self.report('optimize workchain: {}'.format(future.pk))

        return ToContext(optimized=future)

    def set_unitcell(self):
        self.report('set unit cell')

        if 'optimized' in self.ctx:
            self.ctx.final_structure = self.ctx.optimized.out.optimized_structure
            self.out('optimized_data',
                     self.ctx.optimized.out.optimized_structure_data)
        else:
            self.ctx.final_structure = self.inputs.structure

    def create_displacements(self):
        """Create supercells with displacements and primitive cell

        Use phonopy to create the supercells with displacements to
        calculate the force constants by using finite displacements
        methodology

        """

        self.report('create cells for phonon calculation')

        from phonopy import Phonopy

        ph_settings_dict = self.ctx.ph_settings.get_dict()

        # Generate phonopy phonon object
        phonon = Phonopy(
            phonopy_bulk_from_structure(self.ctx.final_structure),
            supercell_matrix=ph_settings_dict['supercell_matrix'],
            primitive_matrix=ph_settings_dict['primitive_matrix'],
            symprec=ph_settings_dict['symmetry_precision'])

        self.ctx.primitive_structure = phonopy_atoms_to_structure(
            phonon.primitive)

        phonon.generate_displacements(distance=ph_settings_dict['distance'])
        cells_with_disp = phonon.supercells_with_displacements
        disp_dataset = {'data_sets': phonon.displacement_dataset,
                        'number_of_displacements': Int(len(cells_with_disp))}
        for i, phonopy_supercell in enumerate(cells_with_disp):
            supercell = phonopy_atoms_to_structure(phonopy_supercell)
            disp_dataset["structure_{}".format(i)] = supercell
        self.ctx.disp_dataset = disp_dataset

    def run_force_and_nac_calculations(self):
        self.report('run force calculations')

        # Forces
        for i in range(self.ctx.disp_dataset['number_of_displacements']):
            label = "structure_{}".format(i)
            supercell = self.ctx.disp_dataset[label]
            builder = generate_inputs(supercell,
                                      self.inputs.es_settings,
                                      calc_type='forces')
            builder.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

        # Born charges
        if bool(self.inputs.is_nac):
            self.report('calculate born charges and dielectric constant')
            builder = generate_inputs(self.ctx.primitive_structure,
                                      self.inputs.es_settings,
                                      calc_type='nac')
            future = self.submit(builder)
            self.report('born_and_epsilon: {}'.format(future.pk))
            self.to_context(**{'born_and_epsilon': future})

    def create_force_sets(self):
        """Build data_sets from forces of supercells with displacments"""

        self.report('create force sets')

        forces = []
        for i in range(self.ctx.disp_dataset['number_of_displacements']):
            out_dict = self.ctx['structure_{}'.format(i)].get_outputs_dict()
            if ('output_forces' in out_dict and
                'final' in out_dict['output_forces'].get_arraynames()):
                forces.append(out_dict['output_forces'].get_array('final'))

        if len(forces) != self.ctx.disp_dataset['number_of_displacements']:
            raise RuntimeError("Forces could not be retrieved.")

        data_sets = self.ctx.disp_dataset['data_sets']
        self.ctx.force_sets = ForceSetsData(data_sets=data_sets)
        self.ctx.force_sets.set_forces(forces)

    def create_nac_params(self):
        self.ctx.nac_data = get_nac_from_data(
            born_charges=self.ctx.born_and_epsilon.out.output_born_charges,
            epsilon=self.ctx.born_and_epsilon.out.output_dielectrics,
            structure=self.ctx.born_and_epsilon.inp.structure)

    def set_default_outputs(self):
        self.out('final_structure', self.ctx.final_structure)
        self.out('force_sets', self.ctx.force_sets)
        if 'nac_data' in self.ctx:
            self.out('nac_data', self.ctx.nac_data)

    def prepare_phonopy_inputs(self):
        self.ctx.band_paths = get_path_using_seekpath(
            self.ctx.primitive_structure,
            band_resolution=30)

    def create_phonopy_builder(self):
        """Generate input parameters for a remote phonopy calculation"""

        self.report("create phonopy process builder")

        ph_settings_dict = self.ctx.ph_settings.get_dict()
        code = Code.get_from_string(ph_settings_dict['code'])
        builder = code.get_builder()

        builder.structure = self.ctx.final_structure
        builder.parameters = self.ctx.ph_settings
        options = ph_settings_dict['options']
        builder.options.resources = options['resources']
        max_wallclock_seconds = options['max_wallclock_seconds']
        builder.options.max_wallclock_seconds = max_wallclock_seconds
        if 'force_constants' in self.ctx:
            builder.force_constants = self.force_constants
        elif 'force_sets' in self.ctx:
            builder.force_sets = self.ctx.force_sets
        if 'nac_data' in self.ctx:
            builder.nac_data = self.ctx.nac_data
        if 'band_paths' in self.ctx:
            builder.bands = self.ctx.band_paths

        if 'force_sets' not in self.ctx and 'force_constants' not in self.ctx:
            raise RuntimeError(
                "Force sets and force constants were not found.")

        self.ctx.phonopy_builder = builder

    def run_phonopy_remote(self):
        if 'code' in self.ctx.ph_settings.get_dict():
            self.report('remote phonopy calculation')
            future = self.submit(self.ctx.phonopy_builder)
            self.report('phonopy calculation: {}'.format(future.pk))
            return ToContext(phonon_properties=future)
        else:
            raise RuntimeError("Code for phonopy is not set.")

    def run_phonopy_in_workchain(self):
        self.report('phonopy calculation in workchain')

        ph_settings_dict = self.ctx.ph_settings.get_dict()

        from phonopy import Phonopy

        phonon = Phonopy(
            phonopy_bulk_from_structure(self.ctx.final_structure),
            supercell_matrix=ph_settings_dict['supercell_matrix'],
            primitive_matrix=ph_settings_dict['primitive_matrix'],
            symprec=ph_settings_dict['symmetry_precision'])

        if 'force_constants' in self.ctx:
            force_constants = self.ctx.force_constants
            phonon.set_force_constants(force_constants.get_data())
        else:
            phonon.set_displacement_dataset(
                self.ctx.force_sets.get_force_sets())
            phonon.produce_force_constants()
            force_constants = ForceConstantsData(
                data=phonon.get_force_constants())

        if 'nac_data' in self.ctx:
            primitive = phonon.get_primitive()
            nac_parameters = self.ctx.nac_data.get_born_parameters_phonopy(
                primitive_cell=primitive.get_cell())
            phonon.set_nac_params(nac_parameters)

        # Normalization factor primitive to unit cell
        normalization_factor = (phonon.unitcell.get_number_of_atoms()
                                / phonon.primitive.get_number_of_atoms())

        # DOS
        phonon.set_mesh(ph_settings_dict['mesh'],
                        is_eigenvectors=True,
                        is_mesh_symmetry=False)
        phonon.set_total_DOS(tetrahedron_method=True)
        phonon.set_partial_DOS(tetrahedron_method=True)

        total_dos = phonon.get_total_DOS()
        partial_dos = phonon.get_partial_DOS()
        dos = PhononDosData(
            frequencies=total_dos[0],
            dos=total_dos[1]*normalization_factor,
            partial_dos=np.array(partial_dos[1]) * normalization_factor,
            atom_labels=np.array(phonon.primitive.get_chemical_symbols()))

        # THERMAL PROPERTIES (per primtive cell)
        phonon.set_thermal_properties()
        t, free_energy, entropy, cv = phonon.get_thermal_properties()

        # Stores thermal properties (per unit cell) data in DB as a workflow
        # result
        thermal_properties = ArrayData()
        thermal_properties.set_array('temperature', t)
        thermal_properties.set_array('free_energy',
                                     free_energy * normalization_factor)
        thermal_properties.set_array('entropy', entropy * normalization_factor)
        thermal_properties.set_array('heat_capacity',
                                     cv * normalization_factor)

        # BAND STRUCTURE
        bands = self.ctx.band_paths
        phonon.set_band_structure(bands.get_bands())
        band_structure = BandStructureData(bands=bands.get_bands(),
                                           labels=bands.get_labels(),
                                           unitcell=bands.get_unitcell())

        band_structure.set_band_structure_phonopy(phonon.get_band_structure())

        self.out('force_constants', force_constants)
        self.out('thermal_properties', thermal_properties)
        self.out('dos', dos)
        self.out('band_structure', band_structure)

        self.report('finish phonon')

    def collect_data(self):
        self.out('force_constants',
                 self.ctx.phonon_properties.out.force_constants)
        self.out('thermal_properties',
                 self.ctx.phonon_properties.out.thermal_properties)
        self.out('dos', self.ctx.phonon_properties.out.dos)
        self.out('band_structure',
                 self.ctx.phonon_properties.out.band_structure)

        self.report('finish phonon')
