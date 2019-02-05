# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction

from aiida.orm import Code, CalculationFactory, load_node, DataFactory, WorkflowFactory
from aiida.work.run import run, submit

from aiida.orm.data.base import Str, Float, Bool, Int

import numpy as np

ForceConstantsData = DataFactory('phonopy.force_constants')
ForceSetsData = DataFactory('phonopy.force_sets')
PhononDosData = DataFactory('phonopy.phonon_dos')
NacData = DataFactory('phonopy.nac')

ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

__testing__ = False


def generate_phono3py_params(structure,
                             parameters,
                             force_sets,
                             nac_data=None,
                             fc2=None,
                             fc3=None,
                             grid_point=None,
                             grid_data=None):
    """
    Generate inputs parameters needed to do a remote phonopy calculation

    :param structure: StructureData Object that constains the crystal structure unit cell
    :param parameters: ParametersData object containing a dictionary with the phonopy input data
    :param force_sets: ForceSetsData object containing the atomic forces and displacement information
    :param nac_data: NacData object containing the dielectric tensor and Born effective charges info
    :param fc2: ForceConstantsData object containing the 2nd order force constants
    :param fc3: ForceConstantsData object containing the 3rd order force constants
    :param grid_point: List containing the grid points to calculate (in distributed calculation)
    :return: Calculation process object, input dictionary
    """

    try:
        code = Code.get_from_string(parameters.dict.code['fc3'])
    except :
        code = Code.get_from_string(parameters.dict.code)

    plugin = code.get_attr('input_plugin')
    PhonopyCalculation = CalculationFactory(plugin)

    # The inputs
    inputs = PhonopyCalculation.process().get_inputs_template()

    # code
    inputs.code = code

    # structure
    inputs.structure = structure

    # Parameters
    if grid_point is not None:
        parameters_dic = parameters.get_dict()
        parameters_dic.update({'grid_point': np.array(grid_point).tolist()})
        parameters = ParameterData(dict=parameters_dic)

    if grid_data is not None:
        inputs.grid_data = grid_data

    inputs.parameters = parameters

    # resources
    inputs._options.resources = parameters.dict.machine['resources']
    inputs._options.max_wallclock_seconds = parameters.dict.machine['max_wallclock_seconds']

    # data_sets & force constants
    if force_sets is not None:
        inputs.data_sets = force_sets
    if fc2 is not None:
        inputs.force_constants = fc2
    if fc3 is not None:
        inputs.force_constants_3 = fc3

    # non-analytical corrections
    if nac_data is not None:
        inputs.nac_data = nac_data

    return PhonopyCalculation.process(), inputs


def get_grid_points(structure, parameters):

    from phono3py.phonon3 import Phono3py
    from phono3py.cui.triplets_info import get_coarse_ir_grid_points
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    # Generate phonopy phonon object
    phono3py = Phono3py(phonopy_atoms_from_structure(structure),
                        supercell_matrix=parameters.dict.supercell,
                        primitive_matrix=parameters.dict.primitive,
                        symprec=parameters.dict.symmetry_tolerance,
                        log_level=1)

    primitive = phono3py.get_primitive()
    (ir_grid_points,
     grid_weights,
     bz_grid_address,
     grid_mapping_table) = get_coarse_ir_grid_points(
         primitive,
         parameters.dict.mesh,
         None,
         None,
         is_kappa_star=True,
         symprec=1e-5)

    return ir_grid_points


class Phono3pyDist(WorkChain):

    @classmethod
    def define(cls, spec):
        super(Phono3pyDist, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("parameters", valid_type=ParameterData)
        # Optional arguments
        spec.input("force_constants", valid_type=ForceConstantsData, required=False)
        spec.input("force_constants_3", valid_type=ForceConstantsData, required=False)
        spec.input("force_sets", valid_type=ForceSetsData)
        spec.input("nac", valid_type=NacData, required=False)
        spec.input("gp_chunks", valid_type=Int, required=False, default=Int(10))

        spec.outline(cls.create_grid_points,
                     cls.calculate_thermal_conductivity,
                     cls.collect_data)

    def create_grid_points(self):

        print('start phono3py distributed (pk={})'.format(self.pid))

        if 'nac' in self.inputs:
            nac_data = self.ctx.inputs.nac
        else:
            nac_data = None

        if __testing__:
            self.ctx._content['gp_0'] = load_node(85570)
            self.ctx._content['gp_1'] = load_node(85572)
            self.ctx._content['gp_2'] = load_node(85574)
            self.ctx._content['gp_3'] = load_node(85576)
            self.ctx._content['gp_4'] = load_node(85578)
            self.ctx._content['gp_5'] = load_node(85580)
            self.ctx.labels = [0, 1, 2, 3, 4, 5]
            return

        grid_points = get_grid_points(self.inputs.structure, self.inputs.parameters)

        def chunks(l, m):
            n = len(l)/m
            for i in range(0, len(l), n):
                yield l[i:i + n]

        calcs = {}
        self.ctx.labels = []
        for label, gp_range in enumerate(chunks(grid_points, int(self.inputs.gp_chunks))):
            JobCalculation, calculation_input = generate_phono3py_params(structure=self.inputs.structure,
                                                                         parameters=self.inputs.parameters,
                                                                         force_sets=self.inputs.force_sets,
                                                                         nac_data=nac_data,
                                                                         grid_point=gp_range)

            future = submit(JobCalculation, **calculation_input)
            print ('gp_{} pk = {}'.format(label, future.pid))
            calcs['gp_{}'.format(label)] = future
            self.ctx.labels.append(label)

        return ToContext(**calcs)

    def calculate_thermal_conductivity(self):

        data_gp = ArrayData()
        for i in self.ctx.labels:
            calc = self.ctx.get('gp_{}'.format(i))
            print ('collect gp_{}'.format(i))
            outputs_dict = calc.get_outputs_dict()

            for key in outputs_dict:
                if key.startswith('kappa') and len(key.split('_')) == 2:
                    num = key.split('_')[1]
                    # print 'num', num
                    for array_name in outputs_dict[key].get_arraynames():
                        array = outputs_dict[key].get_array(array_name)
                        data_gp.set_array(array_name+'_'+num, array)

        if 'nac' in self.inputs:
            nac_data = self.ctx.inputs.nac
        else:
            nac_data = None

        JobCalculation, calculation_input = generate_phono3py_params(structure=self.inputs.structure,
                                                                     parameters=self.inputs.parameters,
                                                                     force_sets=self.inputs.force_sets,
                                                                     nac_data=nac_data,
                                                                     grid_data=data_gp)

        future = submit(JobCalculation, **calculation_input)
        return ToContext(thermal_conductivity=future)

    def collect_data(self):
        self.out('kappa', self.ctx.thermal_conductivity.out.kappa)

if __name__ == '__main__':
    # Silicon structure
    a = 5.404
    cell = [[a, 0, 0],
            [0, a, 0],
            [0, 0, a]]

    symbols=['Si'] * 8
    scaled_positions = [(0.875,  0.875,  0.875),
                        (0.875,  0.375,  0.375),
                        (0.375,  0.875,  0.375),
                        (0.375,  0.375,  0.875),
                        (0.125,  0.125,  0.125),
                        (0.125,  0.625,  0.625),
                        (0.625,  0.125,  0.625),
                        (0.625,  0.625,  0.125)]

    structure = StructureData(cell=cell)
    for i, scaled_position in enumerate(scaled_positions):
        structure.append_atom(position=np.dot(scaled_position, cell).tolist(),
                              symbols=symbols[i])

    # Machine
    machine_dict = {'resources': {'num_machines': 1,
                                  'parallel_env': 'mpi*',
                                  'tot_num_mpiprocs': 16},
                    'max_wallclock_seconds': 3600 * 10,
                    }

    machine = ParameterData(dict=machine_dict)

    # PHONOPY settings
    ph_settings = ParameterData(dict={'supercell': [[2, 0, 0],
                                                    [0, 2, 0],
                                                    [0, 0, 2]],
                                      'primitive': [[0.0, 0.5, 0.5],
                                                    [0.5, 0.0, 0.5],
                                                    [0.5, 0.5, 0.0]],
                                      'distance': 0.01,
                                      'mesh': [20, 20, 20],
                                      'symmetry_tolerance': 1e-5,
                                      'code': 'phono3py@stern_in',
                                      'machine': machine_dict})

    # Chose how to run the calculation
    result = run(Phono3pyDist,
                 structure=structure,
                 parameters=ph_settings,
                 force_sets=load_node(81481),  # load previous data
                 gp_chunks=Int(8)
                 )

    print(result)
