"""General utilities."""

import numpy as np
from aiida.engine import calcfunction
from aiida.plugins import DataFactory
from aiida.orm import Float, Bool, Str, Int, load_node
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import (get_default_displacement_distance,
                                          get_default_physical_units)
from phonopy import Phonopy


Dict = DataFactory('dict')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
StructureData = DataFactory('structure')
BandsData = DataFactory('array.bands')


@calcfunction
def setup_phonopy_calculation(phonon_settings,
                              structure,
                              symmetry_tolerance,
                              run_phonopy,
                              dataset=None):
    """Set up phonopy calculation.

    Valid keys in phonon_settings_info are
        ('supercell_matrix',
         'phonon_supercell_matrix',
         'distance',
         'symmetry_tolerance',
         'number_of_snapshots',
         'random_seed',
         'is_plusminus',
         'is_diagonal',
         'is_trigonal',
         'mesh',
         'fc_calculator').

    Returns
    -------
    dict
        'supercell' : StructureData
            Perfect supercell.
        'supercell_001', 'supercell_002', ... : StructureData
            Supercells with displacements
        'primitive' : StructureData
            Primitive cell.
        'phonon_supercell' : StructureData
            For phono3py. Perfect supercell for harmonic phonon calculation.
        'phonon_supercell_001', 'phonon_supercell_002', ... : StructureData
            For phono3py. Supercells with displacements for harmonic phonon
            calculation.
        'phonon_setting_info' : Dict
            Phonopy setting parameters including those generated in the
            process of displacements creation, e.g., primitive and  sueprcell
            matrix and symmetry information.

    phonon_setting_info contains the following entries:
        'version' : str
            Phonopy version number.
        'supercell_matrix' : array_like
            3x3 integer matrix to generate supercell matrix.
        'phonon_supercell_matrix' : array_like
            3x3 integer matrix to generate fc2 supercell matrix for Phono3py.
        'distance' : float
            Displacement distance.
        'symmetry_tolerance' : float
            Tolerance length used for symmetry finding.
        'displacement_dataset' : dict
            Phonopy.dataset or Phono3py.dataset.
        'primitive_matrix' : array_like
            Phonopy.primitive_matrix.
        'symmetry' : dict
            'number' : Space group number.
            'international' : Space group type.
        'phonon_displacement_dataset' : dict
            Phono3py.phonon_dataset.
        'number_of_snapshots' : int
        'random_seed' : int
        'is_plusminus' : bool
        'is_diagonal' : bool
        'is_trigonal' : bool
        'postprocess_parameters' : dict
            This is given when run_phonopy=True.
            'mesh' : float or list
                Mesh numbers or distance measure of q-point sampling mesh.
            'fc_calculator' : str
                External force constants calculator.

    """
    ph_settings = _get_setting_info(phonon_settings)
    if run_phonopy:
        params = _get_phonopy_postprocess_info(phonon_settings)
        ph_settings['postprocess_parameters'] = params

    ph = _get_phonopy_instance(structure,
                               ph_settings,
                               symmetry_tolerance=symmetry_tolerance.value)
    ph_settings['version'] = ph.version

    if dataset is None:
        supported_keys = (
            'distance', 'is_plusminus', 'is_diagonal', 'is_trigonal',
            'number_of_snapshots', 'random_seed')
        kwargs = {key: ph_settings[key]
                  for key in ph_settings if key in supported_keys}
        ph.generate_displacements(**kwargs)
    else:
        ph.dataset = dataset.get_dict()

    _update_structure_info(ph_settings, ph)
    structures_dict = _generate_phonon_structures(ph)
    return_vals = {'phonon_setting_info': Dict(dict=ph_settings)}
    return_vals.update(structures_dict)

    return return_vals


def _get_phonopy_postprocess_info(phonon_settings):
    """Return phonopy postprocess parameters."""
    valid_keys = ('mesh', 'fc_calculator')
    params = {}
    for key in valid_keys:
        if key in phonon_settings.keys():
            params[key] = phonon_settings[key]

    if 'mesh' not in phonon_settings.keys():
        params['mesh'] = 100.0
    return params


@calcfunction
def generate_phono3py_cells(phonon_settings,
                            structure,
                            symmetry_tolerance,
                            dataset=None):
    """Set up phono3py calculation.

    Returns
    -------
    dict
        'supercell' : StructureData
            Perfect supercell.
        'supercell_001', 'supercell_002', ... : StructureData
            Supercells with displacements
        'primitive' : StructureData
            Primitive cell.
        'phonon_supercell' : StructureData
            For phono3py. Perfect supercell for harmonic phonon calculation.
        'phonon_supercell_001', 'phonon_supercell_002', ... : StructureData
            For phono3py. Supercells with displacements for harmonic phonon
            calculation.
        'phonon_setting_info' : Dict
            Phonopy setting parameters including those generated in the
            process of displacements creation, e.g., primitive and  sueprcell
            matrix and symmetry information.

    phonon_setting_info contains the following entries:
        'supercell_matrix' : array_like
            3x3 integer matrix to generate supercell matrix.
        'phonon_supercell_matrix' : array_like
            3x3 integer matrix to generate fc2 supercell matrix for Phono3py.
        'distance' : float
            Displacement distance.
        'symmetry_tolerance' : float
            Tolerance length used for symmetry finding.
        'displacement_dataset' : dict
            Phonopy.dataset or Phono3py.dataset.
        'primitive_matrix' : array_like
            Phonopy.primitive_matrix.
        'symmetry' : dict
            'number' : Space group number.
            'international' : Space group type.
        'phonon_displacement_dataset' : dict
            Phono3py.phonon_dataset.
        'is_plusminus' : bool
        'is_diagonal' : bool

    """
    ph_settings = _get_setting_info(phonon_settings, code_name='phono3py')

    ph = _get_phono3py_instance(structure, ph_settings, {})
    if dataset is None:
        ph.generate_displacements(distance=ph_settings['distance'])
    else:
        ph.dataset = dataset.get_dict()

    _update_structure_info(ph_settings, ph)
    structures_dict = _generate_phonon_structures(ph)
    return_vals = {'phonon_setting_info': Dict(dict=ph_settings)}
    return_vals.update(structures_dict)

    return return_vals


@calcfunction
def get_force_constants(structure,
                        phonon_settings,
                        force_sets,
                        symmetry_tolerance):
    """Calculate force constants."""
    phonon = _get_phonopy_instance(structure,
                                   phonon_settings,
                                   symmetry_tolerance=symmetry_tolerance.value)
    phonon.dataset = phonon_settings['displacement_dataset']
    phonon.forces = force_sets.get_array('force_sets')

    if 'fc_calculator' in phonon_settings['postprocess_parameters']:
        postprocess_parameters = phonon_settings['postprocess_parameters']
        if postprocess_parameters['fc_calculator'].lower().strip() == 'alm':
            phonon.produce_force_constants(fc_calculator='alm')
    else:
        phonon.produce_force_constants()
    force_constants = ArrayData()
    force_constants.set_array('force_constants', phonon.force_constants)
    force_constants.set_array('p2s_map', phonon.primitive.p2s_map)
    force_constants.label = 'force_constants'

    return force_constants


@calcfunction
def get_phonon_properties(structure,
                          phonon_settings,
                          force_constants,
                          symmetry_tolerance,
                          nac_params=None):
    """Calculate phonon properties."""
    phonon_settings_dict = phonon_settings.get_dict()
    ph = _get_phonopy_instance(structure,
                               phonon_settings_dict,
                               symmetry_tolerance=symmetry_tolerance.value,
                               nac_params=nac_params)
    ph.force_constants = force_constants.get_array('force_constants')
    mesh = phonon_settings['postprocess_parameters']['mesh']

    # Mesh
    total_dos, pdos, thermal_properties = get_mesh_property_data(ph, mesh)

    # Band structure
    bs = _get_bands_data(ph)

    return {'dos': total_dos,
            'pdos': pdos,
            'thermal_properties': thermal_properties,
            'band_structure': bs}


@calcfunction
def get_data_from_node_id(node_id):
    """Collect VASP output from node_id."""
    n = load_node(node_id.value)
    if 'structure' in n.inputs:
        cell = phonopy_atoms_from_structure(n.inputs.structure)
        structure = phonopy_atoms_to_structure(cell)
    else:
        raise RuntimeError("Crystal structure could not be found.")

    if 'born_charges' in n.outputs and 'dielectrics' in n.outputs:
        born = ArrayData()
        born.set_array(
            'born_charges', n.outputs.born_charges.get_array('born_charges'))
        born.label = 'born_charges'
        epsilon = ArrayData()
        epsilon.set_array(
            'epsilon', n.outputs.dielectrics.get_array('epsilon'))
        epsilon.label = 'epsilon'
        return {'born_charges': born, 'dielectrics': epsilon,
                'structure': structure}
    elif 'forces' in n.outputs:
        forces = ArrayData()
        forces.set_array('final', n.outputs.forces.get_array('final'))
        forces.label = 'forces'
        return {'forces': forces, 'structure': structure}
    else:
        raise RuntimeError("Forces or NAC params were not found.")


def compare_structures(cell_ref, cell_calc, symmetry_tolerance):
    """Compare two PhonopyAtoms instances."""
    symprec = symmetry_tolerance.value
    cell_diff = np.subtract(cell_ref.cell, cell_calc.cell)
    if (np.abs(cell_diff) > symprec).any():
        succeeded = Bool(False)
        succeeded.label = "False"
        return succeeded

    positions_ref = [site.position for site in cell_ref.sites]
    positions_calc = [site.position for site in cell_calc.sites]
    diff = np.subtract(positions_ref, positions_calc)
    diff -= np.rint(diff)
    dist = np.sqrt(np.sum(np.dot(diff, cell_ref.cell) ** 2, axis=1))
    if (dist > symprec).any():
        succeeded = Bool(False)
        succeeded.label = "False"
        return succeeded

    succeeded = Bool(True)
    succeeded.label = "True"
    return succeeded


def get_mesh_property_data(ph, mesh):
    """Return total DOS, PDOS, thermal properties."""
    ph.set_mesh(mesh)
    ph.run_total_dos()

    dos = get_total_dos(ph.get_total_dos_dict())

    ph.run_thermal_properties()
    tprops = get_thermal_properties(ph.get_thermal_properties_dict())

    ph.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)
    ph.run_projected_dos()
    pdos = get_projected_dos(ph.get_projected_dos_dict())

    return dos, pdos, tprops


def get_total_dos(total_dos):
    """Return XyData of total DOS."""
    dos = XyData()
    dos.set_x(total_dos['frequency_points'], 'Frequency', 'THz')
    dos.set_y(total_dos['total_dos'], 'Total DOS', '1/THz')
    dos.label = 'Total DOS'
    return dos


def get_projected_dos(projected_dos):
    """Return XyData of PDOS."""
    pdos = XyData()
    pdos_list = [pd for pd in projected_dos['projected_dos']]
    pdos.set_x(projected_dos['frequency_points'], 'Frequency', 'THz')
    pdos.set_y(pdos_list,
               ['Projected DOS', ] * len(pdos_list),
               ['1/THz', ] * len(pdos_list))
    pdos.label = 'Projected DOS'
    return pdos


def get_thermal_properties(thermal_properties):
    """Return XyData of thermal properties."""
    tprops = XyData()
    tprops.set_x(thermal_properties['temperatures'], 'Temperature', 'K')
    tprops.set_y([thermal_properties['free_energy'],
                  thermal_properties['entropy'],
                  thermal_properties['heat_capacity']],
                 ['Helmholtz free energy', 'Entropy', 'Cv'],
                 ['kJ/mol', 'J/K/mol', 'J/K/mol'])
    tprops.label = 'Thermal properties'
    return tprops


def _get_bands_data(ph):
    ph.auto_band_structure()
    labels = [x.replace('$', '').replace('\\', '').replace('mathrm{', '').replace('}', '').upper()
              for x in ph.band_structure.labels]
    frequencies = ph.band_structure.frequencies
    qpoints = ph.band_structure.qpoints
    path_connections = ph.band_structure.path_connections
    label = "%s (%d)" % (ph.symmetry.dataset['international'],
                         ph.symmetry.dataset['number'])

    return get_bands(qpoints, frequencies, labels, path_connections,
                     label=label)


def get_bands(qpoints, frequencies, labels, path_connections, label=None):
    """Return BandsData."""
    qpoints_list = list(qpoints[0])
    frequencies_list = list(frequencies[0])
    labels_list = [(0, labels[0]), ]
    label_index = 1

    for pc, qs, fs in zip(path_connections[:-1], qpoints[1:], frequencies[1:]):
        if labels[label_index] == 'GAMMA' and pc:
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            if label_index < len(labels):
                labels_list.append((len(qpoints_list), labels[label_index]))
            label_index += 1
            qpoints_list += list(qs)
            frequencies_list += list(fs)
        elif pc:
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            label_index += 1
            qpoints_list += list(qs[1:])
            frequencies_list += list(fs[1:])
        else:
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            label_index += 1
            if label_index < len(labels):
                labels_list.append((len(qpoints_list), labels[label_index]))
                label_index += 1
            qpoints_list += list(qs)
            frequencies_list += list(fs)
    labels_list.append((len(qpoints_list) - 1, labels[-1]))

    bs = BandsData()
    bs.set_kpoints(np.array(qpoints_list))
    bs.set_bands(np.array(frequencies_list), units='THz')
    bs.labels = labels_list
    if label is not None:
        bs.label = label

    return bs


def _get_phonopy_instance(structure,
                          phonon_settings_dict,
                          symmetry_tolerance=1e-5,
                          nac_params=None):
    """Create Phonopy instance."""
    phonon = Phonopy(
        phonopy_atoms_from_structure(structure),
        supercell_matrix=phonon_settings_dict['supercell_matrix'],
        primitive_matrix='auto',
        symprec=symmetry_tolerance)
    if nac_params:
        units = get_default_physical_units('vasp')
        factor = units['nac_factor']
        _nac_params = {'born': nac_params.get_array('born_charges'),
                       'dielectric': nac_params.get_array('epsilon'),
                       'factor': factor}
        phonon.nac_params = _nac_params

    return phonon


def _get_phono3py_instance(structure, phonon_settings_dict, params):
    """Create Phono3py instance."""
    from phono3py import Phono3py
    if 'phonon_supercell_matrix' in phonon_settings_dict:
        ph_smat = phonon_settings_dict['phonon_supercell_matrix']
    else:
        ph_smat = None
    ph3py = Phono3py(
        phonopy_atoms_from_structure(structure),
        supercell_matrix=phonon_settings_dict['supercell_matrix'],
        primitive_matrix='auto',
        phonon_supercell_matrix=ph_smat,
        symprec=phonon_settings_dict['symmetry_tolerance'])
    if 'nac_params' in params:
        from phonopy.interface.calculator import get_default_physical_units
        units = get_default_physical_units('vasp')
        factor = units['nac_factor']
        nac_params = {'born': params['nac_params'].get_array('born_charges'),
                      'dielectric': params['nac_params'].get_array('epsilon'),
                      'factor': factor}
        ph3py.nac_params = nac_params

    return ph3py


def phonopy_atoms_to_structure(cell):
    """Convert PhonopyAtoms to StructureData."""
    symbols = cell.symbols
    positions = cell.positions
    structure = StructureData(cell=cell.cell)
    for symbol, position in zip(symbols, positions):
        structure.append_atom(position=position, symbols=symbol)
    return structure


def phonopy_atoms_from_structure(structure):
    """Convert StructureData to PhonopyAtoms."""
    cell = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    return cell


def from_node_id_to_aiida_node_id(node_id):
    """Convert PK or UUID to an AiiDA data type."""
    if type(node_id) is int:
        return Int(node_id)
    elif type(node_id) is str:
        return Str(node_id)
    else:
        raise RuntimeError("%s is not supported in load_node."
                           % type(node_id))


def collect_forces_and_energies(ctx, ctx_supercells, prefix="force_calc"):
    """Collect forces and energies from calculation outputs.

    Parameters
    ----------
    ctx : AttributeDict-like
        AiiDA workchain context.
    ctx_supercells : dict of StructDict
        Supercells. For phono3py, this can be phonon_supercells.
    prefix : str
        Prefix string of dictionary keys of ctx.

    Returns
    -------
    dict
        Forces and energies.

    """
    forces_dict = {}
    for key in ctx_supercells:
        # key: e.g. "supercell_001", "phonon_supercell_001"
        num = key.split('_')[-1]  # e.g. "001"
        calc = ctx["%s_%s" % (prefix, num)]
        if type(calc) is dict:
            calc_dict = calc
        else:
            calc_dict = calc.outputs
        forces_dict["forces_%s" % num] = calc_dict['forces']
        forces_dict["energy_%s" % num] = calc_dict['energy']

    return forces_dict


@calcfunction
def get_force_sets(**forces_dict):
    """Create force sets from supercell forces.

    Parameters
    ----------
    forces_dict : dict
        'forces_001', 'forces_002', ... have to exist.
        'energy_001', 'energy_002', ... are optional.
        'forces_000' and 'energy_000' for perfect supercell are optional.
        The zero-padding length of the numbers can change depending on total
        number of supercell calculations.

    """
    (force_sets,
     energies,
     forces_0_key,
     energy_0_key) = _get_force_set(**forces_dict)

    force_sets_data = ArrayData()
    force_sets_data.set_array('force_sets', force_sets)
    if energies is not None:
        force_sets_data.set_array('energies', energies)
    force_sets_data.label = 'force_sets'
    ret_dict = {'force_sets': force_sets_data}
    if forces_0_key is not None:
        ret_dict['supercell_forces'] = forces_dict[forces_0_key]
    if energy_0_key is not None:
        ret_dict['supercell_energy'] = forces_dict[energy_0_key]

    return ret_dict


def _get_force_set(**forces_dict):
    num_forces = 0
    num_energies = 0
    forces_0_key = None
    energy_0_key = None
    shape = None
    for key in forces_dict:
        value = forces_dict[key]
        if int(key.split('_')[-1]) != 0:
            if 'forces' in key:
                num_forces += 1
                if shape is None:
                    shape = value.get_array('forces').shape
            elif 'energy' in key:
                num_energies += 1
        else:
            if 'forces' in key:
                forces_0_key = key
            elif 'energy' in key:
                energy_0_key = key

    force_sets = np.zeros((num_forces, ) + shape, dtype=float)
    if num_energies > 0:
        energies = np.zeros(num_energies, dtype=float)
    else:
        energies = None
    if forces_0_key is None:
        forces_0 = None
    else:
        forces_0 = forces_dict[forces_0_key].get_array('forces')

    for key in forces_dict:
        value = forces_dict[key]
        num = int(key.split('_')[-1])  # e.g. "001" --> 1
        if 'forces' in key:
            forces = value.get_array('forces')
            if forces_0 is None:
                force_sets[num - 1] = forces
            else:
                force_sets[num - 1] = forces - forces_0
        elif 'energy' in key:
            energies[num - 1] = value.get_array('energy')

    return force_sets, energies, forces_0_key, energy_0_key


def collect_vasp_forces_and_energies(ctx, ctx_supercells, prefix="force_calc"):
    """Collect forces and energies from calculation outputs.

    Parameters
    ----------
    ctx : AttributeDict-like
        AiiDA workchain context.
    ctx_supercells : dict of StructDict
        Supercells. For phono3py, this can be phonon_supercells.
    prefix : str
        Prefix string of dictionary keys of ctx.

    Returns
    -------
    dict
        Forces and energies.

    """
    forces_dict = {}
    for key in ctx_supercells:
        # key: e.g. "supercell_001", "phonon_supercell_001"
        num = key.split('_')[-1]  # e.g. "001"
        calc = ctx["%s_%s" % (prefix, num)]
        if type(calc) is dict:
            calc_dict = calc
        else:
            calc_dict = calc.outputs
        if ('forces' in calc_dict and
            'final' in calc_dict['forces'].get_arraynames()):
            forces_dict["forces_%s" % num] = calc_dict['forces']
        else:
            raise RuntimeError(
                "Forces could not be found in calculation %s." % num)

        if ('misc' in calc_dict and
            'total_energies' in calc_dict['misc'].keys()):  # needs .keys()
            forces_dict["misc_%s" % num] = calc_dict['misc']

    return forces_dict


@calcfunction
def get_vasp_force_sets_dict(**forces_dict):
    """Create force sets from supercell forces."""
    forces = []
    energies = []
    forces_0 = None
    energy_0 = None

    for key in forces_dict:
        num = int(key.split('_')[-1])
        if num == 0:
            continue
        if 'forces' in key:
            forces.append(None)
        elif 'misc' in key:
            energies.append(None)

    for key in forces_dict:
        num = int(key.split('_')[-1])  # e.g. "001" --> 1
        if 'forces' in key:
            forces_ndarray = forces_dict[key].get_array('final')
            if num == 0:
                forces_0 = forces_ndarray
            else:
                forces[num - 1] = forces_ndarray
        elif 'misc' in key:
            for energy_key in ('energy_extrapolated', ):
                if energy_key in forces_dict[key]['total_energies']:
                    energy = forces_dict[key]['total_energies'][energy_key]
                if num == 0:
                    energy_0 = energy
                else:
                    energies[num - 1] = energy
                break

    if forces_0 is not None:
        for forces_ndarray in forces:
            forces_ndarray -= forces_0

    force_sets = ArrayData()
    force_sets.set_array('force_sets', np.array(forces))
    if energies:
        force_sets.set_array('energies', np.array(energies))
    force_sets.label = 'force_sets'
    ret_dict = {'force_sets': force_sets}
    if forces_0 is not None:
        forces_0_array = ArrayData()
        forces_0_array.set_array('forces', forces_0)
        ret_dict['supercell_forces'] = forces_0_array
    if energy_0 is not None:
        ret_dict['supercell_energy'] = Float(energy_0)

    return ret_dict


def _get_setting_info(phonon_settings, code_name='phonopy'):
    """Convert AiiDA inputs to a dict.

    code_name : 'phonopy' or 'phono3py'

    Note
    ----
    Designed to be shared by phonopy and phono3py.

    Returns
    -------
    dict
        'supercell_matrix' : ndarray
            3x3 integer matrix to generate supercell matrix.
        'phonon_supercell_matrix' : ndarray
            3x3 integer matrix to generate fc2 supercell matrix for Phono3py.
        'mesh' : list or float
            Mesh numbers or a distance to represent mesh numbers.
        'distance' : float
            Displacement distance.
        'symmetry_tolerance' : float
            Tolerance length used for symmetry finding.

    """
    ph_settings = {}
    valid_keys = ('supercell_matrix',
                  'phonon_supercell_matrix',
                  'distance',
                  'symmetry_tolerance',
                  'number_of_snapshots',
                  'random_seed',
                  'is_plusminus',
                  'is_diagonal',
                  'is_trigonal')
    for key, value in phonon_settings.get_dict().items():
        if key in valid_keys:
            ph_settings[key] = value
    dim = ph_settings['supercell_matrix']
    ph_settings['supercell_matrix'] = _get_supercell_matrix(dim)
    if 'phonon_supercell_matrix' in ph_settings:
        dim_fc2 = ph_settings['phonon_supercell_matrix']
        ph_settings['phonon_supercell_matrix'] = _get_supercell_matrix(
            dim_fc2, smat_type='phonon_supercell_matrix')
    if 'distance' not in ph_settings:
        if code_name == 'phono3py':
            from phono3py.interface.calculator import (
                get_default_displacement_distance as get_phono3py_ddd)
            distance = get_phono3py_ddd('vasp')
        else:
            distance = get_default_displacement_distance('vasp')
        ph_settings['distance'] = distance

    return ph_settings


def _get_supercell_matrix(dim, smat_type='supercell_matrix'):
    if len(np.ravel(dim)) == 3:
        smat = np.diag(dim)
    else:
        smat = np.array(dim)
    if not np.issubdtype(smat.dtype, np.integer):
        raise TypeError("%s is not integer matrix." % smat_type)
    else:
        return smat.tolist()


def _update_structure_info(ph_settings, ph):
    """Update a phonon_settings dict.

    Parameters
    ----------
    ph_settings : Dict
         Phonopy setting information.
    ph : Phonopy or Phono3py
         A Phonopy or Phono3py instance.

    Returns
    -------
    dict
        'displacement_dataset' : dict
            Phonopy.dataset or Phono3py.dataset.
        'primitive_matrix' : ndarray
            Phonopy.primitive_matrix.
        'symmetry' : dict
            'number' : Space group number.
            'international' : Space group type.
        'phonon_displacement_dataset' : dict
            Phono3py.phonon_dataset.

    """
    ph_settings['displacement_dataset'] = ph.dataset
    ph_settings['primitive_matrix'] = ph.primitive_matrix
    ph_settings['symmetry'] = {
        'number': ph.symmetry.dataset['number'],
        'international': ph.symmetry.dataset['international']}

    if 'phonon_supercell_matrix' in ph.__dir__():
        if ph.phonon_supercell_matrix is not None:
            ph_settings['phonon_displacement_dataset'] = ph.phonon_dataset

    return ph_settings


def _generate_phonon_structures(ph):
    """Generate AiiDA structures of phonon related cells.

    Note
    ----
    Designed to be shared by phonopy and phono3py.
    ph is either an instance of Phonopy or Phono3py.

    Returns
    -------
    dict of StructureData
        'supercell'
            Perfect supercell.
        'supercell_001', 'supercell_002', ...
            Supercells with displacements
        'primitive':
            Primitive cell.
        'phonon_supercell'
            For phono3py. Perfect supercell for harmonic phonon calculation.
        'phonon_supercell_001', 'phonon_supercell_002', ...
            For phono3py. Supercells with displacements for harmonic phonon
            calculation.

    """
    structures_dict = {}

    digits = len(str(len(ph.supercells_with_displacements)))
    for i, scell in enumerate(ph.supercells_with_displacements):
        structure = phonopy_atoms_to_structure(scell)
        label = "supercell_%s" % str(i + 1).zfill(digits)
        structure.label = "%s %s" % (
            structure.get_formula(mode='hill_compact'), label)
        structures_dict[label] = structure

    supercell_structure = phonopy_atoms_to_structure(ph.supercell)
    supercell_structure.label = "%s %s" % (
        supercell_structure.get_formula(mode='hill_compact'), 'supercell')
    structures_dict['supercell'] = supercell_structure

    primitive_structure = phonopy_atoms_to_structure(ph.primitive)
    primitive_structure.label = "%s %s" % (
        primitive_structure.get_formula(mode='hill_compact'), 'primitive cell')
    structures_dict['primitive'] = primitive_structure

    # phono3py
    if 'phonon_supercell_matrix' in ph.__dir__():
        if ph.phonon_supercell_matrix is not None:
            digits = len(str(len(ph.phonon_supercells_with_displacements)))
            for i, scell in enumerate(ph.phonon_supercells_with_displacements):
                structure = phonopy_atoms_to_structure(scell)
                label = "phonon_supercell_%s" % str(i + 1).zfill(digits)
                structure.label = "%s %s" % (
                    structure.get_formula(mode='hill_compact'), label)
                structures_dict[label] = structure
            structure = phonopy_atoms_to_structure(ph.phonon_supercell)
            structure.label = "%s %s" % (
                structure.get_formula(mode='hill_compact'), 'phonon_supercell')
            structures_dict['phonon_supercell'] = structure

    return structures_dict
