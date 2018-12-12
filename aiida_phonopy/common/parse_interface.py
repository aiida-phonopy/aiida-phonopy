from aiida.orm import DataFactory
from aiida.work import workfunction
from aiida.orm.data.base import Int

StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')


@workfunction
def structure_from_trajectory(output_trajectory, pos):
    """Worfunction to extract a structure object from trajectory node and
    keep the provenance

    :param output_trajectory: trajectory node from a QE optimization
    :param pos: number of structure to extract from trajectory node
    :return:

    """
    pos = int(pos)

    # Check maximum
    num_structures = len(output_trajectory.get_array('positions'))
    pos = max(pos, num_structures)

    positions = output_trajectory.get_array('positions')[pos]
    symbols = output_trajectory.get_array('symbols')
    cell = output_trajectory.get_array('cells')[pos]

    structure = StructureData(cell=cell.tolist())
    for i, scaled_position in enumerate(positions):
        structure.append_atom(position=scaled_position.tolist(),
                              symbols=symbols[i])
    return {'structure': structure}


def parse_optimize_calculation(calc):
    """Parse ths information from plugins nodes and set common units

    Stress in kB
    Force in eV/Angstrom

    """

    import numpy as np

    plugin = calc.get_code().get_attr('input_plugin')

    if plugin == 'vasp.vasp':
        forces = calc.out.output_trajectory.get_array('forces')[-1]
        stress = calc.out.output_trajectory.get_array('stress')[-1]

        try:
            structure = calc.out.output_structure
        except:
            structure = structure_from_trajectory(calc.out.output_trajectory,
                                                  Int(-1))['structure']

        energy_wo_entrop = calc.out.output_trajectory.get_array('e_wo_entrp')[-1]
        pressure = np.average(np.diag(stress))
        factor = 0.0006241509125883258  # kBar * A^3 -> eV
        energy = (energy_wo_entrop
                  - structure.get_cell_volume() * pressure * factor)

    elif plugin == 'lammps.optimize':
        forces = calc.out.output_array.get_array('forces')
        stress = calc.out.output_array.get_array('stress')
        structure = calc.out.output_structure
        energy = calc.out.output_parameters.dict.energy

    elif plugin == 'quantumespresso.pw':
        forces = calc.out.output_trajectory.get_array('forces')[-1]
        # GPa to kBar
        stress = calc.out.output_trajectory.get_array('stress')[-1] * 10
        energy = calc.out.output_parameters.dict.energy

        try:
            structure = calc.out.output_structure
        except:
            structure = structure_from_trajectory(
                calc.out.output_trajectory, Int(-1))['structure']

    else:
        return Exception('Not supported plugin')

    output_data = ParameterData(dict={'energy': energy,
                                      'forces': forces.tolist(),
                                      'stress': stress.tolist()})

    return {'output_structure': structure, 'output_data': output_data}
