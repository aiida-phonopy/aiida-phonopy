
def parse_optimize_calculation(calc):

    plugin = calc.get_code().get_attr('input_plugin')

    if plugin == 'vasp.vasp':
        forces = calc.out.output_array.get_array('forces')
        stresses = calc.out.output_array.get_array('stress')

    elif plugin == 'lammps.optimize':
        forces = calc.out.output_array.get_array('forces')
        stresses = calc.out.output_array.get_array('stress')

    elif plugin == 'quantumespresso.pw':
        output_trajectory = calc.output_trajectory
        forces = output_trajectory.get_array('forces')[-1]
        stresses = output_trajectory.get_array('stress')[-1]
    else:
        return Exception('Not supported plugin')

    return {'forces': forces, 'stresses': stresses}