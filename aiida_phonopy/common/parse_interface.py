def parse_optimize_calculation(calc):
    """
    Parse ths information from plugins nodes and set common units
    Stress in kB
    Force in eV/Angstrom
    """

    plugin = calc.get_code().get_attr('input_plugin')

    if plugin == 'vasp.vasp':
        forces = calc.out.output_array.get_array('forces')
        stresses = calc.out.output_array.get_array('stress')

    elif plugin == 'lammps.optimize':
        forces = calc.out.output_array.get_array('forces')
        stresses = calc.out.output_array.get_array('stress')

    elif plugin == 'quantumespresso.pw':
        forces = calc.out.output_trajectory.get_array('forces')[-1]
        stresses = calc.out.output_trajectory.get_array('stress')[-1] * 10
    else:
        return Exception('Not supported plugin')

    return {'forces': forces, 'stresses': stresses}