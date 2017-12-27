from aiida.common.exceptions import InputValidationError
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation


BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')

class Phono3pyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INPUT_DISP_FC3 = 'disp_fc3.yaml'
    _INPUT_FORCES_FC3 = 'FORCES_FC3'

    def _init_internal_params(self):
        super(Phono3pyCalculation, self)._init_internal_params()

        self._default_parser = 'phono3py'

        self._internal_retrieve_list += [self._OUTPUT_DOS,
                                         self._OUTPUT_THERMAL_PROPERTIES,
                                         self._OUTPUT_BAND_STRUCTURE]

        self._calculation_cmd = ['--fc3', '--fc2', '--br']


    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods
        retdict.update(BasePhonopyCalculation._baseclass_use_methods)

        return retdict

    def _create_additional_files(self, tempfolder, inputdict):

        data_sets = inputdict.pop(self.get_linkname('data_sets'), None)
        force_constants = inputdict.pop(self.get_linkname('force_constants'), None)

        if data_sets is None and force_constants is None:
            raise InputValidationError("no force_sets nor force_constants are specified for this calculation")

        parameters_data = inputdict[self.get_linkname('parameters')]
        structure = inputdict[self.get_linkname('structure')]

        if data_sets is not None:
            disp_fc3_filename = tempfolder.get_abs_path(self._INPUT_DISP_FC3)
            disp_f3c_txt = get_disp_fc3_txt(structure, parameters_data, data_sets)
            with open(disp_fc3_filename, 'w') as infile:
                infile.write(disp_f3c_txt)

            forces_txt = get_forces_txt(force_constants)
            forces_filename = tempfolder.get_abs_path(self._INPUT_FORCES_FC3)
            with open(forces_filename, 'w') as infile:
                infile.write(forces_txt)


def get_disp_fc3_txt(structure, parameters_data, dataset):
    from phono3py.file_IO import write_cell_yaml
    from phonopy.structure.cells import get_supercell
    from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure

    import StringIO

    supercell = get_supercell(phonopy_bulk_from_structure(structure),
                              parameters_data.dict.supercell,
                              symprec=parameters_data.dict.symmetry_precision)

    w = StringIO.StringIO()
    w.write("natom: %d\n" %  dataset['natom'])

    num_first = len(dataset['first_atoms'])
    w.write("num_first_displacements: %d\n" %  num_first)
    if 'cutoff_distance' in dataset:
        w.write("cutoff_distance: %f\n" %  dataset['cutoff_distance'])

    num_second = 0
    num_disp_files = 0
    for d1 in dataset['first_atoms']:
        num_disp_files += 1
        num_second += len(d1['second_atoms'])
        for d2 in d1['second_atoms']:
            if 'included' in d2:
                if d2['included']:
                    num_disp_files += 1
            else:
                num_disp_files += 1

    w.write("num_second_displacements: %d\n" %  num_second)
    w.write("num_displacements_created: %d\n" %  num_disp_files)

    w.write("first_atoms:\n")
    count1 = 1
    count2 = num_first + 1
    for disp1 in dataset['first_atoms']:
        disp_cart1 = disp1['displacement']
        w.write("- number: %5d\n" % (disp1['number'] + 1))
        w.write("  displacement:\n")
        w.write("    [%20.16f,%20.16f,%20.16f ] # %05d\n" %
                (disp_cart1[0], disp_cart1[1], disp_cart1[2], count1))
        w.write("  second_atoms:\n")
        count1 += 1

        included = None
        distance = 0.0
        atom2 = -1
        for disp2 in disp1['second_atoms']:
            if atom2 != disp2['number']:
                atom2 = disp2['number']
                if 'included' in disp2:
                    included = disp2['included']
                pair_distance = disp2['pair_distance']
                w.write("  - number: %5d\n" % (atom2 + 1))
                w.write("    distance: %f\n" % pair_distance)
                if included is not None:
                    if included:
                        w.write("    included: %s\n" % "true")
                    else:
                        w.write("    included: %s\n" % "false")
                w.write("    displacements:\n")

            disp_cart2 = disp2['displacement']
            w.write("    - [%20.16f,%20.16f,%20.16f ] # %05d\n" %
                    (disp_cart2[0], disp_cart2[1], disp_cart2[2], count2))
            count2 += 1

    write_cell_yaml(w, supercell)
    w.flush()
    lines = w.read()
    w.close()
    return lines


def get_forces_txt(force_sets):
    import StringIO

    w = StringIO.StringIO()

    from phono3py.file_IO import write_FORCES_FC3

    write_FORCES_FC3(force_sets.get_data_sets3(), force_sets.get_forces3(), fp=w)
    w.flush()
    lines = w.read()
    w.close()
    return lines