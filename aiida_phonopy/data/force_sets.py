from aiida.orm.data.array import ArrayData


class ForceSetsData(ArrayData):
    """
    Store the force constants on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(ForceSetsData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def get_number_of_atoms(self):
        """
        Return the shape of an array (read from the value cached in the
        properties for efficiency reasons).

        :param name: The name of the array.
        """
        return self.get_attr("natom")

    def get_number_of_displacements(self):
        """
        Return the shape of an array (read from the value cached in the
        properties for efficiency reasons).

        :param name: The name of the array.
        """
        return self.get_attr("ndisplacements")

    def get_data_sets(self):
        """
        Return the force constants stored in the node as a numpy array
        """

        natom = self.get_attr("natom")
        ndisplacements = self.get_attr("ndisplacements")

        direction = self.get_array('direction')
        number = self.get_array('number')
        displacement = self.get_array('displacement')

        first_atoms = []
        for i in range(ndisplacements):
            first_atoms.append({'direction': direction[i],
                                'number': number[i],
                                'displacement': displacement[i]})

        return {'natom': natom, 'first_atoms': first_atoms}

    def get_force_sets(self):
        """
        Return the force constants stored in the node as a numpy array
        """

        natom = self.get_attr("natom")
        ndisplacements = self.get_attr("ndisplacements")

        direction = self.get_array('direction')
        number = self.get_array('number')
        displacement = self.get_array('displacement')
        forces = self.get_array('forces')

        first_atoms = []
        for i in range(ndisplacements):
            first_atoms.append({'direction': direction[i],
                                'number': number[i],
                                'forces': forces[i],
                                'displacement': displacement[i]})

        return {'natom': natom, 'first_atoms': first_atoms}

    # {'natom': 64, 'first_atoms': [{'direction': [1, 0, 0], 'number': 0, 'displacement': array([0.01, 0., 0.])}]}

    def set_data_sets(self, data_sets):

        import numpy

        self._set_attr('natom', data_sets['natom'])
        self._set_attr('ndisplacements', len(data_sets['first_atoms']))

        direction = []
        number = []
        displacement = []
        for first_atoms in data_sets['first_atoms']:
            number.append(first_atoms['number'])
            displacement.append(first_atoms['displacement'])
            if 'direction' in first_atoms:
                direction.append(first_atoms['direction'])
            else:
                direction.append([])

        self.set_array('direction', numpy.array(direction))
        self.set_array('number', numpy.array(number))
        self.set_array('displacement', numpy.array(displacement))

    def set_forces(self, forces):

        import numpy
        self.set_array('forces', numpy.array(forces))

    def read_from_phonopy_file(self, filename):
        """
        Read the force constants from a phonopy FORCE_SETS file

        :param filename: FORCE_SETS file name
        """

        from phonopy.file_IO import parse_FORCE_SETS

        data_sets = parse_FORCE_SETS(filename=filename)

        self.set_data_sets(data_sets)
        self.set_forces([displacement['forces'] for displacement in data_sets['first_atoms']])

    # phono3py

    def set_forces3(self, forces):
        import numpy
        self.set_array('forces3', numpy.array(forces))

    def get_forces3(self):
        forces_list = self.get_array('forces3')
        return [forces for forces in forces_list]

    def set_data_sets3(self, data_sets):
        import numpy

        print data_sets['first_atoms'][0].keys()

        self._set_attr('natom', data_sets['natom'])
        self._set_attr('ndisplacements', len(data_sets['first_atoms']))

        direction = []
        number = []
        displacement = []
        pair_distance = []
        ndisplacements_s = []

        direction_f = []
        number_f = []
        displacement_f = []

        for first_atoms in data_sets['first_atoms']:

            direction_f.append(first_atoms['direction'])
            displacement_f.append(first_atoms['displacement'])
            number_f.append(first_atoms['number'])

            direction_s = []
            number_s = []
            displacement_s = []
            pair_distance_s = []
            ndisplacements_s.append(len(first_atoms['second_atoms']))

            for second_atoms in first_atoms['second_atoms']:
                number_s.append(second_atoms['number'])
                displacement_s.append(second_atoms['displacement'])
                direction_s.append(second_atoms['direction'])
                pair_distance_s.append(second_atoms['pair_distance'])
            number.append(number_s)

            displacement.append(displacement_s)
            direction.append(direction_s)
            pair_distance.append(pair_distance_s)

        self.set_array('direction_s', numpy.array(direction))
        self.set_array('number_s', numpy.array(number))
        self.set_array('displacement_s', numpy.array(displacement))
        self.set_array('pair_distance_s', numpy.array(pair_distance))

        self.set_array('direction', numpy.array(direction_f))
        self.set_array('number', numpy.array(number_f))
        self.set_array('displacement', numpy.array(displacement_f))

        self._set_attr('ndisplacements_s', ndisplacements_s)

    def get_data_sets3(self):
        natom = self.get_attr("natom")
        ndisplacements = self.get_attr("ndisplacements")
        ndisplacements_s = self.get_attr("ndisplacements_s")

        direction = self.get_array('direction_s')
        number = self.get_array('number_s')
        displacement = self.get_array('displacement_s')
        pair_distance = self.get_array('pair_distance_s')

        direction_f = self.get_array('direction')
        number_f = self.get_array('number')
        displacement_f = self.get_array('displacement')

        first_atoms = []
        for i, ndisplacements_s in enumerate(ndisplacements_s):
            second_atoms = []
            for j in range(ndisplacements_s):
                second_atoms.append({'direction': direction[i, j],
                                     'number': number[i, j],
                                     'displacement': displacement[i, j],
                                     'pair_distance': pair_distance[i, j]})

            first_atoms.append({'direction': direction_f[i],
                                'displacement': displacement_f[i],
                                'number': number_f[i],
                                'second_atoms': second_atoms})

        return {'natom': natom, 'first_atoms': first_atoms}
