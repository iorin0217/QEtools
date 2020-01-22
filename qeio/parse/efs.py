import sys,math
def out2xtl(f_out):
	# pwout2xtl.py
	# made by mkanzaki@me.com
	# Utility code for Quantum-Espsresso and Vesta.
	# This code converts pw.x output (optimization) to .xtl file,
	# which can be read by Vesta.
	# Assume output file extention as .out 
	# or modify code: base = sys.argv[1].rstrip('.out')
	# pw.x output must be optimizatin (vc-relax or relax) run, 
	# and optimization was successful (i.e., "Final" word exists in the file).
	# Usage:
	# >./pwout2xtl.py test.out
	# Then test.xtl will be produced.
	# Open it from Vesta.
	mat1 = [0.0 for i in range(3)]
	mat2 = [0.0 for i in range(3)]
	mat3 = [0.0 for i in range(3)]
	cell = ['' for i in range(4)]
	# conversion for radian 
	rd = math.pi/180.0
	# conversion factor for atomic unit to angstrom 
	bohr = 0.5291772108
	# reading output file
	# check xtl file name is given or not
	try:
		file1 = open(f_out,'r')
	except IOError:
		print ('Output file open error!')
		exit()
	# Get .out file name
	base = f_out.rstrip('.out')
	# Open output file for .xtl
	f2 = base + '.xtl'
	file2 = open(f2,'w')
	# find a line containing "Final"
	while True:
		ftext = file1.readline()
		if 'Final' in ftext :
			break
		if ftext=="":
			print ('Optimization failed in this output file!')
			file1.close()
			exit()		
	# find cell parameters
	while True:
		ftext = file1.readline()
		if 'CELL_PARAMETERS' in ftext :
			tmp = ftext.strip()
			break
	s = tmp.index('alat=')+5
	t = tmp.index(')')
	alat = float(tmp[s:t])
	# reading cell matrix
	cell[0] = file1.readline().strip()
	list = cell[0].split(' ')
	while list.count('') > 0:
		list.remove('')
	for i in range(0,3):
		mat1[i] = float(list[i])
	cell[1] = file1.readline().strip()
	list = cell[1].split(' ')
	while list.count('') > 0:
		list.remove('')
	for i in range(0,3):
		mat2[i] = float(list[i])
	cell[2] = file1.readline().strip()
	list = cell[2].split(' ')
	while list.count('') > 0:
		list.remove('')
	for i in range(0,3):
		mat3[i] = float(list[i])
	# cell parameters calculation
	cella = math.sqrt(mat1[0]*mat1[0]+mat1[1]*mat1[1]+mat1[2]*mat1[2])
	cellb = math.sqrt(mat2[0]*mat2[0]+mat2[1]*mat2[1]+mat2[2]*mat2[2])
	cellc = math.sqrt(mat3[0]*mat3[0]+mat3[1]*mat3[1]+mat3[2]*mat3[2])
	alpha = math.acos((mat2[0]*mat3[0]+mat2[1]*mat3[1]+mat2[2]*mat3[2])/(cellb*cellc))/rd
	beta = math.acos((mat1[0]*mat3[0]+mat1[1]*mat3[1]+mat1[2]*mat3[2])/(cella*cellc))/rd
	gamma = math.acos((mat1[0]*mat2[0]+mat1[1]*mat2[1]+mat1[2]*mat2[2])/(cella*cellb))/rd
	cella = alat*bohr*cella
	cellb = alat*bohr*cellb
	cellc = alat*bohr*cellc
	# Write title and cell to .xtl file
	file2.write('TITLE ' + 'Produced from pw.x calculation: ' + base + '\n')
	file2.write('CELL\n')
	file2.write(' ' + str(cella) + ' ' + str(cellb) + ' ' + str(cellc) + ' ' + str(alpha) + ' ' + str(beta) + ' ' + str(gamma) + '\n')   
	# find atomic positions
	while True:
		ftext = file1.readline()
		if 'ATOMIC_POSITIONS' in ftext :
			break
	# write atomic positions
	file2.write('SYMMETRY NUMBER 1\n')
	file2.write('SYMMETRY LABEL  P1\n')
	file2.write('ATOMS\n')
	file2.write('NAME         X           Y           Z\n')
	while True:
		ftext = file1.readline()
		if 'End final coordinates' in ftext :
			break
		else:
			file2.write(ftext)
	# write EOF line
	file2.write('EOF\n')
	file1.close()
	file2.close()
	return



"""Sub class of `Data` to handle interatomic force constants produced by the Quantum ESPRESSO q2r.x code."""
from __future__ import absolute_import

import numpy
from six.moves import range

from qe_tools.constants import bohr_to_ang
from aiida.orm import SinglefileData


class ForceConstantsData(SinglefileData):
    """Class to handle interatomic force constants from the Quantum ESPRESSO q2r.x code."""

    def set_file(self, file):
        """Add a file to the node, parse it and set the attributes found.
        :param file: absolute path to the file or a filelike object
        """
        # pylint: disable=redefined-builtin,arguments-differ
        super(ForceConstantsData, self).set_file(file)

        # Parse the force constants file
        dictionary, _, _ = parse_q2r_force_constants_file(self.get_content().splitlines(), also_force_constants=False)

        # Add all other attributes found in the parsed dictionary
        for key, value in dictionary.items():
            self.set_attribute(key, value)

    @property
    def number_of_species(self):
        """Return the number of atom species.
        :return: a scalar
        """
        return self.get_attribute('number_of_species')

    @property
    def number_of_atoms(self):
        """Return the number of atoms.
        :return: a scalar
        """
        return self.get_attribute('number_of_atoms')

    @property
    def cell(self):
        """Return the crystal unit cell where rows are the crystal vectors.
        :return: a 3x3 numpy.array
        """
        return numpy.array(self.get_attribute('cell'))

    @property
    def atom_list(self):
        """Return the list of atoms.
        :return: a list of length-5 tuple (element name, element mass amu_ry, 3 coordinates in cartesian Angstrom)
        """
        return self.get_attribute('atom_list')

    @property
    def has_done_electric_field(self):
        """Return flag to indicate if dielectric tensor and effective charges were computed.
        :return: a boolean
        """
        return self.get_attribute('has_done_electric_field')

    @property
    def dielectric_tensor(self):
        """Return the dielectric tensor matrix.
        :return: a 3x3 tuple
        """
        return self.get_attribute('dielectric_tensor')

    @property
    def effective_charges_eu(self):
        """Return the effective charges for each atom.
        :return: a list of number_of_atoms elements, each being a 3x3 tuple
        """
        return self.get_attribute('effective_charges_eu')

    @property
    def qpoints_mesh(self):
        """Return the number of q-points in each direction.
        :return: a length-3 tuple
        """
        return tuple(self.get_attribute('qpoints_mesh'))


def parse_q2r_force_constants_file(lines, also_force_constants=False):
    """Parse the real-space interatomic force constants file from QE-Q2R.
    :param also_force_constants: True to parse the force constants as well
    :return parsed_data: dictionary with the following keywords:
    - number_of_species: number of atom species ('ntyp' in QE)
    - number_of_atoms: number of atoms ('nat' in QE)
    - cell: unit cell
    - atom_list: list with, for each atom in the cell, a length-5
      tuple of the form (element_name, mass_in_amu_ry, then 3 coordinates in
      cartesian & Angstroms)
    - has_done_electric_field: True if dielectric constants & effective
      charges were computed
    - dielectric_tensor: dielectric constant (3x3 matrix)
    - effective_charges_eu: effective charges (ntyp x 3 x 3 matrix)
    - qpoints_mesh: length-3 tuple with number of qpoints in each dimension
      of the reciprocal lattice
    - force_constants: the real-space force constants: array with 7 indices, of the kind
        C(mi1, mi2, mi3, ji1, ji2, na1, na2) with
        * (mi1, mi2, mi3): the supercell dimensions
        * (ji1, ji2): axis of the displacement of the two atoms (from 1 to 3)
        * (na1, na2): atom numbers in the cell.
    - warnings: a list of warnings
    :return force_constants: the real-space force constants: array with 7 indices, of the kind
        C(mi1, mi2, mi3, ji1, ji2, na1, na2) where:
        * (mi1, mi2, mi3): the supercell dimensions
        * (ji1, ji2): axis of the displacement of the two atoms (from 1 to 3)
        * (na1, na2): atom numbers in the cell.
    """
    # pylint: disable=too-many-statements,too-many-branches,too-many-nested-blocks

    parsed_data = {}
    warnings = []

    try:
        # read first line
        current_line = 0
        first_line = lines[current_line].split()
        ntyp = int(first_line[0])
        nat = int(first_line[1])
        ibrav = int(first_line[2])
        celldm = [float(c) for c in first_line[3:]]
        if len(celldm) != 6:
            warnings.append('Wrong length for celldm')
        if ibrav != 0:
            warnings.append('ibrav ({}) is not 0; q-points path for phonon ' 'dispersion might be wrong'.format(ibrav))
        if any([item != 0 for item in celldm[1:]]):
            warnings.append('celldm[1:] are not all zero; only celldm[0] will ' 'be used')

        parsed_data['number_of_species'] = ntyp
        parsed_data['number_of_atoms'] = nat
        current_line += 1

        # read cell data
        cell = tuple(
            tuple(float(c) * celldm[0] * bohr_to_ang for c in l.split()) for l in lines[current_line:current_line + 3]
        )
        parsed_data['cell'] = cell
        current_line += 3

        # read atom types and masses
        atom_type_list = []
        for ityp in range(ntyp):
            line = lines[current_line].split("'")
            if int(line[0]) == ityp + 1:
                atom_type_list.append(tuple((line[1].strip(), float(line[2]))))
            current_line += 1

        # read each atom coordinates
        atom_list = []
        for _ in range(nat):
            line = [float(c) for c in lines[current_line].split()]
            ityp = int(line[1])
            if 0 < ityp < ntyp + 1:
                line[0] = atom_type_list[ityp - 1][0]  # string with element name
                line[1] = atom_type_list[ityp - 1][1]  # element mass in amu_ry
                # Convert atomic positions (in cartesian) from alat to Angstrom:
                line[2:] = [pos * celldm[0] * bohr_to_ang for pos in line[2:]]
            atom_list.append(tuple(line))
            current_line += 1

        parsed_data['atom_list'] = atom_list

        # read lrigid (flag for dielectric constant and effective charges
        has_done_electric_field = (lines[current_line].split()[0] == 'T')
        parsed_data['has_done_electric_field'] = has_done_electric_field
        current_line += 1

        if has_done_electric_field:
            # read dielectric tensor
            dielectric_tensor = tuple(tuple(float(c) for c in l.split()) for l in lines[current_line:current_line + 3])
            current_line += 3
            effective_charges_eu = []
            for _ in range(nat):
                current_line += 1
                effective_charges_eu.append(
                    tuple(tuple(float(c) for c in l.split()) for l in lines[current_line:current_line + 3])
                )
                current_line += 3

            parsed_data['dielectric_tensor'] = dielectric_tensor
            parsed_data['effective_charges_eu'] = effective_charges_eu

        # read q-points mesh
        qpoints_mesh = tuple(int(c) for c in lines[current_line].split())
        current_line += 1
        parsed_data['qpoints_mesh'] = qpoints_mesh

        force_constants = ()
        if also_force_constants:
            # read force_constants
            force_constants = numpy.zeros(qpoints_mesh + (3, 3, nat, nat), dtype=float)
            for ji1 in range(3):
                for ji2 in range(3):
                    for na1 in range(nat):
                        for na2 in range(nat):

                            indices = tuple([int(c) for c in lines[current_line].split()])
                            current_line += 1
                            if (ji1 + 1, ji2 + 1, na1 + 1, na2 + 1) != indices:
                                raise ValueError('Wrong indices in force constants')

                            for mi3 in range(qpoints_mesh[2]):
                                for mi2 in range(qpoints_mesh[1]):
                                    for mi1 in range(qpoints_mesh[0]):

                                        line = lines[current_line].split()
                                        indices = tuple(int(c) for c in line[:3])

                                        if (mi1 + 1, mi2 + 1, mi3 + 1) != indices:
                                            raise ValueError('Wrong supercell indices in force constants')

                                        force_constants[mi1, mi2, mi3, ji1, ji2, na1, na2] = float(line[3])
                                        current_line += 1

    except (IndexError, ValueError) as exception:
        raise ValueError(str(exception) + '\nForce constants file could not be parsed (incorrect file format)')

    return parsed_data, force_constants, warnings