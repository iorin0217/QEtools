# %%
import spglib
import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import CifWriter
from pymatgen import Structure
from pymatgen import Lattice
from pymatgen import Element
from pymatgen.core.periodic_table import get_el_sp
import seekpath
import json
import os
import re
import itertools
import copy
# %%
with open("/home/CMD35/cmd35stud07/experiments/AgCrSe2_afm/AgCrSe2_R3m.vesta") as f:
    vesta = f.readlines()
# %%
CELLP = [i.strip().split() for i in vesta[int(vesta.index("CELLP\n")) +
                                          1:int(vesta.index("STRUC\n"))] if i.strip().split()[0] != "0"]
STRUC = [j.strip().split() for j in vesta[int(vesta.index("STRUC\n")) +
                                          1:int(vesta.index("THERI 0\n"))] if j.strip().split()[0] != "0"]
VECTR = [k.strip().split()
         for k in vesta[int(vesta.index("VECTR\n"))+1:int(vesta.index("VECTT\n"))]]
a_vesta, b_vesta, c_vesta, alpha_vesta, beta_vesta, gamma_vesta = [
    float(l) for l in CELLP[0]]
species = [STRUC[m][1] for m in range(len(STRUC)) if m % 2 == 0]
frac_coord = [[float(STRUC[n][4]), float(STRUC[n][5]), float(
    STRUC[n][6])] for n in range(len(STRUC)) if n % 2 == 0]
structure = Structure(Lattice.from_parameters(
    a_vesta, b_vesta, c_vesta, alpha_vesta, beta_vesta, gamma_vesta), species, frac_coord)
atomic_numbers = [Element(u).number for u in species]
atoms = set(species)
spin_structure = {}
_atom_count_outer = {atom: 0 for atom in atoms}
_sep = [x for x in range(len(VECTR)) if VECTR[x][0] ==
        "0" and VECTR[x-1][0] != "0"]
_start = 1
# %%
for _end in _sep:
    _atom_count_inner = {atom: 0 for atom in set(species)}
    for z in range(_start, _end):
        _index = int(VECTR[z][0])-1
        if _atom_count_inner[species[_index]] == 0:
            _atom_count_inner[species[_index]] += 1
            # no error when the specie not exsits
            atoms.discard(species[_index])
            # to distinguish same atoms with different spins in seekpath
            atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000
            # to distinguish same atoms with different spins in seekpath
            atoms.add(species[_index] +
                      str(_atom_count_outer[species[_index]]+1))
            spin_structure[species[_index]+str(_atom_count_outer[species[_index]]+1)] = {"frac_coord": [frac_coord[_index]], "vec": [
                float(VECTR[_start-1][1]), float(VECTR[_start-1][2]), float(VECTR[_start-1][3])]}
        else:
            _atom_count_inner[species[_index]] += 1
            atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000
            spin_structure[species[_index]+str(
                _atom_count_outer[species[_index]]+1)]["frac_coord"].append(frac_coord[_index])
    _atom_count_outer[species[_index]] += 1
    _start = _end+2
# %%
skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                    atomic_numbers), with_time_reversal=True, reference_distance=0.025)
dataset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords,
                                       atomic_numbers))
trans_f = dataset["transformation_matrix"]
trans_r = dataset["std_rotation_matrix"]
abc = structure.lattice.matrix
after = skp["primitive_lattice"]
before = dataset["std_lattice"]
# %%
abc
# %%
before
# %%
after
# %%
(before.T @ skp["primitive_transformation_matrix"]).T
# %%
np.dot(after, np.linalg.inv(before))
# %%
(trans_r @ abc.T @ np.linalg.inv(trans_f)).T
# %%
np.linalg.inv(trans_r @ abc.T) @ (
    after.T @ skp["primitive_positions"][3] - trans_r @ abc.T @ structure.frac_coords[3])
# %%
spin_ini = np.array([-2.0, 4.0, -2.0])
# %%
trans_r @ abc.T @ spin_ini
# %%
structure_before = Structure(before, species*3, dataset["std_positions"])
# %%
CifWriter(structure_before).write_file(
    "/home/CMD35/cmd35stud07/experiments/AgCrSe2_afm/AgCrSe2_spg.cif")
# %%
np.linalg.norm(np.dot(skp["primitive_positions"],
                      skp["primitive_lattice"]), axis=1)
structure.get_neighbors_in_shell([0, 0, 0], 13.12, 0.01)
# %%
structure.make_supercell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
abc_super = structure.lattice.matrix
trans_l = abc[0]+abc[1]+abc[2]
super_structure = np.dot(structure.frac_coords, abc_super)-trans_l/2
np.linalg.norm(super_structure)
# %%
structure_cif = CifParser("/home/CMD35/cmd35stud07/experiments/AgCrSe2_afm/AgCrSe2_R3m.cif").get_structures(
    primitive=False)[0]  # when we use mcif, primitive=True
structure_cif.remove_oxidation_states()
atomic_numbers_cif = [
    Element(str(u)).number for u in structure_cif.species]
atomic_numbers_cif = [24, 24, 24, 24, 47, 47,
                      47, 47, 34, 34, 34, 34, 34, 34, 34, 34]
dataset_cif = spglib.get_symmetry_dataset((structure_cif.lattice.matrix, structure_cif.frac_coords,
                                           atomic_numbers_cif))


# %%
structure_cif

# %%
