# %%
import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen import Structure
from pymatgen import Lattice
from pymatgen import Element
from pymatgen.core.periodic_table import get_el_sp
import spglib
import seekpath
import json
import os
import re
import itertools


class Env:
    def __init__(self, structure, variables={}, atoms=[], spin_structure={}):
        """
        main usage
            from_hoge
        input
            structure : pymatgen.Structure
            variables : dict from variables.json
            atoms : ["$atom","$atom", ...]
                distinguish atoms with number like
                    ["Sr1","Sr2","Ru","O","O","O","O"]
            spin_structure: {"$atom": {'frac_coord': [np.array(3),np.array(3), ...], 'vec':np.array(3)}
        """
        # save the input structure
        self.init_crystal, self.init_atoms, self.init_spin = structure, atoms, spin_structure
        # variables
        self.functional = variables["functional"] if variables else "PBE"
        self.reference_distance = variables["reference_distance"] if variables else 0.025
        # elements information
        self.elements, self.lspinorb, self.nspin, self.lda_plus_u, self.ecutwfc, self.ecutrho, self.time_reversal = self._set_elements()
        # structure information
        self.crystal, self.atoms, self.spin, self.bvec, self.bandpath, self.klabel, self.bandpath_linear, self.kticks = self._set_structure()
        # extfields & constraints
        self.extfields = {}
        self.constraints = {}

    def _set_elements(self):
        # setting
        if not self.init_atoms:
            self.init_atoms = [str(u) for u in self.init_crystal.species]
        with open("../settings/elements.json", "r") as f:
            elements_json = json.load(f)
        # pseudopotential & SOC / Hubbard U
        lspinorb = False
        nspin = 1
        lda_plus_u = False
        ecutwfc = 0
        ecutrho = 0
        time_reversal = True
        elements = {}
        for atom in set(self.init_atoms):
            # double atom count exists
            element = re.match(r"\D+", atom).group()
            SOC = elements_json[element]["default"]["SOC"]
            Hubbard = elements_json[element]["default"]["Hubbard"]
            pstype = elements_json[element]["default"]["pstype"]
            if SOC == "fr":
                lspinorb = True
                nspin = 4  # to specify the quantization axis
            if Hubbard:
                lda_plus_u = True
            pseudo = elements_json[element]["pseudopotential"][self.functional][SOC][pstype]
            '''TODO:valence
            elements = {"pseudofile": pseudo["filename"], "valence": pseudo["valence"]
                        "Hubbard": Hubbard, "starting_magnetization": 0}
            '''
            elements[atom] = {"pseudofile": pseudo["filename"],
                              "Hubbard": Hubbard, "starting_magnetization": 0}
            ecutwfc = max(ecutwfc, pseudo["cutoff"])
            ecutrho = max(ecutrho, pseudo["cutoff"] * pseudo["dual"])
        # nspin
        if self.init_spin:
            parallel = False
            for i, j in itertools.permutations(self.init_spin, 2):
                # VESTA VECTR has 1e-5 significatn digits
                parallel = parallel or np.all(
                    np.abs(np.cross(self.init_spin[i]["vec"], self.init_spin[j]["vec"])) < 1e-5)
            if (not lspinorb) and parallel:
                nspin = 2
            else:
                nspin = 4
                time_reversal = False
        return elements, lspinorb, nspin, lda_plus_u, ecutwfc, ecutrho, time_reversal

    def update_elements(self):
        # you can update : functional, pstype, Hubbard U, SOC, spin structure
        # automatically call _set_elements & _set_structure
        # TODO : implement
        pass

    def _set_structure(self):
        # standarize crystal
        _atomic_numbers = [Element(re.findall(r'(\d+|\D+)', i)[0]).number if len(re.findall(r'(\d+|\D+)', i)) == 1 else Element(
            re.findall(r'(\d+|\D+)', i)[0]).number + 1000 * (int(re.findall(r'(\d+|\D+)', i)[1]) - 1) for i in self.init_atoms]
        skp = seekpath.get_explicit_k_path((self.init_crystal.lattice.matrix, self.init_crystal.frac_coords,
                                            _atomic_numbers), with_time_reversal=self.time_reversal, reference_distance=self.reference_distance)
        _avec = skp["primitive_lattice"]
        _duplicated = set(
            [atomnum % 1000 for atomnum in skp["primitive_types"] if atomnum > 1000])
        atoms = [str(get_el_sp(atomnum % 1000)) + f"{atomnum//1000+1}" if atomnum % 1000 in _duplicated else str(
            get_el_sp(atomnum)) for atomnum in skp["primitive_types"]]
        _pos = skp["primitive_positions"]
        crystal = Structure(Lattice(_avec), [re.findall(
            r'(\d+|\D+)', i)[0] for i in atoms], _pos)
        # map spin structure from input to seekpath basis via spglib
        spin = {}
        if self.init_spin:
            spg = spglib.get_symmetry_dataset(
                (self.init_crystal.lattice.matrix, self.init_crystal.frac_coords, _atomic_numbers))
            permutation = np.dot(np.linalg.inv(
                spg["std_lattice"].T @ skp["primitive_transformation_matrix"]), skp["primitive_lattice"].T)
            mapping = np.linalg.inv(
                skp["primitive_lattice"].T) @ spg["std_rotation_matrix"] @ self.init_crystal.lattice.matrix.T @ permutation
            for atom in self.init_spin:
                vec = skp["primitive_lattice"].T @ mapping @ self.init_spin[atom]["vec"] / \
                    np.linalg.norm(skp["primitive_lattice"], axis=1)
                spin.update({atom: {"frac_coord": [
                    mapping @ frac for frac in self.init_spin[atom]["frac_coord"]], "vec": vec}})
                norm = np.linalg.norm(vec)
                if self.nspin == 2:
                    # convert direction to sign
                    if np.abs(vec[2]) > 1e-5:
                        pn_flag = np.sign(vec[2])
                    else:
                        pn_flag = np.sign(vec[1])
                    elements[atom]["starting_magnetization"] = norm * pn_flag
                if self.nspin == 4:
                    elements[atom]["starting_magnetization"] = norm
                    elements[atom].update({
                        "angle1": np.rad2deg(np.arccos(vec[2]/norm))})
                    elements[atom].update({
                        "angle2": np.rad2deg(np.arctan(vec[1]/vec[0]))})
        # reciprocal structure
        bvec = skp["reciprocal_primitive_lattice"]
        bandpath = skp["explicit_kpoints_rel"]
        klabel = skp["point_coords"]
        bandpath_linear = skp["explicit_kpoints_linearcoord"]
        kticks = {}
        for ipath in range(len(skp["path"])):
            _start = skp["explicit_segments"][ipath][0]
            _last = skp["explicit_segments"][ipath-1][1] - 1
            if _start == _last+1:
                kticks[(skp["explicit_kpoints_linearcoord"]
                        [_start]+skp["explicit_kpoints_linearcoord"][_last])/2] = skp["path"][ipath-1][1]+"|"+skp["path"][ipath][0]
            else:
                kticks[skp["explicit_kpoints_linearcoord"]
                       [_start]] = skp["path"][ipath][0]
        kticks[skp["explicit_kpoints_linearcoord"][-1]
               ] = skp["path"][-1][1]  # retrieve the final point
        return crystal, atoms, spin, bvec, bandpath, klabel, bandpath_linear, kticks

    def update_crystal(self):
        # you can update : lattice vectors, atomic positions
        # automatically call _set_structure
        # TODO : implement
        pass

    def set_extfields(self, extfields={"press": 0}):
        # TODO : implement
        pass

    def set_constraints(self, constraints={"symm": False}):
        # TODO : implement
        pass

    def save_cif(self):
        # TODO : implement
        pass

    def save_vesta(self):
        # TODO : implement
        pass

    def save_mcif(self):
        # TODO : implement
        pass

    def save_ase(self):
        # TODO : implement
        pass

    def save_hdf5(self):
        # TODO : implement
        pass

    @classmethod
    def from_cif(cls, ciffile, variables={}, atoms=[], spin_structure={}):
        structure = CifParser(ciffile).get_structures(primitive=False)[0]
        structure.remove_oxidation_states()
        return cls(structure, variables, atoms=atoms, spin_structure=spin_structure)

    @classmethod
    def from_vesta(cls, vestafile, variables={}):
        # read vesta
        spin_structure = {}
        with open(vestafile) as f:
            vesta = f.readlines()
        CELLP = [i.strip().split() for i in vesta[int(vesta.index(
            "CELLP\n"))+1:int(vesta.index("STRUC\n"))] if i.strip().split()[0] != "0"]
        STRUC = [j.strip().split() for j in vesta[int(vesta.index(
            "STRUC\n"))+1:int(vesta.index("THERI 0\n"))] if j.strip().split()[0] != "0"]
        VECTR = [k.strip().split() for k in vesta[int(
            vesta.index("VECTR\n"))+1:int(vesta.index("VECTT\n"))]]
        a_vesta, b_vesta, c_vesta, alpha_vesta, beta_vesta, gamma_vesta = [
            float(l) for l in CELLP[0]]
        species = [STRUC[m][1] for m in range(len(STRUC)) if m % 2 == 0]
        frac_coord = [np.array([float(STRUC[n][4]), float(STRUC[n][5]), float(
            STRUC[n][6])]) for n in range(len(STRUC)) if n % 2 == 0]
        structure = Structure(Lattice.from_parameters(
            a_vesta, b_vesta, c_vesta, alpha_vesta, beta_vesta, gamma_vesta), species, frac_coord)
        atoms = species[:]
        # organize spin_structure
        _atom_count_outer = {atom: 0 for atom in set(atoms)}
        _sep = [x for x in range(len(VECTR)) if VECTR[x]
                [0] == "0" and VECTR[x-1][0] != "0"]
        _start = 1
        for _end in _sep:
            _atom_count_inner = {atom: 0 for atom in set(species)}
            for z in range(_start, _end):
                _index = int(VECTR[z][0])-1
                if _atom_count_inner[species[_index]] == 0:
                    _atom_count_inner[species[_index]] += 1
                    # to distinguish same atom types with different spins in seekpath
                    atoms[_index] = species[_index] + \
                        str(_atom_count_outer[species[_index]]+1)
                    spin_structure[species[_index]+str(_atom_count_outer[species[_index]]+1)] = {"frac_coord": [frac_coord[_index]], "vec": np.array([
                        float(VECTR[_start-1][1]), float(VECTR[_start-1][2]), float(VECTR[_start-1][3])])}
                else:
                    _atom_count_inner[species[_index]] += 1
                    atoms[_index] = species[_index] + \
                        str(_atom_count_outer[species[_index]]+1)
                    spin_structure[species[_index]+str(
                        _atom_count_outer[species[_index]]+1)]["frac_coord"].append(frac_coord[_index])
            _atom_count_outer[species[_index]] += 1
            _start = _end+2
        return cls(structure, variables, atoms=atoms, spin_structure=spin_structure)

    @classmethod
    def from_mcif(cls, mciffile, variables={}):
        # TODO : implement
        pass

    @classmethod
    def from_ase(cls, asestructure, variables={}, atoms=[], spin_structure={}):
        # TODO : implement
        pass

    @classmethod
    def from_xml(cls, xmlfile):
        # TODO : implement
        pass


# %%
cif = "/home/CMD35/cmd35stud07/experiments/Sr2RhO4/Sr2RhO4_mp-757102_conventional_standard.cif"
env = Env.from_cif(cif)

# %%
env.lda_plus_u = False
scf = PW("/home/CMD35/cmd35stud07/experiments/BHO", env, calculation="scf")
nscf = PW("/home/CMD35/cmd35stud07/experiments/BHO",
          env, calculation="nscf")
bands = PW("/home/CMD35/cmd35stud07/experiments/BHO",
           env, calculation="bands")


# %%
