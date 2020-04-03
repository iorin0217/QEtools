# this code will be changed to Env class
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


def create_env(structure_file, variables, extfields={"press": 0}, constraints={"symm": False}):
    # unsupport the occupation != 1 case
    # TODO : conventional seekpath
    # we will deal with mcif: primitive=True
    env = {}
    if os.path.splitext(structure_file)[1][1:] == "cif":
        # read cif
        spin_structure = None
        structure = CifParser(structure_file).get_structures(
            primitive=False)[0]
        structure.remove_oxidation_states()
        atomic_numbers = [
            Element(str(u)).number for u in structure.species]
        atom_types = set([str(v) for v in structure.species])
    elif os.path.splitext(structure_file)[1][1:] == "vesta":
        # read vesta
        spin_structure = {}
        with open(structure_file) as f:
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
        atomic_numbers = [Element(u).number for u in species]
        atom_types = set(species)
        # organize spin_structure
        # frac_coodいらないかも
        _atom_count_outer = {atom: 0 for atom in atom_types}
        _sep = [x for x in range(len(VECTR)) if VECTR[x]
                [0] == "0" and VECTR[x-1][0] != "0"]
        _start = 1
        for _end in _sep:
            _atom_count_inner = {atom: 0 for atom in set(species)}
            for z in range(_start, _end):
                _index = int(VECTR[z][0])-1
                if _atom_count_inner[species[_index]] == 0:
                    _atom_count_inner[species[_index]] += 1
                    # no error if the specie not exsits
                    atom_types.discard(species[_index])
                    # to distinguish same atom_types with different spins in seekpath
                    atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000
                    atom_types.add(
                        species[_index]+str(_atom_count_outer[species[_index]]+1))
                    spin_structure[species[_index]+str(_atom_count_outer[species[_index]]+1)] = {"frac_coord": [frac_coord[_index]], "vec": np.array([
                        float(VECTR[_start-1][1]), float(VECTR[_start-1][2]), float(VECTR[_start-1][3])])}
                else:
                    _atom_count_inner[species[_index]] += 1
                    atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000
                    spin_structure[species[_index]+str(
                        _atom_count_outer[species[_index]]+1)]["frac_coord"].append(frac_coord[_index])
            _atom_count_outer[species[_index]] += 1
            _start = _end+2
    env["raw_crystal"] = structure
    env["raw_spin"] = spin_structure
    # set elemental configuration
    env["lspinorb"] = False
    env["nspin"] = 1
    env["lda_plus_u"] = False
    env["ecutwfc"] = 0
    env["ecutrho"] = 0
    env["time_reversal"] = True
    with open("../../settings/elements.json", "r") as f:
        elements = json.load(f)
    for atom in atom_types:
        element = re.match(r"\D+", atom).group()  # double count exists
        SOC = elements[element]["default"]["SOC"]
        Hubbard = elements[element]["default"]["Hubbard"]
        pstype = elements[element]["default"]["pstype"]
        if SOC == "fr":
            env["lspinorb"] = True
            env["nspin"] = 4  # to specify the quantization axis
        if Hubbard:
            env["lda_plus_u"] = True
        pseudo = elements[element]["pseudopotential"][variables["functional"]][SOC][pstype]
        env[atom] = {"pseudofile": pseudo["filename"], "Hubbard": Hubbard,
                     "starting_magnetization": 0}
        env["ecutwfc"] = max(env["ecutwfc"], pseudo["cutoff"])
        env["ecutrho"] = max(
            env["ecutrho"], pseudo["cutoff"]*pseudo["dual"])
    if spin_structure:
        parallel = False
        for i, j in itertools.permutations(spin_structure, 2):
            parallel = parallel or np.all(np.abs(np.cross(i["vec"],
                                                          j["vec"])) < 1e-5)  # VESTA VECTR has 1e-5 significatn digits
        if (not env["lspinorb"]) and parallel:
            env["nspin"] = 2
        else:
            env["nspin"] = 4
            env["time_reversal"] = False
    # set structural configuration
    skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords, atomic_numbers),
                                       with_time_reversal=env["time_reversal"], reference_distance=variables["reference_distance"])
    env["avec"] = skp["primitive_lattice"]
    env["bvec"] = skp["reciprocal_primitive_lattice"]
    duplicated = set(
        [atomnum % 1000 for atomnum in skp["primitive_types"] if atomnum > 1000])
    env["atoms"] = [str(get_el_sp(atomnum % 1000))+f"{atomnum//1000+1}" if atomnum %
                    1000 in duplicated else str(get_el_sp(atomnum)) for atomnum in skp["primitive_types"]]
    env["pos"] = skp["primitive_positions"]
    # easy to access the structual property
    env["crystal"] = Structure(Lattice(env["avec"]), env["atoms"], env["pos"])
    env["spin"] = None
    if spin_structure:
        # map spin structure from input to seekpath basis
        spg = spglib.get_symmetry_dataset(
            (structure.lattice.matrix, structure.frac_coords, atomic_numbers))
        permutation = np.dot(np.linalg.inv(
            spg["std_lattice"].T @ skp["primitive_transformation_matrix"]), env["avec"].T)
        mapping = np.linalg.inv(
            env["avec"].T) @ spg["std_rotation_matrix"] @ structure.lattice.matrix.T @ permutation
        for atom in spin_structure:
            vec = env["avec"].T @ mapping @ spin_structure[atom]["vec"] / \
                np.linalg.norm(env["avec"], axis=1)
            env["spin"].update({atom: {"frac_coord": [
                               mapping @ frac for frac in spin_structure[atom]["frac_coord"]], "vec": vec}})
            norm = np.linalg.norm(vec)
            if env["nspin"] == 2:
                if np.abs(vec[2]) > 1e-5:
                    pn_flag = np.sign(vec[2])
                else:
                    pn_flag = np.sign(vec[1])
                env[atom]["starting_magnetization"] = norm * pn_flag
            if env["nspin"] == 4:
                env[atom]["starting_magnetization"] = norm
                env[atom].update({
                    "angle1": np.rad2deg(np.arccos(vec[2]/norm))})
                env[atom].update({
                    "angle2": np.rad2deg(np.arctan(vec[1]/vec[0]))})
    # set reciprocal configuration
    env["kpath"] = skp["explicit_kpoints_rel"]
    env["klabel"] = skp["point_coords"]
    env["bandpath"] = skp["explicit_kpoints_linearcoord"]
    env["kticks"] = {}
    for ipath in range(len(skp["path"])):
        _start = skp["explicit_segments"][ipath][0]
        _last = skp["explicit_segments"][ipath-1][1] - 1
        if _start == _last+1:
            env["kticks"][(skp["explicit_kpoints_linearcoord"]
                           [_start]+skp["explicit_kpoints_linearcoord"][_last])/2] = skp["path"][ipath-1][1]+"|"+skp["path"][ipath][0]
        else:
            env["kticks"][skp["explicit_kpoints_linearcoord"]
                          [_start]] = skp["path"][ipath][0]
    env["kticks"][skp["explicit_kpoints_linearcoord"][-1]
                  ] = skp["path"][-1][1]  # retrieve the final point
    _nks = np.round(np.linalg.norm(
        env["bvec"], axis=0) / variables["dk_grid"])
    env["nk"] = np.where(nks < 1, 1, _nks)
    # extfields, constraints
    # hdf5

    return env
