import numpy as np
from pymatgen.io.cif import CifParser
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

def create_env(structure_file, extfields={"press":0}, constrains={"symm":False}):
    env = {}
    if os.path.splitext(structure_file)[1][1:] == "xml":
        pass #parse_xml_to_env
    else:
        if os.path.splitext(structure_file)[1][1:] == "cif":
            spin_calc = False
            structure = CifParser(structure_file).get_structures(primitive=False)[0] # when we use mcif, primitive=True
            structure.remove_oxidation_states()
            atomic_numbers = [Element(str(u)).number for u in structure.species]
            atoms = set([str(v) for v in structure.species])
        elif os.path.splitext(structure_file)[1][1:] == "vesta":
            spin_calc = True
            with open(structure_file) as f:
                vesta = f.readlines()
            CELLP = [i.strip().split() for i in vesta[int(vesta.index("CELLP\n"))+1:int(vesta.index("STRUC\n"))] if i.strip().split()[0]!="0"]
            STRUC = [j.strip().split() for j in vesta[int(vesta.index("STRUC\n"))+1:int(vesta.index("THERI 0\n"))] if j.strip().split()[0]!="0"]
            VECTR = [k.strip().split() for k in vesta[int(vesta.index("VECTR\n"))+1:int(vesta.index("VECTT\n"))]]
            a_vesta,b_vesta,c_vesta,alpha_vesta,beta_vesta,gamma_vesta = [float(l) for l in CELLP[0]]
            species = [STRUC[m][1] for m in range(len(STRUC)) if m%2==0]
            frac_coord = [[float(STRUC[n][4]), float(STRUC[n][5]), float(STRUC[n][6])] for n in range(len(STRUC)) if n%2==0]
            structure = Structure(Lattice.from_parameters(a_vesta,b_vesta,c_vesta,alpha_vesta,beta_vesta,gamma_vesta), species, frac_coord)
            atomic_numbers = [Element(u).number for u in species]
            atoms = set(species)
            spin_structure = {}
            _atom_count_outer = {atom : 0 for atom in atoms}
            _sep = [x for x in range(len(VECTR)) if VECTR[x][0]=="0" and VECTR[x-1][0]!="0"]
            _start = 1
            for _end in _sep:
                _atom_count_inner = {atom : 0 for atom in set(species)}
                for z in range(_start,_end):
                    _index = int(VECTR[z][0])-1
                    if _atom_count_inner[species[_index]] == 0:
                        _atom_count_inner[species[_index]] += 1
                        atoms.discard(species[_index]) # no error when the specie not exsits
                        atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000 # to distinguish same atoms with different spins in seekpath
                        atoms.add(species[_index]+str(_atom_count_outer[species[_index]]+1))
                        spin_structure[species[_index]+str(_atom_count_outer[species[_index]]+1)] = {"frac_coord": [frac_coord[_index]], "vec": [float(VECTR[_start-1][1]), float(VECTR[_start-1][2]), float(VECTR[_start-1][3])]}
                    else:
                        _atom_count_inner[species[_index]] += 1
                        atomic_numbers[_index] += _atom_count_outer[species[_index]]*1000 # to distinguish same atoms with different spins in seekpath
                        spin_structure[species[_index]+str(_atom_count_outer[species[_index]]+1)]["frac_coord"].append(frac_coord[_index])
                _atom_count_outer[species[_index]] += 1
                _start = _end+2
        env["lspinorb"] = False
        env["nspin"] = 1
        env["nbnd"] = 0
        env["ecutwfc"] = 0
        env["ecutrho"] = 0
        env["time_reversal"] = True 
        with open("../../settings/elements.json", "r") as f:
            elements = json.load(f)
        for atom in atoms:
            element = re.match(r"\D+",atom).group()  # double count exists
            SOC = elements[element]["params"]["SOC"]
            if SOC == "fr":
                env["lspinorb"] = True
                env["nspin"] = 4 # to specify the quantization axis
                env["time_reversal"] = False
            pseudo = elements[element]["pseudopotential"]["PBE"][SOC]["ONCV"] # LDA / PAW
            env[atom] = {"pseudo": pseudo, "Hubbard": elements[element]["params"]["Hubbard"], "starting_magnetization": 0}
            env["nbnd"] += pseudo["nwfc"]
            env["ecutwfc"] = max(env["ecutwfc"], pseudo["ecutwfc"])
            env["ecutrho"] = max(env["ecutrho"], pseudo["ecutwfc"]*pseudo["dual"])
        if spin_calc:
            parallel = False
            for i,j in itertools.permutations(spin_structure, 2):
                parallel = parallel or np.all(np.abs(np.cross(np.array(i["vec"]),np.array(j["vec"])))<1e-5) # VESTA VECTR has 1e-5 significatn digits
            if (not env["lspinorb"]) and parallel:
                env["nspin"] = 2
            else:
                env["nspin"] = 4
                env["time_reversal"] = False
        skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords, atomic_numbers), with_time_reversal=env["time_reversal"], reference_distance=0.025) #setting
        env["avec"] = skp["primitive_lattice"]
        env["bvec"] = skp["reciprocal_primitive_lattice"]
        #env["atom"] = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
        #env["nat"] = len(skp["primitive_types"])
        #env[atom]["starting_magnetization"]
        #typ(mag) → ntyp?
        #env["pos"] = skp["primitive_positions"]
        env["nk"] = np.round(np.linalg.norm(bvec, axis=0)/ 0.2) #setting
        #nks
        #kpath
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            final = skp["explicit_segments"][ipath][1] - 1
        print("%5d %8s %10.5f %10.5f %10.5f %8s %10.5f %10.5f %10.5f" % (
            final - start + 1,
            skp["explicit_kpoints_labels"][start],
            skp["explicit_kpoints_rel"][start][0],
            skp["explicit_kpoints_rel"][start][1],
            skp["explicit_kpoints_rel"][start][2],
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]))
        print(len(skp["explicit_kpoints_rel"]), file=f)
        for ik in range(len(skp["explicit_kpoints_rel"])):
            print(" %f %f %f 1.0" % (
                skp["explicit_kpoints_rel"][ik][0],
                skp["explicit_kpoints_rel"][ik][1],
                skp["explicit_kpoints_rel"][ik][2]),
                file=f)

        #press

    return env