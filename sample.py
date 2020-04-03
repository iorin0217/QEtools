'''
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
np.dot(np.linalg.inv(
    before.T @ skp["primitive_transformation_matrix"]), after.T)
# %%
(trans_r @ abc.T @ np.linalg.inv(trans_f)).T
# %%
np.linalg.inv(trans_r @ abc.T) @ (
    after.T @ skp["primitive_positions"][1] - trans_r @ abc.T @ structure.frac_coords[1])
# %%
spin_ini = np.array([0, 0.5, -1])
# %%
np.linalg.inv(after.T) @ trans_r @ abc.T @ spin_ini
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
spin_structure

# %%
atoms

# %%
skp["primitive_positions"][2]

# %%
structure.frac_coords[2]

# %%
skp["primitive_positions"]

# %%
np.linalg.norm(abc.T @ (np.array([2, -4, 2])/np.linalg.norm(abc, axis=1)))


# %%
abc.T @ np.array([1, 1, 1]) @ np.array([0, 0, 1])

# %%
np.array([2, 2, 2])/np.array([1, 1, 1])
'''
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
from pathlib import Path


def create_env(structure_file, variables, extfields={"press": 0}, constraints={"symm": False}):
    # unsupport the occupation != 1 case
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
    with open("/home/CMD35/cmd35stud07/QEtools/settings/elements.json", "r") as f:
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
    env["nk"] = np.round(np.linalg.norm(
        env["bvec"], axis=0) / variables["dk_grid"])
    # extfields, constraints
    # hdf5

    return env

# %%
# https://www.quantum-espresso.org/Doc/INPUT_PW.html


def create_pw_in(path, env, variables, calculation="scf"):
    # useful
    atom_types = list(set(env['atoms']))
    nat = len(env['atoms'])
    # &CONTROL (no "/")
    control = [
        "&CONTROL", f"calculation = '{calculation}'", f"pseudo_dir = '{variables['pseudo_dir']}'"]
    # &SYSTEM (no "/")
    systems = ["&SYSTEM", "ibrav = 0", f"nat = {nat}", f"ntyp = {len(atom_types)}",
               f"ecutwfc = {env['ecutwfc']}", f"ecutrho = {env['ecutrho']}", f"occupations = '{variables['occupations']}''"]
    spin = []
    if env['nspin'] == 2:
        spin = ["nspin = 2"] + \
            [f"starting_magnetization({i+1}) = {env[atom]['starting_magnetization']}" for i,
             atom in enumerate(atom_types) if np.abs(env[atom]['starting_magnetization']) > 1e-3]  # QE has this resolution
    elif env['nspin'] == 4:
        spin = ["noncolin = .true."] + \
            [f"starting_magnetization({i+1}) = {env[atom]['starting_magnetization']}, angle1({i+1}) = {env[atom]['angle1']}, angle2({i+1}) = {env[atom]['angle2']}" for i,
             atom in enumerate(atom_types) if np.abs(env[atom]['starting_magnetization']) > 1e-3]
    soc = []
    if env['lspinorb']:
        soc = ["lspinorb = .true."]
    hubbard = []
    if env['lda_plus_u']:
        hubbard = ["lda_plus_u = .true."] + ["lda_plus_u_kind = 0"] + [f"Hubbard_U({i+1}) = {env[atom]['Hubbard']['U']}" for i, atom in enumerate(
            atom_types) if env[atom]['Hubbard'] and env[atom]['Hubbard'].get('U')] + [f"Hubbard_J0({i+1}) = {env[atom]['Hubbard']['J']}" for i, atom in enumerate(atom_types) if env[atom]['Hubbard'] and env[atom]['Hubbard'].get('J')]
    SSSH = systems + spin + soc + hubbard
    # &ELECTRONS (no "/")
    electrons = ["&ELECTRONS", f"conv_thr = {float(nat)*variables['threshold']}",
                 f"mixing_beta = {variables['mixing_beta']}", f"diagonalization = '{variables['diagonalization']}'", f"diago_full_acc = {variables['diago_full_acc']}"]
    # CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS
    cell_parameters = ["CELL_PARAMETERS angstrom"] + \
        [f"{vec[0]} {vec[1]} {vec[2]}" for vec in env["avec"]]
    atomic_species = ["ATOMIC_SPECIES"] + \
        [f"{atom} -1 {env[atom]['pseudofile']}" for atom in atom_types]
    atomic_positions = ["ATOMIC_POSITIONS crystal"] + [
        f"{env['atoms'][i]} {env['pos'][i][0]} {env['pos'][i][1]} {env['pos'][i][2]}" for i in range(nat)]
    CAA = cell_parameters + atomic_species + atomic_positions

    if calculation == "scf":
        kpoints = ["K_POINTS automatic",
                   f"{int(env['nk'][0])} {int(env['nk'][1])} {int(env['nk'][2])} 0 0 0"]
        scf_in = control + ["/"] + SSSH + ["/"] + \
            electrons + ["/"] + CAA + kpoints
    elif calculation == "vcrelax":
        conv_thr = ["etot_conv_thr = 1.0e-5", "forc_conv_thr = 1.0e-4"]
        ions = ["&IONS", "ion_dynamics = 'bfgs'"]
        cell = ["&CELL", "cell_dynamics = 'bfgs'",
                f"press = {env['press']}"]
        kpoints = ["K_POINTS automatic",
                   f"{int(env['nk'][0])} {int(env['nk'][1])} {int(env['nk'][2])} 0 0 0"]
        vcrelax_in = control + conv_thr + \
            ["/"] + SSSH + ["/"] + electrons + ["/"] + \
            ions + ["/"] + cell + ["/"] + CAA + kpoints
    elif calculation == "nscf":
        kpoints = ["K_POINTS automatic",
                   f"{int(env['nk'][0])*2} {int(env['nk'][1])*2} {int(env['nk'][2])*2} 0 0 0"]
        nscf_in = control + ["/"] + SSSH + \
            ["/"] + electrons + ["/"] + CAA + kpoints
    elif calculation == "bands":
        kpath = ["K_POINTS crystal"] + [f"{len(env['bandpath'])}"] + [
            f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in env['kpath']]
        bands_in = control + ["/"] + SSSH + \
            ["/"] + electrons + ["/"] + CAA + kpath
    tmp = repr('\n')
    eval(
        f"print(*{calculation}_in, sep={tmp}, end={tmp}, file=open('{path}/{calculation}.in','w'))")

# %%


def create_projwfc_in(path, efermi, emin=-5, emax=5, deltae=0.01):
    projwfc_temp = [
        "&PROJWFC", f"emin = {efermi+emin}", f"emax = {efrmi+emax}", f"deltae = {deltae}", "/"]
    print(*projwfc_temp, sep="\n", end="\n",
          file=open(f"{path}/projwfc.in", "w"))

# %%


def create_band_in(path, nspin=1, plot_2d=".false."):
    band_temp = ["&BANDS", f"plot_2d={plot_2d}"]
    if nspin == 1:
        band_temp += ["filband = 'band.out'", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band.in", "w"))
    elif nspin == 2:
        band_temp_cp = band_temp[:]
        band_temp += ["filband = 'band_up.out'", "spin_component = 1", "/"]
        band_temp_cp += ["filband = 'band_dn.out'", "spin_component = 2", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_up.in", "w"))
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_dn.in", "w"))
    elif nspin == 4:
        band_temp += ["filband = 'band_s.out'", "lsigma = .true.", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_s.in", "w"))
# %%


def create_dos_in(path, efermi, emin=-10, emax=10, deltae=0.05):
    dos_temp = ["&DOS", f"emin = {efermi+emin}",
                f"emax = {efrmi+emax}", f"deltae = {deltae}", "/"]
    print(*dos_temp, sep="\n", end="\n", file=open(f"{path}/dos.in", "w"))


# %%
variables = {"reference_distance": 0.025, "dk_grid": 0.2, "occupations": "tetrahedra_opt", "diago_full_acc": ".true.",
             "diagonalization": "david", "mixing_beta": 0.2, "threshold": 1.0e-12, "functional": "PBE", "pseudo_dir": "/home/CMD35/cmd35stud07/QEtools/settings/pseudos"}
# %%
for cif in Path("/home/CMD35/cmd35stud07/experiments/").glob("Sr2*/*.cif"):
    path = cif.parent / "fr"
    env = create_env(cif, variables)
    create_pw_in(path, env, variables, calculation="scf")
    create_pw_in(path, env, variables, calculation="nscf")
    create_pw_in(path, env, variables, calculation="bands")
    create_band_in(path)
# %%


# %%
efermi = 0
create_dos_in(path, efermi)
create_projwfc_in(path, efermi)
