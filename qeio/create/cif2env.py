import numpy as np
import pymatgen
from pymatgen.io.cif import CifParser
from pymatgen.core.periodic_table import get_el_sp
import seekpath
import json
import os

def create_env(cif, spinstructure_file=None, extfields={"press":0}, constrains={"symm":False}):
    env = {}
    structure = CifParser(cif).get_structures(primitive=False)[0] # when we use mcif, primitive=True
    structure.remove_oxidation_states()

    with open("../../settings/elements.json", "r") as f:
        elements = json.load(f)
    env["atom"] = set([str(i) for i in structure.species])
    env["lspinorb"] = False
    env["nspin"] = 1
    env["nbnd"] = 0
    env["ecutwfc"] = 0
    env["ecutrho"] = 0
    env["time_reversal"] = True 
    for atom in env["atom"]:
        SOC = elements[atom]["params"]["SOC"]
        if SOC == "fr":
            env["lspinorb"] = True
            env["nspin"] = 4 # to specify the quantization axis
            env["time_reversal"] = False
        pseudo = elements[atom]["pseudopotential"]["PBE"][SOC]["ONCV"] # LDA / PAW
        env[atom] = {"pseudo": pseudo, "Hubbard": elements[atom]["params"]["Hubbard"], "starting_magnetization": 0}
        env["nbnd"] += pseudo["nwfc"]
        env["ecutwfc"] = max(env["ecutwfc"], pseudo["ecutwfc"])
        env["ecutrho"] = max(env["ecutrho"], pseudo["ecutwfc"]*pseudo["dual"])

    if os.path.splitext(spinstructure_file)[1][1:] == "vesta":
        with open(spinstructure_file) as f:
            spinstructure = f.readlines()
        VECTR = [i.strip().split() for i in spinstructure[int(spinstructure.index("VECTR\n"))+1:int(spinstructure.index("VECTT\n"))] if i.strip().split()[0]!="0"]
        STRUC = [j.strip().split() for j in spinstructure[int(spinstructure.index("STRUC\n"))+1:int(spinstructure.index("THERI 0\n"))] if j.strip().split()[0]!="0"]
        CELLP = [k.strip().split() for k in spinstructure[int(spinstructure.index("CELLP\n"))+1:int(spinstructure.index("STRUC\n"))] if k.strip().split()[0]!="0"]
        #env[atom]["starting_magnetization"]
        #ntyp増やす +1000
        if (not env["lspinorb"]) and parallel:
            env["nspin"] = 2
        else:
            env["nspin"] = 4
            env["time_reversal"] = False
    skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords, [pymatgen.Element(str(spc)).number for spc in structure.species]), with_time_reversal=env["time_reversal"], reference_distance=0.025) #setting
    #env["atom"] = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    #env["nat"] = len(skp["primitive_types"])
    #typ(mag) → ntyp?
    env["avec"] = skp["primitive_lattice"]
    env["bvec"] = skp["reciprocal_primitive_lattice"]
    env["pos"] = skp["primitive_positions"]
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