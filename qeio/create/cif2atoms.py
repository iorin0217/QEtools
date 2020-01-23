import numpy as np
import pymatgen
from pymatgen.core.periodic_table import get_el_sp
import seekpath
import json

def create_atoms(cif,trans=np.eye(3), spins, press, symm):
    atoms = {}
    structure = pymatgen.Structure.from_file(cif)
    # cif ioだと変わる
    structure.remove_oxidation_states()
    skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                        [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                        reference_distance=0.025)
    atoms["avec"] = skp["primitive_lattice"]
    atoms["bvec"] = skp["reciprocal_primitive_lattice"]
    atoms["pos"] = skp["primitive_positions"]
    atoms["nk"] = np.round(np.linalg.norm(bvec, axis=0)/ 0.2)
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
    #nks
    #kpath
    print(len(skp["explicit_kpoints_rel"]), file=f)
    for ik in range(len(skp["explicit_kpoints_rel"])):
        print(" %f %f %f 1.0" % (
            skp["explicit_kpoints_rel"][ik][0],
            skp["explicit_kpoints_rel"][ik][1],
            skp["explicit_kpoints_rel"][ik][2]),
            file=f)
    #ibrav
    #trans
    #press
    atoms["atom"] = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    with open("../../settings/elements.json", "r") as f:
        elements = json.load(f)
    atoms["nat"] = len(skp["primitive_types"])
    #psdir/psfile → nwfc
    #soc, U
    #typ(mag) → ntyp

    return atoms