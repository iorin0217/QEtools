import pymatgen
import json
from ase_espresso import ASEspresso
import pickle
import sys

with open(config.json, "r") as f:
    config = json.load(f)

ciffile = (sys.argv)[1]
structure = pymatgen.Structure.from_file(ciffile)
structure.remove_oxidation_states()
atoms = pymatgen.io.ase.get_atoms(structure)

psdir = config["psdir"]
with open(psdir, "r") as f:
    psdict = json.load(f)

qedir = config["qedir"]
input_data =

calc = ASEspresso(input_data)
atoms.set_calculator(calc)
atoms.get_potential_energy()
fermi_level = calc.get_fermi_level()

input_data['control'].update(
    {'calculation': 'bands', 'restart_mode': 'restart', 'verbosity': 'high'})
calc.set(kpts={ < your Brillouin zone path > }, input_data=input_data)
calc.calculate(atoms)

bs = calc.band_structure()
bs.reference = fermi_energy
bs.plot()
