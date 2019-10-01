import pymatgen
import json
from ase_espresso import ASEspresso
import pickle

ciffile =
structure = pymatgen.Structure.from_file(ciffile)
atoms = pymatgen.io.ase.get_atoms(structure)

psdir =
with open(psdir, "r") as f:
    oncv = json.load(f)

qedir =
input_data = { < your input data > }

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
