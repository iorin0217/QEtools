import sys
from pathlib import Path
import json
sys.path.append(
    str(Path.joinpath(Path(__file__).parent, "../../").resolve()))  # noqa
from qeio.env import Env  # noqa
from qeio.create import PW, PH, Q2R, Matdyn  # noqa

# predefinition
with open(Path.joinpath(Path(__file__).parent, "../../settings/variables.json").resolve(), "r") as f:
    variables = json.load(f)
variables["pstype"] = "US"
#input_file = sys.argv[1]
input_file = "/Users/kurataiori/Master/QExPy/example/Ba2RhO4_ph/Ba2RhO4_experiment.cif"
outpath = str((Path(input_file).parent).resolve())
env = Env.from_cif(input_file, variables)

# workflow
vcrelax = PW(outpath, env, variables, calculation="vcrelax")
# TODO : structure, metal check, backup
if variables["occupations"] == "tetrahedra_opt":
    scf = PW(outpath, env, variables, calculation="scf")
    ph = PH(outpath, env, variables, task="ph")
    elph = PH(outpath, env, variables, task="elph")
    # TODO : alpha2F, matdyn
    job = [vcrelax, scf, ph, elph]
else:
    # TODO : la2F = .true.
    nscf = PW(outpath, env, variables, calculation="nscf")
    elph = PH(outpath, env, variables, task="elph")
    q2r = Q2R(outpath, variables, task="elph")
    matdyn = Matdyn(outpath, env, variables, task="alpha2f")
    # TODO : lambda
    job = [vcrelax, nscf, elph, q2r, matdyn]
