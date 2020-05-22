import sys
from pathlib import Path
import json
sys.path.append(
    str(Path.joinpath(Path(__file__).parent, "../../").resolve()))  # noqa
from jobber import Submit  # noqa
from qeio.env import Env  # noqa
from qeio.create import PW, PH, Q2R, Matdyn  # noqa

# predefinition
with open(Path.joinpath(Path(__file__).parent, "../../settings/variables.json").resolve(), "r") as f:
    variables = json.load(f)
input_file = sys.argv[1]
outpath = str((Path(input_file).parent).resolve())
env = Env.from_cif(input_file, variables)
# workflow
vc_relax = PW(outpath, env, variables, calculation="vc-relax")
# TODO : metal check, backup
if variables["occupation"] == "tetrahedra_opt":
    ph = PH(outpath, env, variables, task="ph")
    elph = PH(outpath, env, variables, task="elph")
    # TODO : alpha2F, matdyn
    job = [vc_relax, ph, elph]
else:
    # TODO : la2F = .true.
    nscf = PW(outpath, env, variables, calculation="nscf")
    elph = PH(outpath, env, variables, task="elph")
    q2r = Q2R(outpath, variables, task="elph")
    matdyn = Matdyn(outpath, env, variables, task="alpha2f")
    # TODO : lambda
    job = [vc_relax, nscf, elph, q2r, matdyn]
submit = Submit(job, outpath)
submit()
