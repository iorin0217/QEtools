import os
import sys
from pathlib import Path
import json
sys.path.append(
    str(Path.joinpath(Path(__file__).parent, "../../").resolve()))  # noqa
from jobber import submit  # noqa
from qeio.env import Env  # noqa
from qeio.create import PW  # noqa

# predefinition
with open(Path.joinpath(Path(__file__).parent, "../../settings/variables.json").resolve(), "r") as f:
    variables = json.load(f)
input_file = sys.argv[1]
outpath = str((Path(input_file).parent).resolve())
env = Env.from_cif(input_file, variables)
# workflow
# pressure
for p in [0, 5, 10]:
    env.extfields["press"] = p
    vc_relax = PW(outpath, env, variables, calculation="vc-relax")
# temperature
# TODO : check imaginary phonon, QHA
submit(job)
