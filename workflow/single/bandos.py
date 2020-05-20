import sys
from pathlib import Path
import json
sys.path.append(
    str(Path.joinpath(Path(__file__).parent, "../../").resolve()))  # noqa
from jobber import submit  # noqa
from qeio.env import Env  # noqa
from qeio.create import PW, Projwfc, Band  # noqa

# predefinition
with open(Path.joinpath(Path(__file__).parent, "../../settings/variables.json").resolve(), "r") as f:
    variables = json.load(f)
input_file = sys.argv[1]
outpath = str((Path(input_file).parent).resolve())
env = Env.from_cif(input_file, variables)
# workflow
scf = PW(outpath, env, variables, calculation="scf")
nscf = PW(outpath, env, variables, calculation="nscf")
# TODO : fermi_velocity / proj
pdos = Projwfc(outpath, variables, task="pdos")
bands = PW(outpath, env, variables, calculation="bands")
band = Band(outpath, env)
fatband = Projwfc(outpath, variables, task="fatband")
# TODO : job class / backup
job = [scf, nscf, pdos, bands, band, fatband]
submit(job)
