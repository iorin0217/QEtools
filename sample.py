from qeio.env import Env
from qeio.create import Projwfc
variables = {}
outpath = ""
cif = "/home/CMD35/cmd35stud07/experiments/Sr2RhO4/Sr2RhO4_mp-757102_conventional_standard.cif"
env = Env.from_cif(cif)
pdos = Projwfc("/home/CMD35/cmd35stud07/QExPy")
