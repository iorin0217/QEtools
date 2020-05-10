from qeio.env import Env
from copy import deepcopy
variables = {}
outpath = ""
cif = "/home/CMD35/cmd35stud07/experiments/Sr2RhO4/Sr2RhO4_mp-757102_conventional_standard.cif"
env = Env.from_cif(cif)
print(env.bvec)
env2 = deepcopy(env)
env2.bvec = 0
print(env.bvec)
print(env2.bvec)
