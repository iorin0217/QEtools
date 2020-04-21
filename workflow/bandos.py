'''
#!/bin/bash -f
#$ -pe smp 8
#$ -cwd
#$ -N job_name
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/pw.x -in scf.in | tee scf.out
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/pw.x -in nscf.in | tee nscf.out
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/fermi_velocity.x -in nscf.in | tee fermi_velocity.out
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/projwfc.x -in pdos.in | tee pdos.out
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/pw.x -in bands.in | tee bands.out
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/bands.x -in band.in
mpirun -np $NSLOTS /home/CMD35/cmd35stud07/QE6.4.1/bin/projwfc.x -in fat.in | tee fat.out
'''
# %%
import sys
from pathlib import Path
sys.path.append(
    str(Path.joinpath(Path(__file__).parent, "../").resolve()))  # noqa
from qeio.env import Env  # noqa
from qeio.create import *  # noqa

variables = {}
outpath = "/home/CMD35/cmd35stud07/QExPy"

cif = "/home/CMD35/cmd35stud07/experiments/Sr2RhO4/Sr2RhO4_mp-757102_conventional_standard.cif"
env = Env.from_cif(cif)

scf = PW(outpath, env, calculation="scf")
nscf = PW(outpath, env, calculation="nscf")
# TODO : fermi_velocity / proj
pdos = Projwfc(outpath, task="pdos")
bands = PW(outpath, env, calculation="bands")
band = Band(outpath, env)
fatband = Projwfc(outpath, task="fatband")
# TODO : job class / backup
job = [scf, nscf, pdos, bands, band, fatband]
