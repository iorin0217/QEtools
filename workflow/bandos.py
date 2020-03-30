variables = {"reference_distance": 0.025, "dk_grid": 0.2, "occupations": "tetrahedra_opt",
             "diagonalization": "cg", "mixing_beta": 0.2, "threshold": 1.0e-10, "functional": "PBE", "pseudo_dir": "/home/CMD35/cmd35stud07/QEtools/settings/pseudos"}
path = "/home/CMD35/cmd35stud07/experiments/Ba2RhO4/fr"
env = create_env(
    "/home/CMD35/cmd35stud07/experiments/Ba2RhO4/Ba2RhO4_experiment.cif", variables)
create_pw_in(path, env, variables, calculation="scf")
create_pw_in(path, env, variables, calculation="nscf")
create_pw_in(path, env, variables, calculation="bands")
create_band_in(path)
# job
# parse
efermi = 0
create_dos_in(path, efermi)
create_projwfc_in(path, efermi)
