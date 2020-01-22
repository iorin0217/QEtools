# https://www.quantum-espresso.org/Doc/INPUT_PW.html

def create_pw_in(path, atoms, calculation="scf", restart_mode="from_scratch"):
    control_noend = ["&CONTROL", f"calculation = '{calculation}'", f"restart_mode = '{restart_mode}'", f"pseudo_dir = '{atoms["psdir"]}'", "/"]
    systems_noend = ["&SYSTEM", "ibrav = 0", f"nat = {atoms['nat']}", f"ntyp = {atoms['ntyp']}", f"ecutwfc = {atoms['ecutwfc']}", f"ecutrho = {atoms['ecutrho']}", "occupations = 'tetrahedra_opt'"]
    atomic_parameters = []
    # spin
    # hubbard
    electrons = ["&ELECTRONS", f"conv_thr = {float(atoms["nat"])*1.0e-10}", "mixing_beta = 0.2", "diagonalization = 'david'", "/"]
    cell_parameters = [" %f %f %f" % (i[0],i[1],i[2]) for i in avec]
    atomic_species = [" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, oncv[str(ityp)]["filename"]) for ityp in typ]
    atomic_positions = [" %s %f %f %f" % (atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]) for iat in range(nat)]
    if calculation == "scf":
        kpoints = ["K_POINTS automatic", f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        scf_in = control_noend + ["/"] + systems_noend + atomic_parameters + ["/"] + electrons + ["CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "vcrelax":
        conv_thr = ["etot_conv_thr = 1.0e-5", "forc_conv_thr = 1.0e-4", "/"]
        ions = ["&IONS", "ion_dynamics = 'bfgs'", "/"]
        cell = ["&CELL", "cell_dynamics = 'bfgs'", f"press = {atoms['press']}", "/"]
        kpoints = ["K_POINTS automatic", f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        vcrelax_in = control_noend + conv_thr + systems_noend + atomic_parameters + ["/"] + electrons + ions + cell + ["CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "nscf":
        nbnd = [f"nbnd = {atoms['nwfc']}", "/"]
        kpoints = ["K_POINTS automatic", f"{int(nk[0])*2} {int(nk[1])*2} {int(nk[2])*2} 0 0 0"]
        nscf_in = control_noend + ["/"] + systems_noend + atomic_parameters + nbnd + electrons + ["CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "bands":
        nbnd = [f"nbnd = {atoms['nwfc']}", "/"]
        #kpath
        nscf_in = control_noend + ["/"] + systems_noend + atomic_parameters + nbnd + electrons + ["CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpath
    print(eval(f"*{calculation}_in, sep='\n', end='\n', file=open(f'{path}/{calculation}.in','w')"))










