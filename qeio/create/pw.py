# https://www.quantum-espresso.org/Doc/INPUT_PW.html

def create_pw_in(path, atoms, calculation="scf", restart_mode="from_scratch"):
    control_noend = ["&CONTROL", f"calculation = '{calculation}'", f"restart_mode = '{restart_mode}'", f"pseudo_dir = '{atoms["psdir"]}'", "/"]
    systems_noend = ["&SYSTEM", "ibrav = 0", f"nat = {atoms['nat']}", f"ntyp = {atoms['ntyp']}", f"ecutwfc = {atoms['ecutwfc']}", f"ecutrho = {atoms['ecutrho']}", "occupations = 'tetrahedra_opt'"]
    spin = []
    if atoms['nspin'] == 2:
        spin = ["nspin = 2"] + [f"starting_magnetization({i}) = {atoms['starting_magnetization'][i]}" for i in atoms['starting_magnetization']]
    elif atoms['nspin'] == 4:
        spin = ["noncolin = .true.", "lspinorb = .true."] + [f"starting_magnetization({i}) = {atoms['starting_magnetization'][i]}, angle1({i}) = {atoms['angle1'][i]}, angle2({i}) = {atoms['angle2'][i]}" for i in atoms['starting_magnetization']]
    hubbard = []
    if atoms['lda_plus_u'] == True:
        hubbard = ["lda_plus_u = .true."] + [f"Hubbard_U({i}) = {atoms['Hubbard_U'][i]}" for i in atoms['Hubbard_U']]
    atomic_parameters = spin + hubbard
    electrons = ["&ELECTRONS", f"conv_thr = {float(atoms["nat"])*1.0e-10}", "mixing_beta = 0.2", "diagonalization = 'david'", "/"]
    cell_parameters = [f"{vec[0]} {vec[1]} {vec[2]}" for vec in avec]
    atomic_species = [f"{ityp} {pymatgen.Element(ityp).atomic_mass} {oncv[str(ityp)]["filename"]}" for ityp in typ]
    atomic_positions = [f"{atom[iat]} {pos[iat][0]} {pos[iat][1]} {pos[iat][2]}" for iat in range(nat)]
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
        kpath = ["K_POINTS crystal"] + [f"{atoms['nks']}"] + [f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in atoms['kpath']]
        nscf_in = control_noend + ["/"] + systems_noend + atomic_parameters + nbnd + electrons + ["CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpath
    print(eval(f"*{calculation}_in, sep='\n', end='\n', file=open(f'{path}/{calculation}.in','w')"))










