# https://www.quantum-espresso.org/Doc/INPUT_PW.html


def create_pw_in(path, env, calculation="scf", restart_mode="from_scratch"):
    # control_noend = ["&CONTROL", f"calculation = '{calculation}'", f"restart_mode = '{restart_mode}'", f"pseudo_dir = '{env["psdir"]}'", "/"]
    systems_noend = ["&SYSTEM", "ibrav = 0", f"nat = {env['nat']}", f"ntyp = {env['ntyp']}",
                     f"ecutwfc = {env['ecutwfc']}", f"ecutrho = {env['ecutrho']}", "occupations = 'tetrahedra_opt'"]

    spin = []
    if env['nspin'] == 2:
        #spin = ["nspin = 2"] + [f"starting_magnetization({i}) = {env['starting_magnetization'][i]}" for i in env['starting_magnetization']]
    elif env['nspin'] == 4:
        #spin = ["noncolin = .true."] + [f"starting_magnetization({i}) = {env['starting_magnetization'][i]}, angle1({i}) = {env['angle1'][i]}, angle2({i}) = {env['angle2'][i]}" for i in env['starting_magnetization']]
    soc = []
    if env['soc'] == True:
        soc = ["lspinorb = .true."]
    hubbard = []
    if env['lda_plus_u'] == True:
        #hubbard = ["lda_plus_u = .true."] + [f"Hubbard_U({i}) = {env['Hubbard_U'][i]}" for i in env['Hubbard_U']]

    electrons = ["&ELECTRONS", f"conv_thr = {float(env["nat"])*1.0e-10}", "mixing_beta = 0.2", "diagonalization = 'david'", "/"]
    #cell_parameters = [f"{vec[0]} {vec[1]} {vec[2]}" for vec in avec]
    # atomic_species = [f"{ityp} {pymatgen.Element(ityp).atomic_mass} {oncv[str(ityp)]["filename"]}" for ityp in typ]
    #atomic_positions = [f"{atom[iat]} {pos[iat][0]} {pos[iat][1]} {pos[iat][2]}" for iat in range(nat)]

    if calculation == "scf":
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        scf_in = control_noend + ["/"] + systems_noend + spin + soc + hubbard + ["/"] + electrons + ["CELL_PARAMETERS angstrom"] + \
            cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + \
            ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "vcrelax":
        conv_thr = ["etot_conv_thr = 1.0e-5", "forc_conv_thr = 1.0e-4", "/"]
        ions = ["&IONS", "ion_dynamics = 'bfgs'", "/"]
        cell = ["&CELL", "cell_dynamics = 'bfgs'",
                f"press = {env['press']}", "/"]
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        vcrelax_in = control_noend + conv_thr + systems_noend + spin + soc + hubbard + ["/"] + electrons + ions + cell + [
            "CELL_PARAMETERS angstrom"] + cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "nscf":
        nbnd = [f"nbnd = {env['nbnd']}", "/"]
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])*2} {int(nk[1])*2} {int(nk[2])*2} 0 0 0"]
        nscf_in = control_noend + ["/"] + systems_noend + spin + soc + hubbard + nbnd + electrons + ["CELL_PARAMETERS angstrom"] + \
            cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + \
            ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpoints
    elif calculation == "bands":
        nbnd = [f"nbnd = {env['nbnd']}", "/"]
        kpath = ["K_POINTS crystal"] + [f"{env['nks']}"] + [
            f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in env['kpath']]
        nscf_in = control_noend + ["/"] + systems_noend + spin + soc + hubbard + nbnd + electrons + ["CELL_PARAMETERS angstrom"] + \
            cell_parameters + ["ATOMIC_SPECIES"] + atomic_species + \
            ["ATOMIC_POSITIONS crystal"] + atomic_positions + kpath

    print(eval(
        f"*{calculation}_in, sep='\n', end='\n', file=open(f'{path}/{calculation}.in','w')"))
