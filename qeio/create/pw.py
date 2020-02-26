# https://www.quantum-espresso.org/Doc/INPUT_PW.html


def create_pw_in(path, env, variables, calculation="scf"):
    # useful
    atom_types = list(set(env['atoms']))
    nat = len(env['atoms'])
    # &CONTROL (no "/")
    control = ["&CONTROL", f"calculation = '{calculation}'", f"pseudo_dir = '{env["psdir"]}'"]
    # &SYSTEM (no "/")
    systems = ["&SYSTEM", "ibrav = 0", f"nat = {nat}", f"ntyp = {len(atom_types)}",
               f"ecutwfc = {env['ecutwfc']}", f"ecutrho = {env['ecutrho']}", f"occupations = {variables['occupations']}"]
    spin = []
    if env['nspin'] == 2:
        spin = ["nspin = 2"] + \
            [f"starting_magnetization({i+1}) = {env[atom]['starting_magnetization']}" for i,
             atom in enumerate(atom_types) if env[atom]['starting_magnetization'] != 0]
    elif env['nspin'] == 4:
        spin = ["noncolin = .true."] + \
            [f"starting_magnetization({i+1}) = {env[atom]['starting_magnetization']}, angle1({i+1}) = {env[atom]['angle1']}, angle2({i+1}) = {env[atom]['angle2']}" for i,
             atom in enumerate(atom_types) if env[atom]['starting_magnetization'] != 0]
    soc = []
    if env['lspinorb']:
        soc = ["lspinorb = .true."]
    hubbard = []
    if env['lda_plus_u']:
        hubbard = ["lda_plus_u = .true."] + ["lda_plus_u_kind = 1"] + [f"Hubbard_U({i+1}) = {env[atom]['Hubbard_U']}" for i,
                                                                       atom in enumerate(atom_types) if env[atom]['Hubbard_U']]
    SSSH = systems + spin + soc + hubbard
    # &ELECTRONS (no "/")
    electrons = ["&ELECTRONS", f"conv_thr = {float(nat)*variables['threshold']}",
                 f"mixing_beta = {variables['mixing_beta']}", f"diagonalization = {variables['diagonalization']}"]
    # CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS
    cell_parameters = ["CELL_PARAMETERS angstrom"] + \
        [f"{vec[0]} {vec[1]} {vec[2]}" for vec in env["avec"]]
    atomic_species = ["ATOMIC_SPECIES"] + \
        [f"{atom} -1 {env[atom]['pseudo']}" for atom in atom_types]
    atomic_positions = ["ATOMIC_POSITIONS crystal"] + [
        f"{atoms[i]} {env['pos'][i][0]} {env['pos'][i][1]} {env['pos'][i][2]}" for i in range(nat)]
    CAA = cell_parameters + atomic_species + atomic_positions

    if calculation == "scf":
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        scf_in = control + ["/"] + SSSH + ["/"] + \
            electrons + ["/"] + CAA + kpoints
    elif calculation == "vcrelax":
        conv_thr = ["etot_conv_thr = 1.0e-5", "forc_conv_thr = 1.0e-4"]
        ions = ["&IONS", "ion_dynamics = 'bfgs'"]
        cell = ["&CELL", "cell_dynamics = 'bfgs'",
                f"press = {env['press']}"]
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])} {int(nk[1])} {int(nk[2])} 0 0 0"]
        vcrelax_in = control + conv_thr + \
            ["/"] + SSSH + ["/"] + electrons + ["/"] + \
            ions + ["/"] + cell + ["/"] + CAA + kpoints
    elif calculation == "nscf":
        nbnd = [f"nbnd = {env['nbnd']}"]
        kpoints = ["K_POINTS automatic",
                   f"{int(nk[0])*2} {int(nk[1])*2} {int(nk[2])*2} 0 0 0"]
        nscf_in = control + ["/"] + SSSH + nbnd + \
            ["/"] + electrons + ["/"] + CAA + kpoints
    elif calculation == "bands":
        nbnd = [f"nbnd = {env['nbnd']}"]
        kpath = ["K_POINTS crystal"] + [f"{len(env['bandpath'])}"] + [
            f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in env['kpath']]
        nscf_in = control + ["/"] + SSSH + nbnd + \
            ["/"] + electrons + ["/"] + CAA + kpath

    print(eval(
        f"*{calculation}_in, sep='\n', end='\n', file=open(f'{path}/{calculation}.in','w')"))
