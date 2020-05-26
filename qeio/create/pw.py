# https://www.quantum-espresso.org/Doc/INPUT_PW.html
import numpy as np


class PW:
    def __init__(self, outpath, env, variables={}, calculation="scf"):
        '''
        main usage
            scf = PW(outpath, env, variables, calculation="scf")
                $calculation.in is created
                sub scf.command
        input
            outpath : path
            env : Env
            variables : dict from variables.json
            calculation = scf/nscf/bands/vcrelax
        '''
        # &CONTROL basic
        self.calculation = calculation
        self.pseudo_dir = variables['pseudo_dir'] if variables else './'
        # TODO : restart_mode
        # &SYSTEM basic
        self.ecutwfc = env.ecutwfc
        self.ecutrho = env.ecutrho
        self.degauss = variables['degauss'] if variables else 0.01
        # TODO : ibrav neq 0
        '''
        self.nbnd = sum(
            list([env.elements[atom]["valence"][1] * 2 + 1 for atom in env.atoms]))
        '''
        self.nbnd = sum(list(
            [sum(list([valence[1] * 2 + 1 for valence in env.elements[atom]["valence"]])) for atom in env.atoms]))
        # &ELECTRONS basic
        self.conv_thr = len(
            env.atoms) * variables['threshold'] if variables else len(env.atoms) * 1.0e-12
        self.mixing_beta = variables['mixing_beta'] if variables else 0.2
        self.diagonalization = variables['diagonalization'] if variables else 'david'
        self.diago_full_acc = variables['diago_full_acc'] if variables else '.true.'
        # CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS are determined by only reading env
        # calculation type dependent
        self.nk = variables["nk"] if variables else [8, 8, 8]
        self.occupations = variables['occupations'] if variables else 'tetrahedra_opt'
        self.etot_conv_thr = variables['etot_conv_thr'] if variables else 1.0e-5
        self.forc_conv_thr = variables['forc_conv_thr'] if variables else 1.0e-4
        self.ion_dynamics = variables['ion_dynamics'] if variables else 'bfgs'
        self.cell_dynamics = variables['cell_dynamics'] if variables else 'bfgs'
        # create input file
        self._save_in(outpath, env)
        # run command
        self.command = [["pw.x",
                         f" -in {self.calculation}.in | tee {self.calculation}.out"]]

    def _save_in(self, outpath, env):
        _atom_types, _nat = list(set(env.atoms)), len(
            env.atoms)  # useful for iteration
        # &CONTROL
        _calc = "vc-relax" if self.calculation == "vcrelax" else self.calculation
        control = [
            "&CONTROL", f"calculation = '{_calc}'", f"pseudo_dir = '{self.pseudo_dir}'"]
        # &SYSTEM
        _systems = ["&SYSTEM", "ibrav = 0", f"nat = {_nat}", f"ntyp = {len(_atom_types)}",
                    f"ecutwfc = {self.ecutwfc}", f"ecutrho = {self.ecutrho}"]
        _spin = []
        if env.nspin == 2:
            _spin = ["nspin = 2"] + [f"starting_magnetization({i+1}) = {env.elements[atom]['starting_magnetization']}" for i,
                                     atom in enumerate(_atom_types) if np.abs(env.elements[atom]['starting_magnetization']) > 1e-3]  # QE has this resolution
        elif env.nspin == 4:
            _spin = ["noncolin = .true."] + \
                [f"starting_magnetization({i+1}) = {env.elements[atom]['starting_magnetization']}, angle1({i+1}) = {env.elements[atom]['angle1']}, angle2({i+1}) = {env.elements[atom]['angle2']}" for i,
                 atom in enumerate(_atom_types) if np.abs(env.elements[atom]['starting_magnetization']) > 1e-3]
        _soc = []
        if env.lspinorb:
            _soc = ["lspinorb = .true."]
        _hubbard = []
        if env.lda_plus_u:
            _hubbard = ["lda_plus_u = .true."] + ["lda_plus_u_kind = 1"] + [f"Hubbard_U({i+1}) = {env.elements[atom]['Hubbard']['U']}" for i, atom in enumerate(
                _atom_types) if env.elements[atom]['Hubbard'] and env.elements[atom]['Hubbard'].get('U')] + [f"Hubbard_J0({i+1}) = {env.elements[atom]['Hubbard']['J']}" for i, atom in enumerate(_atom_types) if env.elements[atom]['Hubbard'] and env.elements[atom]['Hubbard'].get('J')]
        SSSH = _systems + _spin + _soc + _hubbard
        OCC = [f"occupations = '{self.occupations}'"] if self.occupations == "tetrahedra_opt" else [
            f"occupations = 'smearing'", f"degauss = {self.degauss}"]
        # &ELECTRONS
        electrons = ["&ELECTRONS", f"conv_thr = {self.conv_thr}", f"mixing_beta = {self.mixing_beta}",
                     f"diagonalization = '{self.diagonalization}'", f"diago_full_acc = {self.diago_full_acc}"]
        # CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS
        _cell_parameters = ["CELL_PARAMETERS angstrom"] + \
            [f"{vec[0]} {vec[1]} {vec[2]}" for vec in env.crystal.lattice.matrix]
        _atomic_species = ["ATOMIC_SPECIES"] + \
            [f"{atom} -1 {env.elements[atom]['pseudofile']}" for atom in _atom_types]
        _atomic_positions = ["ATOMIC_POSITIONS crystal"] + [
            f"{env.atoms[i]} {env.crystal.frac_coords[i][0]} {env.crystal.frac_coords[i][1]} {env.crystal.frac_coords[i][2]}" for i in range(_nat)]
        CAA = _cell_parameters + _atomic_species + _atomic_positions
        # calculation type
        if self.calculation == "scf":
            kpoints = ["K_POINTS automatic",
                       f"{int(self.nk[0])} {int(self.nk[1])} {int(self.nk[2])} 0 0 0"]
            scf_in = control + ["/"] + SSSH + OCC + ["/"] + \
                electrons + ["/"] + CAA + kpoints
        elif self.calculation == "vcrelax":
            conv_thr = [
                f"etot_conv_thr = {self.etot_conv_thr}", f"forc_conv_thr = {self.forc_conv_thr}"]
            ions = ["&IONS", f"ion_dynamics = {self.ion_dynamics}"]
            cell = ["&CELL", f"cell_dynamics = {self.cell_dynamics}",
                    f"press = {env.extfields['press']}"]
            kpoints = ["K_POINTS automatic",
                       f"{int(self.nk[0])} {int(self.nk[1])} {int(self.nk[2])} 0 0 0"]
            vcrelax_in = control + conv_thr + \
                ["/"] + SSSH + [f"occupations = 'smearing'", f"degauss = {self.degauss}", "/"] + electrons + ["/"] + \
                ions + ["/"] + cell + ["/"] + CAA + kpoints
        elif self.calculation == "nscf":
            # TODO : nbnd
            kpoints = ["K_POINTS automatic",
                       f"{int(self.nk[0])*2} {int(self.nk[1])*2} {int(self.nk[2])*2} 0 0 0"]
            nscf_in = control + ["/"] + SSSH + \
                [f"nbnd = {self.nbnd}"] + OCC + ["/"] + \
                electrons + ["/"] + CAA + kpoints
        elif self.calculation == "bands":
            # TODO : nbnd
            kpath = ["K_POINTS crystal"] + [f"{len(env.bandpath)}"] + [
                f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in env.bandpath]
            bands_in = control + ["/"] + SSSH + \
                [f"nbnd = {self.nbnd}", "/"] + electrons + ["/"] + CAA + kpath
        nline = repr('\n')
        eval(
            f"print(*{self.calculation}_in, sep={nline}, end={nline}, file=open('{outpath}/{self.calculation}.in','w'))")
