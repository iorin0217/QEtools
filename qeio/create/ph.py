# https://www.quantum-espresso.org/Doc/INPUT_PH.html


class PH:
    def __init__(self, outpath, env, variables={}, task="ph"):
        '''
        main usage
            scf = PH(outpath, env, variables, task="ph")
                $task.in is created
                sub ph.command
        input
            outpath : path
            env : Env
            variables : dict from variables.json
            task = ph/elph
        '''
        self.task = task
        # check the calc type
        self.functional = variables['functional'] if variables else 'PBE'
        self.pstype = variables['pstype'] if variables else 'US'
        self.lspinorb = env.lspinrob
        self.lda_plus_u = env.lda_plus_u
        # TODO : validate calc type
        # &INPUTPH basic
        self.nq = variables["nq"] if variables else [4, 4, 4]
        self.occupations = variables['occupations'] if variables else 'tetrahedra_opt'
        # for elph & tetrahedra_opt
        self.nk = variables["nk"] if variables else [8, 8, 8]
        # TODO : dynamical setting
        self.nfreq = variables["nfreq"] if variables else 500
        # create input file
        self._save_in(outpath, env)
        # run command
        self.command = [["ph.x",
                         f" -in {self.task}.in | tee {self.task}.out"]]

    def _save_in(self, outpath, env):
        # calculation type
        base = ["&INPUTPH", "fildyn = 'ph.dyn'", "fildvscf = 'dv'", "fildrho = 'drho'",
                "ldisp = .true.", f"nq1 = {self.nq[0]}", f"nq2 = {self.nq[1]}", f"nq3 = {self.nq[2]}"]
        if self.task == "ph":
            if self.occupations == "tetrahedra_opt":
                ph_in = ["phonon"] + base + ["lshift_q = .true.", "/"]
            else:
                ph_in = ["phonon"] + base
        elif self.task == "elph":
            if self.occupations == "tetrahedra_opt":
                elph_in = ["electron-phonon"] + base + ["lshift_q = .true.", "electron_phonon = 'lambda_tetra",
                                                        f"nk1 = {self.nk[0] * 2}", f"nk2 = {self.nk[1] * 2}", f"nk3 = {self.nk[2] * 2}", "/"] + ["&INPUTa2F", f"nfreq = {self.nfreq}", "/"]
            else:
                elph_in = ["electron-phonon"] + base + \
                    ["electron_phonon = 'interpolated'", "/"]
        nline = repr('\n')
        eval(
            f"print(*{self.task}_in, sep={nline}, end={nline}, file=open('{outpath}/{self.task}.in','w'))")
