# https://gitlab.com/QEF/q-e/-/blob/develop/PHonon/PH/matdyn.f90
# asr = 'crystal' only
# loto_2d, fd is not implemented
# TODO : fldyn


class Matdyn:
    def __init__(self, outpath, env, variables={}, task="phband"):
        '''
        main usage
            phband = Matdyn(outpath, variables, task="phband")
                phband.in is created
                sub phband.command
        input
            outpath : path
            variables : dict from variables.json
            task = phband/phdos/alpha2f
        '''
        self.task = task
        self.nk = variables["nk"]
        # create input file
        self._save_in(env, outpath)
        # run command
        self.command = [["matdyn.x", f" -in {task}.in | tee {task}.out"]]

    def _save_in(self, env, outpath):
        # calculation type
        base = ["&INPUT", "asr = crystal", "flfrc = 'ph.ifc'"]
        if self.task == "phband":
            kpath = [f"{len(env.bandpath)}"] + \
                [f"{kcoord[0]} {kcoord[1]} {kcoord[2]} 1.0" for kcoord in env.bandpath]
            phband_in = base + ["q_in_band_form = .true.", "/"] + kpath
        elif self.task == "phdos":
            phdos_in = base + ["dos = .true.", f"nq1 = {self.nq[0] * 2}",
                               f"nq2 = {self.nq[1] * 2}", f"nq3 = {self.nq[2] * 2}", "/"]
        elif self.task == "alpha2f":
            alpha2f_in = base + ["dos = .true.", "la2F = .true.", f"nq1 = {self.nq[0] * 2}",
                                 f"nq2 = {self.nq[1] * 2}", f"nq3 = {self.nq[2] * 2}", "/"]
        nline = repr('\n')
        eval(
            f"print(*{self.task}_in, sep={nline}, end={nline}, file=open('{outpath}/{self.task}.in','w'))")
