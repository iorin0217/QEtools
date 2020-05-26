# https://gitlab.com/QEF/q-e/blob/develop/PHonon/PH/q2r.f90
# zasr = 'crystal' only
# loto_2d is not implemented


class Q2R:
    def __init__(self, outpath, variables={}, task="ph"):
        '''
        main usage
            q2r = Q2R(outpath, variables, task)
                q2r.in is created
                sub q2r.command
        input
            outpath : path
            variables : dict from variables.json
            task : ph/elph
        '''
        self.task = task
        self.occupations = variables['occupations'] if variables else 'tetrahedra_opt'
        # TODO : validation
        # create input file
        self._save_in(outpath)
        # run command
        self.command = [["q2r.x", " -in q2r.in | tee q2r.out"]]

    def _save_in(self, outpath):
        # task type
        base = ["&INPUT", "fildyn = 'ph.dyn'",
                "flfrc = 'ph.ifc'", "zasr = 'crystal'"]
        if self.task == "ph":
            q2r_in = base + ["/"]
        elif self.task == "elph":
            q2r_in = base + ["la2F = .true.", "/"]
        print(*q2r_in, sep='\n', end='\n', file=open(f'{outpath}/q2r.in', 'w'))
