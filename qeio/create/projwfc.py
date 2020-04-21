# https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html


class Projwfc:
    def __init__(self, outpath, variables={}, task="pdos"):
        '''
        main usage
            pdos = Projwfc(outpath)
                pdos.in is created
                sub pdos.command
        input
            outpath : path
            variables : dict from variables.json
            task : pdos/fatband
        '''
        self.task = task
        self.deltae = variables['deltae'] if variables else 0.01
        # create input file
        self._save_in(outpath)
        # run command
        self.command = [f"projwfc.x -in {self.task}.in | tee {self.task}.out"]

    def _save_in(self, outpath):
        projwfc_temp = ["&PROJWFC"]
        pdos_in = projwfc_temp + [f"deltae = {self.deltae}", "/"]
        fatband_in = projwfc_temp + ["lsym=.false.", "/"]
        nline = repr('\n')
        eval(
            f"print(*{self.task}_in, sep={nline}, end={nline}, file=open('{outpath}/{self.task}.in','w'))")
