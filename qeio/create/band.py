# https://www.quantum-espresso.org/Doc/INPUT_BANDS.html


class Band:
    def __init__(self, outpath, env):
        '''
        main usage
            band = Band(outpath, env)
                band*.in is created
                sub band.command
        input
            outpath : path
            env : Env
        '''
        if env.nspin == 1:
            self.fnames = ["band"]
        elif env.nspin == 2:
            self.fnames = ["band_up", "band_dn"]
        elif env.nspin == 4:
            self.fnames = ["band", "band_sx", "band_sy", "band_sz"]
        # create input file
        self._save_in(outpath)
        # run command
        self.command = [["bands.x",
                         f" -in {fname}.in | tee {fname}.out"] for fname in self.fnames]

    def _save_in(self, outpath):
        band_temp = ["&BANDS"]
        band_in = band_temp + ["filband = 'band.out'", "/"]
        band_up_in = band_temp + \
            ["filband = 'band_up.out'", "spin_component = 1", "/"]
        band_dn_in = band_temp + \
            ["filband = 'band_dn.out'", "spin_component = 2", "/"]
        band_sx_in = band_temp + \
            ["filband = 'band_sx.out'", "lsigma(1) = .true.", "/"]
        band_sy_in = band_temp + \
            ["filband = 'band_sy.out'", "lsigma(1) = .true.", "/"]
        band_sz_in = band_temp + \
            ["filband = 'band_sz.out'", "lsigma(1) = .true.", "/"]
        nline = repr('\n')
        for fname in self.fnames:
            eval(
                f"print(*{fname}_in, sep={nline}, end={nline}, file=open('{outpath}/{fname}.in','w'))")
