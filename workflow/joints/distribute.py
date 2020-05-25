class Distribute():
    def __init__(self, outpath):
        self.commands = ["update_params.py", f"{outpath}/"]
        self.task = task

    def __call__(self, calcs):
        if self.task == "structure":
            env =
            variables =
        for calc in calcs:
            calc(env, variables)


if __name__ == "__main__":
    '''
    python update_params.py job1.pkl
        job1.pkl : (update, [nscf, ph])
    '''
    import pickle
    import sys
    pkl = pickle.load(open(sys.argv[1], 'rb'))
    pkl[0](pkl[1])
