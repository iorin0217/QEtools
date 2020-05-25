class Update():
    def __init__(self, outpath):
        self.commands = ["update.py", f"{outpath}/"]
        self.task = task

    def __call__(self, calcs):
        if self.task == "structure":
            env =
            variables =
        for calc in calcs:
            calc(env, variables)


if __name__ == "__main__":
    '''
    python update.py job1.pkl
        job1.pkl : (update, [nscf, ph])
    '''
    import pickle
    import sys
    pkl = pickle.load(open(sys.argv[1], 'rb'))
    pkl[0](pkl[1])
