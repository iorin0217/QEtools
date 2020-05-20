# https://gitlab.com/QEF/q-e/blob/63daaaeec919bb4c7efb3a5cd3c8312e53993b13/PHonon/PH/lambda.f90


class Lambda:
    def __init__(self, outpath):
        '''
        main usage
            lambda_ = Lambda(outpath, variables)
                lambda.in is created
                sub lambda.command
        input
            outpath : path
        '''
        # TODO : search phonon output in outpath and read grids & weights
        self._save_in(outpath)
        # run command
        self.command = [["lambda.x", " -in lambda.in | tee lambda.out"]]

    def _save_in(self, outpath):
        lambda_in =
        print(*lambda_in, sep='\n', end='\n',
              file=open(f'{outpath}/lambda.in', 'w'))
