from pathlib import Path
import subprocess
import json
# TODO : job dependent


class Submit():
    '''
    main usage
        submit = Submit(job, workdir)
        submit()
    input
        job : [(scf1, scf2), nscf]
            {scf1, scf2} : set of array job
            [scf, nscf] : list of dependent job
            scf : PW
        workdir : path
    job type is read from "server" in paths.json
        - qsub
        - local
    job command is read from instance like scf.command
    '''

    def __init__(self, job, workdir):
        self.workdir = workdir
        with open(Path.joinpath(Path(__file__).parent, "../settings/paths.json").resolve(), "r") as f:
            paths = json.load(f)
        self.QExPy = paths["QExPy"]
        self.QE = paths["QE"]
        self.mpirun = paths["mpirun"]
        self.python = paths["python"]
        self.server = paths["server"]
        self.jobchain = []
        if self.server == "docker":
            self.sh = "./"
            self._docker_jobchain(job)
        elif self.server == "RB":
            self.sh = "qsub"
            self._RB_jobchain(job)

    def __call__(self):
        self.sh
        self.jobchain
        print("job submitted")

    def _commander(self, exe):
        _commands = []
        for task in exe:
            if task.command[0].split(".")[-1] == "x":
                _commands.append(
                    f"{self.mpirun} - np {self.np} {self.QE}/{task.command[0]} -npool {self.npool} {task.commdand[1]}")
            elif task.command[0].split(".")[-1] == "py":
                _commands.append(
                    f"{self.python} {self.QExPy}/{task.command[1]}")
        return _commands

    def _RB_jobchain(self, job):
        with open(f"{self.workdir}/RBconfig.json", "r") as f:
            RBconfig = json.load(f)
        header1 = ["#!/bin/bash -f", f"#PBS -q {RBconfig['qname']}", f"#PBS -l select=1:ncpus={self.np}:mpiprocs={self.np}:ompthreads=1",
                   f"#PBS -l walltime={RBconfig['walltime']}", f"#PBS -W group_list={RBconfig['group']}", "#PBS -m abe", f"#PBS -M {RBconfig['mail']}"]
        header2 = ["cd $PBS_O_WORKDIR"] + \
            [f"module load {m}" for m in RBconfig['module']]
        num = 0
        _commands_dep = []
        for exes in job:
            if isinstance(exes, set):
                if _commands_dep:
                    print(*(header1 + _commands_dep), sep="\n", end="\n",
                          file=open(f"{self.workdir}/job_{num}.sh", "w"))  # write upto here
                    self.jobchain.append(f"job_{num}")
                    num += 1
                    _commands_dep = []
                _commands_array = []  # array job
                for exe in exes:
                    _commands_array.extend(self._commander(exe))
                print(*(header1 + header_array + header2 + _commands_array), sep="\n",
                      end="\n", file=open(f"{self.workdir}/job_{num}.sh", "w"))
                self.jobchain.append(f"job_{num}")  # TODO : array header
                num += 1
            else:
                _commands_dep.extend(self._commander(exes))
        if _commands_dep:
            print(*(header1 + header2 + _commands_dep), sep="\n", end="\n",
                  file=open(f"{self.workdir}/job_{num}.sh", "w"))
            self.jobchain.append(f"job_{num}")

    def _docker_jobchain(self, job):
        header = ["!/bin/bash"]
        self.jobchain = ["job_0"]  # array job is serialized
        _commands_all = []
        for exes in job:
            if isinstance(exes, set):
                for exe in exes:
                    _commands_all.extend(self._commander(exe))
            else:
                _commands_all.extend(self._commander(exes))
        print(*(header + _commands_all), sep="\n", end="\n",
              file=open(f"{self.workdir}/job_0.sh", "w"))
