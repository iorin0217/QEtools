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
        job : [scf, nscf]
            scf : PW
        workdir : path
    job type is read from "server" in paths.json
        - qsub
        - local
    job command is read from instance like scf.command
    '''

    def __init__(self, job, workdir, jobname):
        self.workdir = workdir
        self.jobname = jobname
        with open(Path.joinpath(Path(__file__).parent, "../settings/paths.json").resolve(), "r") as f:
            paths = json.load(f)
        self.QE = paths["QE"]
        self.mpi = paths["mpi"]
        self.python = paths["python"]
        self.server = paths["server"]
        if self.server == "local":
            self.sh = "./"
            self.header = ["!/bin/bash"]
        elif self.server == "OBCX":
            self.sh = "qsub"
            self.header = ["#!/bin/bash -f",
                           f"#$ -pe smp {self.np}", "#$ -cwd", "#$ -N job_name"]
        self.commands = []
        for tasks in job:
            for task in tasks:
                if task.command[0].split(".")[-1] == "x":
                    self.commands.append(
                        f"{self.mpi} - np {self.np} {self.QE}/{task.command[0]} -npool {self.npool} {task.commdand[1]}")
                elif task.command[0].split(".")[-1] == "py":
                    self.commands.append(f"{self.python} {task.command[1]}")
        self.jobsh = self.header + self.commands
        print(*self.jobsh, sep="\n", end="\n",
              file=open(f"{self.workdir}/job.sh"))

    def __call__(self):
        # cd self.submitdir
        subprocess.run(f"{self.sh} {self.workdir}/job.sh")
        print("job submitted")
