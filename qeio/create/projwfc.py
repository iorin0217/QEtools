# https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
def create_projwfc_in(path, variables):
    projwfc_temp = ["&PROJWFC"]
    pdos_in = projwfc_temp+[f"deltae = {variables['deltae']}", "/"]
    fatband_in = projwfc_temp+["lsym=.false.", "/"]
    print(*pdos_in, sep="\n", end="\n", file=open(f"{path}/pdos.in", "w"))
    print(*fatband_in, sep="\n", end="\n",
          file=open(f"{path}/fatbands.in", "w"))
