# https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
def create_projwfc_in(path, deltae=0.01):
    projwfc_temp = ["&PROJWFC", f"deltae = {deltae}", "/"]
    print(*projwfc_temp, sep="\n", end="\n",
          file=open(f"{path}/projwfc.in", "w"))
