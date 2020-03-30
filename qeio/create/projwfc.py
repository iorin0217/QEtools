# https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
def create_projwfc_in(path, efermi, emin=-5, emax=5, deltae=0.01):
    projwfc_temp = ["&PROJWFC", f"emin = {efermi+emin}",
                    f"emax = {efrmi+emax}", f"deltae = {deltae}", "/"]
    print(*projwfc_temp, sep="\n", end="\n",
          file=open(f"{path}/projwfc.in", "w"))
