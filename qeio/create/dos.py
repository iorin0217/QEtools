# https://www.quantum-espresso.org/Doc/INPUT_DOS.html
def create_dos_in(path, efermi, emin = -10, emax = 10, deltae = 0.05):
    dos_temp = ["&DOS", f"emin = {efermi+emin}", f"{efrmi+emax}", f"deltae = {deltae}", "/"]
    print(*dos_temp, sep="\n", end="\n", file=open(f"{path}/dos.in","w"))