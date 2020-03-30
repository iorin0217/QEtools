# https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
def create_band_in(path, nspin=1, plot_2d=".false."):
    band_temp = ["&BANDS", f"plot_2d={plot_2d}"]
    if nspin == 1:
        band_temp += ["filband = 'band.out'", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band.in", "w"))
    elif nspin == 2:
        band_temp_cp = band_temp[:]
        band_temp += ["filband = 'band_up.out'", "spin_component = 1", "/"]
        band_temp_cp += ["filband = 'band_dn.out'", "spin_component = 2", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_up.in", "w"))
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_dn.in", "w"))
    elif nspin == 4:
        band_temp += ["filband = 'band_s.out'", "lsigma = .true.", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_s.in", "w"))
