# https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
def create_band_in(path, nspin=1, plot_2d=".false."):
    band_temp = ["&BANDS", f"plot_2d={plot_2d}"]
    if nspin == 1:
        band_temp += ["filband = 'band.out'", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band.in", "w"))
    elif nspin == 2:
        band_temp_dn = band_temp[:]
        band_temp += ["filband = 'band_up.out'", "spin_component = 1", "/"]
        band_temp_dn += ["filband = 'band_dn.out'", "spin_component = 2", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_up.in", "w"))
        print(*band_temp_dn, sep="\n", end="\n",
              file=open(f"{path}/band_dn.in", "w"))
    elif nspin == 4:
        band_temp_y = band_temp[:]
        band_temp_z = band_temp[:]
        band_temp += ["filband = 'band_s.out'", "lsigma(1) = .true.", "/"]
        band_temp_y += ["filband = 'band_s.out'", "lsigma(2) = .true.", "/"]
        band_temp_z += ["filband = 'band_s.out'", "lsigma(3) = .true.", "/"]
        print(*band_temp, sep="\n", end="\n",
              file=open(f"{path}/band_sx.in", "w"))
        print(*band_temp_y, sep="\n", end="\n",
              file=open(f"{path}/band_sy.in", "w"))
        print(*band_temp_z, sep="\n", end="\n",
              file=open(f"{path}/band_sz.in", "w"))
