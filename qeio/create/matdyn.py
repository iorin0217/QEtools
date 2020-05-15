# https://gitlab.com/QEF/q-e/-/blob/develop/PHonon/PH/matdyn.f90

&input
asr = 'simple',
dos = .true.,
amass(1) = 28.0855,
flfrc = 'si.fc',
fldos = 'si.phdos',
nk1 = 50, nk2 = 50, nk3 = 50,
/
# lambda
   &input
        asr = 'simple'
        flfrc = 'al.fc', flfrq = 'al.freq', la2F = .true., dos = .true.
        fldos = 'phonon.dos', nk1 = 16, nk2 = 16, nk3 = 16
    /
