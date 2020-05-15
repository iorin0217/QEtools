# https://gitlab.com/QEF/q-e/blob/develop/PHonon/PH/q2r.f90

&input
fildyn = 'si.dyn',
zasr = 'simple',
flfrc = 'si.fc',
/
# lambda
   &input
        zasr = 'simple',  fildyn = 'al.dyn', flfrc = 'al.fc', la2F = .true.,
    /
# tetra?
&INPUT
 fildyn = 'matdyn'
   la2f = .true.
 lshift_q = .true.
   zasr = 'crystal'
  flfrc = 'ifc.dat'
/
