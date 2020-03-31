mpirun - np 1 fermi_velocity.x - in nscf.in
fermisurfer vfermi.frmsf
mpirun - np 1 ~/bin/fermi_proj.x - in pdos.in
fermisurfer proj.frmsf
