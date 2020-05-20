# Some useful constants
from numpy import pi

tpi = 2.0 * pi
fpi = 4.0 * pi

C_SI = 2.99792458E+8                # m sec^-1
H_PLANCK_SI = 6.62606896E-34        # J s
K_BOLTZMANN_SI = 1.3806504E-23      # J K^-1
HARTREE_SI = 4.35974394E-18         # J
BOHR_RADIUS_SI = 0.52917720859E-10  # m
RYDBERG_SI = HARTREE_SI/2.0         # J
AU_SEC = H_PLANCK_SI/tpi/HARTREE_SI
AU_PS = AU_SEC * 1.0E+12
AU_TERAHERTZ = AU_PS
AU_GPA = HARTREE_SI / BOHR_RADIUS_SI ** 3 / 1.0E+9
K_BOLTZMANN_RY = K_BOLTZMANN_SI / RYDBERG_SI
RY_TO_THZ = 1.0 / AU_TERAHERTZ / fpi
RY_TO_GHZ = RY_TO_THZ * 1000.0
RY_TO_CMM1 = 1.0E+10 * RY_TO_THZ / C_SI
RY_KBAR = 10.0 * AU_GPA / 2.0

kb1 = 1.0 / K_BOLTZMANN_RY / RY_TO_CMM1  # inverse Boltzmann constant in cm^{
    -1
}/K
ev_to_ry = 0.073498618

# -*- coding: utf-8 -*-
"""
Physical or mathematical constants.
Since every code has its own conversion units, this module defines what
QE understands as for an eV or other quantities.
Whenever possible, we try to use the constants defined in
:py:mod:aiida.common.constants: , but if some constants are slightly different
among different codes (e.g., different standard definition), we define
the constants in this file.
"""

# These have been put here from the one of QE, taken directly from
# those in aiida.common.constants
bohr_to_ang = 0.52917720859
ang_to_m = 1.e-10
bohr_si = bohr_to_ang * ang_to_m
ry_to_ev = 13.6056917253
ry_si = 4.35974394 / 2. * 10**(-18)
hartree_to_ev = ry_to_ev * 2.
timeau_to_sec = 2.418884326155573e-17
invcm_to_THz = 0.0299792458

# From the definition of Quantum ESPRESSO, conversion from atomic mass
# units to Rydberg units:
#  REAL(DP), PARAMETER : : AMU_SI           = 1.660538782E-27_DP  ! Kg
#  REAL(DP), PARAMETER : : ELECTRONMASS_SI  = 9.10938215E-31_DP   ! Kg
#  REAL(DP), PARAMETER : : AMU_AU           = AMU_SI / ELECTRONMASS_SI
#  REAL(DP), PARAMETER : : AMU_RY           = AMU_AU / 2.0_DP
amu_Ry = 911.4442421323

prefix
outdir
lsym
no_overlap

ngauss
degauss
pawproj = .false.
filpdos
filproj
lwrite_overlaps = .false.
lbinary_data = .false.
kresolveddos = .false.
tdosinboxes = .false.
n_proj_boxes = 1

CONSTRAINTS
OCCUPATIONS
ATOMIC_FORCES

tstress
tprnfor
etot_conv_thr
forc_conv_thr

tefield
~: extfield
gate


assume_isolated

It is the convergence threshold for the force to an ionic minimization.
 Any force for all elements must be less than this value(3.8E-4 Ry/Bohr=0.01 eV/Å).




def convert_constraints(atoms):
    '''
    Convert some of ase's constraints to pw.x constraints for pw.x internal
    relaxation returns constraints which are simply expressed as setting
    force components as first list and other contraints that are
    implemented in espresso as second list
    '''

    if atoms.constraints:
        n = len(atoms)
        if n == 0:
            return [], []
        forcefilter = []
        otherconstr = []
        for c in atoms.constraints:
            if isinstance(c, constraints.FixAtoms):
                if len(forcefilter) == 0:
                    forcefilter = np.ones((n, 3), np.int)
                forcefilter[c.index] = [0, 0, 0]
            elif isinstance(c, constraints.FixCartesian):
                if len(forcefilter) == 0:
                    forcefilter = np.ones((n, 3), np.int)
                forcefilter[c.a] = c.mask
            elif isinstance(c, constraints.FixBondLengths):
                for d in c.constraints:
                    otherconstr.append("'distance' %d %d" %
                                       (d.indices[0]+1, d.indices[1]+1))
            elif isinstance(c, constraints.FixBondLength):
                otherconstr.append("'distance' %d %d" %
                                   (c.indices[0]+1, c.indices[1]+1))
            elif isinstance(c, constraints.FixInternals):
            # we ignore the epsilon in FixInternals because there can only be one global
            # epsilon be defined in espresso for all constraints
                for d in c.constraints:
                    if isinstance(d, constraints.FixInternals.FixBondLengthAlt):
                        otherconstr.append("'distance' %d %d %s" % (d.indices[0]+1, d.indices[1]+1, num2str(d.bond)))
                    elif isinstance(d, constraints.FixInternals.FixAngle):
                        otherconstr.append("'planar_angle' %d %d %d %s" % (d.indices[0]+1, d.indices[1]+1, d.indices[2]+1,num2str(np.arccos(d.angle)*180./np.pi)))
                    elif isinstance(d, constraints.FixInternals.FixDihedral):
                        otherconstr.append("'torsional_angle' %d %d %d %d %s" % (d.indices[0]+1, d.indices[1]+1, d.indices[2]+1,d.indices[3]+1,num2str(np.arccos(d.angle)*180./np.pi)))
                    else:
                        raise NotImplementedError('constraint {} from FixInternals not implemented\n'
                            'consider ase-based relaxation with this constraint instead'.format(d.__name__))
            else:
                raise NotImplementedError('constraint {} not implemented\n'
                    'consider ase-based relaxation with this constraint instead'.format(c.__name__))
        return forcefilter, otherconstr
    else:
        return [], []
