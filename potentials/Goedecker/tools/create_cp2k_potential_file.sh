#!/bin/sh
#
  potential_file=$(pwd)/../cp2k/GTH_POTENTIALS
  cat <<*** >$potential_file
################################################################################
#
# Potential data base file for CP2K (Quickstep)
#
################################################################################
#
# Pseudopotentials of Goedecker, Teter and Hutter (GTH)
# -----------------------------------------------------
#
# History:    - Creation (12.12.1999, Matthias Krack)
#             - Electronic configurations added (11.05.2000,MK)
#             - GTH-PP for first-row transition metal added (18.03.2003,MK)
#             - Automatic update (16.12.2003,MK)
#             - PBE GTH-PPs for the Lanthanides added (30.11.2012,MK)
#             - Last update ($(date +%d.%m.%y),MK)
#
# Literature: - S. Goedecker, M. Teter, and J. Hutter,
#               Phys. Rev. B 54, 1703 (1996)
#             - C. Hartwigsen, S. Goedecker, and J. Hutter,
#               Phys. Rev. B 58, 3641 (1998)
#             - M. Krack,
#               Theor. Chem. Acc. 114, 145 (2005)
#
# GTH-potential format:
#
# Element symbol  Name of the potential  Alias names
# n_elec(s)  n_elec(p)  n_elec(d)  ...
# r_loc   nexp_ppl        cexp_ppl(1) ... cexp_ppl(nexp_ppl)
# nprj
# r(1)    nprj_ppnl(1)    ((hprj_ppnl(1,i,j),j=i,nprj_ppnl(1)),i=1,nprj_ppnl(1))
# r(2)    nprj_ppnl(2)    ((hprj_ppnl(2,i,j),j=i,nprj_ppnl(2)),i=1,nprj_ppnl(2))
#  .       .               .
#  .       .               .
#  .       .               .
# r(nprj) nprj_ppnl(nprj) ((hprj_ppnl(nprj,i,j),j=i,nprj_ppnl(nprj)),
#                                               i=1,nprj_ppnl(nprj))
#
# n_elec   : Number of electrons for each angular momentum quantum number
#            (electronic configuration -> s p d ...)
# r_loc    : Radius for the local part defined by the Gaussian function
#            exponent alpha_erf
# nexp_ppl : Number of the local pseudopotential functions
# cexp_ppl : Coefficients of the local pseudopotential functions
# nprj     : Number of the non-local projectors => nprj = SIZE(nprj_ppnl(:))
# r        : Radius of the non-local part for angular momentum quantum number l
#            defined by the Gaussian function exponents alpha_prj_ppnl
# nprj_ppnl: Number of the non-local projectors for the angular momentum
#            quantum number l
# hprj_ppnl: Coefficients of the non-local projector functions
#
***
  for xcfun in blyp bp hcth120 hcth407 pade pbe pbesol olyp; do
    cd ../cp2k/$xcfun
    XCFUN=$(echo $xcfun | tr [:lower:] [:upper:])
    cat <<*** >>$potential_file
################################################################################
#
# $XCFUN functional
#
################################################################################
#
***
    for e in H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar\
             K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se\
             Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn\
             Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy\
             Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb\
             Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf\
             Es Fm Md No Lr; do
      for f in $(ls ${e}-q* 2>/dev/null); do
        if [[ -s $f ]]; then
          cat $f >>$potential_file
          echo "#" >>$potential_file
        fi
      done
    done
    cd - >/dev/null
  done
  echo New CP2K/Quickstep PP database file created
