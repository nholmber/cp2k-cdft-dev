SHELL = /bin/bash
# PGI 5.0
#FC90     = pgf90 -Mfree -c
#FC77     = pgf90 -Mfixed -c
#FFLAGS   = -fastsse -r8 -tp p7
#LFLAGS   = -lacml
#LD       = pgf90 $(FFLAGS)
# Intel Fortran compiler 9.0 and Cluster MKL 8.0
#FC90     = ifort -FR -c
#FC77     = ifort -FI -c
#LD       = ifort $(FFLAGS) -i-static
#FFLAGS   = -O3 -r8
#MKLPATH  = /apps64/intel/cmkl/8.1/lib/32
#LFLAGS   = $(FFLAGS) -i-static
#LIBS     = $(MKLPATH)/libmkl_lapack.a\
#           $(MKLPATH)/libmkl_ia32.a\
#           $(MKLPATH)/libguide.a\
#           -lpthread
#G95
FC90     = g95 -ffree-form
FC77     = g95 -ffixed-form 
LD       = g95
FFLAGS   = -c
# -g -ftrace=full -fbounds-check
LFLAGS   =
LIBS     = -L/users/famo/lib/Linux_HAMMER32SSE2_2 -llapack_LINUX -lf77blas -lcblas -latlas
# NAG compiler
#FC90     = f95 -c
#FC77     = f95 -c
#LD       = f95
#FFLAGS   = -O0 -g -C=all -mismatch_all -f77 -nan
#FFLAGS   = -O3 -mismatch_all -f77
#LFLAGS   = $(FFLAGS) -L$(HOME)/lib/32
#LIBS     = -llapack_nag -lblas_nag
#
OBJECTS  = basic_data_types.o mathlib.o \
mod_atom.o mod_energies.o mod_pspot.o mod_rint.o mod_upd.o mod_xcfcn.o \
atomcode.o \
c2a.o calcnn.o chargeden.o coulomb.o calc_energy.o \
denmat.o diag.o diag_ks.o diag_oneeh.o \
evxc.o exp_grad.o exp_opt.o  \
expdiis.o \
fac.o fio.o fitin.o functionals.o \
gam.o gdiis90.o gint.o \
hcth.o hermite_mat.o hermite_matd.o input.o \
kbpp.o kbppfit.o kinetic_int.o kinetic_intd.o ks_matrix.o \
lbfgs90.o lda_scf.o linesearch.o lsfbtr.o \
nuc_att.o nuc_attd.o numd_wfn.o \
overlap.o overlapd.o \
ppinit.o pseudopot.o pseudopotd.o  \
rint.o \
test_grad.o test_grad_direct.o \
show.o sint.o trace.o \
update.o updalpha.o \
v_add.o v_add_d.o wfn_ortho.o
#
#
all : 
	make OPTBAS
#
OPTBAS:	$(OBJECTS)
	$(LD) -o optbas $(OBJECTS) $(LFLAGS) $(LIBS)
#
%.o : %.f90
	$(FC90) $(FFLAGS) $*.f90
#
%.o : %.f
	$(FC77) $(FFLAGS) $*.f
#
clean :
	rm $(OBJECTS) *.mod
#
#
#
# DEPENDENCES BELOW THIS LINE WILL BE DELETED BY NEXT RUN OF MKF90DEP
atomcode.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 mod_xcfcn.f90 
c2a.o: mod_pspot.f90 
calc_energy.o: mod_atom.f90 mod_energies.f90 mod_pspot.f90 mod_rint.f90 
calcnn.o: mod_atom.f90 
chargeden.o: mod_atom.f90 
coulomb.o: mod_atom.f90 
denmat.o: mod_atom.f90 
diag.o: mathlib.f90
diag_ks.o: mod_atom.f90 
diag_oneeh.o: mod_atom.f90 
evxc.o: mod_atom.f90 
exp_grad.o: mod_atom.f90 mod_pspot.f90 
exp_opt.o: mod_atom.f90 mod_upd.f90 
expdiis.o: mod_atom.f90 mod_upd.f90 
functionals.o: mod_xcfcn.f90 
hermite_mat.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 
hermite_matd.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 
input.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 mod_xcfcn.f90 
kbpp.o: mod_atom.f90 mod_pspot.f90 
kbppfit.o: mod_atom.f90 mod_pspot.f90 
kinetic_int.o: mod_atom.f90 
kinetic_intd.o: mod_atom.f90 
ks_matrix.o: mod_atom.f90 mod_rint.f90 
lda_scf.o: mod_atom.f90 mod_energies.f90 mod_pspot.f90 mod_rint.f90 
mod_pspot.o: mod_atom.f90 
mod_upd.o: mod_atom.f90 
nuc_att.o: mod_atom.f90 
nuc_attd.o: mod_atom.f90 
numd_wfn.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 
overlap.o: mod_atom.f90 
overlapd.o: mod_atom.f90 
ppinit.o: mod_pspot.f90 
pseudopot.o: mod_atom.f90 mod_pspot.f90 
pseudopotd.o: mod_atom.f90 mod_pspot.f90 
show.o: mod_atom.f90 mod_pspot.f90 mod_rint.f90 mod_xcfcn.f90 
test_grad.o: mod_atom.f90 mod_energies.f90 mod_upd.f90 
test_grad_direct.o: mod_atom.f90 mod_energies.f90 mod_upd.f90 
trace.o: mod_atom.f90 
updalpha.o: mod_atom.f90 mod_upd.f90 
update.o: mod_atom.f90 
v_add.o: mod_atom.f90 
wfn_ortho.o: mod_atom.f90 
