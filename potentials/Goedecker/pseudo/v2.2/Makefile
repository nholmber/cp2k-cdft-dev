# Targets:

LIBRARY		= libout.a
PROGRAM		= pseudo.x
TARGET		= $(PROGRAM)

#   Compilers, Programs and Flags
#   for Intel Fortran Compiler (Pentium4)
#F90C		= ifc -FR
#FC		= ifc -FI
#FFLAGS		= -O3 -tpp7
#   Linker/Loader
#LD		= ifc
#LFLAGS		= $(FFLAGS) -Vaxlib -L/usr/local/intel/mkl60/lib/32/
#LIBS		= -lmkl_lapack -lmkl_ia32 -lguide

#
#   for PGI Fortran Compiler
F90C		= pgf90 -Mfree
FC		= pgf90 -Mfixed
FFLAGS		= -fast
#   Linker/Loader
LD		= pgf90
LFLAGS		= $(FFLAGS)
LIBS		= -llapack -lblas
#
#   for ibm/aix:
#   Fortran Compiler
#F90C		= xlf95
#FC		= xlf
#FFLAGS		= -O -q32 -qarch=auto -qmaxmem=-1 -qtune=auto
#   Linker/Loader
#LD		= xlf95
#LFLAGS		= $(FFLAGS) -L/users/krack/lib/32 -L/apps/mass/lib
#LIBS		= -lmass -llapack3 -lessl
# to link libaries static: -bnso -bI:/usr/lib/syscalls.exp

#   for SGI/IRIX 6.3 (O2):
#   Fortran Compiler
#FC		= f77
#FFLAGS		= -O2 -mips4 -n32 -r10000
#   Linker/Loader
#LD		= f77
#LFLAGS		= $(FFLAGS)
#LIBS		= -lcomplib.sgimath -lfastm

#   for DEC
#   Fortran Compiler
#FC		= f90 
#   Linker/Loader
#LD              = f90 -ldxml

#   for cray-td3
#FC              = f90 
#FFLAGS          = -dp -O2  
# Linker/Loader
#LD              = f90
#LDFLAGS         = -l /opt/ctl/craylibs/craylibs/libblas.a  


# Shell Command
#SHELL		= /bin/ksh

#   Sources Headers and Objects 
#   if the librarys lapack & blas are not available
#   add the files sub_lapack.f and sub_blas.f


SCRS	= pseudo.f \
	crtvh.f \
	penalty.f \
	amoeba.f \
	gatom.f \
	wave.f \
	wave2.f \
	wave3.f \
	detnp.f \
	resid.f \
	etot.f \
	pj2test.f \
	xpown.f \
	gamma.f \
	ppack.f \
	functionals.f\
	radgrid.f\
	zero.f\
	xcfunction.f\
	zbrent.f\
	ggaenergy_15.f
#        sub_lapack.f\
#        sub_blas.f

OBJS	= pseudo.o \
	crtvh.o \
	penalty.o \
	amoeba.o \
	gatom.o \
	wave.o \
	wave2.o \
	wave3.o \
	detnp.o \
	resid.o \
	etot.o \
	pj2test.o \
	xpown.o \
	gamma.o \
	ppack.o \
	functionals.o\
	radgrid.o\
	zero.o\
	xcfunction.o\
	zbrent.o\
	ggaenergy_15.o\
	hcth.o
#        sub_lapack.o\
#        sub_blas.o

#-----------------------------------------

$(PROGRAM):	$(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
		@echo "done"

hcth.o: hcth.f
	$(F90C) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
		rm -f *.o

realclean:
		rm -f *.o $(PROGRAM)

