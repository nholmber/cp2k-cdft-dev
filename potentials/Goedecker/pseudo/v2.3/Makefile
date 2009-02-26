# Targets:

LIBRARY		= libout.a
PROGRAM		= pseudo
TARGET		= $(PROGRAM)

F90C		= ifort -free
FC		= ifort -fixed
FFLAGS		= -g 
# -warn all 
LD		= ifort
LFLAGS		= -Vaxlib -L/opt/intel/mkl/10.0.2.018/lib/em64t -Wl,-rpath,/opt/intel/mkl/10.0.2.018/lib/em64t
#LFLAGS		= -Vaxlib -L/opt/intel/mkl/10.0.2.018/lib/32 -Wl,-rpath,/opt/intel/mkl/10.0.2.018/lib/32
LIBS		= -lmkl_sequential -lmkl_lapack -lmkl_em64t

#FC = g95 -ffixed-form
#FFLAGS =
#LD = g95
#LDFLAGS = $(FFLAGS) 
#LIBS = -L/users/famo/lib/Linux_HAMMER32SSE2_2 -llapack_LINUX -lptf77blas -lptcblas -latlas -lpthread
##   Compilers, Programs and Flags

#   for ibm/aix:
#   Fortran Compiler
#FC		= xlf
#FFLAGS		= -O2 -qarch=com -qmaxmem=-1
## POWER2
##FFLAGS		= -O2 -qarch=pwr2 -w
##   Linker/Loader
#LD		= xlf
#LFLAGS		= $(FFLAGS) -L/home/pool30/qs/lib
#LIBS		= -lmass -llapack-opt -lessl
## to link libaries static: -bnso -bI:/usr/lib/syscalls.exp

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
        hcth.o\
        xc_b97.o
#        sub_lapack.o\
#        sub_blas.o

#-----------------------------------------

$(PROGRAM):	$(SCRS) $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
		@echo "done"

functionals.o : xc_b97.o

%.o: %.f90
	$(F90C) $(FFLAGS) -c $<


#%.o:%.f
#	$(FC) -c -$(FFLAGS) %<
clean:
		rm -f *.o
