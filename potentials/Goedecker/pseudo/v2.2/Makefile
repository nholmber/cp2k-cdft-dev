# Targets:

LIBRARY		= libout.a
PROGRAM		= pseudo.x
TARGET		= $(PROGRAM)

FC = gfortran
F90C = gfortran -ffree-form
FFLAGS = -O2 -g
LD = gfortran
LDFLAGS = -static
LIBS = -llapack -lblas

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
#       sub_lapack.f\
#       sub_blas.f

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
#       sub_lapack.o\
#       sub_blas.o

$(PROGRAM):	$(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
		@echo "done"

functionals.o : xc_b97.o

%.o: %.f90
	$(F90C) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
		rm -f *.o *.mod

distclean:
		rm -f *.o *.mod $(PROGRAM)
