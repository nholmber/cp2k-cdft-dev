SHELL = /bin/sh
GFITHOME  := $(PWD)/..
PROG =	gfit.exe
MAKEOBJ      := $(shell $(GFITHOME)/bin/makeobj)
MAKEDEP      := $(shell $(GFITHOME)/bin/makedep)
####################################################################
# Objects List
####################################################################
include objlist
####################################################################

CC = gcc
CFLAGS = -O2
##
## xlf compiler family
##
#F90 = xlf95_r
#F90FLAGS = -O3  -qsuffix=f=f90 -qfree=f90 -qstrict 
##
## nag compiler family
##
F90 = f95
F90FLAGS = -C=all -mtrace=size -nan -colour -free 


all: 	$(PROG)
	@echo "GFIT compilation done!"

install: 
	 @$(MAKE) all
	 @mv gfit.exe $(GFITHOME)/bin/
	 @ls -l $(GFITHOME)/bin/$(PROG)
	 @echo "Executable has been moved in "$(GFITHOME)/bin/" directory!"
	 @echo "Happy computing!"

$(PROG): $(OBJS)
	 $(F90) -o $@ $(OBJS) $(LIBS)

clean:
	@rm  -f $(PROG) $(OBJS) *.mod *.d
	@rm  -f objlist dependencies
	@rm  -f $(GFITHOME)/bin/$(PROG) 

.SUFFIXES:
.SUFFIXES: .o .f90 

.f90.o:
	$(F90) $(F90FLAGS) -c $<

####################################################################
# Module Dependencies
####################################################################
include dependencies
####################################################################

