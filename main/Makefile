################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq rmpmk
.PHONY: all


# main programs and required modules 

MAIN = cylinder # halfstrip infstrip

RANDOM = ranlxs ranlxd # gauss
# expo rootexpo gaussalt symexpo symrootexpo

INTEGRATION = integrate methods gaussmom pigreco

SIMULATIONS = harmosc metropolis

RUNGEKUTTA = ode2

CLUSTER = jackknife

START = start utils

EXTRAS = extras zeros

NUMEROV = numerov

MODULES = $(RANDOM) $(EXTRAS) $(START) 
# $(CLUSTER) $(INTEGRATION) $(SIMULATIONS) $(RUNGEKUTTA) $(NUMEROV)


# search path for modules

MDIR = ../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/integration:$(MDIR)/simulations:$(MDIR)/cluster:$(MDIR)/plot:$(MDIR)/rungekutta:$(MDIR)/numerov


# additional include directories

INCPATH = /u/shared/modules/libraries ../include  /usr/include /usr/include/gsl  # /usr/include/igraph  /usr/local/include/igraph


# additional libraries to be included 

LIBS = m  gsl  # igraph

LIBPATH = 


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -lm -lgsl -lgslcblas -O -g -DSSE4 -Wall # -std=gnu99  #-lcblas
 

############################## do not change ###################################

SHELL=/bin/bash

CC=mpicc #gcc

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory and executables only remain

rmpmk:
	@ -rm -rf *.d *.o .tmp


# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
