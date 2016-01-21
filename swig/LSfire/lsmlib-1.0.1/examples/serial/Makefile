##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for serial example codes
##

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

LIB_DIRS     = -L$(LSMLIB_LIB_DIR)

PROGS        = grid_management_demo

SUBDIRS = curvature_example         \
          fast_marching_method      \
          reinitialization_example  \

all:  $(PROGS)
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

grid_management_demo: grid_management_demo.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $< $(LSMLIB_LIBS) -o $@

clean:
	rm -f *.o 
	rm -f *.dat
	rm -f $(PROGS)
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

		
spotless:  clean
	rm -f *.log
	rm -f *.ascii
	rm -f *.binary
