##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for toolbox functions
##

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

LIB_DIRS     = -L$(LSMLIB_LIB_DIR)

PROGS        = calculus_toolbox_demo

SUBDIRS = fast_marching_method      \
          boundary_conditions

all:  $(PROGS)
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

calculus_toolbox_demo:  calculus_toolbox_demo.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

clean:
	rm -f *.o 
	rm -f $(PROGS)
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done
		
spotless:  clean
	rm -f *.log
