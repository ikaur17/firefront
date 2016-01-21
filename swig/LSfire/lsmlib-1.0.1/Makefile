##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: top-level makefile for level set method library
##              (builds archive files)
##

SRC_DIR = .
VPATH = .
BUILD_DIR = .
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


all:  library 

directories:
	if [ ! -d "include" ]; then mkdir include; fi
	if [ ! -d "lib" ]; then mkdir lib; fi

includes:  directories
	cp -f config/LSMLIB_config.h include
	cd src; make $@ || exit 1


library:  directories includes 
	cd src; make $@ || exit 1
	make lsm_toolbox || exit 1
	make lsm_serial || exit 1
	if [ ! -z "" ]; then                             \
	  make lsm_parallel || exit 1;                                  \
	fi
	if [ ! -z "" ]; then                                  \
	  make matlab || exit 1;                                        \
	fi


lsm_toolbox:
	find src/toolbox -name "*.o" > objs_file.tmp
	# remove test code
	sed -e '/test/d' objs_file.tmp > objs_file_no_test.tmp 
	mv -f -f objs_file_no_test.tmp objs_file.tmp
	cat objs_file.tmp | xargs ar -ru $(BUILD_DIR)/lib/liblsm_toolbox.a 
	ranlib lib/liblsm_toolbox.a
	rm -f objs_file.tmp 


lsm_serial:
	find src/serial -name "*.o" > objs_file.tmp
	# remove test code
	sed -e '/test/d' objs_file.tmp > objs_file_no_test.tmp 
	mv -f -f objs_file_no_test.tmp objs_file.tmp
	cat objs_file.tmp | xargs ar -ru $(BUILD_DIR)/lib/liblsm_serial.a 
	ranlib lib/liblsm_serial.a
	rm -f objs_file.tmp 


lsm_parallel:
	find src/parallel -name "*.o" > objs_file.tmp
	# remove test code
	sed -e '/test/d' objs_file.tmp > objs_file_no_test.tmp 
	mv -f -f objs_file_no_test.tmp objs_file.tmp
	cat objs_file.tmp | xargs ar -ru $(BUILD_DIR)/lib/liblsm_parallel.a 
	ranlib lib/liblsm_parallel.a
	rm -f objs_file.tmp 


INSTALL_DIR = /usr/local

install:
	install -d -m 755 $(INSTALL_DIR)/lib
	install -d -m 755 $(INSTALL_DIR)/include
	install -d -m 755 $(INSTALL_DIR)/share/LSMLIB
	install -c -m 644 $(SRC_DIR)/LICENSE $(INSTALL_DIR)/share/LSMLIB
	install -d -m 755 $(INSTALL_DIR)/share/LSMLIB/config
	install -c -m 644 config.status $(INSTALL_DIR)/share/LSMLIB/config
	install -c -m 644 config/Makefile.config $(INSTALL_DIR)/share/LSMLIB/config
	sed 's/BUILD_DIR/LSMLIB_DIR/' $(INSTALL_DIR)/share/LSMLIB/config/Makefile.config > Makefile.config.install;                                   \
	mv -f Makefile.config.install $(INSTALL_DIR)/share/LSMLIB/config/Makefile.config;                                                                \
	for i in include/*; do                                          \
	  install -c -m 644 $$i $(INSTALL_DIR)/include;               \
	done
	if test -f lib/liblsm_toolbox.a; then                           \
	  for i in lib/*.a; do                                          \
	    install -c -m 644 $$i $(INSTALL_DIR)/lib;                 \
	    ranlib $(INSTALL_DIR)/$$i;                                \
	  done                                                          \
	fi
	if test -f lib/liblsm_toolbox.so; then                          \
	  for i in lib/*.so; do                                         \
	    install -c -m 644 $$i $(INSTALL_DIR)/lib;                 \
	    ranlib $(INSTALL_DIR)/$$i;                                \
	  done                                                          \
	fi
	if test -d $(SRC_DIR)/doc/lsmlib-dox; then                            \
  	cp -f -r $(SRC_DIR)/doc/lsmlib-dox $(INSTALL_DIR)/share/LSMLIB/doc;  \
  fi
	cp -f -r $(BUILD_DIR)/examples $(INSTALL_DIR)/share/LSMLIB/
	find $(INSTALL_DIR) -name Makefile > examples_makefiles.tmp
	for i in `cat examples_makefiles.tmp`; do                       \
	  awk '$$1=="BUILD_DIR" && $$2=="=" {$$3 = "$(INSTALL_DIR)"} {print}' $$i > $$i.install;           \
	  mv -f $$i.install $$i;                                         \
	  awk '$$1=="include" && $$2=="$$(BUILD_DIR)/config/Makefile.config" {$$2 = "$(INSTALL_DIR)/share/LSMLIB/config/Makefile.config"} {print}' $$i > $$i.install;           \
	  mv -f $$i.install $$i;                                         \
	  sed 's/BUILD_DIR/LSMLIB_DIR/' $$i > $$i.install;            \
	  mv -f $$i.install $$i;                                         \
	done
	rm -f examples_makefiles.tmp


dox:
	cd $(SRC_DIR)/doc/doxygen; doxygen Doxyfile

cleandox:
	if [ -d "doc/lsmlib-dox" ]; then rm -f -r doc/lsmlib-dox; fi


mex:  matlab

matlab:
	if [ ! -z "" ]; then                                  \
	  cd src/matlab; make library || exit 1;                        \
	else                                                              \
	  echo "";                                                        \
	  echo "LSMLIB configured without MATLAB...unable to generate";   \
	  echo "MEX files";                                               \
	  echo "";                                                        \
	fi

clean:
	cd src; make $@ || exit 1
	rm -f *.o
	rm -f *.tmp

spotless: clean cleandox
	cd src; make $@ || exit 1
	cd examples; make $@ || exit 1
	if [ -d "include" ]; then rm -f -rf include/*; rmdir include; fi
	if [ -d "lib" ]; then rm -f -rf lib/*; rmdir lib; fi

