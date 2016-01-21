# makefile by Undy
flag_omp     = yes
flag_omp_for = yes
gprof = no
o = 3

SHELL = /bin/sh

DIR_SRC = src
DIR_OBJ = build
DIR_LSMLIB_H = ./lsmlib-1.0.1/include
DIR_LSMLIB_A = ./lsmlib-1.0.1/lib
DIR_FIRELIB  = ./firelib-1.0.4
CC = gcc

FLAG_O = -O$(o)

ifeq ($(gprof),no)
FLAG_GPROF =
else
FLAG_O = -O0
FLAG_GPROF = -pg
endif

ifeq ($(flag_omp),yes)
FLAG_OMP     = -fopenmp
FLAG_OMP_FOR = -cpp -fopenmp
else
FLAG_OMP      = -Wno-pragmas
FLAG_OMP_FOR = -cpp 
endif


FFLAGS     = $(FLAG_O) $(FLAG_GPROF) -Wno-unused-result -Wno-format $(FLAG_OMP) -I$(DIR_LSMLIB_H) -I$(DIR_FIRELIB)
FFLAGS_FOR = $(FLAG_O) $(FLAG_GPROF) $(FLAG_OMP_FOR) -fdefault-real-8 -fdefault-double-8 #-g
#LFLAGS     = $(FFLAGS) -L$(DIR_LSMLIB_A) -llsm_serial -llsm_toolbox -lm
LFLAGS     = $(FLAG_O) $(FLAG_GPROF) $(FLAG_OMP) -L$(DIR_LSMLIB_A) -llsm_serial -llsm_toolbox -lm

#SOURCES    = $(DIR_FIRELIB)/fireLib.c $(DIR_SRC)/LSfire+.c $(DIR_SRC)/parameters.h $(DIR_SRC)/input_fire.c $(DIR_SRC)/input_adiff.c $(DIR_SRC)/libsanity.c
HEADERS    = $(DIR_SRC)/LSfire+.h
OBJECTS    = $(DIR_OBJ)/fireLib.o $(DIR_OBJ)/LSfire+.o $(DIR_OBJ)/input_fire.o $(DIR_OBJ)/input_adiff.o $(DIR_OBJ)/libsanity.o $(DIR_OBJ)/lsm2lib.o
LINKER     = @$(CC) $(OBJECTS) $(LFLAGS) -o $@
EXECUTABLE = LSfire+

#$(OBJECTS): $(SOURCES) $(HEADERS)
#	$(CC) -c -o $@ $< $(CFLAGS)
	
$(EXECUTABLE) : mkdir-build $(OBJECTS) 
	@echo "...linking   :" $(OBJECTS) $(LFLAGS)
	@$(CC) $(OBJECTS) $(LFLAGS) -o $@ 
	@echo "...executable:" $(EXECUTABLE)
	
mkdir-build:
	@if [ ! -d ./build ] ; then mkdir build ; fi

$(DIR_OBJ)/LSfire+.o : $(DIR_SRC)/LSfire+.c $(DIR_SRC)/LSfire+.h $(DIR_SRC)/parameters.h $(DIR_SRC)/libsanity.c $(DIR_SRC)/input_fire.c $(DIR_SRC)/input_adiff.c $(DIR_SRC)/lsm2lib.f90 $(DIR_FIRELIB)/fireLib.c
	@echo "...compiling :" $<
	@$(CC) $(FFLAGS) -c -std=c99 $< -o $@

$(DIR_OBJ)/lsm2lib.o : $(DIR_SRC)/lsm2lib.f90
	@echo "...compiling :" $<
	@$(CC) $(FFLAGS_FOR) -Wno-c-binding-type -c $< -o $@
	@mv lsm2.mod $(DIR_OBJ)/lsm2.mod
		
$(DIR_OBJ)/input_fire.o : $(DIR_SRC)/input_fire.c $(DIR_SRC)/LSfire+.h $(DIR_SRC)/parameters.h
	@echo "...compiling :" $<
	@$(CC) $(FFLAGS) -c $< -o $@
	
$(DIR_OBJ)/input_adiff.o : $(DIR_SRC)/input_adiff.c $(DIR_SRC)/LSfire+.h $(DIR_SRC)/parameters.h
	@echo "...compiling :" $<
	@$(CC) $(FFLAGS) -c $< -o $@
	
$(DIR_OBJ)/libsanity.o : $(DIR_SRC)/libsanity.c $(DIR_SRC)/LSfire+.h $(DIR_SRC)/parameters.h
	@echo "...compiling" : $<
	@$(CC) $(FFLAGS) -c $< -o $@
		
$(DIR_OBJ)/fireLib.o : $(DIR_FIRELIB)/fireLib.c
	@echo ."..compiling" : $<
	@$(CC) $(FFLAGS) -std=c99 -Wno-format -c $< -o $@

run: $(EXECUTABLE)
	@echo "\n...running executable "$(EXECUTABLE)"...\n"
	@./$(EXECUTABLE)
	
.PHONY : clean cleanobj cleanexe cleantxt cleandat cleanLSout

clean : cleanobj cleanexe

cleanall : cleanobj cleanexe cleanLSout
	
cleanLSout :
	@./clean LSout
	
cleantxt :
	@./clean txt

cleandat :
	@./clean dat
	
cleanobj :
	@./clean obj
#@i=0; for FILE in $(OBJECTS) ; do if [ -f $$FILE ] ; then rm $$FILE ; i=$$((i+1)) ; REMOVED=$$REMOVED" "$$FILE; fi ; done ; if [ "$$i" -gt "0" ] ; then echo "...removed"$$REMOVED ; fi
	
cleanexe :
	@if [ -f $(EXECUTABLE) ] ; then rm $(EXECUTABLE) ; echo "...deleting exe: "$(EXECUTABLE) ; else echo "...deleting exe: (None)" ; fi
#@./clean exe

spotless :
	if [ -d ./build ] ; then rm -f ./build/* ; rmdir ./build ; rm $(EXECUTABLE) ; fi
