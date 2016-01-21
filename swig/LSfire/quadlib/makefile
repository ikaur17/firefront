# makefile by Undy

SHELL = /bin/sh
CC = gfortran
FLAGS = #-Wall # -g


#$(OBJECTS): $(SOURCES) $(HEADERS)
#	$(CC) -c -o $@ $< $(CFLAGS)
	
compile: gen_hermite_rule.exe gen_laguerre_rule.exe

gen_hermite_rule.exe: gen_hermite_rule.o 
	@echo "...linking   :" $< 
	@$(CC) $(FLAGS) $< -o $@ 
	@echo "...executable: "$@

gen_laguerre_rule.exe: gen_laguerre_rule.o 
	@echo "...linking   :" $< 
	@$(CC) $(FLAGS) $< -o $@ 
	@echo "...executable: "$@

	
gen_hermite_rule.o : gen_hermite_rule.f90
	@echo "...compiling :" $<
	@$(CC) $(FLAGS) -c $< -o $@
	
gen_laguerre_rule.o : gen_laguerre_rule.f90
	@echo "...compiling :" $<
	@$(CC) $(FLAGS) -c $< -o $@
	
.PHONY : clean cleanobj cleanexe cleantxt

run: gen_hermite_rule.exe gen_laguerre_rule.exe
	@./gen_hermite_rule.exe  10 0.0 0.0 1.0 quad_hermite_10
	@./gen_hermite_rule.exe  20 0.0 0.0 1.0 quad_hermite_20
	@./gen_hermite_rule.exe  30 0.0 0.0 1.0 quad_hermite_30
	@./gen_hermite_rule.exe  40 0.0 0.0 1.0 quad_hermite_40
	@./gen_hermite_rule.exe  50 0.0 0.0 1.0 quad_hermite_50
	@./gen_hermite_rule.exe  60 0.0 0.0 1.0 quad_hermite_60
	@./gen_hermite_rule.exe  70 0.0 0.0 1.0 quad_hermite_70
	@./gen_hermite_rule.exe  80 0.0 0.0 1.0 quad_hermite_80
	@./gen_hermite_rule.exe  90 0.0 0.0 1.0 quad_hermite_90
	@./gen_hermite_rule.exe 100 0.0 0.0 1.0 quad_hermite_100
	@./gen_hermite_rule.exe 120 0.0 0.0 1.0 quad_hermite_120
	@./gen_hermite_rule.exe 140 0.0 0.0 1.0 quad_hermite_140
	@./gen_hermite_rule.exe 160 0.0 0.0 1.0 quad_hermite_160
	@./gen_hermite_rule.exe 180 0.0 0.0 1.0 quad_hermite_180
	@./gen_hermite_rule.exe 200 0.0 0.0 1.0 quad_hermite_200
	@./gen_laguerre_rule.exe  10 0.0 0.0 1.0 quad_laguerre_10
	@./gen_laguerre_rule.exe  20 0.0 0.0 1.0 quad_laguerre_20
	@./gen_laguerre_rule.exe  30 0.0 0.0 1.0 quad_laguerre_30
	@./gen_laguerre_rule.exe  40 0.0 0.0 1.0 quad_laguerre_40
	@./gen_laguerre_rule.exe  50 0.0 0.0 1.0 quad_laguerre_50
	@./gen_laguerre_rule.exe  60 0.0 0.0 1.0 quad_laguerre_60
	@./gen_laguerre_rule.exe  70 0.0 0.0 1.0 quad_laguerre_70
	@./gen_laguerre_rule.exe  80 0.0 0.0 1.0 quad_laguerre_80
	@./gen_laguerre_rule.exe  90 0.0 0.0 1.0 quad_laguerre_90
	@./gen_laguerre_rule.exe 100 0.0 0.0 1.0 quad_laguerre_100
	
clean :
	rm -f *.o *.exe *.txt
