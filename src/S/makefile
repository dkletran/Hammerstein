
# CCFLAG=-c -g
CCFLAG=-c -g  -O4 -fpermissive -Wno-deprecated
LIBFLAG=-ru

INCLUDEUTIL=../include/
INCLUDEIOS=/usr/include/c++/4.3/backward
OBJ=../../obj/

LIBUTIL=../../bin/

help:  
	@echo      "make simul		: simulations - main.cc"
	@echo      "make lib		: créer le bibliothèaque libclasse.a"
	@echo      "make test_estsp		: test  estimation polyspectre-test_estsp.cc"
	@echo      "make test_std		: test  steepest descet algorithme- test_std.cc"
	@echo      "make test_lma		: test  LMA algorithme test_lma.cc"

simul :$(LIBUTIL)/simul
	@echo "Done"
	$(LIBUTIL)/simul

lib :   $(LIBUTIL)/libclasse.a
	@echo "Done"


test_estsp: $(LIBUTIL)/test_estsp
	@echo "Done"
	$(LIBUTIL)/test_estsp

test_lma: $(LIBUTIL)/test_lma
	@echo "Done"
	$(LIBUTIL)/test_lma

test_std: $(LIBUTIL)/test_std
	@echo "Done"
	$(LIBUTIL)/test_std

OBJGENERAL=$(OBJ)/Volume.o \
		$(OBJ)/Image.o \
		$(OBJ)/Vecteur.o\
		$(OBJ)/Moment3_1D.o\
		$(OBJ)/Moment4_1D.o\
		$(OBJ)/steepest_descent.o\
		$(OBJ)/EstimationMoment.o\
		$(OBJ)/EstimationCumulant.o\
		$(OBJ)/Noise.o\
		$(OBJ)/FiltrageSOVM.o\
		$(OBJ)/GenereSignal.o\
		$(OBJ)/fft.o\
		$(OBJ)/EstimationSpectre.o\
		$(OBJ)/ConjugateGradient.o\
		$(OBJ)/LM_algo.o\

$(OBJ)/%.o: %.cc
	g++  $(CCFLAG) -o $@ $?  -I$(INCLUDEUTIL) -I$(INCLUDEIOS)



$(LIBUTIL)/libclasse.a : $(OBJGENERAL)
	ar $(LIBFLAG) $@  $?

$(LIBUTIL)/test_estsp: $(OBJ)/test_estsp.o $(LIBUTIL)/libclasse.a
	g++  -o $@  $< -L$(LIBUTIL) -l classe

$(LIBUTIL)/test_lma: $(OBJ)/test_lma.o $(LIBUTIL)/libclasse.a
	g++  -o $@  $< -L$(LIBUTIL) -l classe

$(LIBUTIL)/test_std: $(OBJ)/test_std.o $(LIBUTIL)/libclasse.a
	g++  -o $@  $< -L$(LIBUTIL) -l classe

$(LIBUTIL)/simul: $(OBJ)/main.o $(LIBUTIL)/libclasse.a
	g++  -o $@  $< -L$(LIBUTIL) -l classe
