#FC = gfortran
FC = ifort
FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
FFLAGS_PAR= -parallel -par_threshold75 -par_report3
FFLAGS_HYD= -O3 -ftz -align -xHost
FFLAGS_HYD_PAR= -O3 -ftz -align -xHost -parallel -par_threshold75 -par_report3
FFLAGS_DEB= -debug -C -traceback
FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
FFLAGS_MIO= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xB -ipo
FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
FFLAGS_LAP= -O3 -ftz -align -xAVX -no-prec-div 
FFLAGS_GF= -O3
LDLAGS= 
OMP_IF= -openmp -DOMP
#OMP_GF= -fopenmp -DOMP 
#LDLAGS= --static --disable-shared --enable-static
#FFLAGS= $(FFLAGS_DEB)
#FFLAGS= $(FFLAGS_GF_DEB)
FFLAGS= -O3 -ftz -align all ${OMP_IF} #$(FFLAGS_HYD_PAR) -fpp #$(FFLAGS_OMP_IF) -pg #-check bounds 
CPP= -fpp
MPI= -DMPI
#LIBS = -L/usr/lib -llapack  -lblas 
#LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIBS =  $(MKL_LIBS_MT) ${MKL_SCA_LIBS}
#LIBS = -L/usr/lib -llapack  -lblas -lfftw3
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group#
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
#LIBSF = -lm -L/usr/lib/x86_64-linux-gnu -lfftw3
LIB_CURIE = $(MKL_LIBS_MT) ${MKL_SCA_LIBS}
#INC = -I/usr/include/ 
#INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/

BINPATH=../bin
WAVET=../../WaveT/src/

OBJ_TDPLAS= constants.o global_tdplas.o pedra_friends.o readio_medium.o MathTools.o BEM_medium.o scf.o td_contmed.o QM_coupling.o main_tdplas.o
OBJ_FREQ= constants.o global_tdplas.o pedra_friends.o readio_medium.o MathTools.o BEM_medium.o scf.o td_contmed.o main_freq.o 

# rules

all: $(BINPATH)  $(BINPATH)/freq_bem.x $(BINPATH)/tdplas.x
	
$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/tdplas.x: $(OBJ_TDPLAS) 
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_TDPLAS) $(LIBS) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/freq_bem.x: $(OBJ_FREQ) 
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_FREQ) $(LIBS) $(LIBSF) $(LIB_CURIE)

clean:
	rm -f *.o *.mod *.il *.optrpt $(BINPATH)/*.x

constants.o: constants.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

global_tdplas.o: global_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

pedra_friends.o: pedra_friends.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

readio_medium.o: readio_medium.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

MathTools.o: MathTools.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

BEM_medium.o: BEM_medium.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

QM_coupling.o: QM_coupling.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

scf.o: scf.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

td_contmed.o: td_contmed.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

main_tdplas.o: main_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

main_freq.o: main_freq.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

global_wavet.o: $(WAVET)/global_wavet.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

# dependencies

global_tdplas.o: constants.o

pedra_friends.o: constants.o global_tdplas.o

readio_medium.o: constants.o global_tdplas.o pedra_friends.o

BEM_medium.o: constants.o global_tdplas.o pedra_friends.o readio_medium.o

QM_coupling.o: constants.o global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

scf.o: constants.o global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

td_contmed.o: constants.o global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o scf.o

main_tdplas.o: readio_medium.o BEM_medium.o

main_freq.o: readio_medium.o td_contmed.o


