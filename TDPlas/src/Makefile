# compilers

FC = gfortran
#FC = ifort

# useful predefined compiler flags for different architectures and compilers
# they are not necessarily used in the compiling unless you refer to them in FFLAGS

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

# open mp flags
# comment first line and uncommment the second for ifort

OMP= -fopenmp -DOMP
#OMP= -openmp -DOMP

# mpi flags (not used here)

MPI= -DMPI

# other flags - to clean up?

LDLAGS=
#LDLAGS= --static --disable-shared --enable-static

# actual flags used in the compilation
# comment the first and uncomment the second if compiler is ifort

FFLAGS= $(FFLAGS_GF) ${OMP} #-heap-arrays 1024
#FFLAGS= -O3 -ftz -align all ${OMP_IF} -heap-arrays 1024 #-check bounds #$(FFLAGS_HYD_PAR) -fpp #$(FFLAGS_OMP_IF) -pg #-check bounds

# preprocessor flag
# comment the first and uncomment the second if compiler is ifort

CPP=-cpp
#CPP= -fpp

# libraries to be linked

LIBS = -L/usr/lib -llapack -lblas -lpthread
#LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
#LIBS = -L/usr/lib -llapack  -lblas -lfftw3
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group#
#INC =  -I/local/home/corni/emanuele/fftw-3.3.4/include/

# executable directory

BINPATH=../bin

# tdplas library file

TDPLASLIB= libtdplas.a

# QM code to be interfaced with TDPlas

QMCODE= WaveT
QMCODE_OK= $(shell echo $(QMCODE) | tr A-Z a-z)
INTERFACE_MODULE= interface_$(QMCODE_OK)

# objects to be linked with executables

OBJ= tdplas.o td_contmed.o BEM_medium.o MathTools.o readio_medium.o pedra_friends.o interface_qmcode.o constants.o

# compilation rules

all: $(BINPATH) $(TDPLASLIB) $(BINPATH)/freq_bem.x $(BINPATH)/tdplas.x

$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/tdplas.x: main_tdplas.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/freq_bem.x: main_freq.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

clean:
	rm -f *.o *.mod *.il *.optrpt *.a $(BINPATH)/*.x

lib: $(TDPLASLIB)

$(TDPLASLIB): $(OBJ)
	ar rcs $@ $^
	ranlib $@

constants.o: constants.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

interface_qmcode.o: $(INTERFACE_MODULE).f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

pedra_friends.o: pedra_friends.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

readio_medium.o: readio_medium.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

MathTools.o: MathTools.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

BEM_medium.o: BEM_medium.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

td_contmed.o: td_contmed.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

tdplas.o: tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

main_tdplas.o: main_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

main_freq.o: main_freq.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

# dependencies

interface_qmcode.o: constants.o

pedra_friends.o: constants.o interface_qmcode.o

readio_medium.o: constants.o interface_qmcode.o pedra_friends.o

MathTools.o: constants.o interface_qmcode.o pedra_friends.o readio_medium.o

BEM_medium.o: constants.o interface_qmcode.o pedra_friends.o readio_medium.o MathTools.o

td_contmed.o: constants.o interface_qmcode.o pedra_friends.o readio_medium.o MathTools.o BEM_medium.o

tdplas.o: constants.o interface_qmcode.o pedra_friends.o readio_medium.o MathTools.o BEM_medium.o td_contmed.o

main_tdplas.o: tdplas.o

main_freq.o: tdplas.o
