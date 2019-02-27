FC = gfortran

# open mp flags
# comment first line and uncommment the second for ifort

OMP= -fopenmp -DOMP
#OMP= -openmp -DOMP

# mpi flags

MPI=

# other flags - to clean up?

LDLAGS=

# actual flags used in the compilation
FFLAGS= -O3 ${OMP} #$(FFLAGS_HYD_PAR) -fpp #$(OMP) -pg

# preprocessor flag
# comment the first and uncomment the second if compiler is ifort

CPP= -cpp

# libraries to be linked

LIBS = -L/usr/lib -llapack -lblas -lfftw3
INC = -I/usr/include/
LIBSPAR= -L/opt/scali/lib64 -lfmpi -lmpi
INCPAR= -I/opt/scali/include64 -I/local/home/corni/emanuele/fftw-3.3.4/include/

# executable, tools and TDPlas directories

BINPATH=../bin
TDPLAS=../../TDPlas/src/
TOOLSPATH=tools

# tdplas library file

TDPLASLIB= $(TDPLAS)/libtdplas.a

# tdplas include statement (to fetch only tdplas.mod)

TDPLASINC= -I$(TDPLAS)

# objects

OBJ= constants.o interface_tdplas.o random.o readio.o spectra.o dissipation.o propagate.o main.o
OBJ_SPECTRA= constants.o pedra_friends.o interface_tdplas.o global_tdplas.o readio.o readio_medium.o MathTools.o BEM_medium.o QM_coupling.o scf.o td_contmed.o spectra.o main_spectra.o
OBJ_VIB= constants.o readio.o $(TOOLSPATH)/vib.o $(TOOLSPATH)/vib_tool.o
OBJ_COH= constants.o $(TOOLSPATH)/coherence_tool.o
OBJ_FFT= constants.o $(TOOLSPATH)/fft_tool.o
OBJ_POP_COH= constants.o $(TOOLSPATH)/pop_coh_tool.o
OBJ_NUCL= constants.o $(TOOLSPATH)/nucl_wavep_tool.o

# compilation and linking rules

serial: $(BINPATH) $(BINPATH)/WaveT-serial.x $(BINPATH)/make_spectra.x $(BINPATH)/make_vib.x $(BINPATH)/make_coherence.x $(BINPATH)/make_2dfft.x $(BINPATH)/make_pop_coh.x $(BINPATH)/make_nucl_wavep.x

parallel: $(BINPATH) $(BINPATH)/WaveT-parallel.x

$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/make_spectra.x: main_spectra.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ) $(TDPLASLIB) $(LIB_CURIE) $(LIBS) $(LIBSF)

$(BINPATH)/WaveT-serial.x: main.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ) $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/WaveT-parallel.x: main.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(INCPAR) $(TDPLASLIB) $(LIBS) $(LIBSPAR) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/make_vib.x: $(OBJ_VIB)
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_VIB)

$(BINPATH)/make_coherence.x: $(OBJ_COH)
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_COH)

$(BINPATH)/make_2dfft.x: $(OBJ_FFT)
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_FFT) $(LIBS) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/make_pop_coh.x: $(OBJ_POP_COH)
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_POP_COH)

$(BINPATH)/make_nucl_wavep.x: $(OBJ_NUCL)
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_NUCL)

clean:
	rm -f *.o *.mod *.il *.optrpt $(BINPATH)/*.x $(TOOLSPATH)/*.o $(TOOLSPATH)/*.mod $(TOOLSPATH)/*.il $(TOOLSPATH)/*.optrpt

clean_obj:
	rm -f *.o *.mod *.il *.optrpt $(TOOLSPATH)/*.o $(TOOLSPATH)/*.mod $(TOOLSPATH)/*.il $(TOOLSPATH)/*.optrpt

constants.o: constants.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

interface_tdplas.o: interface_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) $(TDPLASINC) -c $< -o $@

random.o: random.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readio.o: readio.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

spectra.o: spectra.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

main_spectra.o: main_spectra.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

dissipation.o: dissipation.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

propagate.o: propagate.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/coherence_tool.o: $(TOOLSPATH)/coherence_tool.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/fft_tool.o: $(TOOLSPATH)/fft_tool.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

$(TOOLSPATH)/vib.o: $(TOOLSPATH)/vib.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/vib_tool.o: $(TOOLSPATH)/vib_tool.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/pop_coh_tool.o: $(TOOLSPATH)/pop_coh_tool.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/nucl_wavep_tool.o: $(TOOLSPATH)/nucl_wavep_tool.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

main.o: main.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

# dependencies
# all tdplas dependencies are contained in the static linking to the TDPlas library

random.o:

constants.o:

interface_tdplas.o: constants.o

readio.o: interface_tdplas.o random.o

spectra.o: readio.o

dissipation.o: random.o constants.o interface_tdplas.o readio.o

propagate.o: spectra.o dissipation.o

main_spectra.o: interface_tdplas.o readio.o spectra.o

main.o: constants.o interface_tdplas.o readio.o spectra.o dissipation.o propagate.o

vib.o: readio.o

vib_tool.o: readio.o vib.o

coherence_tool.o:

fft_tool.o:

pop_coh_tool.o:

nucl_wavep_tool.o:
