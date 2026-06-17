#FC = mpiifort
FC = gfortran
# open mp flags
 # comment first line and uncommment the second for ifort

 OMP= -fopenmp -DOMP
 #OMP= -qopenmp -DOMP

 # mpi flags

 #MPI= -DMPI

 # other flags - to clean up?

 #LDLAGS=

 # actual flags used in the compilation
 #FFLAGS= -O3 ${OMP} #$(FFLAGS_HYD_PAR) -fpp #$(OMP) -pg
 ###########FFLAGS= -Og -g ${OMP} #$(FFLAGS_HYD_PAR) -fpp #$(OMP) -pg
 FFLAGS= -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace ${OMP} 
 #FFLAGS= -O3 -ftz -align all ${OMP} ${MPI} -heap-arrays 1024 -check bounds 

 # preprocessor flag
 # comment the first and uncomment the second if compiler is ifort

 CPP= -cpp


#LIBS = -qmkl -lfftw3 -lm
#INC = -I/home/coccia/Libs/fftw-3.3.8
#INC = -I/home/coccia/Codes/WaveT_nondipole/WaveT/WaveT/src
INC = -I/usr/include/x86_64-linux-gnu
LIBS = -L/usr/lib/x86_64-linux-gnu -llapack -lblas -lfftw3
LIBSPAR= #-L/opt/scali/lib64 -lfmpi -lmpi
INCPAR= #-I/opt/scali/include64 -I/local/home/corni/emanuele/fftw-3.3.4/include/

# executable and tools directories

BINPATH=../bin
TOOLSPATH=tools

# objects

WAVETOBJ = random.o constants.o readio.o spectra.o  WTMathTools.o interface_classic.o scf.o QM_coupling.o dissipation.o initialise.o propagate.o
OBJ= $(WAVETOBJ) main.o
OBJ_SPECTRA= $(WAVETOBJ) main_spectra.o
OBJ_VIB= random.o constants.o readio.o $(TOOLSPATH)/vib.o $(TOOLSPATH)/vib_tool.o
OBJ_COH= constants.o $(TOOLSPATH)/coherence_tool.o
OBJ_FFT= constants.o $(TOOLSPATH)/fft_tool.o
OBJ_POP_COH= constants.o $(TOOLSPATH)/pop_coh_tool.o
OBJ_NUCL= constants.o $(TOOLSPATH)/nucl_wavep_tool.o

# tdplas binding flag
# write no if you want to compile wavet without any link to tdplas

TDPLASLINK= YES 
TDPLASLIB=
TDPLASINC=
TDPLASPP=

#ifeq ($(TDPLASLINK),YES)

# tdplas directory

TDPLAS=../../TDPlas/src/

# tdplas library file

TDPLASLIB= $(TDPLAS)/libtdplas.a

# tdplas include statement (to fetch only tdplas.mod)

TDPLASINC= -I$(TDPLAS)

# tdplas preprocessor flag

TDPLASPP= -DTDPLAS

#endif

# compilation and linking rules

serial: $(BINPATH) scf.o $(BINPATH)/WaveT-serial.x $(BINPATH)/make_spectra.x $(BINPATH)/make_vib.x $(BINPATH)/make_coherence.x $(BINPATH)/make_2dfft.x $(BINPATH)/make_pop_coh.x $(BINPATH)/make_nucl_wavep.x

parallel: $(BINPATH) $(BINPATH)/WaveT-parallel.x

$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/make_spectra.x: main_spectra.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ_SPECTRA) $(TDPLASLIB) $(LIB_CURIE) $(LIBS) $(LIBSF)

$(BINPATH)/WaveT-serial.x: main.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ) $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

$(BINPATH)/WaveT-parallel.x: main.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $(OBJ) $(INCPAR) $(TDPLASLIB) $(LIBS) $(LIBSPAR) $(LIBSF) $(LIB_CURIE)

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

random.o: random.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

constants.o: constants.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

interface_classic.o: interface_classic.f90
	$(FC) $(CPP) $(TDPLASPP) $(FFLAGS) $(INC) $(TDPLASINC) -c $< -o $@

readio.o: readio.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

spectra.o: spectra.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

QM_coupling.o: QM_coupling.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

scf.o: scf.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

initialise.o: initialise.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

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

main_spectra.o: main_spectra.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

WTMathTools.o: WTMathTools.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

# dependencies

random.o:

constants.o:

WTMathTools.o: readio.o constants.o

interface_classic.o: constants.o WTMathTools.o

readio.o: random.o

spectra.o: readio.o

QM_coupling.o: constants.o readio.o interface_classic.o initialise.o

scf.o: constants.o readio.o interface_classic.o initialise.o

initialise.o: constants.o readio.o 

dissipation.o: random.o constants.o readio.o interface_classic.o

propagate.o: random.o spectra.o dissipation.o interface_classic.o initialise.o

main_spectra.o: readio.o spectra.o interface_classic.o

main.o: constants.o readio.o spectra.o initialise.o dissipation.o propagate.o interface_classic.o

vib.o: readio.o

vib_tool.o: readio.o vib.o

coherence_tool.o:

fft_tool.o:

pop_coh_tool.o:

nucl_wavep_tool.o:
