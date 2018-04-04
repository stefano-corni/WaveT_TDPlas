#FC = gfortran
#FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
#FFLAGS_PAR= -parallel -par_threshold75 -par_report3
#FFLAGS_HYD= -O3 -ftz -align -xHost
#FFLAGS_HYD_PAR= -O3 -ftz -align -xHost -parallel -par_threshold75 -par_report3
#FFLAGS_DEB= -debug -C -traceback
#FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
#FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
#FFLAGS_MIO= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xB -ipo
#FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
#FFLAGS_LAP= -O3 -ftz -align -xAVX -no-prec-div 
#FFLAGS_GF= -O3
#LDLAGS= 
#FFLAGS= $(FFLAGS_GF_DEB)
#LIBS =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/

#FC = gfortran
FC = ifort
FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
FFLAGS_PAR= -parallel -par_threshold75 -par_report3
#FFLAGS_DEB= -debug -C -g
FFLAGS_DEB= -debug -traceback -C -g -diag-disable  warn -parallel -par_threshold75 -par_report3
FFLAGS_DEB_PAR= -g -diag-disable  warn -parallel -par_threshold75 -par_report3
FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
FFLAGS_MIO= -O3 -ftz -align -xAVX
FFLAGS_HYD= -O3 -ftz -align -xHost
FFLAGS_HYD_PAR= -O3 -ftz -align -xHost -parallel -par_threshold75 -par_report3 -check bounds 
FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
FFLAGS_GF= -O3
LDLAGS=
FFLAGS= $(FFLAGS_HYD_PAR)
#LIBS = -L/usr/lib -llapack  -lblas 
LIBS =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/



BINPATH=../bin
TDPLAS=../../TDPlas/src/
TOOLSPATH=tools

OBJ= constants.o global_wavet.o global_tdplas.o random.o readio.o readio_medium.o MathTools.o BEM_medium.o QM_coupling.o scf.o td_contmed.o spectra.o dissipation.o propagate.o main.o pedra_friends.o  
OBJ_SPECTRA= constants.o pedra_friends.o global_wavet.o global_tdplas.o readio.o readio_medium.o MathTools.o BEM_medium.o QM_coupling.o scf.o td_contmed.o spectra.o main_spectra.o 
OBJ_VIB= constants.o readio.o vib.o main_vib.o
OBJ_COH= constants.o $(TOOLSPATH)/coherence_tool.o 
OBJ_FFT= constants.o $(TOOLSPATH)/fft_tool.o
OBJ_NUCL= constants.o $(TOOLSPATH)/nucl_wavep_tool.o


# rules

all: $(BINPATH) $(BINPATH)/WaveT.x $(BINPATH)/make_spectra.x $(BINPATH)/make_vib.x $(BINPATH)/make_coherence.x $(BINPATH)/make_2dfft.x $(BINPATH)/make_nucl_wavep.x
	
$(BINPATH):
	mkdir $(BINPATH)
	
$(BINPATH)/make_spectra.x: $(OBJ_SPECTRA) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_SPECTRA) $(LIBS) $(LIBSF)

$(BINPATH)/WaveT.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS) $(LIBSF)

$(BINPATH)/make_vib.x: $(OBJ_VIB)
	$(FC) $(FFLAGS) -o $@ $(OBJ_VIB)

$(BINPATH)/make_coherence.x: $(OBJ_COH)
	$(FC) $(FFLAGS) -o $@ $(OBJ_COH)

$(BINPATH)/make_2dfft.x: $(OBJ_FFT)
	$(FC) $(FFLAGS) -o $@ $(OBJ_FFT) $(LIBS) $(LIBSF)
$(BINPATH)/make_nucl_wavep.x: $(OBJ_NUCL)
	$(FC) $(FFLAGS) -o $@ $(OBJ_NUCL)
clean:
	rm -f *.o *.mod *.il *.optrpt $(BINPATH)/*.x $(TOOLSPATH)/*.o $(TOOLSPATH)/*.mod $(TOOLSPATH)/*.il $(TOOLSPATH)/*.optrpt

constants.o: constants.f90
	$(FC) $(FFLAGS) -c $< -o $@

global_wavet.o: global_wavet.f90
	$(FC) $(FFLAGS) -c $< -o $@

random.o: random.f90
	$(FC) $(FFLAGS) -c $< -o $@

readio.o: readio.f90
	$(FC) $(FFLAGS) -c $< -o $@

spectra.o: spectra.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

main_spectra.o: main_spectra.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

dissipation.o: dissipation.f90
	$(FC) $(FFLAGS) -c $< -o $@

propagate.o: propagate.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/coherence_tool.o: $(TOOLSPATH)/coherence_tool.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(TOOLSPATH)/fft_tool.o: $(TOOLSPATH)/fft_tool.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

$(TOOLSPATH)/nucl_wavep_tool.o: $(TOOLSPATH)/nucl_wavep_tool.f90
	$(FC) $(FFLAGS) -c $< -o $@

vib.o: vib.f90
	$(FC) $(FFLAGS) -c $< -o $@

main_vib.o: main_vib.f90
	$(FC) $(FFLAGS) -c $< -o $@

main.o: main.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

MathTools.o: $(TDPLAS)/MathTools.f90
	$(FC) $(FFLAGS) -c $< -o $@

td_contmed.o: $(TDPLAS)/td_contmed.f90
	$(FC) $(FFLAGS) -c $< -o $@

global_tdplas.o: $(TDPLAS)/global_tdplas.f90
	$(FC) $(FFLAGS) -c $< -o $@

pedra_friends.o: $(TDPLAS)/pedra_friends.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

readio_medium.o: $(TDPLAS)/readio_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

BEM_medium.o: $(TDPLAS)/BEM_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

QM_coupling.o: $(TDPLAS)/QM_coupling.f90
	$(FC) $(FFLAGS) -c $< -o $@

scf.o: $(TDPLAS)/scf.f90
	$(FC) $(FFLAGS) -c $< -o $@


# dependencies

random.o: 

readio.o: global_wavet.o

global_wavet.o: td_contmed.o

spectra.o: global_wavet.o readio.o pedra_friends.o readio_medium.o

main_spectra.o: readio.o readio_medium.o spectra.o

dissipation.o: global_wavet.o readio.o random.o

propagate.o: global_wavet.o random.o readio.o spectra.o readio_medium.o pedra_friends.o td_contmed.o

vib.o: readio.o

main_vib.o: readio.o vib.o

coherence_tool.o: 

fft_tool.o: 

nucl_wavep_tool.o:

main.o: readio.o spectra.o dissipation.o propagate.o readio_medium.o QM_coupling.o  

global_tdplas.o: 

pedra_friends.o: global_tdplas.o

readio_medium.o: global_tdplas.o pedra_friends.o

BEM_medium.o: global_tdplas.o pedra_friends.o readio_medium.o

QM_coupling.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

scf.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

td_contmed.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o scf.o
