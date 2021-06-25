# compilers
#FC = ifort
FC = gfortran
#FC = mpiifort
# useful predefined compiler flags for different architectures and compilers
# they are not necessarily used in the compiling unless you refer to them in FFLAGS
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

OMP= -fopenmp -DOMP
#OMP= -qopenmp -DOMP

# mpi flags (not used here)
MPI= -DMPI

LDLAGS=
#LDLAGS= --static --disable-shared --enable-static

# actual flags used in the compilation
#FFLAGS= $(FFLAGS_GF) ${OMP} #-heap-arrays 1024
#FFLAGS= -O3 -ftz -align all ${MPI} ${OMP} -heap-arrays 1024 #-check bounds #$(FFLAGS_HYD_PAR) -fpp #$(FFLAGS_OMP_IF) -pg #-check bounds
#FFLAGS= -O -g ${OMP} #-heap-arrays 1024
#FFLAGS= -O3 -ftz -align all ${OMP} -heap-arrays 1024 #-check bounds #$(FFLAGS_HYD_PAR) -fpp #$(FFLAGS_OMP_IF) -pg #-check bounds
#FFLAGS= -Og -g -fcheck-array-temporaries ${OMP} -fbounds-check #debug flag list to check coherence of array dimension 
FFLAGS= -Og -g ${OMP} #-heap-arrays 1024  

# preprocessor flag
# comment the first and uncomment the second if compiler is ifort
CPP=-cpp
#CPP= -fpp

# libraries to be linked - paths to be set by the user
#LIBS = -L/usr/lib -llapack -lblas -lpthread #FOR PADOVA
LIBS = -L/usr/lib -llapack  -lblas -lfftw3 #MARTA PC
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group#FOR HYDRA
INC =  -I/usr/local/include
# executable directory
BINPATH=../../bin
# tdplas library file

TDPLASLIB= libtdplas.a

# QM code to be interfaced with TDPlas

QMCODE= WaveT
QMCODE_OK= $(shell echo $(QMCODE) | tr A-Z a-z)
#INTERFACE_MODULE= interface_$(QMCODE_OK)

INTERFACE_MODULE = global_quantum

# objects to be linked with executables

OBJ= tdplas.o constants.o auxiliary_functions.o dictionary_entry.o tdplas_readio_dictionary.o input_keyword_type.o lower_and_check_namelist.o write_header_out_tdplas.o readio_input_namelist.o readio_medium.o init.o global_quantum.o global_tdplas.o debye_epsilon.o drudel_epsilon.o readfile_epsilon.o do_epsilon.o cavity_types.o sphe_surface.o pedra_friends.o read_pot.o quantum_coupling_modes.o 


# compilation rules

all: $(BINPATH) $(TDPLASLIB) $(BINPATH)/freq_bem.x $(BINPATH)/tdplas.x $(BINPATH)/print_eps.x

$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/tdplas.x: main_tdplas.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) 

$(BINPATH)/freq_bem.x: main_freq.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) 

$(BINPATH)/eps_gold.x: main_eps.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

clean:
	rm -f *.o *.mod *.il *.optrpt *.a $(BINPATH)/*.x

lib: $(TDPLASLIB)

$(TDPLASLIB): $(OBJ)
	ar rcs $@ $^
	ranlib $@


constants.o: constants.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

global_quantum.o: $(INTERFACE_MODULE).f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

auxiliary_functions.o: auxiliary_functions.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

dictionary_entry.o: dictionary_entry.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

input_keyword_type.o: input_keyword_type.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

cavity_types.o: cavity_types.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

pedra_friends.o: pedra_friends.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

sphe_surface.o: sphe_surface.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

debye_epsilon.o: debye_epsilon.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

drudel_epsilon.o: drudel_epsilon.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readfile_epsilon.o: readfile_epsilon.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

do_epsilon.o: do_epsilon.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

quantum_coupling_modes.o: quantum_coupling_modes.f90
	 $(FC) $(CPP) $(FFLAGS) -c $< -o $@

read_pot.o: read_pot.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

tdplas_readio_dictionary.o: tdplas_readio_dictionary.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

lower_and_check_namelist.o: lower_and_check_namelist.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readio_input_namelist.o: readio_input_namelist.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readio_medium.o: readio_medium.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

global_tdplas.o: global_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

init.o: init.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

write_header_out_tdplas.o: write_header_out_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

#MathTools.o: MathTools.f90
#	$(FC) $(CPP) $(FFLAGS) -c $< -o $@
#
#BEM_medium.o: BEM_medium.f90
#	$(FC) $(CPP) $(FFLAGS) -c $< -o $@
#
#td_contmed.o: td_contmed.f90
#	$(FC) $(CPP) $(FFLAGS) -c $< -o $@
#
#tdplas.o: tdplas.f90
#	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

tdplas.o: tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

main_tdplas.o: main_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

main_freq.o: main_freq.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

main_eps.o: main_eps.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

# dependencies

auxiliary_functions.o: constants.o


dictionary_entry.o: constants.o

tdplas_readio_dictionary.o: constants.o dictionary_entry.o

input_keyword_type.o: constants.o

lower_and_check_namelist.o: tdplas_readio_dictionary.o dictionary_entry.o 

write_header_out_tdplas.o: global_tdplas.o sphe_surface.o pedra_friends.o drudel_epsilon.o debye_epsilon.o readfile_epsilon.o do_epsilon.o quantum_coupling_modes.o 


readio_input_namelist.o: constants.o input_keyword_type.o lower_and_check_namelist.o

readio_medium.o: constants.o input_keyword_type.o readio_input_namelist.o tdplas_readio_dictionary.o global_quantum.o


init.o: constants.o input_keyword_type.o dictionary_entry.o tdplas_readio_dictionary.o global_tdplas.o 


global_quantum.o: constants.o

global_tdplas.o: constants.o auxiliary_functions.o read_pot.o  sphe_surface.o pedra_friends.o quantum_coupling_modes.o readfile_epsilon.o drudel_epsilon.o debye_epsilon.o do_epsilon.o



debye_epsilon.o: constants.o

drudel_epsilon.o: constants.o

readfile_epsilon.o: constants.o

do_epsilon.o: constants.o readfile_epsilon.o drudel_epsilon.o debye_epsilon.o



cavity_types.o: constants.o

sphe_surface.o: constants.o auxiliary_functions.o

pedra_friends.o: constants.o cavity_types.o


read_pot.o: constants.o global_quantum.o

quantum_coupling_modes.o: constants.o


#MathTools.o: constants.o global_quantum.o pedra_friends.o global_tdplas.o

#BEM_medium.o: constants.o global_quantum.o pedra_friends.o sphe_surface.o global_tdplas.o MathTools.o eps_module.o quantum_coupling_modes.o debye_epsilon.o drudel_epsilon.o inout_epsilon.o cavity_types.o

#td_contmed.o: constants.o global_quantum.o pedra_friends.o sphe_surface.o global_tdplas.o MathTools.o BEM_medium.o inout_epsilon.o drudel_epsilon.o debye_epsilon.o 

tdplas.o: constants.o auxiliary_functions.o dictionary_entry.o tdplas_readio_dictionary.o input_keyword_type.o lower_and_check_namelist.o write_header_out_tdplas.o readio_input_namelist.o readio_medium.o init.o global_quantum.o global_tdplas.o debye_epsilon.o drudel_epsilon.o readfile_epsilon.o do_epsilon.o cavity_types.o sphe_surface.o pedra_friends.o read_pot.o quantum_coupling_modes.o  

main_tdplas.o: tdplas.o

main_freq.o: tdplas.o

main_eps.o: tdplas.o
