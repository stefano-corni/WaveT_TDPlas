# compilers
#FC = ifort
FC = gfortran
#FC = mpiifort

OMP= -fopenmp -DOMP

# mpi flags (not used here)
MPI= -DMPI

LDLAGS=
#Debug flags
#FFLAGS= -Og -g ${OMP} -fPIC -m64 #-heap-arrays 1024  
#Production flags
FFLAGS= -O3 ${OMP} -fPIC -m64 #-heap-arrays 1024  

# preprocessor flag
# comment the first and uncomment the second if compiler is ifort
CPP=-cpp
#CPP= -fpp

# libraries to be linked - paths to be set by the user
LIBS = -L/usr/lib -llapack -lblas -lpthread 
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group#FOR HYDRA
INC =  -I/usr/local/include
# executable directory
BINPATH=../bin
# tdplas library file

TDPLASLIB= libtdplas.a

# QM code to be interfaced with TDPlas = wt, ocpy
QMCODE= wt
QMCODE_OK= $(shell echo $(QMCODE) | tr A-Z a-z)
INTERFACE_MODULE= interface_global_quantum_$(QMCODE_OK)


# objects to be linked with executables
 
OBJ = tdplas.o tdplas_constants.o auxiliary_functions.o cavity_types.o dictionary_entry.o readio_tdplas_mod.o read_inputfile_tdplas.o check_global.o lower_and_check_namelist.o init.o global_quantum.o global_tdplas.o  debye_epsilon.o drudel_epsilon.o readfile_freq.o readfile_epsilon.o dielectric_function.o  sphe_surface.o pedra_friends.o user_input_type.o user_input_read_namelist.o user_input_check.o user_input_and_flags_dictionary.o write_header_out_tdplas.o MathTools.o BEM_medium.o td_contmed.o 


# compilation rules

all: $(BINPATH) $(TDPLASLIB) $(BINPATH)/freq_bem.x $(BINPATH)/tdplas.x $(BINPATH)/main_eps.x

$(BINPATH):
	mkdir $(BINPATH)

$(BINPATH)/tdplas.x: main_tdplas.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) 

$(BINPATH)/freq_bem.x: main_freq.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) 

$(BINPATH)/main_eps.x: main_eps.o
	$(FC) $(CPP) $(FFLAGS) -o $@ $< $(TDPLASLIB) $(LIBS) $(LIBSF) $(LIB_CURIE)

clean:
	rm -f *.o *.mod *.il *.optrpt *.a $(BINPATH)/*.x

lib: $(TDPLASLIB)

$(TDPLASLIB): $(OBJ)
	ar rcs $@ $^
	ranlib $@


tdplas_constants.o: tdplas_constants.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

global_quantum.o: $(INTERFACE_MODULE).f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

auxiliary_functions.o: auxiliary_functions.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

cavity_types.o: cavity_types.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

dictionary_entry.o: dictionary_entry.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

read_inputfile_tdplas.o: read_inputfile_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readio_tdplas_mod.o: readio_tdplas_mod.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

check_global.o: check_global.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

lower_and_check_namelist.o: lower_and_check_namelist.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

init.o: init.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

global_tdplas.o: global_tdplas.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

readfile_freq.o: readfile_freq.f90
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

dielectric_function.o: dielectric_function.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

user_input_and_flags_dictionary.o: user_input_and_flags_dictionary.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

user_input_check.o: user_input_check.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

user_input_read_namelist.o: user_input_read_namelist.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

user_input_type.o: user_input_type.f90
	$(FC) $(CPP) $(FFLAGS) -c $< -o $@

write_header_out_tdplas.o: write_header_out_tdplas.f90
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

main_eps.o: main_eps.f90
	$(FC) $(CPP) $(FFLAGS) $(INC) -c $< -o $@

# dependencies

auxiliary_functions.o: tdplas_constants.o


dictionary_entry.o: tdplas_constants.o

user_input_and_flags_dictionary.o: tdplas_constants.o dictionary_entry.o

user_input_type.o: tdplas_constants.o

readfile_freq.o: tdplas_constants.o tdplas_constants.o

user_input_check.o: tdplas_constants.o auxiliary_functions.o user_input_type.o

user_input_read_namelist.o: tdplas_constants.o user_input_type.o lower_and_check_namelist.o

lower_and_check_namelist.o: user_input_and_flags_dictionary.o dictionary_entry.o 

write_header_out_tdplas.o: global_tdplas.o sphe_surface.o pedra_friends.o drudel_epsilon.o debye_epsilon.o readfile_epsilon.o  

read_inputfile_tdplas.o: tdplas_constants.o user_input_type.o user_input_read_namelist.o  

readio_tdplas_mod.o: tdplas_constants.o user_input_type.o read_inputfile_tdplas.o user_input_check.o init.o check_global.o write_header_out_tdplas.o user_input_and_flags_dictionary.o



init.o: tdplas_constants.o user_input_type.o dictionary_entry.o user_input_and_flags_dictionary.o global_tdplas.o pedra_friends.o 


global_quantum.o: tdplas_constants.o user_input_type.o read_inputfile_tdplas.o user_input_check.o init.o check_global.o write_header_out_tdplas.o user_input_and_flags_dictionary.o global_tdplas.o readfile_freq.o

global_tdplas.o: tdplas_constants.o auxiliary_functions.o sphe_surface.o readfile_epsilon.o drudel_epsilon.o debye_epsilon.o dielectric_function.o 

check_global.o: tdplas_constants.o auxiliary_functions.o global_tdplas.o pedra_friends.o 

debye_epsilon.o: tdplas_constants.o

drudel_epsilon.o: tdplas_constants.o

readfile_epsilon.o: tdplas_constants.o

dielectric_function.o: tdplas_constants.o auxiliary_functions.o readfile_epsilon.o drudel_epsilon.o debye_epsilon.o  



cavity_types.o: tdplas_constants.o

sphe_surface.o: tdplas_constants.o auxiliary_functions.o

pedra_friends.o: tdplas_constants.o cavity_types.o global_tdplas.o




MathTools.o: tdplas_constants.o global_quantum.o pedra_friends.o global_tdplas.o readfile_freq.o 

BEM_medium.o: tdplas_constants.o global_quantum.o readfile_freq.o pedra_friends.o sphe_surface.o global_tdplas.o MathTools.o  debye_epsilon.o drudel_epsilon.o readfile_epsilon.o dielectric_function.o cavity_types.o

td_contmed.o: tdplas_constants.o global_quantum.o pedra_friends.o sphe_surface.o global_tdplas.o MathTools.o BEM_medium.o readfile_epsilon.o drudel_epsilon.o debye_epsilon.o readfile_freq.o 

tdplas.o: tdplas_constants.o auxiliary_functions.o cavity_types.o dictionary_entry.o read_inputfile_tdplas.o readio_tdplas_mod.o check_global.o lower_and_check_namelist.o init.o global_quantum.o global_tdplas.o  debye_epsilon.o drudel_epsilon.o readfile_epsilon.o dielectric_function.o  sphe_surface.o pedra_friends.o user_input_type.o user_input_read_namelist.o user_input_check.o user_input_and_flags_dictionary.o write_header_out_tdplas.o MathTools.o BEM_medium.o td_contmed.o  

main_tdplas.o: tdplas.o

main_freq.o: tdplas.o

main_eps.o: tdplas.o
