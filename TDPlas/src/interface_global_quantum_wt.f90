        module global_quantum
            use tdplas_constants
            use user_input_type
            use read_inputfile_tdplas 
            use user_input_check
            use init
            use check_global
            use write_header_out_tdplas
            use user_input_and_flags_dictionary
            use global_tdplas
            use readfile_freq

#ifdef MPI
#ifndef SCALI
            use mpi
#endif
#endif

            implicit none

#ifdef MPI
#ifdef SCALI
            include 'mpif.h'
#endif
#endif

            ! block of global variables to be supplied by WaveT
            character(flg)            :: quantum_Ffld                             ! shape of impulse                            #test pot
            character(flg)            :: quantum_Fbin                             ! binary output                               # prop_mdm per stampare
            character(flg)            :: quantum_Fopt                             ! OMP-optimized matrix/vector multiplication  # non serve perchè ci guarda da solo
            real(dbl)                 :: quantum_dt                               ! time step
            integer(i4b)              :: quantum_n_ci, quantum_n_ci_read          ! number of CIS states
            real(dbl), allocatable    :: quantum_e_ci(:)                          ! CIS energies                                #no
            complex(cmp), allocatable :: quantum_c_i(:)                          ! CIS coefficients                             #no
            real(dbl), allocatable    :: quantum_mut(:,:,:)                       ! CIS transition dipoles
            real(dbl)                 :: quantum_mol_cc(3)                        ! molecule center                              #
            real(dbl)                 :: quantum_fmax(3,10),quantum_omega(10)     ! field amplitude and frequency                #boh ref_mu
            !real(dbl)                 :: quantum_tdelay(10),quantum_pshift(10)    ! time delay and phase shift
            integer(i4b)              :: quantum_n_out,quantum_n_f                ! auxiliaries for output
            integer(i4b)              :: quantum_nthr = 1                            ! number of threads
            !character(1)              :: quantum_medium_res                       !restart for medium
            integer(i4b)              :: quantum_n_res                            ! frequency for restart





            public quantum_init, readio_and_init_tdplas_for_wt

            contains


            subroutine quantum_init(dt, mol_cc, n_ci, n_ci_read, mut, e_ci, c_i, &
                                   fmax, omega, Ffld, n_out, n_f, &
                                   Fbin,  Fopt, n_res)
!------------------------------------------------------------------------
! @brief Set global variables for medium
!
! @date Created: S. Pipolo
! Modified: E. Coccia 28/11/17
!------------------------------------------------------------------------
                implicit none

                real(dbl)     , intent(in)  :: dt                                    ! time step
                integer(i4b)  , intent(in)  :: n_ci, n_ci_read                     ! number of CIS states
                real(dbl)     , intent(in)  :: e_ci(:)
                complex(cmp)  , intent(in)  :: c_i(:)
                real(dbl)     , intent(in)  :: mut(:,:,:)               ! CIS transition dipoles
                real(dbl)     , intent(in)  :: mol_cc(3)                ! molecule center
                real(dbl)     , intent(in)  :: fmax(3,10),omega(10)     ! field amplitude and frequency
                character(3)  , intent(in)  :: Ffld                     ! shape of impulse
                character(3)  , intent(in)  :: Fbin                     ! binary output
                character(3)  , intent(in)  :: Fopt                     ! matrix/vector multiplication
                integer(i4b)  , intent(in)  :: n_out,n_f                ! auxiliaries for output
                integer(i4b)  , intent(in)  :: n_res                    ! frequency for restart

                quantum_dt=dt
                quantum_mol_cc=mol_cc
                quantum_n_ci=n_ci
                quantum_n_ci_read = n_ci_read
                allocate(quantum_e_ci(quantum_n_ci))
                quantum_e_ci = e_ci
                quantum_c_i = c_i
                allocate(quantum_mut(3,quantum_n_ci,quantum_n_ci))
                quantum_mut=mut
                quantum_fmax=fmax
                quantum_omega=omega
                quantum_Ffld=Ffld
                quantum_n_out=n_out
                quantum_n_f=n_f
                quantum_Fbin=Fbin
                quantum_Fopt=Fopt
                quantum_n_res=n_res

                return

            end subroutine


             subroutine readio_and_init_tdplas_for_wt(nthr)
                character(flg)  :: calculation_exe = "propagation"
                type(tdplas_user_input)  :: user_input
                integer :: nthr

                call set_default_input_tdplas_for_wt(user_input)
                call dict_tdplas_init                                          !init dictionary to convert from user friendly values to internal

                call read_input_tdplas_for_propagation(user_input)

                call check_tdplas_input_for_wt(user_input)
                call init_tdplas(calculation_exe, nthr, user_input)
                call check_global_var
                call read_gau_out_medium(quantum_n_ci)
                call write_out(calculation_exe)
            end subroutine






        subroutine set_default_input_tdplas_for_wt(user_input)
            type(tdplas_user_input)  :: user_input

           user_input%test_type          = "non"
           user_input%debug_type         = "non"
           user_input%out_level          = "low"

           user_input%medium_type       =  "nanop"
           user_input%medium_init0      =  "frozen"
           user_input%medium_pol        =  "charge"
           user_input%bem_type          =  "diagonal"
           user_input%bem_read_write    =  "read"
           user_input%local_field       =  "non"
           user_input%normalization     =  "non"       

           user_input%surface_type               = "mesh"
           user_input%input_mesh                 = "non"
           user_input%particles_number           = 1
           user_input%spheres_number             = 0
           user_input%spheroids_number           = 0
           user_input%object_shape               = "non"
           user_input%sphere_position_x(nsmax)   = zero
           user_input%sphere_position_y(nsmax)   = zero
           user_input%sphere_position_z(nsmax)   = zero
           user_input%sphere_radius(nsmax)       = zero
           user_input%spheroid_position_x(nsmax) = zero
           user_input%spheroid_position_y(nsmax) = zero
           user_input%spheroid_position_z(nsmax) = zero
           user_input%spheroid_radius(nsmax)     = zero
           user_input%spheroid_axis_x(nsmax)     = zero
           user_input%spheroid_axis_y(nsmax)     = zero
           user_input%spheroid_axis_z(nsmax)     = zero
           user_input%inversion                  = "non"
           user_input%find_spheres               = "yes"

           user_input%epsilon_omega     = "non"
           user_input%eps_A             = -1.0
           user_input%eps_gm            = -1.0
           user_input%eps_w0            = -1.0
           user_input%f_vel             = -1.0
           user_input%eps_0             = -1.0
           user_input%eps_d             = -1.0
           user_input%tau_deb           = -1.0
           user_input%propagation_pole  =  "velocity-verlet"

           user_input%n_omega           = -1
           user_input%omega_ini         = -1.0
           user_input%omega_end         = -1.0

           user_input%propagation_software = "wavet"
           user_input%propagation_type     = "charge-ief"
           user_input%interaction_stride  = 1
           user_input%mix_coef            = 0.2
           user_input%max_cycles          = 600
           user_input%threshold           = 10
           user_input%interaction_type    = "pcm"
           user_input%interaction_init    = "nsc"
           user_input%medium_relax        = "non"
           user_input%medium_restart      = "non"

           user_input%gamess             = "non"
           user_input%print_lf_matrix    = "non"

           user_input%n_print_charges  = 0
           user_input%charge_mopac     = "non"

           user_input%n_modes  = 10
           user_input%quantum_calculation = "non"

           user_input%pert_type ="field"
           user_input%direction(1)         = 0.0
           user_input%direction(2)         = 0.0
           user_input%direction(3)         = 1.0
           user_input%eet ="non"
        end subroutine



        subroutine check_tdplas_input_for_wt(user_input)
            type(tdplas_user_input) ::  user_input

            if((user_input%medium_type.ne."nanop").and.(user_input%medium_type.ne."sol")) call wavet_error("medium_type")
            if(user_input%propagation_software.ne."wavet") call wavet_error("propagation_software")


            !if(user_input%propagation_type.eq."dipole") then
            !  if(user_input%epsilon_omega.eq."drude-lorentz") call tdplas_couple_error("propagation_type", "epsilon_omega")
            !endif
            call check_surface(user_input)
            call check_spheres_or_spheroids(user_input)
            call check_eps(user_input)
            end subroutine


        end module
