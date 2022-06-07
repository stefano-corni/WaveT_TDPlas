        module global_quantum
            use tdplas_constants
            use user_input_type
            use read_inputfile_tdplas 
            use user_input_check
            use init
            use check_global
            use write_header_out_tdplas
            use user_input_and_flags_dictionary
            use readfile_freq

            implicit none



            character(flg)            :: quantum_Ffld  = "non"
            character(flg)            :: quantum_Fbin  = "non"
            character(flg)            :: quantum_Fopt  = "non"
            integer(i4b)              :: quantum_n_ci_read = 0
            real(dbl), allocatable    :: quantum_e_ci(:)
            complex(cmp), allocatable :: quantum_c_i(:)
            real(dbl), allocatable    :: quantum_mut(:,:,:)
            real(dbl)                 :: quantum_mol_cc(3)
            real(dbl)                 :: quantum_fmax(3,10),quantum_omega(1)
            integer(i4b)              :: quantum_n_out = 0
            integer(i4b)              :: quantum_n_f   = 0
            integer(i4b)              :: quantum_n_res = 0





            real(dbl)                 :: quantum_dt                               ! time step
            integer(i4b)              :: quantum_n_ci          ! number of CIS states                  #

            public quantum_init, readio_and_init_tdplas_for_ocpy

            contains

            subroutine quantum_init(dt, &
                                   vts, n_tessere, n_ci)
!------------------------------------------------------------------------
! @brief Set global variables for medium
!
! @date Created: S. Pipolo
! Modified: E. Coccia 28/11/17
!------------------------------------------------------------------------
                implicit none

                real(dbl)     , intent(in)  :: dt                          ! time step
                integer(i4b)  , intent(in)  :: n_ci, n_tessere             ! number of CIS states
                real(dbl)     , intent(in)  :: vts(:,:,:)

                quantum_dt=dt
                quantum_n_ci=n_ci
                allocate (quantum_vts(n_tessere, n_ci, n_ci))
                quantum_vts = vts
                return

            end subroutine



            subroutine readio_and_init_tdplas_for_ocpy(nthr)

                character(flg)  :: calculation_exe = "propagation"
                type(tdplas_user_input)  :: user_input
                integer :: nthr
                call set_default_input_tdplas_for_ocpy(user_input)
                call dict_tdplas_init                                          !init dictionary to convert from user friendly values to internal

                call read_input_tdplas_for_propagation(user_input)

                call check_tdplas_input_for_ocpy(user_input)
                call init_tdplas(calculation_exe, nthr, user_input)

                call check_global_var
                call write_out(calculation_exe)

            end subroutine


            subroutine set_default_input_tdplas_for_ocpy(user_input)
            type(tdplas_user_input)  :: user_input

           user_input%test_type          = "non"
           user_input%debug_type         = "non"
           user_input%out_level          = "low"

           user_input%medium_type       =  "nanop"
           user_input%medium_init0      =  "frozen"
           user_input%medium_pol        =  "charge"
           user_input%bem_type          =  "diagonal"
           user_input%bem_read_write    =  "read"
           user_input%local_field         = "local"
           

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
 

           user_input%propagation_software = "ocpy"
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



            subroutine check_tdplas_input_for_ocpy(user_input)
                type(tdplas_user_input) ::  user_input

                if(user_input%test_type.ne."non") call ocpy_error("test_type")
                if(user_input%debug_type.ne."non") call ocpy_error("debug_type")
                if(user_input%out_level.ne."low") call ocpy_error("out_level")
                if((user_input%medium_type.ne."nanop").and.(user_input%medium_type.ne."sol")) call ocpy_error("medium_type")
                if(user_input%medium_pol.ne."charge") call ocpy_error("medium_pol")
                if(user_input%bem_read_write.ne."read") call ocpy_error("bem_read_write")
                if(user_input%surface_type.ne."mesh") call ocpy_error("surface_type")
                if(user_input%inversion.ne."non") call ocpy_error("inversion")
                if(user_input%n_omega.ne.-1) call ocpy_error("n_omega")
                if(user_input%omega_ini.ne.-1.0) call ocpy_error("omega_ini")
                if(user_input%omega_end.ne.-1.0) call ocpy_error("omega_end)")
                if(user_input%propagation_software.ne."ocpy") call ocpy_error("propagation_software")
                if(user_input%propagation_type.ne."charge-ief") call ocpy_error("propagation_type")
                if(user_input%interaction_stride.ne.1) call ocpy_error("interaction_stride")
                if(user_input%interaction_type.ne."pcm") call ocpy_error("interaction_type")
                if(user_input%interaction_init.ne."nsc") call ocpy_error("interaction_init")
                if(user_input%medium_relax.ne."non") call ocpy_error("medium_relax")
                if(user_input%medium_restart.ne."non") call ocpy_error("medium_restart")
                if(user_input%gamess.ne."non") call ocpy_error("gamess")
                if(user_input%print_lf_matrix.ne."non") call ocpy_error("print_lf_matrix")
                if(user_input%charge_mopac.ne."non") call ocpy_error("charge_mopac")
                if(user_input%quantum_calculation.ne."non") call ocpy_error("quantum_calculation")

                call check_tdplas_keys_for_ocpy(user_input)
                call check_surface(user_input)
                call check_spheres_or_spheroids(user_input)
                call check_eps(user_input)

            end subroutine

            subroutine check_tdplas_keys_for_ocpy(user_input)
                    type(tdplas_user_input) ::  user_input

                    if(user_input%bem_type.eq."diagonal") then
                            if(user_input%epsilon_omega.eq."general") call tdplas_couple_error("bem_type", "epsilon_omega")
                    elseif(user_input%bem_type.eq."standard") then
                            if(user_input%epsilon_omega.ne."general") call tdplas_couple_error("bem_type", "epsilon_omega")
                    endif
            end subroutine



        end module
