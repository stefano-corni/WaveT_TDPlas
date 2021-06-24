      module user_input_read_namelist

        use tdplas_constants
        use user_input_type
        use lower_and_check_namelist


        implicit none

        save



        private
        public read_eps_function_nml, read_medium_nml, read_propagate_nml, read_surface_nml, read_out_nml, &
               read_system_nml, read_do_eps_nml, read_quantum_coupling_nml, &
               read_ext_pert_nml

        contains

        subroutine read_eps_function_nml(user_input)

            type(tdplas_user_input) :: user_input

            real(dbl) :: eps_0, eps_d, eps_A, eps_gm, eps_w0, f_vel, tau_deb
            character(flg) :: epsilon_omega, propagation_pole

            namelist /eps_function/epsilon_omega,        &
                                   eps_0,                &
                                   eps_d,                &
                                   eps_A,                &
                                   eps_gm,               &
                                   eps_w0,               &
                                   f_vel,                &
                                   tau_deb,              &
                                   propagation_pole



            epsilon_omega  =    user_input%epsilon_omega
            eps_0          =    user_input%eps_0
            eps_d          =    user_input%eps_d
            eps_A          =    user_input%eps_A
            eps_gm         =    user_input%eps_gm
            eps_w0         =    user_input%eps_w0
            f_vel          =    user_input%f_vel
            tau_deb        =    user_input%tau_deb
            propagation_pole = user_input%propagation_pole


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=eps_function)
            call lower_and_check_allowed_values_eps_nml(epsilon_omega,propagation_pole)


             user_input%epsilon_omega        = epsilon_omega
             user_input%eps_0                = eps_0
             user_input%eps_d                = eps_d
             user_input%eps_A                = eps_A
             user_input%eps_gm               = eps_gm
             user_input%eps_w0               = eps_w0
             user_input%f_vel                = f_vel
             user_input%tau_deb              = tau_deb
             user_input%propagation_pole     = propagation_pole


        end subroutine


        subroutine read_do_eps_nml(user_input)

            type(tdplas_user_input) :: user_input

            real(dbl) :: n_omega, omega_ini, omega_end
            character(flg) :: epsilon_omega


            namelist /do_eps/n_omega,              &
                             omega_ini,            &
                             omega_end

            n_omega     =   user_input%n_omega
            omega_ini   =   user_input%omega_ini
            omega_end   =   user_input%omega_end



            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=do_eps)
             user_input%n_omega              = n_omega
             user_input%omega_ini            = omega_ini
             user_input%omega_end            = omega_end

        end subroutine

        subroutine read_ext_pert_nml(user_input)

            type(tdplas_user_input) :: user_input

            real(dbl) :: direction(3)
            character(flg)    :: pert_type,eet,pl_type
            character(flg)    :: print_field, print_surface_charges
            real(dbl)         :: gamma_vacuum
            real(dbl)         :: dipole_vacuum(3)
            integer(i4b)      :: n_excited
            integer(i4b)      :: nstate
            real(dbl)         :: quantum_mol_cc(3),pl_omega_abs,pl_omega_emi


            namelist /ext_pert/pert_type,              &
                             direction,gamma_vacuum,   &
                             eet,n_excited, nstate,    &
                             quantum_mol_cc, &
                             dipole_vacuum, pl_type,   &
                             pl_omega_abs,pl_omega_emi,&
                             print_field,print_surface_charges

            pert_type      =   user_input%pert_type
            direction      =   user_input%direction
            eet            =   user_input%eet
            gamma_vacuum   =   user_input%gamma_vacuum
            dipole_vacuum  =   user_input%dipole_vacuum
            n_excited      =   user_input%n_excited
            nstate         =   user_input%nstate
            quantum_mol_cc =   user_input%quantum_mol_cc
            pl_type        =   user_input%pl_type
            pl_omega_abs   =   user_input%pl_omega_abs
            pl_omega_emi   =   user_input%pl_omega_emi
            print_field    =   user_input%print_field
            print_surface_charges  =   user_input%print_surface_charges               



            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=ext_pert)
            call lower_and_check_allowed_values_ext_pert_nml(pert_type,eet)
            user_input%pert_type            = pert_type
            user_input%direction            = direction
            user_input%eet                  = eet        
            user_input%gamma_vacuum         = gamma_vacuum
            user_input%dipole_vacuum        = dipole_vacuum
            user_input%n_excited            = n_excited
            user_input%nstate               = nstate
            user_input%quantum_mol_cc       = quantum_mol_cc
            user_input%pl_omega_abs         = pl_omega_abs
            user_input%pl_omega_emi         = pl_omega_emi
            user_input%pl_type              = pl_type
            user_input%print_field          = print_field
            user_input%print_surface_charges= print_surface_charges 

        end subroutine



        subroutine read_system_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: test_type, debug_type, out_level


            namelist /system/test_type,                  &
                            debug_type,                  &
                            out_level

            test_type = user_input%test_type
            debug_type = user_input%debug_type
            out_level = user_input%out_level



            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=system)
            call lower_and_check_allowed_values_system_nml(test_type, debug_type, out_level)


             user_input%test_type        = test_type
             user_input%debug_type       = debug_type
             user_input%out_level        = out_level

        end subroutine


        subroutine read_out_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: gamess, print_lf_matrix



             namelist /out_matrix/gamess,                 &
                                  print_lf_matrix



            gamess                =       user_input%gamess
            print_lf_matrix       =       user_input%print_lf_matrix


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=out_matrix)
            call lower_and_check_allowed_values_out_nml(gamess, print_lf_matrix)


              user_input%gamess            =  gamess
              user_input%print_lf_matrix   =  print_lf_matrix


        end subroutine


        subroutine read_propagate_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: propagation_software,propagation_type,&
                              interaction_init,interaction_type,medium_relax,medium_restart
            integer(i4b)     :: interaction_stride, max_cycles
            real(dbl)        :: mix_coef, threshold

            namelist /propagate/propagation_software,    &
                                propagation_type,        &
                                interaction_stride,      &
                                interaction_init,        &
                                mix_coef,                &
                                max_cycles,              &
                                threshold,               &
                                interaction_type,        &
                                medium_relax,            &
                                medium_restart


            propagation_software  =   user_input%propagation_software
            propagation_type      =   user_input%propagation_type
            interaction_stride    =   user_input%interaction_stride
            mix_coef              =   user_input%mix_coef
            threshold             =   user_input%threshold
            max_cycles            =   user_input% max_cycles
            interaction_init      =   user_input%interaction_init
            interaction_type      =   user_input%interaction_type
            medium_relax          =   user_input%medium_relax
            medium_restart        =   user_input%medium_restart


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=propagate)
            call lower_and_check_allowed_values_propagate_nml(propagation_software,   &
                                                              propagation_type,       &
                                                              interaction_type,       &
                                                              interaction_init,       &
                                                              medium_relax,           &
                                                              medium_restart )




            user_input%propagation_software   = propagation_software
            user_input%propagation_type       = propagation_type
            user_input%interaction_stride     = interaction_stride
            user_input%mix_coef               = mix_coef
            user_input%threshold              = threshold
            user_input%max_cycles             = max_cycles
            user_input%interaction_init       = interaction_init
            user_input%interaction_type       = interaction_type
            user_input%medium_relax           = medium_relax
            user_input%medium_restart         = medium_restart

        end subroutine


        subroutine read_medium_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: medium_type, medium_init0, medium_pol, bem_type, bem_read_write 
            character(flg) :: local_field, normalization, bem_symmetric

            namelist /medium/medium_type,                 &
                             medium_init0,                &
                             medium_pol,                  &
                             bem_type,                    &
                             bem_read_write,              &
                             local_field,                 &
                             normalization,               &
                             bem_symmetric


            medium_type     =  user_input%medium_type
            medium_init0    =  user_input%medium_init0
            medium_pol      =  user_input%medium_pol
            bem_type        =  user_input%bem_type
            bem_read_write  =  user_input%bem_read_write
            local_field     =  user_input%local_field
            normalization   =  user_input%normalization
            bem_symmetric   =  user_input%bem_symmetric

            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=medium)
            call lower_and_check_allowed_values_medium_nml(medium_type, medium_init0, medium_pol, bem_type, &
                                                      bem_read_write, local_field)

            user_input%medium_type    =     medium_type
            user_input%medium_init0   =     medium_init0
            user_input%medium_pol     =     medium_pol
            user_input%bem_type       =     bem_type
            user_input%bem_read_write =     bem_read_write
            user_input%local_field    =     local_field
            user_input%normalization  =     normalization
            user_input%bem_symmetric  =     bem_symmetric

        end subroutine



        subroutine read_quantum_coupling_nml(user_input)

            type(tdplas_user_input) :: user_input
            integer(i4b)   :: n_modes
            integer(i4b)   :: qm_modes(nmmax)   
            integer(i4b)   :: print_charges(nmmax)   
            character(flg) :: quantum_calculation
            character(flg) :: charge_mopac
            integer(i4b)   :: n_print_charges

            namelist /quantum_coupling/n_modes, qm_modes, & 
                                    quantum_calculation, &
                                    n_print_charges,   &
                                    print_charges,   &
                                    charge_mopac

            n_modes    =   user_input%n_modes
            qm_modes   =   user_input%qm_modes
            quantum_calculation      =   user_input%quantum_calculation
            charge_mopac      =   user_input%charge_mopac
            n_print_charges   =   user_input%n_print_charges
            print_charges   =   user_input%print_charges

            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=quantum_coupling)
            call lower_and_check_allowed_values_quantum_coupling_nml(charge_mopac,&
                     quantum_calculation)

            user_input%n_modes    = n_modes  
            user_input%qm_modes   = qm_modes 
            user_input%quantum_calculation  =  quantum_calculation  
            user_input%charge_mopac    = charge_mopac
            user_input%n_print_charges = n_print_charges
            user_input%print_charges = print_charges

        end subroutine



        subroutine read_surface_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: surface_type, input_mesh, object_shape, inversion, find_spheres
            real(dbl)      :: sphere_position_x(nsmax), sphere_position_y(nsmax), sphere_position_z(nsmax), &
                              sphere_radius(nsmax), &
                              spheroid_position_x(nsmax), spheroid_position_y(nsmax), spheroid_position_z(nsmax), &
                              spheroid_radius(nsmax), &
                              spheroid_axis_x(nsmax), spheroid_axis_y(nsmax), spheroid_axis_z(nsmax)
            integer(i4b)   :: spheres_number, spheroids_number, particles_number


            namelist /surface/surface_type,                &
                              input_mesh,                  &
                              particles_number,            &
                              object_shape,                &
                              spheres_number,              &
                              spheroids_number,            &
                              sphere_position_x,           &
                              sphere_position_y,           &
                              sphere_position_z,           &
                              spheroid_axis_x,             &
                              spheroid_axis_y,             &
                              spheroid_axis_z,             &
                              spheroid_position_x,         &
                              spheroid_position_y,         &
                              spheroid_position_z,         &
                              sphere_radius,               &
                              spheroid_radius,             &
                              inversion,                   &
                              find_spheres


            surface_type         =       user_input%surface_type
            input_mesh           =       user_input%input_mesh
            particles_number     =       user_input%particles_number
            object_shape         =       user_input%object_shape
            spheres_number       =       user_input%spheres_number
            spheroids_number     =       user_input%spheroids_number
            sphere_position_x    =       user_input%sphere_position_x
            sphere_position_y    =       user_input%sphere_position_y
            sphere_position_z    =       user_input%sphere_position_z
            spheroid_axis_x      =       user_input%spheroid_axis_x
            spheroid_axis_y      =       user_input%spheroid_axis_y
            spheroid_axis_z      =       user_input%spheroid_axis_z
            spheroid_position_x  =       user_input%spheroid_position_x
            spheroid_position_y  =       user_input%spheroid_position_y
            spheroid_position_z  =       user_input%spheroid_position_z
            sphere_radius        =       user_input%sphere_radius
            spheroid_radius      =       user_input%spheroid_radius
            inversion            =       user_input%inversion
            find_spheres         =       user_input%find_spheres




            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=surface)
            call lower_and_check_allowed_values_surface_nml(surface_type, input_mesh, object_shape, inversion, find_spheres)

            user_input%surface_type             = surface_type
            user_input%input_mesh               = input_mesh
            user_input%object_shape             = object_shape
            user_input%particles_number         = particles_number
            user_input%spheres_number           = spheres_number
            user_input%spheroids_number         = spheroids_number
            user_input%sphere_position_x        = sphere_position_x
            user_input%sphere_position_y        = sphere_position_y
            user_input%sphere_position_z        = sphere_position_z
            user_input%spheroid_axis_x          = spheroid_axis_x
            user_input%spheroid_axis_y          = spheroid_axis_y
            user_input%spheroid_axis_z          = spheroid_axis_z
            user_input%spheroid_position_x      = spheroid_position_x
            user_input%spheroid_position_y      = spheroid_position_y
            user_input%spheroid_position_z      = spheroid_position_z
            user_input%sphere_radius            = sphere_radius
            user_input%spheroid_radius          = spheroid_radius
            user_input%inversion                = inversion
            user_input%find_spheres             = find_spheres

        end subroutine

 end module


