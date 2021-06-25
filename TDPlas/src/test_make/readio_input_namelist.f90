      module readio_input_namelist

        use constants
        use input_keyword_type
        use lower_and_check_namelist


        implicit none

        save



        private
        public read_eps_function_nml, read_medium_nml, read_propagate_nml, read_surface_nml, read_out_nml, &
               read_system_nml, read_print_charges_nml

        contains


        subroutine read_eps_function_nml(user_input)
            
            type(tdplas_user_input) :: user_input
            
            real(dbl) :: eps_0, eps_d, eps_A, eps_gm, eps_w0, f_vel, tau_deb, n_omega, omega_ini, omega_end
            character(flg) :: epsilon_omega
            
            namelist /eps_function/epsilon_omega,        &
                                   eps_0,                &
                                   eps_d,                &
                                   eps_A,                &
                                   eps_gm,               &
                                   eps_w0,               &
                                   f_vel,                &
                                   tau_deb,              &
                                   n_omega,              &
                                   omega_ini,            &
                                   omega_end            
            
            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=eps_function)
            call lower_and_check_allowed_values_eps_nml(epsilon_omega)
            
            
             user_input%epsilon_omega        = epsilon_omega
             user_input%eps_0                = eps_0
             user_input%eps_d                = eps_d
             user_input%eps_A                = eps_A
             user_input%eps_gm               = eps_gm
             user_input%eps_w0               = eps_w0
             user_input%f_vel                = f_vel
             user_input%tau_deb              = tau_deb
             user_input%n_omega              = n_omega
             user_input%omega_ini            = omega_ini
             user_input%omega_end            = omega_end
            
        end subroutine

        
        subroutine read_system_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: test_type, debug_type, out_level 


            namelist /system/test_type,                  &
                            debug_type,                  &
                            out_level

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
                              interaction_init,interaction_type,local_field,medium_relax,medium_restart 
            integer(i4b)     :: interaction_stride

            namelist /propagate/propagation_software,    &
                                propagation_type,        &
                                interaction_stride,      &
                                interaction_init,        &
                                interaction_type,        &
                                local_field,             &
                                medium_relax,            &
                                medium_restart


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=propagate)
            call lower_and_check_allowed_values_propagate_nml(propagation_software,   &
                                                              propagation_type,       &
                                                              interaction_init,       &
                                                              local_field,            &
                                                              medium_relax,           &
                                                              medium_restart )




            user_input%propagation_software   = propagation_software 
            user_input%propagation_type       = propagation_type
            user_input%interaction_stride     = interaction_stride
            user_input%interaction_init       = interaction_init 
            user_input%interaction_type       = interaction_type
            user_input%local_field            = local_field
            user_input%medium_relax           = medium_relax
            user_input%medium_restart         = medium_restart

        end subroutine 


        subroutine read_medium_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: medium_type, medium_init0, medium_pol, bem_type, bem_read_write

            namelist /medium/medium_type,                 &
                             medium_init0,                &
                             medium_pol,                  &
                             bem_type,                    &
                             bem_read_write



            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=medium)
            call lower_and_check_allowed_values_medium_nml(medium_type, medium_init0, medium_pol, bem_type, bem_read_write)

            user_input%medium_type    =     medium_type
            user_input%medium_init0   =     medium_init0
            user_input%medium_pol     =     medium_pol
            user_input%bem_type       =     bem_type
            user_input%bem_read_write =     bem_read_write

        end subroutine 


        subroutine read_print_charges_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: charge_mopac
            integer(i4b)   :: n_print_charges

            namelist /print_charges/n_print_charges,   &
                                    charge_mopac


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=print_charges)
            call lower_and_check_allowed_values_print_charges_nml(charge_mopac)


            user_input%charge_mopac    = charge_mopac
            user_input%n_print_charges = n_print_charges

        end subroutine 



        subroutine read_surface_nml(user_input)

            type(tdplas_user_input) :: user_input

            character(flg) :: input_surface, create_cavity, spheres_shape, inversion 
            real(dbl)      :: sphere_position_x, sphere_position_y, sphere_position_z, sphere_radius, &
                              spheroid_position_x, spheroid_position_y, spheroid_position_z, spheroid_radius, &
                              spheroid_axis_x, spheroid_axis_y, spheroid_axis_z
            integer(i4b)   :: spheres_number, spheroids_number


            namelist /surface/input_surface,               &
                              create_cavity,               &
                              spheres_shape,               &
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
                              inversion


            open(888,file="namelist_tdplas.inp")
            rewind(888)
            read(888, nml=surface)
            call lower_and_check_allowed_values_surface_nml(input_surface, create_cavity, spheres_shape, inversion)



            user_input%input_surface            = input_surface 
            user_input%create_cavity            = create_cavity
            user_input%spheres_shape            = spheres_shape
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

        end subroutine

 end module


