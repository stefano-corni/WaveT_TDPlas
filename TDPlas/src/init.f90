      module init
        use tdplas_constants
        use user_input_type

        use dictionary_entry
        use user_input_and_flags_dictionary

        use global_tdplas
        use pedra_friends
        implicit none

        save


        private
        public init_tdplas, init_parallel

        contains



        subroutine init_tdplas(calculation_exe, nthr, user_input)

            type(tdplas_user_input) ::  user_input
            character(flg) :: calculation_exe
            integer :: nthr

            call init_Fglobal(calculation_exe, user_input)
            call init_parallel(nthr)
            !surface
            select case(global_surf_Fsurf)
                case('sphe')
                    call sph_surf_init(d_entry_convert_to_internal(dict_Fshape, user_input%object_shape), &
                                       user_input%spheres_number,       &
                                       user_input%spheroids_number,     &
                                       user_input%sphere_position_x,    &
                                       user_input%sphere_position_y,    &
                                       user_input%sphere_position_z,    &
                                       user_input%sphere_radius,        &
                                       user_input%spheroid_position_x,  &
                                       user_input%spheroid_position_y,  &
                                       user_input%spheroid_position_z,  &
                                       user_input%spheroid_radius,      &
                                       user_input%spheroid_axis_x,      &
                                       user_input%spheroid_axis_y,      &
                                       user_input%spheroid_axis_z)

                case('cav')
                    call pedra_surf_init(d_entry_convert_to_internal(dict_Fcav, user_input%input_mesh),    &
                                    d_entry_convert_to_internal(dict_Finv, user_input%inversion),        &
                                    d_entry_convert_to_internal(dict_Ffind, user_input%find_spheres),        &
                                    user_input%particles_number,  &
                                    user_input%spheres_number,    &
                                    user_input%sphere_position_x, &
                                    user_input%sphere_position_y, &
                                    user_input%sphere_position_z, &
                                    user_input%sphere_radius, &     
                                    global_medium_Fmdm)


                end select

           !epsilon

                select case(global_eps_Feps)
                   case("drl")
                       call drudel_eps_init(user_input%eps_A, user_input%eps_gm, user_input%eps_w0, user_input%f_vel)
                   case("deb")
                       call debye_eps_init(user_input%eps_d, user_input%tau_deb, user_input%eps_0)
                   case("gen")
                       if(user_input%eps_d.eq.-1.) user_input%eps_d = 1.
                       if(user_input%eps_0.eq.-1) user_input%eps_0 = 1000.
                       call readf_eps_init(user_input%eps_d, user_input%eps_0)
                end select


                if((calculation_exe.eq."frequency").or.(calculation_exe.eq."epsilon")) then
                    call dielectric_func_init(user_input%n_omega, user_input%omega_ini, user_input%omega_end)
                endif
        end subroutine

        subroutine init_parallel(nthr)
                integer :: nthr
                if ((sph_surf_nsph.gt.ntst).or.(pedra_surf_n_tessere.gt.ntst)) then
                    call global_calculation_init('omp', nthr)
                else
                    call global_calculation_init('non', nthr)
                endif

        end subroutine




        subroutine init_Fglobal(calculation_exe, user_input)


            type(tdplas_user_input) ::  user_input
            character(flg) :: calculation_exe

                !init are in global_tdplas take as argument user input values converted from user friendly to internal
                !all variables have default values so all init can be called
        

            call global_system_init(calculation_exe,                                                    &
                                    d_entry_convert_to_internal(dict_Ftest,   user_input%test_type),    &
                                    d_entry_convert_to_internal(dict_Fdeb,    user_input%debug_type),   &
                                    d_entry_convert_to_internal(dict_Fwrite,  user_input%out_level))




            call global_eps_init(d_entry_convert_to_internal(dict_Feps, user_input%epsilon_omega),       &
                                 d_entry_convert_to_internal(dict_typ_prop, user_input%propagation_pole))




            call global_medium_init(d_entry_convert_to_internal(dict_Fmdm,      user_input%medium_type),   &
                                    d_entry_convert_to_internal(dict_Finit,     user_input%medium_init0),  &
                                    d_entry_convert_to_internal(dict_Fpol,      user_input%medium_pol),    &
                                    d_entry_convert_to_internal(dict_Fbem,      user_input%bem_type),      &
                                    d_entry_convert_to_internal(dict_Floc,      user_input%local_field),   &
                                    d_entry_convert_to_internal(dict_read_write,user_input%bem_read_write),&
                                    d_entry_convert_to_internal(dict_Fnorm,     user_input%normalization), &
                                    d_entry_convert_to_internal(dict_bem_sym,   user_input%bem_symmetric))




            call global_surface_init(d_entry_convert_to_internal(dict_Fsurf, user_input%surface_type))

            call global_out_init(d_entry_convert_to_internal(dict_Fgamess,   user_input%gamess),            &
                                 d_entry_convert_to_internal(dict_Fprint_lf_matrix, user_input%print_lf_matrix))



            call global_propagation_init(d_entry_convert_to_internal(dict_Fprop, user_input%propagation_type),            &
                                         user_input%interaction_stride,                                                   &
                                         user_input%mix_coef,                                                             &
                                         user_input%threshold,                                                            &
                                         user_input%max_cycles,                                                           &
                                         d_entry_convert_to_internal(dict_Fint,         user_input%interaction_type),     &
                                         d_entry_convert_to_internal(dict_Finit_int,    user_input%interaction_init),     &
                                         d_entry_convert_to_internal(dict_Fmdm_relax,   user_input%medium_relax),         &
                                         d_entry_convert_to_internal(dict_Fmdm_res,     user_input%medium_restart))


            call global_quantum_coupling_init(user_input%n_modes, user_input%qm_modes, user_input%n_print_charges, &
                                                                                       user_input%print_charges,   &
                                                 d_entry_convert_to_internal(dict_Fmop, user_input%charge_mopac),  &
                                                 d_entry_convert_to_internal(dict_Fqcp, user_input%quantum_calculation))
            call global_ext_pert_init(user_input%pert_type, user_input%direction, user_input%eet, user_input%gamma_vacuum, &
                                    user_input%n_excited,user_input%nstate,user_input%quantum_mol_cc, &
                                    user_input%dipole_vacuum, user_input%pl_type,user_input%pl_omega_abs,user_input%pl_omega_emi,&
                                    user_input%print_field,user_input%print_surface_charges)


        end subroutine






 end module init



