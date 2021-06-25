      module init
        use constants
        use input_keyword_type 

        use dictionary_entry
        use tdplas_readio_dictionary

        use global_tdplas

        implicit none

        save


        private
        public init_and_check_tdplas 

        contains




    subroutine init_and_check_tdplas(calculation_exe, user_input)


        type(tdplas_user_input) ::  user_input
        character(flg) :: calculation_exe



                !init are in global_tdplas take as argument user input values converted from user friendly to internal
                !all variables have default values so all init can be called


                call global_system_init(calculation_exe,                                                    &
                                        d_entry_convert_to_internal(dict_Ftest,   user_input%test_type),    &
                                        d_entry_convert_to_internal(dict_Fdeb,    user_input%debug_type),   &
                                        d_entry_convert_to_internal(dict_Fwrite,  user_input%out_level))




                
                call global_eps_init(d_entry_convert_to_internal(dict_Feps, user_input%epsilon_omega),               &
                                     user_input%eps_d,       &
                                     user_input%tau_deb,     &
                                     user_input%eps_0,       &
                                     user_input%eps_A,       &
                                     user_input%eps_gm,      &
                                     user_input%eps_w0,      &
                                     user_input%f_vel,       &
                                     user_input%n_omega,     &
                                     user_input%omega_ini,   &
                                     user_input%omega_end)




                call global_medium_init(d_entry_convert_to_internal(dict_Fmdm,     user_input%medium_type),   &
                                        d_entry_convert_to_internal(dict_Finit,    user_input%medium_init0),  &
                                        d_entry_convert_to_internal(dict_Fpol,     user_input%medium_pol),    &
                                        d_entry_convert_to_internal(dict_Fbem,     user_input%bem_type),      &
                                        d_entry_convert_to_internal(dict_FinitBEM, user_input%bem_read_write),&
                                        d_entry_convert_to_internal(dict_Floc,     user_input%local_field))




                call global_surface_init(d_entry_convert_to_internal(dict_Fsurf, user_input%input_surface),    &
                                  d_entry_convert_to_internal(dict_Fcav,         user_input%create_cavity),    &
                                  d_entry_convert_to_internal(dict_Finv,         user_input%inversion),        &
                                  d_entry_convert_to_internal(dict_Fshape,       user_input%spheres_shape),    &
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


                call global_out_init(d_entry_convert_to_internal(dict_Fgamess,   user_input%gamess),            &
                              d_entry_convert_to_internal(dict_Fprint_lf_matrix, user_input%print_lf_matrix))




                call global_propagation_init(d_entry_convert_to_internal(dict_Fprop, user_input%propagation_type),     &
                                          user_input%interaction_stride,                                                   &
                                          d_entry_convert_to_internal(dict_Fint,         user_input%interaction_type),     &
                                          d_entry_convert_to_internal(dict_Finit_int,    user_input%interaction_init),     &
                                          d_entry_convert_to_internal(dict_Fmdm_relax,   user_input%medium_relax),         &
                                          d_entry_convert_to_internal(dict_Fmdm_res,     user_input%medium_restart))


                call global_print_charges_init(d_entry_convert_to_internal(dict_Fmop, user_input%charge_mopac),        &
                                               user_input%n_print_charges)





                if (((sph_surf_nsph.gt.ntst).or.(pedra_surf_n_tessere.gt.ntst)).and.(quantum_nthr.gt.1)) then
                    call global_calculation_init('omp')
                else
                    call global_calculation_init('non')
                endif

        end subroutine



 end module init



