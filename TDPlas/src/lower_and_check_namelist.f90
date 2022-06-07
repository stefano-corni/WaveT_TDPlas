module lower_and_check_namelist

        use dictionary_entry
        use user_input_and_flags_dictionary

        implicit none

        save


        private

        public lower_and_check_allowed_values_system_nml,  &
               lower_and_check_allowed_values_surface_nml, &
               lower_and_check_allowed_values_medium_nml,  &
               lower_and_check_allowed_values_propagate_nml, &
               lower_and_check_allowed_values_eps_nml, &
               lower_and_check_allowed_values_out_nml, &
               lower_and_check_allowed_values_quantum_coupling_nml, &
               lower_and_check_allowed_values_ext_pert_nml

        contains



    subroutine lower_and_check_allowed_values_system_nml(test_type, &
                                                        debug_type, &
                                                        out_level)

        character(flg)  :: test_type
        character(flg)  :: debug_type
        character(flg)  :: out_level

        call To_lower(test_type)
        call To_lower(debug_type)
        call To_lower(out_level)

        call d_entry_check_allowed_val(dict_Ftest, test_type)
        call d_entry_check_allowed_val(dict_Fdeb, debug_type)
        call d_entry_check_allowed_val(dict_Fwrite, out_level)

    end subroutine

    subroutine lower_and_check_allowed_values_surface_nml(surface_type, &
                                                          input_mesh, &
                                                          object_shape, &
                                                          inversion,    &
                                                          find_spheres)

        character(flg)   :: surface_type
        character(flg)   :: object_shape
        character(flg)   :: input_mesh
        character(flg)   :: inversion
        character(flg)   :: find_spheres

        call To_lower(surface_type)
        call To_lower(input_mesh)
        call To_lower(object_shape)
        call To_lower(inversion)
        call To_lower(find_spheres)

        call d_entry_check_allowed_val(dict_Fsurf, surface_type)
        call d_entry_check_allowed_val(dict_Fcav, input_mesh)
        call d_entry_check_allowed_val(dict_Fshape, object_shape)
        call d_entry_check_allowed_val(dict_Finv, inversion)
        call d_entry_check_allowed_val(dict_Ffind, find_spheres)
    end subroutine

    subroutine lower_and_check_allowed_values_medium_nml(medium_type,  &
                                                         medium_init0, &
                                                         medium_pol,   &
                                                         bem_type,     &
                                                         bem_read_write,&
                                                         local_field)
                                                         
        character(flg)  :: medium_type
        character(flg)  :: medium_init0
        character(flg)  :: medium_pol
        character(flg)  :: bem_type
        character(flg)  :: bem_read_write
        character(flg)  :: local_field

        call To_lower(medium_type)
        call To_lower(medium_init0)
        call To_lower(medium_pol)
        call To_lower(bem_type)
        call To_lower(bem_read_write)
        call To_lower(local_field)

        call d_entry_check_allowed_val(dict_Fmdm, medium_type)
        call d_entry_check_allowed_val(dict_Finit, medium_init0)
        call d_entry_check_allowed_val(dict_Fpol, medium_pol)
        call d_entry_check_allowed_val(dict_Fbem, bem_type)
        call d_entry_check_allowed_val(dict_read_write, bem_read_write)
        call d_entry_check_allowed_val(dict_Floc, local_field)

    end subroutine

    subroutine lower_and_check_allowed_values_eps_nml(epsilon_omega,propagation_pole)

        character(flg)    :: epsilon_omega,propagation_pole

        call To_lower(epsilon_omega)
        call To_lower(propagation_pole)

        call d_entry_check_allowed_val(dict_Feps, epsilon_omega)
        call d_entry_check_allowed_val(dict_typ_prop, propagation_pole)
    end subroutine


    subroutine lower_and_check_allowed_values_propagate_nml(propagation_software, &
                                                            propagation_type,     &
                                                            interaction_type, &
                                                            interaction_init, &
                                                            medium_relax,     &
                                                            medium_restart)


        character(flg)   :: propagation_software
        character(flg)   :: propagation_type


        character(flg)   :: interaction_type
        character(flg)   :: interaction_init
        character(flg)   :: medium_relax
        character(flg)   :: medium_restart


        call To_lower(propagation_software)
        call To_lower(propagation_type)
        call To_lower(interaction_type)
        call To_lower(interaction_init)
        call To_lower(medium_relax)
        call To_lower(medium_restart)


        call d_entry_check_allowed_val(dict_Fsoft, propagation_software)
        call d_entry_check_allowed_val(dict_Fprop, propagation_type)
        call d_entry_check_allowed_val(dict_Fint, interaction_type)
        call d_entry_check_allowed_val(dict_Finit_int, interaction_init)
        call d_entry_check_allowed_val(dict_Fmdm_relax, medium_relax)
        call d_entry_check_allowed_val(dict_Fmdm_res, medium_restart)

    end subroutine

    subroutine lower_and_check_allowed_values_out_nml(gamess, &
                                                      print_lf_matrix)
        character(flg)   :: gamess
        character(flg)   :: print_lf_matrix


        call To_lower(gamess)
        call To_lower(print_lf_matrix)

        call d_entry_check_allowed_val(dict_Fgamess, gamess)
        call d_entry_check_allowed_val(dict_Fprint_lf_matrix, print_lf_matrix)
    end subroutine

    subroutine lower_and_check_allowed_values_quantum_coupling_nml(charge_mopac, &
                                                quantum_calculation)
        character(flg)  :: charge_mopac
        character(flg)  :: quantum_calculation
        call To_lower(charge_mopac)
        call To_lower(quantum_calculation)
        call d_entry_check_allowed_val(dict_Fmop, charge_mopac)
        call d_entry_check_allowed_val(dict_Fqcp, quantum_calculation)

    end subroutine

    subroutine lower_and_check_allowed_values_ext_pert_nml(pert_type,eet)
        character(flg)  :: pert_type    
        character(flg)  :: eet
        call To_lower(pert_type)
        call To_lower(eet)
        call d_entry_check_allowed_val(dict_Ftyp, pert_type)
        call d_entry_check_allowed_val(dict_Feet, eet)

    end subroutine

   !MR http://rosettacode.org/wiki/String_case#Fortran
   subroutine To_lower(str)
     character(*), intent(in out) :: str
     integer :: i
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do
   end subroutine To_Lower


end module lower_and_check_namelist



