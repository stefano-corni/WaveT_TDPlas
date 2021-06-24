!> The module reads the input of the medium.
!> It initialize all input namelist variables to default values (when default exists)
!> Checks on typos on keys values are done
!> Internal variables are initialized and stored in global_tdplas or in the different modules
!> WARNING: internal initialization order is important as in global_tdplas checks are performed on calculation consistency:
!> 1 - system_init, medium_init, eps_init,_surface_init (any order)
!> 2 - out_init, print_charges_init, propagation_init/freq_calculation_init (any order)


      module readio_medium
        use constants
        use input_keyword_type 

        use readio_input_namelist 

        use tdplas_readio_dictionary


        use global_quantum, only: quantum_Ffld   !serve per i test



        implicit none

        save


        private

        public read_namelist

        contains



        subroutine read_namelist(calculation_exe, user_input)

            character(flg)  :: calculation_exe
            type(tdplas_user_input), intent(out) :: user_input

            call dict_tdplas_init                                          !init dictionary to convert from user friendly values to internal


            select case(calculation_exe)
            !set calculation dependent default
                case("propagation")
                    user_input%propagation_type = "charge-ief"
                    call read_propagation_input(user_input)
                case("tdplas")
                    user_input%create_cavity    = "built"
                    call read_tdplas_input(user_input)
                case("frequency")
                    call read_freq_input(user_input)
                case("epsilon")
                    user_input%n_omega = 5000
                    user_input%omega_ini = 0.01
                    user_input%omega_end = 0.35
                    call read_epsilon_input(user_input)


            end select

            call set_default_test(user_input)  !WARNING: modifies things without check because test!

        end subroutine read_namelist

        subroutine read_propagation_input(user_input)
       
            type(tdplas_user_input) :: user_input


            call read_system_nml(user_input)
            call read_medium_nml(user_input)
            call read_surface_nml(user_input)
            call read_eps_function_nml(user_input)
            call read_propagate_nml(user_input)
            call read_print_charges_nml(user_input)
            call read_out_nml(user_input)

        end subroutine

        subroutine read_tdplas_input(user_input)
       
            type(tdplas_user_input) :: user_input


            call read_system_nml(user_input)
            call read_medium_nml(user_input)
            call read_surface_nml(user_input)
            call read_eps_function_nml(user_input)
            call read_print_charges_nml(user_input)
            call read_out_nml(user_input)

        end subroutine

        subroutine read_epsilon_input(user_input)

            type(tdplas_user_input) :: user_input

            call read_system_nml(user_input)
            call read_eps_function_nml(user_input)
            call read_out_nml(user_input)

        end subroutine

        subroutine read_freq_input(user_input)

            type(tdplas_user_input) :: user_input


            call read_system_nml(user_input)
            call read_medium_nml(user_input)
            call read_surface_nml(user_input)
            call read_eps_function_nml(user_input)
            call read_out_nml(user_input)

        end subroutine


        subroutine set_default_test(user_input)
            type(tdplas_user_input) :: user_input

            character (flg), dimension(2) :: test_type_loc    = [Character(len=20) :: "n-l", "s-l"]
            character (flg), dimension(3) :: test_type_nonloc    = [Character(len=20) :: "n-r", "s-r", "qmt"]

        !set local_field. No warning for inconsistency because test!
            if (any(user_input%test_type.eq.test_type_loc)) then
                user_input%local_field = "local"
                quantum_Ffld = 'snd'   ! se non è stato chiamatol'init di interface_qmcode questa var non esiste
            elseif(any(user_input%test_type.eq.test_type_nonloc)) then
                user_input%local_field = "non"
            end if
        end subroutine


        



 end module



