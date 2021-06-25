!> The module reads the input of the medium.
!> It initialize all input namelist variables to default values (when default exists)
!> Checks on typos on keys values are done
!> Internal variables are initialized and stored in global_tdplas or in the different modules

      module read_tdplas
        use tdplas_constants
        use user_input_type

        use user_input_read_namelist

        use global_quantum, only: quantum_Ffld   !serve per i test

        implicit none

        save


        private

        public read_input_tdplas_internal, read_input_tdplas_for_propagation

        contains

        subroutine read_input_tdplas_internal(calculation_exe, user_input)
            character(flg)  :: calculation_exe
            type(tdplas_user_input), intent(out) :: user_input

            call read_system_nml(user_input)
            call read_out_nml(user_input)

            select case(calculation_exe)
                case("tdplas")
                    user_input%create_cavity    = "built"
                    call read_medium_nml(user_input)
                    call read_surface_nml(user_input)
                    call read_eps_function_nml(user_input)
                    call read_print_charges_nml(user_input)
                    call read_quantum_coupling_nml(user_input)

                case("frequency")
                    call read_medium_nml(user_input)
                    call read_surface_nml(user_input)
                    call read_eps_function_nml(user_input)
                    call read_do_eps_nml(user_input)

                case("epsilon")
                    user_input%n_omega = 5000
                    user_input%omega_ini = 0.01
                    user_input%omega_end = 0.35

                    call read_eps_function_nml(user_input)
                    call read_do_eps_nml(user_input)

            end select

            call set_default_test(user_input)  !WARNING: modifies things without check because test!
        end subroutine

        subroutine read_input_tdplas_for_propagation(user_input)
            type(tdplas_user_input), intent(out) :: user_input


            call read_system_nml(user_input)
            call read_medium_nml(user_input)
            call read_surface_nml(user_input)
            call read_eps_function_nml(user_input)
            call read_propagate_nml(user_input)
            call read_out_nml(user_input)


            call set_default_test(user_input)  !WARNING: modifies things without check because test!
        end subroutine






 !--------------------------set default test, modifies things without checking----------------!

       subroutine set_default_test(user_input)
            type(tdplas_user_input) :: user_input

            character (flg), dimension(2) :: test_type_loc    = [Character(len=20) :: "n-l", "s-l"]
            character (flg), dimension(3) :: test_type_nonloc    = [Character(len=20) :: "n-r", "s-r", "qmt"]

        !set local_field. No warning for inconsistency because test!
            if (any(user_input%test_type.eq.test_type_loc)) then
                user_input%local_field = "local"
                !quantum_Ffld = 'snd'   ! se non è stato chiamatol'init di interface_qmcode questa var non esiste
            elseif(any(user_input%test_type.eq.test_type_nonloc)) then
                user_input%local_field = "non"
            end if
        end subroutine


 end module



