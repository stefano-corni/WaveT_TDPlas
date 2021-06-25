      module readio_tdplas_mod
        use tdplas_constants
        use user_input_type
        use read_inputfile_tdplas
        use user_input_check
        use init
        use check_global
        use write_header_out_tdplas
        use user_input_and_flags_dictionary
#ifdef MPI
            use mpi
#endif


        implicit none

        save


        private

        public readio_tdplas

        contains



        subroutine readio_tdplas(calculation_exe)

            character(flg)  :: calculation_exe
            type(tdplas_user_input)  :: user_input
            integer :: nthr

#ifdef OMP
                nthr=omp_get_max_threads()
#endif
#ifndef OMP
                nthr=1
#endif
            call set_default_input_tdplas(user_input)
            call dict_tdplas_init

            call read_input_tdplas(calculation_exe, user_input)
            call check_tdplas_user_input(user_input)
            call init_tdplas(calculation_exe, nthr, user_input)
            call check_global_var
            call write_out(calculation_exe)

        end subroutine



        subroutine check_tdplas_user_input(user_input)
            type(tdplas_user_input) ::  user_input

            call check_surface(user_input)
            call check_spheres_or_spheroids(user_input)
            call check_eps(user_input)

        end subroutine




 end module






