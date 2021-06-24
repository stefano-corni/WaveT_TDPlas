      module auxiliary_functions
          use tdplas_constants 


          contains

        subroutine mpi_error(error_line1, error_line2, error_line3)
        character(*), intent(in)  :: error_line1, error_line2, error_line3
        character(200)             :: string_tot

        string_tot = trim(error_line1) // " "  // trim(error_line2) // " "  // trim(error_line3)
        write(6,*) string_tot
#ifdef MPI
               call mpi_finalize(tp_ierr_mpi)
#endif
            stop
        end subroutine


        subroutine ocpy_error(key)
          character(*), intent(in)  :: key
          write(6,*) key," value is wrong, not all tdplas values are allowed in ocpy calculations"
#ifdef MPI
               call mpi_finalize(tp_ierr_mpi)
#endif
          stop

        end subroutine


        subroutine wavet_error(key)
          character(*), intent(in)  :: key
          write(6,*) "ERROR: ", key," value is wrong, not all tdplas values are allowed in wavet calculations"
#ifdef MPI
               call mpi_finalize(tp_ierr_mpi)
#endif
          stop

        end subroutine


        subroutine tdplas_couple_error(key1,key2)
          character(*), intent(in)  :: key1, key2
          write(6,*) "ERROR: choiches for ", key1," value and ", key2, " value are not allowed together"
#ifdef MPI
               call mpi_finalize(tp_ierr_mpi)
#endif
          stop

        end subroutine




      end module

