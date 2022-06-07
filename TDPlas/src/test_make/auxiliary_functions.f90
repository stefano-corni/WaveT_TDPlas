      module auxiliary_functions
          use constants


          contains

        subroutine mpi_error(error_line1, error_line2, error_line3)
        character(*), intent(in)  :: error_line1, error_line2, error_line3
        character(200)             :: string_tot

        string_tot = trim(error_line1) // " "  // trim(error_line2) // " "  // trim(error_line3)
        write(6,*), string_tot 
#ifdef MPI
               call mpi_finalize(ierr_mpi)
#endif
            stop
        end subroutine




!        logical function optimized_matmul_tessere(n_tessere, size_thr, nthreads) result(logical_output)
!
!            integer(i4b) n_tessere, size_thr, nthreads

!            if (n_tessere.gt.size_thr.and.nthreads.gt.1) then
!                logical_output = .true
!            else
!                logical_output = .false
!                endif

!       end subroutine






      end module

