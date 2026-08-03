      program main_tdplas
      use tdplas
      use BEM_medium
      use readio_tdplas_mod
      implicit none
      integer :: st,current,rate
!
!     read in the input parameter for the present evolution
!      call system_clock(st,rate)
!
!      type(tdplas_user_input) user_input

      character(flg) :: calculation = "tdplas"


      call system_clock(st,rate)

      call readio_tdplas(calculation) 

      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)


!     Silvio 02/08/2026 added the following, cannot use do_BEM_prop if
!                       do_BEM is not called
      call do_BEM
!     diagonalise matrix
      call do_BEM_prop
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
         
      stop
      end
