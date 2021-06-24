      program main_tdplas
      use tdplas
      implicit none
      integer :: st,current,rate
!
!     read in the input parameter for the present evolution
!      call system_clock(st,rate)
!
      type(tdplas_user_input) user_input
!
      character(flg) :: calculation = "tdplas"
      call read_namelist(calculation, user_input)
      call init_and_check_tdplas(calculation, user_input)
      call write_out(calculation)
! 
!      call system_clock(current)
!      write(6,'("Done reading input, took", &
!            F10.3,"s")') real(current-st)/real(rate)
!
!!
!!     diagonalise matrix
!      call do_BEM_prop
!      call system_clock(current)
!      write(6,'("Done , total elapsed time", &
!            F10.3,"s")') real(current-st)/real(rate)
!!         
!      stop
      end
