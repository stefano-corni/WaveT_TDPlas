      program tdcis
      use tdplas
      implicit none
      integer :: st,current,rate

!      type(tdplas_user_input) user_input
      character(flg) :: calculation = "frequency" 
      call system_clock(st,rate)

!      call read_input_tdplas(calculation, user_input)
!      call check_input(user_input)
!      call init_tdplas(calculation, user_input)
!      call check_global_var
!      call write_out(calculation)

      call readio_tdplas(calculation)

      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)

!!     diagonalise matrix
      call do_BEM_freq 
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
         

      stop

      end program tdcis
