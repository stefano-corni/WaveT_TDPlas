      program tdplas
      use readio_medium
      use BEM_medium
      implicit none
      integer :: st,current,rate
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      call read_medium_tdplas
      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)
!
!     diagonalise matrix
      call do_BEM_prop
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      stop
      end
