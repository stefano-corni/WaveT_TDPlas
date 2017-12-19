      program tdcis
      use readio_medium
      use BEM_medium
      implicit none
      integer :: st,current,rate
      !integer(4) :: n_omega,i
      integer :: i 
      real(8),allocatable :: omega_list(:)
      !real(8) :: omega_ini,omega_end
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      call read_medium_freq
      !read(5,*) n_omega,omega_ini,omega_end
      allocate(omega_list(n_omega))
      do i=1,n_omega
       omega_list(i)=(omega_end-omega_ini)/(n_omega-1)*(i-1)+omega_ini
       !write(6,*) 'omega',omega_list(i)
      enddo
      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)
!
!     diagonalise matrix
      call do_BEM_freq(omega_list,n_omega) 
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      deallocate(omega_list)

      stop

      end program tdcis
