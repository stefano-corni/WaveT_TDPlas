      program main_eps
      use tdplas
      implicit none
      integer :: st,current,rate
      integer :: i 
      real(dbl),allocatable :: omega_list(:)

!     read in the input parameters
      call system_clock(st,rate)
      call read_medium_eps

!     printing eps function in file
      open(1,file='eps.inp')
      open(2,file='real_imag_eps.inp')
      write(1,*) n_omega
      write(2,*) n_omega
      allocate(omega_list(n_omega))
      do i=1,n_omega
       omega_list(i)=(omega_end-omega_ini)/(n_omega-1)*(i-1)+omega_ini
       omega(1) = omega_list(i)
       select case( Feps )
       case('deb')
        ! debye eps
        call do_eps_deb
       case('drl')
        ! drude-lorentz eps
        call do_eps_drl
       case('gen')
        ! for now gold case 
        ! extra case should be place here selecting possible material
        eps = eps_gold(omega(1))
       end select
       write(1,*) omega(1), eps
       write(2,*) omega(1), real(eps,dbl), dimag(eps)
      enddo
      close(1)
      close(2)
      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)


      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      deallocate(omega_list)


      end program main_eps


