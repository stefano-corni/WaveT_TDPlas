      program main_eps
      use tdplas
      implicit none
      integer :: st,current,rate
      integer :: i 
      complex(cmp) :: eps

!     read in the input parameters
!      call system_clock(st,rate)
!      call read_medium_eps
!
!      open(1,file="eps.inp", status="new")
!      open(2,file="real_imag_eps.inp", status="new")
!      write(1,*) inout_eps_n_omega
!      write(2,*) inout_eps_n_omega
!      do i=1,inout_eps_n_omega
!          select case( global_eps_Feps )
!              case('deb')
!              ! debye eps
!                  eps = debye_do_eps(inout_eps_omegas_i(i))
!              case('drl')
!              ! drude-lorentz eps
!                  eps = drudel_do_eps(inout_eps_omegas_i(i))
!              case('print_gen')
!              ! for now gold case
!              ! extra case should be place here selecting possible material
!                  eps = eps_gold(omega(1))
!          end select
!          write(1,*) omega(1), eps
!          write(2,*) omega(1), real(eps,dbl), dimag(eps)
!      enddo
!      close(1)
!      close(2)
!      write(6,'("Done reading input, took", &
!            F10.3,"s")') real(current-st)/real(rate)
!
!
!      call system_clock(current)
!      write(6,'("Done , total elapsed time", &
!            F10.3,"s")') real(current-st)/real(rate)
!!         
!

      end program main_eps


