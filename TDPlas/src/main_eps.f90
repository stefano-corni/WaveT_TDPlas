      program tdplas
      use readio_medium
      use BEM_medium
      use constants
      implicit none

      integer(i4b) :: i
      real(dbl) :: omega, omega_start, omega_fin, delta_omega

!     read in the input parameter for the present evolution
      call read_medium_tdplas

!     data
      omega_start = zero
      omega_fin  = 100.0d0/27.211d0
      delta_omega = 0.1d0/27.211d0
      npts = int( ( omega_fin - omega_start ) / delta_omega )

!     printing eps function in file
      open(1,file='eps.inp')
      open(2,file='real_eps.inp')
      write(1,*) npts
      write(2,*) npts
      do i=1, npts
       omega = ( omega_fin - omega_start ) * i/npts + omega_start
       call do_eps_drl(omega)
       write(1,*) omega, eps
       write(2,*) omega, real(eps)
      enddo
      close(1)
      close(2)

      end
