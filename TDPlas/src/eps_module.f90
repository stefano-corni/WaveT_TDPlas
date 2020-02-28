      module eps_module

      use constants
      implicit none

      private
      public :: eps_gold

      contains

       ! gold dielectric function from Etchegoin et al., J. Chem. Phys. 125, 164705 (2006)
       ! --(errata) J. Chem. Phys. 127, 189901 (2007)-- 
       ! fitting data from Johnson and Christy, Phys. Rev. B 6, 4370 (1972).
       complex(cmp) function eps_gold(the_omega)

        implicit none

        real(dbl), intent(in) :: the_omega            ! in au
        real(dbl)             :: the_lambda           ! in nm

        ! parameters from Etchegoin et al., J. Chem. Phys. 125, 164705 (2006).
        !real(dbl), parameter :: phi_1 = -0.25d0*pi
        !real(dbl), parameter :: phi_2 = -0.25d0*pi
        !real(dbl), parameter :: eps_inf = 1.53d0
        !real(dbl), parameter :: lambda_p = 145.0d0   ! in nm
        !real(dbl), parameter :: gamma_p = 17000.0d0  ! in nm
        !real(dbl), parameter :: a_1 = 0.94d0
        !real(dbl), parameter :: lambda_1 = 468.0d0   ! in nm
        !real(dbl), parameter :: gamma_1 = 2300.0d0   ! in nm
        !real(dbl), parameter :: a_2 = 1.36d0
        !real(dbl), parameter :: lambda_2 = 331.0d0   ! in nm
        !real(dbl), parameter :: gamma_2 = 940.0d0    ! in nm

        ! parameters from (errata) Etchegoin et al., J. Chem. Phys. 127, 189901 (2007).
        !real(dbl), parameter :: phi_1 = -0.25d0*pi
        !real(dbl), parameter :: phi_2 = -0.25d0*pi
        real(dbl), parameter :: eps_inf  = 1.54d0
        real(dbl), parameter :: lambda_p = 143.0d0    ! in nm
        real(dbl), parameter :: gamma_p  = 14500.0d0  ! in nm
        !real(dbl), parameter :: a_1 = 1.27d0
        !real(dbl), parameter :: lambda_1 = 470.0d0   ! in nm
        !real(dbl), parameter :: gamma_1 = 1900.0d0   ! in nm
        !real(dbl), parameter :: a_2 = 1.1d0
        !real(dbl), parameter :: lambda_2 = 325.0d0   ! in nm
        !real(dbl), parameter :: gamma_2 = 1060.0d0   ! in nm
  
        real(dbl), parameter :: a_i(2)      = (/ 1.27d0  , 1.1d0    /)
        real(dbl), parameter :: lambda_i(2) = (/ 470.0d0 , 325.0d0  /)   ! in nm
        real(dbl), parameter :: gamma_i(2)  = (/ 1900.0d0, 1060.0d0 /)   ! in nm

        ! inverse of fine structure constant
        real(dbl), parameter :: alpham1 = 137.035999139d0

        ! conversion factor
        real(dbl), parameter :: bohr_to_nm = TOANGS * 0.1d0

        ! extra auxiliar -> exp(-i phi)
        complex(cmp), parameter :: zz = dcmplx(dsqrt(two)*pt5,-dsqrt(two)*pt5)
        complex(cmp), parameter :: zc = dconjg(zz)

        integer :: i
 
        ! initialization
        the_lambda = 0.0d0
        eps_gold = 0.0d0

        ! from energy to wavelength (in atomic units)
        the_lambda = twp * alpham1 / the_omega

        ! converting wavelength to atomic units
        the_lambda = the_lambda * bohr_to_nm

        ! drude-lorentz term
        eps_gold = eps_inf * onec - one/(lambda_p*lambda_p)/(onec/(the_lambda*the_lambda)+ui/(gamma_p*the_lambda))

        ! for interband transitions
        do i = 1, 2
         eps_gold = eps_gold + a_i(i)/lambda_i(i) * ( zz/(onec/lambda_i(i)-onec/the_lambda-ui/gamma_i(i)) + &
                                                      zc/(onec/lambda_i(i)+onec/the_lambda+ui/gamma_i(i))   )
        end do
                 
       end function eps_gold

      end module eps_module
