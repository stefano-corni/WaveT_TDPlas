   module do_epsilon

         use constants
         use drudel_epsilon
         use debye_epsilon
         use readfile_epsilon
         implicit none

       ! variables read from eps.inp in the case of the general
       ! dielectric function case (i.e., eps_omega = 'gen')


            integer(i4b)              :: do_eps_n_omega

            real(dbl)                 :: do_eps_omega_ini       !< initial value of $\omega$ for spectrum
            real(dbl)                 :: do_eps_omega_end       !< final value of $\omega$ for spectrum

            real(dbl), allocatable    :: do_eps_omegas(:)     !< sampling frequencies for the complex dielectric function
            complex(cmp), allocatable :: do_eps_epsilons(:)        !< complex dielectric function values for the sampling frequencies



      public do_eps_n_omega, do_eps_omega_ini, do_eps_omega_end, do_eps_omegas, do_eps_epsilons, &
             do_eps_init, do_eps, &
             do_eps_fromfile, do_eps_drl, do_eps_deb, do_eps_gold 

      contains




            subroutine do_eps_init(n_omega, omega_ini, omega_end)

                integer   :: n_omega               !< number of points of the discretized spectrum
                real(dbl) :: omega_ini             !< initial value of $\omega$ for spectrum
                real(dbl) :: omega_end             !< final value of $\omega$ for spectrum

                integer   :: i

                do_eps_n_omega   = n_omega
                do_eps_omega_ini = omega_ini
                do_eps_omega_end = omega_end


                allocate(do_eps_omegas(do_eps_n_omega))
                do i=1,n_omega
                   do_eps_omegas(i)=(do_eps_omega_end-do_eps_omega_ini)/&
                                    (do_eps_n_omega-1)*(i-1)+do_eps_omega_ini
                enddo
            end subroutine



            subroutine do_eps(Feps)
                integer :: i
                character(flg)  :: Feps

                allocate(do_eps_epsilons(do_eps_n_omega))
                select case(Feps)
                    case('deb')
                        do i=1,do_eps_n_omega
                            do_eps_epsilons(i) = do_eps_deb(do_eps_omegas(i))
                        enddo
                    case('drl')
                        do i=1,do_eps_n_omega
                            do_eps_epsilons(i) = do_eps_drl(do_eps_omegas(i))
                        enddo
                    case('read')
                        do i=1,do_eps_n_omega
                            do_eps_epsilons(i) = do_eps_fromfile(do_eps_omegas(i))
                        enddo
                    case('gold')
                        do i=1,do_eps_n_omega
                            do_eps_epsilons(i) = do_eps_gold(do_eps_omegas(i))
                        enddo
                end select
            end subroutine do_eps



            complex(cmp) function do_eps_deb(omega)
!------------------------------------------------------------------------
! @brief Compute deb cmplx eps(\omega) and (3*eps(\omega))/(2*eps(\omega)+1)
!
! @date Created: S. Pipolo
! Modified: M.Rosa
!------------------------------------------------------------------------
                real(dbl), intent(in) :: omega

                do_eps_deb = dcmplx(debye_eps_d,zero)+dcmplx(debye_eps_0-debye_eps_d,zero)/ &
                             dcmplx(one,-omega*debye_eps_tau)
            end function do_eps_deb




            complex(cmp) function do_eps_drl(omega)
!------------------------------------------------------------------------
! @brief Compute drl cmplx eps(\omega) and (eps(\omega)-1)/(eps(\omega)+2)
!
! @date Created: S. Pipolo
! Modified: M.Rosa
!------------------------------------------------------------------------
                real(dbl), intent(in) :: omega
                !drudel_eps_gm=drudel_eps_gm+drudel_eps_f_vel/pedra_surf_spheres(1)%r
                do_eps_drl=dcmplx(drudel_eps_A,zero)/dcmplx(drudel_eps_w0**2-omega**2,-omega*drudel_eps_gm)
                do_eps_drl=do_eps_drl+onec

            end function do_eps_drl


          complex(cmp) function do_eps_fromfile(omega)
!------------------------------------------------------------------------------
! @brief Compute gen cmplx eps(\omega) from points through linear interpolation
!
! @date Created: G. Gil
! Modified: M.Rosa
!------------------------------------------------------------------------------
                real(dbl), intent(in) :: omega
                integer(i4b) :: min, max, half

                !questo possiamo metterlo nel check dell'init di freq
                !if((omega.lt.readf_eps_omega_ini).or.(omega.gt.readf_eps_omega_end)) then
                !    call mpi_error("ERROR: you want to calculate the dielectric function value eps_i corresponding to "&
                !                   " a frequency omega_i which is outside the frequency interval available")
                !end if

               ! bisection search of the right frequency interval
                min = 1
                max = readf_eps_n_omega

                do while( min.le.max-1 )
                    half=(min+max)/2
                    if (omega.ge.readf_eps_omegas(half)) then
                        min=half
                    else
                        max=half
                    endif
                enddo

                ! linear interpolation in the right frequency interval
                do_eps_fromfile = (readf_eps_epsilons(max)-readf_eps_epsilons(min))&
                                   /(readf_eps_omegas(max)&
                                   -readf_eps_omegas(min))*(readf_eps_omegas(min))&
                                   +readf_eps_epsilons(min)

            end function do_eps_fromfile




       complex(cmp) function do_eps_gold(omega)

        implicit none

        real(dbl), intent(in) :: omega            ! in au
        real(dbl)             :: lambda           ! in nm

        ! parameters from (errata) Etchegoin et al., J. Chem. Phys. 127, 189901 (2007).
        real(dbl), parameter :: eps_inf  = 1.54d0
        real(dbl), parameter :: lambda_p = 143.0d0    ! in nm
        real(dbl), parameter :: gamma_p  = 14500.0d0  ! in nm
        real(dbl), parameter :: a_i(2)      = (/ 1.27d0  , 1.1d0    /)
        real(dbl), parameter :: lambda_i(2) = (/ 470.0d0 , 325.0d0  /)   ! in nm
        real(dbl), parameter :: gamma_i(2)  = (/ 1900.0d0, 1060.0d0 /)   ! in nm
        ! extra auxiliar -> exp(-i phi)
        complex(cmp), parameter :: zz = dcmplx(dsqrt(two)*pt5,-dsqrt(two)*pt5)
        complex(cmp), parameter :: zc = dconjg(zz)

        integer :: i

        ! initialization
        lambda = 0.0d0
        do_eps_gold = 0.0d0

        ! from energy to wavelength (in atomic units)
        lambda = twp * alpham1 / omega

        ! converting wavelength to atomic units
        lambda = lambda * bohr_to_nm

        ! drude-lorentz term
        do_eps_gold = eps_inf * onec - one/(lambda_p*lambda_p)/(onec/(lambda*lambda)+ui/(gamma_p*lambda))

        ! for interband transitions
        do i = 1, 2
         do_eps_gold = do_eps_gold + a_i(i)/lambda_i(i) * ( zz/(onec/lambda_i(i)-onec/lambda-ui/gamma_i(i)) + &
                                                      zc/(onec/lambda_i(i)+onec/lambda+ui/gamma_i(i))   )
        end do

       end function do_eps_gold


end module





