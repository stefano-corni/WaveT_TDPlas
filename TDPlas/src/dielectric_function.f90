   module dielectric_function
         use auxiliary_functions
         use tdplas_constants
         use drudel_epsilon
         use debye_epsilon
         use readfile_epsilon
         implicit none

       ! variables read from eps.inp in the case of the general
       ! dielectric function case (i.e., eps_omega = 'gen')


            integer(i4b)              :: dielectric_func_n_omega
            real(dbl)                 :: dielectric_func_omega_ini       !< initial value of $\omega$ for spectrum
            real(dbl)                 :: dielectric_func_omega_end       !< final value of $\omega$ for spectrum

            real(dbl), allocatable    :: dielectric_func_omegas(:)     !< sampling frequencies for the complex dielectric function
            complex(cmp), allocatable :: dielectric_func_epsilons(:)        !< complex dielectric function values for the sampling frequencies



      public dielectric_func_n_omega, dielectric_func_omega_ini, dielectric_func_omega_end,&
             dielectric_func_omegas, dielectric_func_epsilons, &
             dielectric_func_init, dielectric_func_do_eps, &
             dielectric_func_eps_fromfile, dielectric_func_eps_drl, &
             dielectric_func_eps_deb, dielectric_func_eps_gold

      contains




            subroutine dielectric_func_init(n_omega, omega_ini, omega_end)

                integer   :: n_omega               !< number of points of the discretized spectrum
                real(dbl) :: omega_ini             !< initial value of $\omega$ for spectrum
                real(dbl) :: omega_end             !< final value of $\omega$ for spectrum

                integer   :: i

                dielectric_func_n_omega   = n_omega
                dielectric_func_omega_ini = omega_ini
                dielectric_func_omega_end = omega_end


                allocate(dielectric_func_omegas(dielectric_func_n_omega))
                do i=1,n_omega
                   dielectric_func_omegas(i)=(dielectric_func_omega_end-dielectric_func_omega_ini)/&
                                    (dielectric_func_n_omega-1)*(i-1)+dielectric_func_omega_ini
                enddo
            end subroutine



            subroutine dielectric_func_do_eps(Feps)
                integer :: i
                character(flg)  :: Feps

                allocate(dielectric_func_epsilons(dielectric_func_n_omega))
                select case(Feps)
                    case('deb')
                        do i=1,dielectric_func_n_omega
                            dielectric_func_epsilons(i) = dielectric_func_eps_deb(dielectric_func_omegas(i))
                        enddo
                    case('drl')
                        do i=1,dielectric_func_n_omega
                            dielectric_func_epsilons(i) = dielectric_func_eps_drl(dielectric_func_omegas(i))
                        enddo
                    case('gen')
                        do i=1,dielectric_func_n_omega
                            dielectric_func_epsilons(i) = dielectric_func_eps_fromfile(dielectric_func_omegas(i))
                        enddo
                    case('gold')
                        do i=1,dielectric_func_n_omega
                            dielectric_func_epsilons(i) = dielectric_func_eps_gold(dielectric_func_omegas(i))
                        enddo
                end select
            end subroutine dielectric_func_do_eps



            complex(cmp) function dielectric_func_eps_deb(omega)
!------------------------------------------------------------------------
! @brief Compute deb cmplx eps(\omega) and (3*eps(\omega))/(2*eps(\omega)+1)
!
! @date Created: S. Pipolo
! Modified: M.Rosa
!------------------------------------------------------------------------
                real(dbl), intent(in) :: omega

                dielectric_func_eps_deb = dcmplx(debye_eps_d,zero)+dcmplx(debye_eps_0-debye_eps_d,zero)/ &
                             dcmplx(one,-omega*debye_eps_tau)
            end function dielectric_func_eps_deb




            complex(cmp) function dielectric_func_eps_drl(omega)
!------------------------------------------------------------------------
! @brief Compute drl cmplx eps(\omega) and (eps(\omega)-1)/(eps(\omega)+2)
!
! @date Created: S. Pipolo
! Modified: M.Rosa
!------------------------------------------------------------------------
                real(dbl), intent(in) :: omega
                !drudel_eps_gm=drudel_eps_gm+drudel_eps_f_vel/pedra_surf_spheres(1)%r
                dielectric_func_eps_drl=dcmplx(drudel_eps_A,zero)/dcmplx(drudel_eps_w0**2-omega**2,-omega*drudel_eps_gm)
                dielectric_func_eps_drl=dielectric_func_eps_drl+onec

            end function dielectric_func_eps_drl


          complex(cmp) function dielectric_func_eps_fromfile(omega)
!------------------------------------------------------------------------------
! @brief Compute gen cmplx eps(\omega) from points through linear interpolation
!
! @date Created: G. Gil
! Modified: M.Rosa
!------------------------------------------------------------------------------
                real(dbl), intent(in) :: omega
                integer(i4b) :: min, max, half

                !questo andrebbe da un'altra parte
                if((omega.lt.readf_eps_omega_ini).or.(omega.gt.readf_eps_omega_end)) then
                    call mpi_error("ERROR: you want to calculate the dielectric function value eps_i corresponding to ", &
                                   " a frequency omega_i which is outside the frequency interval available", " ")
                end if

                ! bisection search of the right frequency interval
                min = 1
                max = readf_eps_n_omega
                half = (min+max)/2
                do while( min.le.half-1 )
                    if (omega.gt.readf_eps_omegas(half)) then
                        min=half
                    else
                        max=half
                    endif
                    half=(min+max)/2
                enddo


 !               min = 1
 !               max = 1
 !
 !               do while(omega.ge.readf_eps_omegas(max))
 !                   max = max+1
 !               end do
 !               min = max-1

                if (omega.eq.readf_eps_omegas(min)) then
                    dielectric_func_eps_fromfile = readf_eps_epsilons(min)
                else
                    dielectric_func_eps_fromfile = (readf_eps_epsilons(max)-readf_eps_epsilons(min))&
                                   /(readf_eps_omegas(max)-readf_eps_omegas(min))&
                                   *(omega-readf_eps_omegas(min))+readf_eps_epsilons(min)
                endif




            end function dielectric_func_eps_fromfile






       complex(cmp) function dielectric_func_eps_gold(omega)

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
        dielectric_func_eps_gold = 0.0d0

        ! from energy to wavelength (in atomic units)
        lambda = twp * alpham1 / omega

        ! converting wavelength to atomic units
        lambda = lambda * bohr_to_nm

        ! drude-lorentz term
        dielectric_func_eps_gold = eps_inf * onec - one/(lambda_p*lambda_p)/(onec/(lambda*lambda)+ui/(gamma_p*lambda))

        ! for interband transitions
        do i = 1, 2
         dielectric_func_eps_gold = dielectric_func_eps_gold + a_i(i)/lambda_i(i) *&
                                                    ( zz/(onec/lambda_i(i)-onec/lambda-ui/gamma_i(i)) + &
                                                      zc/(onec/lambda_i(i)+onec/lambda+ui/gamma_i(i))   )
        end do

       end function dielectric_func_eps_gold


end module





