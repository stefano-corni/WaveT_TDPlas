    module readfile_epsilon
         use tdplas_constants

         implicit none

       ! variables read from eps.inp in the case of the general
       ! dielectric function case (i.e., eps_omega = 'gen')


            integer(i4b)              :: readf_eps_n_omega
            real(dbl)                 :: readf_eps_omega_ini       !< initial value of $\omega$ for spectrum
            real(dbl)                 :: readf_eps_omega_end       !< final value of $\omega$ for spectrum
            real(dbl)                 :: readf_eps_0  = 1000.
            real(dbl)                 :: readf_eps_d  = 1.

            real(dbl), allocatable    :: readf_eps_omegas(:)     !< sampling frequencies for the complex dielectric function
            complex(cmp), allocatable :: readf_eps_epsilons(:)        !< complex dielectric function values for the sampling frequencies



            real(dbl), allocatable    :: readf_eps_re_deps(:)    !< real part of the derivative of the dielectric function at the sampling frequencies
            real(dbl), allocatable    :: readf_eps_im_deps(:)    !< imaginary part of the derivative of the dielectric function at the sampling frequencies
            real(dbl), allocatable    :: readf_eps_func(:)         !< first-order Taylor expansion for eps around frequency of the pole
            real(dbl), allocatable    :: readf_eps_dfunc(:)        !< first derivative of func_eps


      public readf_eps_n_omega, readf_eps_omega_ini, readf_eps_omega_end, readf_eps_omegas, readf_eps_epsilons, &
             readf_eps_re_deps, readf_eps_im_deps, readf_eps_func, readf_eps_dfunc, &
             readf_eps_init


      contains


            subroutine readf_eps_init(eps_d, eps_0)

                    real(dbl)       :: eps_d            !<  $\omega \rightarrow \infty$ limits of $\epsilon(\omega)$
                    real(dbl)       :: eps_0


                    real(dbl)::a,b,c,eps_real,eps_imag
                    integer(i4b)::i,j



                    readf_eps_d   =  eps_d
                    readf_eps_0   =  eps_0

                    open(1,file="eps.inp")
                    read(1,*) readf_eps_n_omega
                    allocate(readf_eps_omegas(readf_eps_n_omega))
                    allocate(readf_eps_epsilons(readf_eps_n_omega))

                    allocate(readf_eps_re_deps(readf_eps_n_omega))
                    allocate(readf_eps_im_deps(readf_eps_n_omega))
                    allocate(readf_eps_func(readf_eps_n_omega))
                    allocate(readf_eps_dfunc(readf_eps_n_omega))
                    do i=1, readf_eps_n_omega
                        read(1,*) readf_eps_omegas(i), eps_real, eps_imag
                        if (i.eq.1) then
                            readf_eps_omega_ini = readf_eps_omegas(i)
                        elseif (i.eq.readf_eps_n_omega) then
                            readf_eps_omega_end = readf_eps_omegas(i)
                        endif
                        readf_eps_epsilons(i)=cmplx(eps_real,eps_imag)
                    enddo
                    call fivepts_stencil(real(readf_eps_epsilons(:),dbl),readf_eps_omegas(:),readf_eps_re_deps(:))
                    call fivepts_stencil(aimag(readf_eps_epsilons(:)),readf_eps_omegas(:),readf_eps_im_deps(:))
                    ! assumption for the first derivative
                    readf_eps_re_deps(1) = readf_eps_re_deps(2)
                    readf_eps_im_deps(1) = readf_eps_im_deps(2)
                    readf_eps_func(1)=real(readf_eps_epsilons(1),dbl)
                    open(3,file='func_eps.inp')
                    open(4,file='re_deps.inp')
                    write(3,*) readf_eps_omegas(1), readf_eps_func(1)
                    write(4,*) readf_eps_omegas(1), readf_eps_re_deps(1)
                    do i=2, readf_eps_n_omega
                        readf_eps_func(i) = real(readf_eps_epsilons(i),dbl)+&
                                            aimag(readf_eps_epsilons(i))*(readf_eps_im_deps(i)/readf_eps_re_deps(i))
                        write(3,*) readf_eps_epsilons(i), readf_eps_func(i)
                        write(4,*) readf_eps_epsilons(i), readf_eps_re_deps(i)
                    enddo
                    call fivepts_stencil(readf_eps_func(:),readf_eps_omegas(:),readf_eps_dfunc(:))
                    ! assumption for the first derivative
                    readf_eps_dfunc(1) = readf_eps_dfunc(2)
                    close(3)
                    close(4)
                    close(1)


            end subroutine

       subroutine fivepts_stencil(func,var,dfunc)
         implicit none
         real(dbl), intent(in) :: func(:)
         real(dbl), intent(in) :: var(:)
         real(dbl), intent(out) :: dfunc(:)
         integer(i4b) :: i
         real(dbl), parameter :: par1 = real(8.0d0,dbl)
         real(dbl), parameter :: par2 = real(12.0d0,dbl)
         ! missing first point derivative, which needs to be assumed outside this subroutine
         ! forward finite differences where there are no enough points to compute 5pts stencil
         dfunc(2) = (func(2)-func(1))/(var(2)-var(1))
         ! computing 5pts stencil
         do i=3, readf_eps_n_omega-2
           !dfunc(i) = (func(i)-func(i-1))/(var(i)-var(i-1))
           dfunc(i) = (-func(i+2)+par1*func(i+1)-par1*func(i-1)+func(i-2))/(var(i)-var(i-1))/par2
         enddo
         ! forward finite differences where there are no enough points to compute 5pts stencil
         dfunc(readf_eps_n_omega-1) = (func(readf_eps_n_omega-1)-func(readf_eps_n_omega-2))&
                                       /(var(readf_eps_n_omega-1)-var(readf_eps_n_omega-2))
         dfunc(readf_eps_n_omega) = (func(readf_eps_n_omega)-func(readf_eps_n_omega-1))&
                                       /(var(readf_eps_n_omega)-var(readf_eps_n_omega-1))
       end subroutine fivepts_stencil


    end module readfile_epsilon



