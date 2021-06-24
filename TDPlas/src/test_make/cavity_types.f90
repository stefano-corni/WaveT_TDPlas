        Module cavity_types
            use constants
            implicit none

            type tessera
             real(dbl) :: z
             real(dbl) :: phi
             real(dbl) :: fz
             real(dbl) :: dfdz
             real(dbl) :: dz
             real(dbl) :: area
            end type
!
            type tess_pcm
             real(dbl) :: x
             real(dbl) :: y
             real(dbl) :: z
             real(dbl) :: area
             real(dbl) :: n(3)
             real(dbl) :: rsfe
            end type
!
            type sfera
             real(dbl) :: x
             real(dbl) :: y
             real(dbl) :: z
             real(dbl) :: r
            end type

            type poles_t
              real(dbl), allocatable    :: omega_p(:)          !< real part of the poles of the diagonal Kerne
              real(dbl), allocatable    :: gamma_p(:)          !< imaginary part of the poles of the diagonal
              complex(cmp), allocatable :: eps_omega_p(:)      !< complex dielectric function valued on the re
              real(dbl), allocatable    :: re_deps_domega_p(:) !< real part of the derivative of the complex d
              real(dbl), allocatable    :: im_deps_domega_p(:) !< real part of the derivative of the complex d
              real(dbl), allocatable    :: A_coeff_p(:)        !< numerator in DL expansion
            end type

            save 
            public tessera, tess_pcm, sfera, poles_t

            contains

      end module
