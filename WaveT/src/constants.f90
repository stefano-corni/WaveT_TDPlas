      module constants    
      implicit none
      INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)
      INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
      INTEGER, PARAMETER :: i4b = selected_int_kind(9)
      INTEGER, PARAMETER :: cmp = dbl
      INTEGER, PARAMETER :: flg = 10 
      INTEGER, PARAMETER :: nvibmax = 1000
!
      ! constants, conversion factors, and number literals
      real(dbl), parameter :: pi=3.141592653589793D+00
      real(dbl), parameter :: ev_to_au=0.0367493
      real(dbl), parameter :: debye_to_au=0.393456
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
      real(dbl), parameter :: clight=1.37036d2
      real(dbl), parameter :: au_to_cm=219474.63137
      real(dbl), parameter :: cm_to_au=1.d0/au_to_cm
      real(dbl), parameter :: zero=0.d0
      real(dbl), parameter :: one=1.d0
      real(dbl), parameter :: two=2.d0
      real(dbl), parameter :: twp=two*pi
      real(dbl), parameter :: pt5=0.5d0
      real(dbl), parameter :: three=3.d0
      real(dbl), parameter :: plus_inf=HUGE(zero)
      real(dbl), parameter :: zeromin = tiny(zero) 

      integer(i4b), parameter :: one_i=1
      double complex, parameter :: zeroc=(zero,zero)      
      double complex, parameter :: onec=(one,zero)                
      double complex, parameter :: twoc=(two,zero)                
      double complex, parameter :: ui=(zero,one)

      end module constants 
