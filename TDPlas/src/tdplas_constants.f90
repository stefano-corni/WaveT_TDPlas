      module tdplas_constants
      implicit none
      INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)
      INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
      INTEGER, PARAMETER :: i4b = selected_int_kind(9)
      INTEGER, PARAMETER :: cmp = dbl
      INTEGER, PARAMETER :: flg = 20
      INTEGER, PARAMETER :: usr_flg = 40
      INTEGER, PARAMETER :: path = 500
      INTEGER, PARAMETER :: nvibmax = 1000
      INTEGER, PARAMETER :: npulsemax = 10
      integer, parameter :: nsmax     =  10   !< Maximum number of speres/spheroids to read from input!
      integer, parameter :: nmmax     =  10   !< Maximum number of plasmon modes to read from input!
      integer(i4b),parameter :: nts_max=100
      integer(i4b),parameter :: letter_max=100
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

      ! inverse of fine structure constant
      real(dbl), parameter :: alpham1 = 137.035999139d0
      ! conversion factor
      real(dbl), parameter :: bohr_to_nm = TOANGS * 0.1d0



      integer(i4b), parameter :: one_i=1
      complex(cmp), parameter :: zeroc=(zero,zero)
      complex(cmp), parameter :: onec=(one,zero)
      complex(cmp), parameter :: twoc=(two,zero)
      complex(cmp), parameter :: ui=(zero,one)




      !Variables for MPI calculations
      integer(i4b)              :: tp_myrank
      integer(i4b)              :: tp_nproc
      integer(i4b)              :: tp_ierr_mpi
      !MR 16112019 moved here from readio
      integer(i4b)              :: ntst  = 150          !< Threshold value for optimized loops, hard coded

      end module tdplas_constants

