        module global_quantum 
            use constants


#ifdef MPI
#ifndef SCALI
            use mpi
#endif
#endif

            implicit none

#ifdef MPI
#ifdef SCALI
            include 'mpif.h'
#endif
#endif

            ! block of global variables to be supplied by WaveT
            character(flg)            :: quantum_Ffld                             ! shape of impulse
            character(flg)            :: quantum_Fbin                             ! binary output
            character(flg)            :: quantum_Fopt                             ! OMP-optimized matrix/vector multiplication
            real(dbl)                 :: quantum_dt                               ! time step
            integer(i4b)              :: quantum_n_ci, quantum_n_ci_read          ! number of CIS states
            real(dbl), allocatable    :: quantum_c_i(:), quantum_e_ci(:)          ! CIS coefficients and energies
            real(dbl), allocatable    :: quantum_mut(:,:,:)                       ! CIS transition dipoles
            real(dbl)                 :: quantum_mol_cc(3)                        ! molecule center
            real(dbl)                 :: quantum_fmax(3,10),quantum_omega(10)     ! field amplitude and frequency
            !real(dbl)                 :: quantum_tdelay(10),quantum_pshift(10)    ! time delay and phase shift
            integer(i4b)              :: quantum_n_out,quantum_n_f                ! auxiliaries for output
            integer(i4b)              :: quantum_nthr                             ! number of threads
            !character(1)              :: quantum_medium_res                       !restart for medium
            integer(i4b)              :: quantum_n_res                            ! frequency for restart


            contains

            subroutine quantum_init(dt, mol_cc, n_ci, n_ci_read, mut, e_ci, c_i, &
                                   fmax, omega, Ffld, n_out, n_f, &
                                   Fbin,  Fopt, nthr, n_res)
!------------------------------------------------------------------------
! @brief Set global variables for medium
!
! @date Created: S. Pipolo
! Modified: E. Coccia 28/11/17
!------------------------------------------------------------------------
                implicit none

                real(dbl)     , intent(in)  :: dt                                    ! time step
                integer(i4b)  , intent(in)  :: n_ci, n_ci_read                     ! number of CIS states
                real(dbl)     , intent(in)  :: e_ci(:), c_i(:)
                real(dbl)     , intent(in)  :: mut(:,:,:)               ! CIS transition dipoles
                real(dbl)     , intent(in)  :: mol_cc(3)                ! molecule center
                real(dbl)     , intent(in)  :: fmax(3,10),omega(10)     ! field amplitude and frequency
                character(3)  , intent(in)  :: Ffld                     ! shape of impulse
                character(3)  , intent(in)  :: Fbin                     ! binary output
                character(3)  , intent(in)  :: Fopt                     ! matrix/vector multiplication
                integer(i4b)  , intent(in)  :: n_out,n_f                ! auxiliaries for output
                integer(i4b)  , intent(in)  :: nthr                     ! number of threads
                integer(i4b)  , intent(in)  :: n_res                    ! frequency for restart

                quantum_dt=dt
                quantum_mol_cc=mol_cc
                quantum_n_ci=n_ci
                quantum_n_ci_read = n_ci_read
                allocate(quantum_e_ci(quantum_n_ci))
                quantum_e_ci = e_ci                
                quantum_c_i = c_i
                allocate(quantum_mut(3,quantum_n_ci,quantum_n_ci))
                quantum_mut=mut
                quantum_fmax=fmax
                quantum_omega=omega
                quantum_Ffld=Ffld
                quantum_n_out=n_out
                quantum_n_f=n_f
                quantum_Fbin=Fbin
                quantum_Fopt=Fopt
                quantum_nthr=nthr
                quantum_n_res=n_res

                return

            end subroutine

        end module
