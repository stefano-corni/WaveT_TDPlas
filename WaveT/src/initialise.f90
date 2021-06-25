!------------------------------------------------------------------------------
!        TDPLAS - initialise 
!------------------------------------------------------------------------------
! MODULE        : initialise 
! DATE          : 11 Jul 2020
! REVISION      : V 0.00
!> @authors 
!> S.Pipolo   
!
! DESCRIPTION:
!> Module for initializing Hilbert space used in the propagation. Can
!> perform both SCF or Quantum Coupling with the environment.
!------------------------------------------------------------------------------
      Module initialise    
      use constants    
      use interface_tdplas
      use readio       
      use scf            
      use QM_coupling
      use, intrinsic :: iso_c_binding
#ifdef OMP
      use omp_lib
#endif
#ifdef MPI
      use mpi
#endif
      implicit none
                                                      !> This description comes first.
      !character(flg) :: FQBEM                        !< Flag driving the QM calculation mode
      real(dbl), allocatable :: energies(:)           !<Energies of the states     
      real(dbl), allocatable :: trans_dipoles(:,:,:)  !<Transition dipoles between states
      complex(cmp), allocatable :: coeff0(:)          !<Initial coefficients 
      integer(i4b) :: nstates                         !<Dimension of Hilbert space

      save
      private
      public init_propagation, & ! subroutines
             energies, trans_dipoles, nstates, coeff0   ! variables   
!
!
      contains
!
!------------------------------------------------------------------------
!>    @brief Driver routine of initialise. 
!>    @date Created: 11 Jul 2020 
!>    @author S.Pipolo
!>    @note 
!----------------------------------------------------------------------------
      subroutine init_propagation
       implicit none 
       integer :: ici
#ifndef MPI
       myrank=0
#endif
       call init_initialise  
       if (this_Finit_int.eq."scf") then
         !> SCF initialisation 
         !call do_scf(nstates,energies,trans_dipoles)
       elseif (this_Finit_int.eq."qmt") then
         !> Quantum Coupling initialisation 
         call do_QM_coupling(nstates,energies,trans_dipoles)
       else
         energies=e_ci
         trans_dipoles=mut
         !> Input initialisation as in ci_*.inp 
       endif
       ! The following lines need to be moified in order to initialize
       ! the system in plexciton states greater than n_ci
       do ici=1,n_ci
         coeff0(ici)=c_i(ici)
       enddo
       do ici=n_ci+1,nstates
         coeff0(ici)=zeroc
       enddo
       !> Deallocate matrices: nothing to deallocate, quantities still
       !employed in scf subroutines                                   
       !call fin_initialise  

      return
      end subroutine init_propagation
!
!
!------------------------------------------------------------------------
!     @brief Init routine of initialise   
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm_dim,Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
      subroutine init_initialise
       implicit none
       ! The charge mode w=0 is counted in this_nmodes for testing purposes
       if (this_Finit_int.eq."qmt") call do_BEM_quant_in_wavet
       if (this_Finit_int.eq."scf") then
         write(6,*) "System initialised with self-consistent procedure"
         nstates=n_ci 
       elseif (this_Finit_int.eq."qmt") then
         write(6,*) "System initialised with Quantum Coupling"
         if(global_sys_Ftest.eq."qmt") then 
           this_nmodes=3
           this_qmmodes(1)=2 
           this_qmmodes(2)=3 
           this_qmmodes(3)=4 
         endif
         nstates=n_ci*(this_nmodes+1)
       else
         write(6,*) "System initialised as in input files"
         nstates=n_ci 
         !stop
       endif
       write(6,*) "Hilbert Space with ",nstates, " states."
       allocate(energies(nstates))
       allocate(coeff0(nstates))
       allocate(trans_dipoles(3,nstates,nstates))
       coeff0=zeroc
       energies=zero
       trans_dipoles=zero
      return
      end subroutine init_initialise
!
!
!------------------------------------------------------------------------
!     @brief Finalise routine of initialise   
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
      subroutine fin_initialise
      implicit none
       deallocate(c_i)
      return
      end subroutine fin_initialise

      end module
