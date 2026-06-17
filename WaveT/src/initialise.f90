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
      use readio       
      use interface_classic
      use QM_coupling      
      use scf              
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
      real(dbl), allocatable :: energies(:)           !< Energies of the states     
      real(dbl), allocatable :: trans_dipoles(:,:,:)  !< Transition dipoles between states
      real(dbl), allocatable :: trans_mag(:,:,:)      !< Mag. Trans. dipoles between states - MM - test
      complex(cmp), allocatable :: coeff0(:)          !< Initial coefficients 
      integer(i4b) :: nstates                         !< Dimension of Hilbert space

      save
      private
      public init_Hspace, & ! subroutines
             energies, trans_dipoles, nstates, coeff0, &   ! variables   
             trans_mag ! MM
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
      subroutine init_Hspace(f0)
       implicit none 
       real(dbl), intent(in) :: f0(3)      !< Initial field 
       integer :: ici
#ifndef MPI
       myrank=0
#endif
       call init_initialise  
       if (Fmdm.ne."vac") then
         ! SP: better define initial external field, this is the linearly polarized one
         write(6,*) "Initializing the environment"
         call init_environment(c_i,fmax(:,1))
       endif
       if (this_Finit_int.eq."sce") then
         !> SCF initialisation 
         call do_scf(nstates,energies,trans_dipoles,f0)
       elseif (this_Finit_int.eq."qmt") then
         !> Quantum Coupling initialisation 
         call do_QM_coupling(nstates,energies,trans_dipoles,f0)
       else
         !> Input initialisation as in ci_*.inp 
         energies=e_ci
         trans_dipoles=mut
         if (Fmag.eq.'mag') then
            trans_mag=lt !MM
         endif
       endif
       ! The following lines need to be modified in order to initialize
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
      end subroutine init_Hspace
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
       if (Fmdm.ne."vac") then
         if (this_Finit_int.eq."sce") then
           write(6,*)"System initialised with self-consistent procedure"
           nstates=n_ci 
         elseif (this_Finit_int.eq."nsc") then
           write(6,*)"System initialised as from input files"
           nstates=n_ci 
         elseif (this_Finit_int.eq."qmt") then
           write(6,*) "System initialised in plexciton states"
           if(global_sys_Ftest.eq."qmt") then 
             nmodes=3
             qmmodes(1)=2 
             qmmodes(2)=3 
             qmmodes(3)=4 
           endif
           nstates=n_ci*nmodes
         else
           write(6,*)"WARNING: No initialisation specified, "
           write(6,*)"  using input file initialisation. "
           nstates=n_ci 
         endif
       ! here one may add a scf initialisation with a static electric
       ! field
       else
         !write(6,*) "System initialised as in input files"
         nstates=n_ci 
         !stop
       endif
       write(6,*) "Hilbert Space has ",nstates, " states."
       allocate(energies(nstates))
       allocate(coeff0(nstates))
       allocate(trans_dipoles(3,nstates,nstates))
       if (Fmag.eq.'mag') then
          allocate(trans_mag(3,nstates,nstates)) !MM
       endif
       coeff0=zeroc
       energies=zero
       trans_dipoles=zero
       if (Fmag.eq.'mag') then !MM
          trans_mag=zero
       endif
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
