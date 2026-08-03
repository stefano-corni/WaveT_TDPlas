module interface_quantum
      use tdplas_constants
      use global_tdplas
      use MathTools
      use pedra_friends
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

      !character(flg) :: 
      !real(dbl), allocatable :: 
      !integer(i4b) :: 
      real(dbl), allocatable :: mut(:,:,:),e_ci(:)
      real(dbl), allocatable :: quantum_vts(:,:,:),quantum_vtsn(:)
      integer(i4b) :: n_ci,nts
     
      public  read_molecule_file

      contains
  

!-----------------------------------------------------------------------------
! @brief read ci_mut.inp and ci_energy.inp and ci_pot.inp
!        when requested for a pl calculation
!
! @date Created: 07/05/2021 G. Dall'Osto
! Modified: Silvio Pipolo 01/09/26
!-----------------------------------------------------------------------------
      subroutine read_molecule_file
         character(4) :: junk
         integer(i4b) :: i, j,its,nts
         real(dbl)  :: scr


         n_ci=global_ext_pert_n_ci+1
         nts=pedra_surf_n_tessere

         open(7,file="ci_mut.inp",status="old")
       allocate (mut(3,n_ci,n_ci))
       do i=1,n_ci
        if (i.le.n_ci) then
         read(7,*)junk,junk,junk,junk,mut(1,1,i),mut(2,1,i),mut(3,1,i)
         mut(:,i,1)=mut(:,1,i)
        else
         read(7,*)
        endif
       enddo
       do i=2,n_ci
         do j=2,i
          if (i.le.n_ci.and.j.le.n_ci) then
           read(7,*)junk,junk,junk,junk,mut(1,i,j),mut(2,i,j),mut(3,i,j)
           mut(:,j,i)=mut(:,i,j)
          else
           read(7,*)
          endif
         enddo
       enddo
       close(7)

       open(7,file="ci_energy.inp",status="old")
       allocate (e_ci(n_ci))
       e_ci(1)=0.d0
       write (6,*) "Excitation energies read from input, in Hartree"
       do i=2,n_ci
         read(7,*) junk,junk,junk,e_ci(i)
         e_ci(i)=e_ci(i)*ev_to_au
         write (6,*) e_ci(i)
       enddo
       close(7)
       ! Reading potential from read_gau_out_medium in WaveT
       open(7,file="ci_pot.inp",status="old")
       read(7,*) nts
       allocate (quantum_vts(nts,n_ci,n_ci))
       allocate (quantum_vtsn(nts))
       quantum_vts=zero
       quantum_vtsn=zero
       ! V00
       read(7,*)
       do its=1,nts
        read(7,*) quantum_vts(its,1,1),scr,quantum_vtsn(its)
       enddo
       !all the others
10     read(7,*,end=20) i,j
       i=i+1
       j=j+1
       if (i.le.n_ci.and.j.le.n_ci) then
        do its=1,nts
         read(7,*) quantum_vts(its,i,j)
         quantum_vts(its,j,i)=quantum_vts(its,i,j)
        enddo
       else
        do its=1,nts
         read(7,*)
        enddo
       endif
       goto 10
20     close(7)
       do i=1,n_ci
        do its=1,nts
         quantum_vts(its,i,i)=quantum_vts(its,i,i)+quantum_vtsn(its)
        enddo
       enddo
       write (6,*) "Done reading in potentials from ci_pot.inp"

       return
       end subroutine read_molecule_file


!-----------------------------------------------------------------------------
! @brief get potential for calculations                  
! @date Created: 01/09/26 Silvio Pipolo 
! Modified: 
!-----------------------------------------------------------------------------
      subroutine get_potential(istate,estate,pot)
       integer(i4b), intent(IN)  :: istate,estate
       real(dbl)   , intent(OUT) :: pot(nts)
        
       pot(:)=quantum_vts(:,istate,estate)

       return
      end subroutine get_potential
 

!-----------------------------------------------------------------------------
! @brief get dipole for calculations                  
! @date Created: 01/09/26 Silvio Pipolo 
! Modified: 
!-----------------------------------------------------------------------------
      subroutine get_dipole(istate,estate,mu)
       integer(i4b), intent(IN)  :: istate,estate
       real(dbl)   , intent(OUT) :: mu(3)
        
       mu(:)=mut(:,istate,estate)

       return
      end subroutine get_dipole    


!-----------------------------------------------------------------------------
! @brief get energy for calculations                  
! @date Created: 01/09/26 Silvio Pipolo 
! Modified: 
!-----------------------------------------------------------------------------
      subroutine get_energy(istate,ene)
       integer(i4b), intent(IN)  :: istate
       real(dbl)   , intent(OUT) :: ene
        
       ene=e_ci(istate)

       return
      end subroutine get_energy
 
end module interface_quantum


