        module readfile_freq
               use tdplas_constants
               use global_tdplas


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

            real(dbl), allocatable :: quantum_vts(:,:,:)       !<transition potentials on tesserae from cis
            real(dbl), allocatable :: quantum_vtsn(:)          !<nuclear potential on tesserae
            real(dbl)  :: tomega,mu_trans(3)

            public  read_molecule_file,read_gau_out_medium

            contains

!-----------------------------------------------------------------------------
! @brief read ci_mut.inp and ci_energy.inp when requested for a pl calculation
!
! @date Created: 07/05/2021 G. Dall'Osto
! Modified:
!-----------------------------------------------------------------------------

      subroutine read_molecule_file()
         real(dbl), allocatable :: mut(:,:,:),e_ci(:)
         character(4) :: junk
         integer(i4b) :: i, j
         integer(i4b) :: n_ci, nstate


         n_ci=global_ext_pert_n_ci+1
         nstate=global_ext_pert_nstate+1

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
       mu_trans=mut(:,1,nstate)
       write(*,*) "Transition dipole_moment considered", mu_trans
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
       tomega = e_ci(nstate)
       write (6,*) "Energy of the chosen state", tomega
       close(7)

       end subroutine read_molecule_file


!------------------------------------------------------------------------
! @brief Read transition potentials on tesserae
!
! @date Created: S. Pipolo
! Modified: E. Coccia
! Modified: S.Corni (27/06/2020): now the state pair is read from ci_pot,
!           we do not assume upper or lower triangular. Should work
!           for current gamess version as well
!------------------------------------------------------------------------
      subroutine read_gau_out_medium(quantum_n_ci)
       integer(i4b), intent(in)  :: quantum_n_ci
       integer(i4b) :: i,j,its,nts
       real(dbl)  :: scr

       open(7,file="ci_pot.inp",status="old")
       read(7,*) nts
       allocate (quantum_vts(nts,quantum_n_ci,quantum_n_ci))
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
       if (i.le.quantum_n_ci.and.j.le.quantum_n_ci) then
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
       do i=1,quantum_n_ci
        do its=1,nts
         quantum_vts(its,i,i)=quantum_vts(its,i,i)+quantum_vtsn(its)
        enddo
       enddo
       write (6,*) "Done reading in potentials from ci_pot.inp"


       return

      end subroutine read_gau_out_medium


!         subroutine read_gau_out_medium()
!
!              integer(i4b)  :: nts
!              integer(i4b) :: i,j, its
!              real(dbl)    :: scr
!
!
!              if((global_medium_FinitBEM.eq."rea").and.(global_prop_Fprop(1:3).eq."chr")) then
!                    open(7,file="ci_pot.inp",status="old")
!                    read(7,*) nts
!                    allocate (quantum_vts(nts, quantum_n_ci, quantum_n_ci))
!                    allocate (quantum_vtsn(nts))
!                    ! V00
!                    read(7,*)
!                    do its=1, nts
!                        read(7,*) quantum_vts(its,1,1),scr, quantum_vtsn(its)
!                        quantum_vts(its,1,1)=quantum_vts(its,1,1)+quantum_vtsn(its)
!                    enddo
!                    !V0j
!                    do j=2,quantum_n_ci_read
!                        read(7,*)
!                        if (j.le.quantum_n_ci) then
!                            do its=1,nts
!                                read(7,*) quantum_vts(its,1,j)
!                            enddo
!                            quantum_vts(:,j,1)=quantum_vts(:,1,j)
!                        else
!                            do its=1,nts
!                                read(7,*)
!                            enddo
!                        endif
!                    enddo
!                    !Vij
!                    do i=2,quantum_n_ci_read
!                        do j=2,i
!                            read(7,*)
!                            if (i.le.quantum_n_ci.and.j.le.quantum_n_ci) then
!                                do its=1,nts
!                                    read(7,*) quantum_vts(its,i,j)
!                                enddo
!                                quantum_vts(:,j,i)=quantum_vts(:,i,j)
!                            else
!                                do its=1,nts
!                                    read(7,*)
!                                enddo
!                            endif
!                        enddo
!                        ! add nuclear potential
!                        if (i.le.quantum_n_ci) then
!                            do its=1,nts
!                                quantum_vts(its,i,i)=quantum_vts(its,i,i)+quantum_vtsn(its)
!                            enddo
!                        endif
!                    enddo
!
!                    close(7)
!            endif
!
!            return
!
!
!
!
!
!
!
!      end subroutine



      
      end module
