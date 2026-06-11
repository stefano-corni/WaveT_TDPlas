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

            real(dbl)  :: tomega,mu_trans(3)

            public  read_molecule_file

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







      
      end module
