      Module scf            
      use constants
      use readio    
      use interface_tdplas
      use, intrinsic :: iso_c_binding

#ifdef MPI
      use mpi
#endif

      implicit none

      real(dbl), allocatable :: Htot(:,:)    !< Hamiltonian matrix in SCF cycle
      real(dbl), allocatable :: eigt_c(:,:)  !< Eigenvectors of Htot at current cycle
      real(dbl), allocatable :: eigv_c(:)    !< Eigenvalues of Htot at current cycle
      real(dbl), allocatable :: eigt_cp(:,:) !< Eigenvectors of Htot at previous cycle
      real(dbl), allocatable :: eigv_cp(:)   !< Eigenvalues of Htot at previous cycle
      real(dbl) :: maxv                      !< Max values if eigenvectors differences wrt previous cycle                   
      real(dbl) :: maxe                      !< Max values if eigenvalues  differences wrt previous cycle       
      ! Working arrays
      real(dbl) :: mu(3)                    !< Temporary array containing dipole in SCF cycle
      real(dbl), allocatable :: pot(:)      !< Temporary array containing potential in SCF cycle
      real(dbl), allocatable :: c_c(:)      !< Temporary array containing coefficients in old basis 

      save
      private
      public do_scf
!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  DRIVER  ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!------------------------------------------------------------------------
! @brief SCF friver routine 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_scf(c_prev)

       implicit none
       complex(cmp), intent(INOUT) :: c_prev(:)  !< basis state coefficients
       integer(i4b) :: ncyc=1                   !< cycle number 
       logical :: docycle=.true.                !< choice on continue cycling
       real(dbl) :: thre,thrv                   !< thresholds
       real(dbl) :: e_scf, e_ini                !< GS energies
       real(dbl) :: fld(3)                      !  field from charges
       integer(i4b):: max_p(1)   
       integer(i4b):: its 

#ifndef MPI
       myrank=0
#endif

       if (myrank.eq.0) then
          write(6,*) "SCF Cycle"
          write(6,*) "Max_cycle ", this_ncycmax
       endif
       thrv=10**(-this_thrshld+2)
       thre=10**(-this_thrshld)
       if (myrank.eq.0) write(6,*) "Threshold ", thrv,thre
       call init_scf ! Initialize/allocate
       if (Fmdm.ne."vac") call init_environment_scf(c_c,f)
       ! scf cycle
       do while (docycle.and.ncyc.le.this_ncycmax) 
         ! Build the diagonal part of the Hamiltonian 
         call do_htot_ene
         if (Fmdm.ne."vac") call add_interaction_h(Htot)
         ! Diagonalize Hamiltonian           
         eigt_c=Htot
         call diag_mat_in_wavet(eigt_c,eigv_c,n_ci)       
         ! Update charges or field with new coefficients 
         call do_c_oldbasis
         call update_energies(e_scf,e_ini)
         if (Fmdm.ne."vac") call update_environment_scf(c_c,f)
         ! Check convergence                                
         if (ncyc.gt.2) then 
           call check_conv(maxe,maxv,n_ci)       
!           if (maxe.le.thre.and.maxv.le.thrv) docycle=.false.         
! SC 24/4/2016: check convergence only on egeinvalues:
!               in case of degeneracy the variation of eigenvector
!               can be erratic
           if (maxe.le.thre) docycle=.false.         
           if (myrank.eq.0) write(6,*) "cycle ", ncyc, e_scf, e_ini
         endif
         eigt_cp=eigt_c
         eigv_cp=eigv_c
         ncyc=ncyc+1 
         if (myrank.eq.0) then
            write(6,*) "Max Diff on Eigenvalue ", maxe
            write(6,*) "Max Diff on Eigenvector ", maxv
         endif
       enddo
       if (myrank.eq.0) write(6,*) "SCF Done"
       ! Write-out integrals/properties in the new basis 
       if (Fmdm.ne."vac") call out_environment_scf(eigt_c)
       if (myrank.eq.0) then
          call out_dipoles
          call out_energies
       endif
       !  find the new eigenvector that is most similar to the old one
       c_c=abs(matmul(c_i,eigt_c))
       max_p=maxloc(c_c)
       if (myrank.eq.0) write(6,*) 'maxloc',max_p(1)
       c_i=0.d0
       c_i(max_p(1))=1.d0
       c_prev=c_i

       return

      end subroutine do_scf

!------------------------------------------------------------------------
! @brief Init/allocation SCF 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_scf

       allocate(eigv_c(n_ci),eigt_c(n_ci,n_ci))
       allocate(eigv_cp(n_ci),eigt_cp(n_ci,n_ci))
       allocate(Htot(n_ci,n_ci))
       allocate(c_c(n_ci))

       return

      end subroutine init_scf

!------------------------------------------------------------------------
! @brief Finalize/deallocation SCF 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine finalize_scf

       deallocate(eigv_c,eigt_c)
       deallocate(eigv_cp,eigt_cp)
       deallocate(Htot)
       deallocate(c_c)

       return

      end subroutine finalize_scf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  SCF  ROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!------------------------------------------------------------------------
! @brief Write the new coefficients in the old basis 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_c_oldbasis

       implicit none

       integer(i4b):: max_p(1),i    
! SP 12/07/17: avoiding use of automatic arrays, especially in cycles   
       !real(dbl) :: c_c(n_ci)

#ifndef MPI
       myrank=0
#endif

       ! find the new eigenvector that is most similar to the old one
       c_c=abs(matmul(c_i,eigt_c))
       max_p=maxloc(c_c)
       if(this_Fwrite.eq."high") then 
          if (myrank.eq.0) write(6,*) 'maxloc',max_p(1)
       endif
       c_c=0.d0
       c_c(max_p(1))=1.d0
       ! This is the new state on the basis of the old states  
       c_c=matmul(eigt_c,c_c) 
       ! write the state
       if(this_Fwrite.eq."high") then
         if (myrank.eq.0) then
             write(6,*) "State on the basis of original states"
         endif
         do i=1,n_ci
          if (myrank.eq.0) write(6,*) i, c_c(i)
         enddo
         write(6,*)
       endif

       return

      end subroutine do_c_oldbasis


!------------------------------------------------------------------------
! @brief Compute Hamiltonian with charges or fields 
!
! @date Created: S. Pipolo
! Modified: S.Corni 
!------------------------------------------------------------------------
      subroutine do_htot_ene
       integer(4)::i,j,k 
        do j=1,n_ci
          Htot(j,j)=Htot(j,j)+e_ci(j)
          if(this_Fwrite.eq."high") write(6,*) j,Htot(j,j)
        enddo
       return

      end subroutine do_htot_ene


!------------------------------------------------------------------------
! @brief Compute dipoles from coefficients should probably go in a
!        MathTools module
!
! @date Created: S. Pipolo
! Modified:         
!------------------------------------------------------------------------
      subroutine do_dip_from_coeff(c,m)
       complex(cmp), intent(in) :: c(n_ci) !> (1:n_ci) - molecular wavefunction coefficients
       complex(cmp), intent(out):: m(3) !> (1:n_ci)    - molecular dipole                   
       integer(4)::i, 
       do i=1,3
         mu(i)=dot_product(c,matmul(mut(i,:,:),c))
       enddo
       return
      end subroutine do_htot_ene


!------------------------------------------------------------------------
! @brief Check convergence 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine check_conv(mxe,mxv,Mdim)

       real(dbl),intent(inout) :: mxe,mxv
       integer(i4b),intent(in) :: Mdim 
       integer(i4b) :: i,j
       real(dbl):: diff               

       mxv=zero                
       mxe=zero          
!       write(6,*)
       do i=1,Mdim   
         diff=sqrt((eigv_c(i)-eigv_cp(i))**2)
         if(diff.gt.mxe) mxe=diff
!         write (6,*) i,eigt_c(i,:)
         do j=1,Mdim   
! SC 24/4/2016: changed below, otherwise a change of sign result in non-convergence
           diff=abs(eigt_c(j,i)**2-eigt_cp(j,i)**2)
           if(diff.gt.mxv) mxv=diff
         enddo
       enddo

       return

      end subroutine check_conv


!------------------------------------------------------------------------
! @brief Define total energy 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_energies(e_scf,e_ini)

       implicit none

       integer(4) :: i
       real(8) :: e_scf,e_ini

       e_scf=0.d0
       e_ini=0.d0

       do i=1,n_ci
        e_scf=e_scf+abs(c_i(i))*abs(c_i(i))*eigv_c(i)
        e_ini=e_ini+abs(c_i(i))*abs(c_i(i))*e_ci(i)
       enddo

       return

      end subroutine update_energies
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! OUTPUT/TRANSFORMATION ROUTINES !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!------------------------------------------------------------------------
! @brief Transform to the new basis and write out dipole integrals (mut)  
!
! @date Created: S. Pipolo
! Modified: L. Biancorosso 9/23
!------------------------------------------------------------------------
      subroutine out_dipoles

       implicit none

       integer(i4b) :: its,i,j

#ifndef MPI
       myrank=0
#endif

       do its=1,3
        mut(its,:,:)=matmul(mut(its,:,:),eigt_c)
        mut(its,:,:)=matmul(transpose(eigt_c),mut(its,:,:))
       enddo

       if (Fmag.eq.'mag') then
          do its=1,3
             lt(its,:,:)=matmul(lt(its,:,:),eigt_c)
             lt(its,:,:)=matmul(transpose(eigt_c),lt(its,:,:))
          enddo
       endif


! write out the scf dipoles
       open(unit=7,file="ci_mut_scf.inp",status="unknown", &
           form="formatted")
       do i=1,n_ci
        write(7,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0, 'and',i-1,mut(1,1,i),mut(2,1,i),mut(3,1,i)
       enddo
       do i=2,n_ci
         do j=2,i   
           write(7,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1, 'and',i-1,mut(1,i,j),mut(2,i,j),mut(3,i,j)
         enddo
       enddo
       close(unit=7)

       if (Fmag.eq.'mag') then
           open(unit=8,file='ci_lt_scf.inp',status='unknown', &
           form="formatted")
           do i=1,n_ci
              write(8,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0,'and',i-1,lt(1,1,i),lt(2,1,i),lt(3,1,i)
           enddo
           do i=2,n_ci
              do j=2,i
               write(8,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1,'and',i-1,lt(1,i,j),lt(2,i,j),lt(3,i,j)
              enddo
           enddo
           close(unit=8)
       endif

         
       if (myrank.eq.0) write(6,*) "Written out the SCF dipoles"

       return 

      end subroutine out_dipoles      

!------------------------------------------------------------------------
! @brief Write out new energies and reset the zero of energy
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_energies

       implicit none

       integer(i4b) :: i

#ifndef MPI
       myrank=0
#endif

       e_ci=eigv_c
       open(unit=7,file="ci_energy_scf.inp",status="unknown", &
           form="formatted")
       do i=2,n_ci
        e_ci(i)=e_ci(i)-e_ci(1)
        write(7,'(A,I6,X,A,f15.8)') 'Root',i-1,':',e_ci(i)/ev_to_au
       enddo
       e_ci(1)=0.d0
       close(unit=7)
       if (myrank.eq.0) then
       write(6,*) "Written out the SCF energies,", &
               " GS has been given zero energy!"
       endif

       return 

      end subroutine out_energies      

      end module
