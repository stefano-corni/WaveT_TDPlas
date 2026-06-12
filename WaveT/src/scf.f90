      Module scf            
      use constants
      use readio    
      use interface_classic
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
      complex(cmp), allocatable :: c_old(:)  !< Temporary array containing coefficients in old basis 
      complex(cmp), allocatable :: c_new(:)  !< Temporary array containing coefficients in old basis 
       real(dbl) :: e_scf, e_ini                !< GS energies

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
      subroutine do_scf

       implicit none
       integer(i4b) :: ncyc=1                   !< cycle number 
       logical :: docycle=.true.                !< choice on continue cycling
       real(dbl) :: thre,thrv                   !< thresholds
       real(dbl) :: f(3)                        !< external field this should be fixed        
       integer(i4b):: max_p(1),i   

#ifndef MPI
       myrank=0
#endif

       if (myrank.eq.0) then
          write(6,*) "SCF Cycle"
          write(6,*) "Max_cycle ", this_ncycmax
       endif
       thrv=10**(-this_thrshld+2)
       thre=10**(-this_thrshld)
       if (myrank.eq.0) write(6,*) "Thresholds ", thrv,thre
       ! Initialize/allocate
       if (myrank.eq.0) write(6,*) "Initialising SCF "
       call init_scf 
       if (Fmdm.ne."vac") call update_environment_scf(c_i,f)
       ! scf cycle
       if (myrank.eq.0) then 
         write(6,*) "Starting SCF Cycle"
         write(6,*) "Cycle, e_scf, e_ini, Max_Diff_Eigenval, " 
         write(6,*) "     Max_Diff_Eigenvec "
       endif
       do while (docycle.and.ncyc.le.this_ncycmax) 
         ! Build the diagonal part of the Hamiltonian 
         Htot(:,:)=zero
         call do_htot_ene
         if (Fmdm.ne."vac") call do_interaction(Htot)
         ! Diagonalize Hamiltonian           
         eigt_c=Htot
         call diag_mat_in_wavet(eigt_c,eigv_c,n_ci)       
         ! Transform the new state on the old basis      
         call do_c_oldbasis
         ! compute scf and initial energies
         call update_energies
         ! Update charges or field with new coefficients, 
         ! Test
         !call transform_dipoles
         !if(Fmdm.ne."vac") call transform_environment_scf(eigt_c)
         !if (Fmdm.ne."vac") call update_environment_scf(c_new,f)
         ! Test
         if (Fmdm.ne."vac") call update_environment_scf(c_old,f)
         ! Check convergence                                
         if (ncyc.gt.2) then 
           call check_conv(maxe,maxv,n_ci)       
           ! If (maxe.le.thre.and.maxv.le.thrv) docycle=.false.         
           ! SC 24/4/2016: check convergence only on egeinvalues:
           ! in case of degeneracy the variation of eigenvector
           ! can be erratic
           if (maxe.le.thre) docycle=.false.         
         endif
         if (myrank.eq.0) write(6,*) ncyc, e_scf, e_ini, maxe, maxv
         eigt_cp=eigt_c
         eigv_cp=eigv_c
         ncyc=ncyc+1 
       enddo
       if (myrank.eq.0) write(6,*) "SCF Done"
       ! Write-out integrals/properties in the new basis 
       call transform_dipoles
       if(Fmdm.ne."vac") call transform_environment_scf(eigt_c)
       c_i=c_new
       e_ci=eigv_c
       if (myrank.eq.0) then
          if(Fmdm.ne."vac") call out_environment_scf
          call out_dipoles
          call out_energies
       endif
       ! Update the initial coefficients      

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
       allocate(c_old(n_ci))
       allocate(c_new(n_ci))
       eigv_c=e_ci
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
       deallocate(c_old)
       deallocate(c_new)

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
       real(dbl), allocatable :: c_tmp(:)     !< coefficients

#ifndef MPI
       myrank=0
#endif
       allocate(c_tmp(n_ci))
       c_tmp(:)=real(c_i(:))
       ! find the new eigenvector that is most similar to the old one
       c_tmp=abs(matmul(c_tmp,eigt_c))
       ! c_tmp here contains the coefficient of the old occupied state
       ! for each of the new states
       max_p=maxloc(c_tmp)
       ! max_p is the position of the maximum value in c_tmp
       if(this_Fwrite.eq."high") then 
          if (myrank.eq.0) write(6,*) 'maxloc',max_p(1)
       endif
       if (max_p(1).ne.1) stop
       c_tmp=0.d0
       c_tmp(max_p(1))=1.d0
       ! c_tmp now has value equal 1 only for the new eigenvector that
       ! is most similar to the old one
       do i=1,n_ci
          c_new(i)=complex(c_tmp(i),0.d0)
       enddo
       ! c_tmp below is the new state on the basis of the old states  
       c_tmp=matmul(eigt_c,c_tmp)
       do i=1,n_ci
          c_old(i)=complex(c_tmp(i),0.d0)
       enddo
       ! write the state
       if(this_Fwrite.eq."high") then
         if (myrank.eq.0) then
             write(6,*) "State on the basis of original states"
         endif
         do i=1,n_ci
          if (myrank.eq.0) write(6,*) i, c_old(i)
         enddo
         write(6,*)
       endif
       deallocate(c_tmp)
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
          !Htot(j,j)=Htot(j,j)+eigv_c(j)
          !Htot(j,j)=Htot(j,j)+e_ci(j)
          Htot(j,j)=+e_ci(j)
          if(this_Fwrite.eq."high") write(6,*) j,Htot(j,j)
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

       real(dbl),intent(out) :: mxe,mxv
       integer(i4b),intent(in) :: Mdim 
       integer(i4b) :: i,j
       real(dbl):: diff               

       mxv=zero                
       mxe=zero          
       do i=1,Mdim   
         diff=sqrt((eigv_c(i)-eigv_cp(i))**2)
         if(diff.gt.mxe) mxe=diff
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
      subroutine update_energies

       implicit none
       integer(4) :: i

       e_scf=zero
       e_ini=zero

       do i=1,n_ci
        e_scf=e_scf+abs(c_new(i))*abs(c_new(i))*eigv_c(i)
        e_ini=e_ini+abs(c_new(i))*abs(c_new(i))*e_ci(i)
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
! @brief Transform dipole integrals to the new basis  
!
! @date Created: S. Pipolo
! Modified: L. Biancorosso 9/23
!------------------------------------------------------------------------
      subroutine transform_dipoles

       implicit none

       integer(i4b) :: its,i,j

#ifndef MPI
       myrank=0
#endif
       do i=1,3
        mut(i,:,:)=matmul(mut(i,:,:),eigt_c)
        mut(i,:,:)=matmul(transpose(eigt_c),mut(i,:,:))
       enddo
       if (Fmag.eq.'mag') then
          do i=1,3
             lt(i,:,:)=matmul(lt(i,:,:),eigt_c)
             lt(i,:,:)=matmul(transpose(eigt_c),lt(i,:,:))
          enddo
       endif
       return 
      end subroutine transform_dipoles      

!------------------------------------------------------------------------
! @brief write out dipole integrals (mut)  
!
! @date Created: S. Pipolo
! Modified: L. Biancorosso 9/23
!------------------------------------------------------------------------
      subroutine out_dipoles

       implicit none

       integer(i4b) :: i,j

#ifndef MPI
       myrank=0
#endif

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
