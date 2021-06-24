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
      subroutine do_scf(q_or_f,c_prev)

       implicit none
       real(dbl), intent(INOUT):: q_or_f(:)     !< charges or field  
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
       ! scf cycle
       do while (docycle.and.ncyc.le.this_ncycmax) 
         ! Build the Hamiltonian
         if(this_Fint.eq."ons") then
           if(this_Fprop.eq."chr-ief".or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") then 
             call do_field_from_charges_in_wavet(q_or_f,fld)
             call do_htot(fld)
           endif
         else
           call do_htot(q_or_f)
         endif
         ! Diagonalize Hamiltonian           
         eigt_c=Htot
         call diag_mat_in_wavet(eigt_c,eigv_c,n_ci)       
         ! Update charges or field with new coefficients 
         call do_c_oldbasis
         if(this_Fprop.eq."dip") call do_field(q_or_f)
         if(this_Fprop.eq."chr-ief".or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") then
            call do_charges(q_or_f)
         endif
         call do_energies(e_scf,e_ini)
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
       if (this_Fprop.eq.'chr-ief'.or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") then
       ! transform vts to the SCF state basis
        do its=1,this_nts_act
         this_vts(its,:,:)=matmul(this_vts(its,:,:),eigt_c)
         this_vts(its,:,:)=matmul(transpose(eigt_c),this_vts(its,:,:))
        enddo
         if (myrank.eq.0) then
            call out_charges(q_or_f)
            call out_vts
         endif
       endif
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
       if(this_Fprop.eq."chr-ief".or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") allocate(pot(this_nts_act))
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
       if(allocated(pot))deallocate(pot)
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
! @brief Compute field from dipole 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_field(f)

       implicit none 

       real(dbl), intent(OUT):: f(3)     
       integer(i4b)::i    

       do i=1,3
         mu(i)=dot_product(c_c,matmul(mut(i,:,:),c_c))
       enddo
       ! SP 04/0717 matmul for general spheroid orientation
       f=(1.-this_mix_coef)*f+this_mix_coef*matmul(this_mat_f0,mu)

       return

      end subroutine do_field

!------------------------------------------------------------------------
! @brief Compute charges from potential 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_charges(q)

       implicit none 

       real(dbl), intent(OUT):: q(this_nts_act)     
       integer(i4b)::i    

       do i=1,this_nts_act
         pot(i)=dot_product(c_c,matmul(this_vts(i,:,:),c_c))
       enddo 

       q=(1.-this_mix_coef)*q+this_mix_coef*matmul(this_BEM_Q0,pot)
! SC 12/8/2016: apparently for NP, charge compensation is needed
       if (Fmdm.eq.'cnan'.or.Fmdm.eq.'qnan') q=q-sum(q)/this_nts_act

       return

      end subroutine do_charges

!------------------------------------------------------------------------
! @brief Compute Hamiltonian with charges or fields 
!
! @date Created: S. Pipolo
! Modified: S.Corni 
!------------------------------------------------------------------------
      subroutine do_htot(q_or_f)

       real(dbl), intent(IN):: q_or_f(:)     
       integer(4)::i,j,k 
       if (this_Fint.eq."pcm") then
          do j=1,n_ci
             do k=j,n_ci   
                 Htot(k,j)=dot_product(this_vts(:,k,j),q_or_f(:)-this_q0(:))
                 Htot(j,k)=Htot(k,j)
             enddo
             Htot(j,j)=Htot(j,j)+e_ci(j)
             if(this_Fwrite.eq."high") write(6,*) j,Htot(j,j)
          enddo
       else
          do j=1,n_ci
             do k=j,n_ci
           Htot(k,j)=dot_product(mut(:,k,j),q_or_f(:)-this_fr_0(:))
           Htot(j,k)=Htot(k,j)
         enddo
         Htot(j,j)=Htot(j,j)+e_ci(j)
         if(this_Fwrite.eq."high") write(6,*) j,Htot(j,j)
        enddo
       endif

       return

      end subroutine do_htot


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
      subroutine do_energies(e_scf,e_ini)

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

      end subroutine do_energies
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! OUTPUT/TRANSFORMATION ROUTINES !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!------------------------------------------------------------------------
! @brief Write out the charges in the charges0_scf.dat file 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_charges(q)

       implicit none

       real(dbl), intent(IN):: q(this_nts_act)     
       integer(i4b) its

#ifndef MPI
       myrank=0
#endif

       open(unit=7,file="charges0_scf.inp",status="unknown", &
            form="formatted")
         write (7,*) this_nts_act
         do its=1,this_nts_act
          write (7,'(E22.8,F22.10)') q(its)
         enddo
       close(unit=7)
       if (myrank.eq.0) write(6,*) "Written out the SCF charges"
       ! SP 23/10/16: update the q0 vector to have a consistent correction if a
       !              propagation is performed after the SCF cycle
       this_q0(:)=q(:)

       return 

      end subroutine out_charges     


!------------------------------------------------------------------------
! @brief Transform to the new basis and write out potential integrals on
! tesserae (vts) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_vts

       implicit none

       integer(i4b) :: its,i,j

#ifndef MPI
       myrank=0
#endif




       open(unit=7,file="ci_pot_scf.inp",status="unknown", &
          form="formatted")
       write (7,*) this_nts_act
       write (7,*) "V0  check Vnuc"
       do its=1,this_nts_act
        write (7,*) this_vts(its,1,1),0.d0,this_vtsn(its)
       enddo
       do j=2,n_ci
         write(7,*) 0,j-1
         do its=1,this_nts_act
          write(7,*) this_vts(its,1,j)
         enddo
       enddo
       !Vij
       do i=2,n_ci
          do j=i,n_ci   
             write(7,*) i-1,j-1
             do its=1,this_nts_act
                write(7,*) this_vts(its,i,j)             
             enddo
          enddo
       enddo
       close(unit=7)
       if (myrank.eq.0) write(6,*) "Written out the SCF potentials"

       return 

      end subroutine out_vts      

!------------------------------------------------------------------------
! @brief Transform to the new basis and write out dipole integrals (mut)  
!
! @date Created: S. Pipolo
! Modified:
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
