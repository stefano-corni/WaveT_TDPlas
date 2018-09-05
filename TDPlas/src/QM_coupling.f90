      Module QM_coupling    
      use constants    
      use global_tdplas
      use readio_medium
      use pedra_friends
      use MathTools 
      use BEM_medium
      use, intrinsic :: iso_c_binding
#ifdef OMP
      use omp_lib
#endif

#ifdef MPI
#ifndef SCALI
      use mpi
#endif
#endif

      implicit none

      real(dbl), allocatable :: occ(:)       !<Occupations                      
      real(dbl), allocatable :: Hqm(:,:)     !<QM-coupling super-matrix   
      real(dbl), allocatable :: Hqm_evt(:,:) !<eigenvalues of QM-coupling super-matrix 
      real(dbl), allocatable :: Hqm_evl(:)   !<eigenvectors of QM-coupling super-matrix 
      integer(i4b) :: Hqm_dim  !<dimension of of QM-coupling super-matrix 
      integer(i4b) :: nmodes  ! quantum plasmonic modes                  

      save
      private
      public do_QM_coupling, & ! subroutines
             Hqm_evt,Hqm_evl ! variables   
!
!------------------------------------------------------------------------
! @brief Module for molecule-environment QM coupling.      
! @param 
!------------------------------------------------------------------------
!
      contains
!
      subroutine do_QM_coupling                         
!------------------------------------------------------------------------
!     @brief Driver routine of QM_coupling  
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param  
!----------------------------------------------------------------------------
       ! allocate matrices and initialize                                   

#ifndef MPI
       myrank=0
#endif

       call init_QM_coupling 
       if (myrank.eq.0) write(6,*) "QM_coupling correcty initialized"
       ! Build Hamiltonian Super-matrix
       call do_matrix
       if (myrank.eq.0) write(6,*) "Super matrix has been built"
       if(Ftest.eq."qmt") then
         call do_vts_from_dip
         if (myrank.eq.0) call test_QM_coupling
       elseif(FQBEM(1:4)=='diag') then 
         ! Diagonalize Super-matrix        
         Hqm_evt=Hqm
         call diag_mat(Hqm_evt,Hqm_evl,Hqm_dim)
         if (myrank.eq.0)write(6,*) "Super matrix has been diagonalized"
         ! Print Output                    
         if (myrank.eq.0) call out_QM_coupling
       else
         if (myrank.eq.0)write(6,*)"Calculation not implemented yet"
#ifdef MPI
       call mpi_finalize(ierr_mpi)
#endif
         stop
       endif
       ! deallocate                                    
       call fin_QM_coupling 
      return
      end subroutine
!
!
      subroutine init_QM_coupling                         
!------------------------------------------------------------------------
!     @brief Init routine of QM_coupling  
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm_dim,Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
       call do_BEM_quant
       if(FQBEM(1:8)=='diag-all') then ! couple with all modes but only one occupied
         ! Mode 1 is the charge mode w=0
         nmodes=nts_act
         Hqm_dim=n_ci*(nmodes+1)
       else
         write(6,*) "FQBEM=",FQBEM," not implemented yet"
#ifdef MPI 
       call mpi_finalize(ierr_mpi)
#endif
         stop
       endif
       allocate(occ(nts_act))
       if(FQBEM(1:8)=='diag-all') occ=1.d0
       allocate(Hqm(Hqm_dim,Hqm_dim))
       Hqm(:,:)=0.d0
       allocate(Hqm_evt(Hqm_dim,Hqm_dim))
       allocate(Hqm_evl(Hqm_dim))
      return
      end subroutine
!
!
      subroutine fin_QM_coupling                         
!------------------------------------------------------------------------
!     @brief Finalize routine of QM_coupling  
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
       call deallocate_BEM_public
       deallocate(Hqm,Hqm_evt,Hqm_evl)
       deallocate(occ)
      return
      end subroutine
!     
      subroutine do_matrix                       
!------------------------------------------------------------------------
!     @brief Build Hamiltonian Super-matrix
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
       real(dbl):: omega_p  !< mode frequency 
       real(dbl):: we  !< energy factor in coupling 
       real(dbl):: gFi !< molecule-semiclassical_field coupling 
       real(dbl), allocatable:: dp(:) !< \f$ \vec{s}\cdot\vec{F} \f$
       integer(4)::i,j,k,p,s !< indices    
       !
       if(FQBEM(1:4)=='prop') allocate(dp(nts_act))
       ! Build the diagonal superblocs:
       ! H11

       if(FQBEM(1:4)=='prop') then ! propagation_semiclassical
         ! Introduces the coupling with the field for propagation
         do j=1,n_ci
           do k=1,n_ci
             Hqm(k,j)=-dot_product(mut(:,k,j),fmax(:,1))
           enddo
         enddo

         do i=1,nmodes 
           dp(i)=cts_act(i)%x*fmax(1,1)+cts_act(i)%y*fmax(2,1)+       &
                                      cts_act(i)%z*fmax(3,1) 
         enddo 
       endif

       ! CI energies, diagonal 
       do j=1,n_ci
         Hqm(j,j)=Hqm(j,j)+e_ci(j)
       enddo

       gFi=0.d0
       do i=2,nmodes   
         omega_p=sqrt(BEM_W2(i)) 
         we=sqrt((omega_p**2-eps_w0**2)/(two*omega_p))
         ! Introduces the coupling with the field for propagation
         if(FQBEM(1:4)=='prop') gFi=-dot_product(BEM_Modes(i,:),dp(:))*we
         do j=1,n_ci
           p=(i-1)*n_ci+j
           do k=j,n_ci
             s=(i-1)*n_ci+k
             ! Hii 
             Hqm(s,p)=Hqm(k,j)
             Hqm(p,s)=Hqm(s,p)
             ! H1i checked indices j,s simmetrize k,p in H1i block
             Hqm(k,p)=dot_product(BEM_Modes(i,:),vts(:,k,j))*we
             Hqm(j,s)=Hqm(k,p)
             ! Hi1 checked indices p,k simmetrize s,j in Hi1 block
             Hqm(s,j)=Hqm(j,s)
             Hqm(p,k)=Hqm(s,j)
           enddo
           !Hqm(p,p)=Hqm(p,p)+omega_p*(occ(i)+pt5)    
           Hqm(p,p)=Hqm(p,p)+omega_p*occ(i)    
           Hqm(p,j)=Hqm(p,j)+gFi
           Hqm(j,p)=Hqm(j,p)+gFi
         enddo
       enddo

       if(FQBEM.eq."pr_sc") deallocate(dp) 
      return
      end subroutine
!
!
      subroutine out_QM_coupling                         
!------------------------------------------------------------------------
!     @brief Writes output of QM_coupling  
!    
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm_evl  
!----------------------------------------------------------------------------
       integer(i4b) :: i,j   
       character(len=32) :: my_fmt
       open(7,file="Hqm.mat",status="unknown")
       open(8,file="Hqm.ene",status="unknown")
       write(8,*) "Energies: "
       write(my_fmt,'(a,i0,a)') "(",Hqm_dim,"F10.6)"
       write(7,*) "Super-matrix: ", my_fmt
       do i=1,Hqm_dim   
         write(7,my_fmt) (Hqm(i,j), j=1,Hqm_dim)
       enddo
       do i=1,Hqm_dim
         if(i.le.n_ci) then
           write(8,"(i0,3F10.6)") i,Hqm_evl(i),e_ci(i),sqrt(BEM_W2(i))
         else
           write(8,"(i0,F10.6)") i, Hqm_evl(i)
         endif
       enddo
       close(7)
       close(8) 
      return
      end subroutine
!
      subroutine test_QM_coupling                         
!------------------------------------------------------------------------
!     @brief Writes the coupling terms to compare with the dipole approximation
!    
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm_evl  
!----------------------------------------------------------------------------
       real(dbl):: omega_p,we,g,g_ref,r,mud                
       integer(i4b) :: i,j,k   
       real(dbl), allocatable :: sp(:)               
       real(dbl), allocatable :: tot(:,:),ref(:,:)               

#ifndef MPI
       myrank=0
#endif

       allocate(sp(3),tot(n_ci,n_ci),ref(n_ci,n_ci))
       open(7,file="g.mat",status="unknown")
       write(7,*) "# Test for dipolar-mode couplings" 
       r=sqrt(sfe_act(1)%x**2+sfe_act(1)%y**2+sfe_act(1)%z**2)
       sp(1)=sfe_act(1)%x 
       sp(2)=sfe_act(1)%y 
       sp(3)=sfe_act(1)%z 
       write(7,*) "# Sphere radius distance (bohr) and position"
       write(7,"(5F10.4)") cts_act(1)%rsfe,r,sp(1),sp(2),sp(3)
       write(7,*) "# g=dot_product(BEM_Modes(p,:),vts(:,i,j))*we" 
       write(7,*) "# g_ref=mud*sqrt(2*omega_p*cts_act(1)%rsfe^3)/(r^3)"
       write(7,*) "#" 
       write(7,*) "#  p    i    j            g                 g_ref" 
       tot=0.0d0
       do i=1,4        
         omega_p=sqrt(BEM_W2(i)) 
         we=sqrt((omega_p**2-eps_w0**2)/(two*omega_p))
         do j=1,n_ci
           do k=j,n_ci
             mud=dot_product(mut(:,k,j),sp(:))/r
             g=dot_product(BEM_Modes(i,:),vts(:,k,j))*we
             tot(k,j)=tot(k,j)+g*g
             !g_ref: Garcia-Vidal PRL 112, 253601 (2014)
             !g_ref=sqrt(2*mut(i-1,k,j)**2*omega_p*cts_act(1)%rsfe**3)/(r**3)
             g_ref=mud*sqrt(2*sqrt(eps_A/3)*cts_act(1)%rsfe**3)/(r**3)
             ref(k,j)=g_ref
             write(7,"(3i5, 3F20.12)") i,j,k,g,g_ref
           enddo
         enddo
         write(7,*) "" 
       enddo
       do j=1,n_ci
         do k=j,n_ci
           write(7,"(3i5, 3F20.12)") 0,j,k,sqrt(tot(k,j)),ref(k,j)
         enddo
       enddo
       write(7,*) ""
       write(7,*) "# Dipolar resonance frequancy (a.u.)"
       write(7,*) "#  p          omega_p            sqrt(A/3)" 
       do i=2,4        
         write(7,"(i5, 2F20.12)")i, sqrt(BEM_W2(i)), sqrt(eps_A/3)
       enddo
       close(7)
       if (myrank.eq.0) then 
          write(6,*) "Test for dipolar-mode couplings...DONE" 
          write(6,*) "  Results in the g.mat file. " 
       endif
       deallocate(sp,tot,ref)
      return
      end subroutine
!
!
      end module
