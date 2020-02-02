      Module QM_coupling    
      use constants
      use interface_tdplas
      use readio       
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
                                               !> This description comes first.
      character(flg) :: FQBEM                  !< Flag driving the QM calculation mode
      real(dbl), allocatable :: omega_p(:)     !<Plasmon energies                                        
      real(dbl), allocatable :: we(:)          !<Plasmon coupling energy terms in g                      
      real(dbl), allocatable :: g(:,:,:)       !<Plexcitons coupling terms 
      real(dbl), allocatable :: occ(:)         !<Plasmon modes occupations set to 1                      
      real(dbl), allocatable :: Hqm(:,:)       !<QM-coupling matrix \f$ \mathcal{H}_{\text{QM}} \f$ 
      real(dbl), allocatable :: Hqm_int(:,:)   !<Plexcitons-classic_field interaction hamiltonian \f$ \mathcal{H}_{\text{int}} \f$
      real(dbl), allocatable :: plexd(:,:,:)   !<Plexcitons dipole integrals \f$ \boldsymbol{\mu}_{rs}\oplus \mathbf{g}_{Fp} \f$
      real(dbl), allocatable :: Hqm_evt(:,:)   !<Eigenvalues of \mathcal{H}_{\text{QM}} or \mathcal{H}_{\text{SC}} if static field (public)
      real(dbl), allocatable :: Hqm_evl(:)     !<Eigenvectors of \mathcal{H}_{\text{QM}} or \mathcal{H}_{\text{SC}} if static field (public)
      real(dbl), allocatable :: qg(:,:)        !<Charges associated to each mode    
      integer(i4b) :: Hqm_dim  !<Dimension of \mathcal{H}_{\text{QM}} matrix 
      integer(i4b) :: nmodes   !<Quantum plasmonic modes to couple to the molecule                  

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

       implicit none 

#ifndef MPI
       myrank=0
#endif

       call init_QM_coupling 
       if (myrank.eq.0) write(6,*) "QM_coupling correcty initialized"
       call export_mdm_qmcoup

       ! Build Hamiltonian Super-matrix
       !if (allocated(this_vts)) deallocate(this_vts) !solves seg fault due to vts dimensioned as ci_pot.ini
       !allocate (this_vts(this_nts_act,n_ci,n_ci))
       call do_vts_from_dip_in_wavet
       call do_matrix
       if (myrank.eq.0) write(6,*) "Super matrix has been built"
       if(FQBEM(1:4)=='prop') then 
         ! Diagonalize Super-matrix        
         Hqm_evt=Hqm
         call diag_mat_in_wavet(Hqm_evt,Hqm_evl,Hqm_dim)
         if (myrank.eq.0)write(6,*) "Super matrix has been diagonalized"
         ! Print Output                    
       endif
       if(Ftest.eq."qmt") then
         call test_QM_coupling
       endif
         if (myrank.eq.0) call out_QM_coupling
#ifdef MPI
       call mpi_finalize(ierr_mpi)
#endif
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

       implicit none
       FQBEM='diag-all' !enforces use of correct QM_coupling flag
       call do_BEM_quant_in_wavet
!       if(FQBEM(1:8)=='diag-all') then ! couple with all modes but only one occupied
         ! Mode 1 is the charge mode w=0

       nmodes=this_nts_act
       Hqm_dim=n_ci*(this_nts_act+1) !maybe the +1 is not needed?
       
!       else
!         write(6,*) "FQBEM=",FQBEM," not implemented yet"
!#ifdef MPI 
!       call mpi_finalize(ierr_mpi)
!#endif
!         stop
!       endif
       allocate(g(nmodes,n_ci,n_ci))
       allocate(we(nmodes))
       allocate(omega_p(nmodes))
       allocate(qg(nmodes,this_nts_act))
       allocate(occ(nmodes))
       occ=1.d0
       allocate(Hqm(Hqm_dim,Hqm_dim))
       Hqm(:,:)=0.d0
       allocate(Hqm_evt(Hqm_dim,Hqm_dim))
       allocate(Hqm_evl(Hqm_dim))
       write(*,*) "QM initialization done"
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
       implicit none
       call deallocate_BEM_public_in_wavet
       deallocate(we,omega_p,g,qg)
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
       implicit none
!       real(dbl):: omega_p  !< mode frequency 
       real(dbl):: gFi !< molecule-semiclassical_field coupling 
       real(dbl), allocatable:: dp(:) !< \f$ \vec{s}\cdot\vec{F} \f$
       integer(4)::i,j,k,p,s !< indices    
       !
       if(FQBEM(1:4)=='prop') allocate(dp(this_nts_act))
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
           dp(i)=this_cts_act(i)%x*fmax(1,1)+this_cts_act(i)%y*fmax(2,1)+       &
                                      this_cts_act(i)%z*fmax(3,1) 
         enddo 
       endif

       ! CI energies, diagonal 
       do j=1,n_ci
         Hqm(j,j)=Hqm(j,j)+e_ci(j)
       enddo

       gFi=0.d0
       do i=2,nmodes   
         omega_p(i)=sqrt(this_BEM_W2(i)) 
         we(i)=sqrt((omega_p(i)**2-this_eps_w0**2)/(two*omega_p(i)))
         ! Introduces the coupling with the field for propagation
         if(FQBEM(1:4)=='prop') gFi=-dot_product(this_BEM_Modes(i,:),dp(:))*we(i)
         do j=1,n_ci
           p=(i-1)*n_ci+j
           do k=j,n_ci
             s=(i-1)*n_ci+k
             ! Hii 
             Hqm(s,p)=Hqm(k,j)
             Hqm(p,s)=Hqm(s,p)
             ! H1i checked indices j,s simmetrize k,p in H1i block
             Hqm(k,p)=dot_product(this_BEM_Modes(i,:),this_vts(:,k,j))*we(i)
             Hqm(j,s)=Hqm(k,p)
             ! Hi1 checked indices p,k simmetrize s,j in Hi1 block
             Hqm(s,j)=Hqm(j,s)
             Hqm(p,k)=Hqm(s,j)
           enddo
           !Hqm(p,p)=Hqm(p,p)+omega_p*(occ(i)+pt5)    
           Hqm(p,p)=Hqm(p,p)+omega_p(i)*occ(i)    
           Hqm(p,j)=Hqm(p,j)+gFi
           Hqm(j,p)=Hqm(j,p)+gFi
         enddo
       enddo

       if(allocated(dp)) deallocate(dp) 
      return
      end subroutine
      

!
!
!------------------------------------------------------------------------
!>    @brief Writes output of QM_coupling  
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!>    @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine out_QM_coupling                         
       integer(i4b) :: i,j   
       character(len=32) :: my_fmt
       open(7,file="Hqm_matrix.dat",status="unknown")
       open(8,file="Hqm_energies.dat",status="unknown")
       write(8,*) "Energies: "
       write(my_fmt,'(a,i0,a)') "(",Hqm_dim,"F10.6)"
       write(7,*) "Super-matrix: ", my_fmt
       do i=1,Hqm_dim   
         write(7,my_fmt) (Hqm(i,j), j=1,Hqm_dim)
       enddo
       close(7)
       close(8) 
       call out_gcharges_in_wavet
      return
      end subroutine
!
!------------------------------------------------------------------------
!>     @brief computes the molecule-environment quantum couplig elements "g"
!>     @date Created: 25 Sep 2019
!>     @author S.Pipolo
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine do_gcharges  
       integer(i4b) :: i   
       omega_p(1)=zero
       we(1)=zero
       qg(1,:)=zero
       do i=2,nmodes  
         omega_p(i)=sqrt(this_BEM_W2(i)) 
         we(i)=sqrt((omega_p(i)**2-eps_w0**2)/(two*omega_p(i)))
         qg(i,:)=this_BEM_Modes(i,:)*we(i)
       enddo
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>     @brief computes the molecule-environment quantum couplig elements "g"
!>     @date Created: 25 Sep 2019
!>     @author S.Pipolo
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine do_couplings
       integer(i4b) :: i,j,k   
       real(dbl) :: tmp
#ifndef MPI
       myrank=0
#endif
       do i=1,nmodes   
         do j=1,n_ci
           do k=j,n_ci
             g(i,k,j)=dot_product(qg(i,:),this_vts(:,k,j))
           enddo
         enddo
       enddo
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>     @brief Writes out the coupling factors to compare with the dipole approximation
!>     @date Created: 02 May 2017
!>     @author S.Pipolo
!>     @modified J.Fregoni 18 December 2019 (moved print of .pqr to
!out_gcharges)
!>      and MOPAC interface
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine test_QM_coupling                         
       real(dbl):: r,d,mud,wl,gref,gloc
       character(54):: my_fmt
       integer(i4b) :: i,j,k,pmax,p
       real(dbl), allocatable :: sp(:)               
       real(dbl), allocatable :: tot(:,:),ref(:,:)               

#ifndef MPI
       myrank=0
#endif

       allocate(sp(3),tot(n_ci,n_ci),ref(n_ci,n_ci))
       call do_couplings
       call do_gcharges
       open(37,file="g.dat",status="unknown")
       write(37,*) "# Test for dipolar-mode couplings" 
       r=sqrt(this_cts_act(1)%x**2+this_cts_act(1)%y**2+this_cts_act(1)%z**2)
       sp(1)=this_cts_act(1)%x 
       sp(2)=this_cts_act(1)%y 
       sp(3)=this_cts_act(1)%z 
       write(37,*) "# Sphere radius distance (bohr) and position"
       write(37,"(5F10.4)") this_cts_act(1)%rsfe,r,sp(1),sp(2),sp(3)
       write(37,*) "# g=dot_product(this_BEM_Modes(p,:),vts(:,i,j))*we" 
       write(37,*) "# g_ref=mud*sqrt(2*omega_p*cts_act(1)%rsfe^3)/(r^3)"
       write(37,*) "#" 
       write(37,*) "#  p    i    j            g                 g_ref" 
       tot=0.0d0
       do i=2,nmodes 
         omega_p(i)=sqrt(this_BEM_W2(i)) 
         we(i)=sqrt((omega_p(i)**2-this_eps_w0**2)/(two*omega_p(i)))
         do j=1,n_ci
           do k=j,n_ci
             mud=dot_product(mut(:,k,j),sp(:))/r
             gloc=dot_product(this_BEM_Modes(i,:),this_vts(:,k,j))*we(i)
             tot(k,j)=tot(k,j)+gloc*gloc
             gref=mud*sqrt(2*sqrt(this_eps_A/3)*this_cts_act(1)%rsfe**3)/(r**3)
             ref(k,j)=gref
             write(37,"(3i5, 3F20.12)") i,j,k,sqrt(tot(k,j)),ref(k,j)
           enddo
         enddo
         write(37,*) "" 
       enddo
       write(37,*) ""
       write(37,*) "# Dipolar resonance frequency (a.u.)"
       write(37,*) "#  p          omega_p            sqrt(A/3)" 
       do i=2,nmodes        
         write(37,"(i5, 2E20.12)")i, sqrt(BEM_W2(i)), sqrt(eps_A/3)
       enddo
       close(37)
       if (myrank.eq.0) then 
          write(6,*) "Test for dipolar-mode couplings...DONE" 
          write(6,*) "  Results in the g.dat file. " 
       endif
       deallocate(sp,tot,ref)
      return
      end subroutine
end module
