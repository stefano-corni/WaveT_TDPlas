!------------------------------------------------------------------------------
!        TDPLAS - QM COUPLING
!------------------------------------------------------------------------------
! MODULE        : QM_coupling
! DATE          : 02 May 2017
! REVISION      : V 0.00
!> @authors 
!> S.Pipolo   
!
! DESCRIPTION:
!> Module for molecule-environment QM coupling. 
!> The basis for a quantum treatment of the Molecule-NanoParticle Hamiltonian 
!> \f$ \mathcal{H}_{MP} \f$ is the generalized eigenmodes decomposition for the 
!> \f$\mathbf{DA}\f$ matrix recently provided in ref. \cite corni2015JPCA .
!> \f{align}{
!>       \mathbf{S}^{\text{-1/2}}\mathbf{DA}\mathbf{S}^{\text{1/2}}=\mathbf{T}\boldsymbol{\Lambda}\mathbf{T}^{\dagger} \label{test}
!> \f}
!> In equation \f$\ref{test}\f$ (reference problem) matrices \f$\mathbf{D}\f$, \f$\mathbf{A}\f$ and \f$\mathbf{S}\f$, are built from the environment's surface discretization. More details on the BEM approach and the eigenmode decomposition of the continuum-environment response matrix implemented in TDPlas can be found in the module \ref bem_medium. Testing the reference to a variable \ref bem_medium.poles_t.
!>  and  what if I continue here?
!>    \f{align}{\nonumber
!>    \left[\mathbf{H}_{\text{MP}}\right]_{rs,p}=-\sqrt{\frac{\omega_p^2-\omega_0^2}{2\omega_p}} \left[ \mathbf{T}^{\dagger}\mathbf{S}^{\text{-1/2}}\right]_p~\mathbf{V}_{rs}
!>    \f}  
!> Array:
!>           \f{array}{{ccc c ccc}
!>             \ddots && \f$\mathbf{H}_{\text{MF}}\f$&\hspace{0.3cm}&\ddots && \f$\mathbf{H}_{\text{MP}}\f$ \\
!>             &\mathbf{H}^0_{\text{M}}+\mathbf{H}^0_{\text{P}}&&\hspace{0.3cm}&& ~~~~\mathbf{H}_{\text{PF}}~~~~ &\\
!>             \mathbf{H}_{\text{MF}}&&\ddots &\hspace{0.3cm}& \mathbf{H}_{\text{MP}}&&\ddots\\                   
!>           \f} 
!------------------------------------------------------------------------------
      Module QM_coupling    
      use constants    
      use interface_tdplas
!      use global_tdplas
      use readio       
!      use pedra_friends
!      use MathTools 
!      use BEM_medium
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
             Hqm_evt,Hqm_evl   ! variables   
!
      contains
!
!
!------------------------------------------------------------------------
!>    @brief Driver routine of QM_coupling. 
!>    @date Created: 02 May 2017 
!>    @author S.Pipolo
!>    @note this routine should change to accomodate quantum external fields
!----------------------------------------------------------------------------
      subroutine do_QM_coupling                         
#ifndef MPI
       myrank=0
#endif
       !> Allocate and initialize matrices
       call init_QM_coupling 
       if (myrank.eq.0) write(6,*) "QM_coupling correcty initialized"
       !> if debugging performs the dipolar test on the spherical couplings and exit
       if(Ftest.eq."qmt") then
         if (allocated(vts)) deallocate(vts)
         allocate (vts(nts_act,n_ci,n_ci))
         call do_vts_from_dip
       endif
       if (myrank.eq.0) write(6,*) "Integrals from dipoles computed"
       !> Build Plexcitons coplings terms "g"    
       call do_couplings      
       if (myrank.eq.0) write(6,*) "couplings computed"
       !> Testing against dipolar model of Garcia-Vidal PRL 112, 253601 (2014)
       if(Ftest.eq."qmt".and.myrank.eq.0) call test_QM_coupling
       !> Build Plexcitons matrix: do_Hqm_matrix 
       call do_Hqm_matrix
       if (myrank.eq.0) write(6,*) "Plexcitons matrix built"
       if (FQBEM(1:4)=='prop') then
         !> Diagonalize Plexciton matrix              
         Hqm_evt=Hqm 
         call diag_mat(Hqm_evt,Hqm_evl,Hqm_dim)
         if (myrank.eq.0)write(6,*) "Plexcitons matrix diagonalized"
         !> Print Energies and eigenstates.            
         if (myrank.eq.0) call out_QM_coupling
       endif
       !> Prepare perturbation integrals if external field is present: do_plexd_matrix
       if(mdl(fmax(:,1)).gt.0.) call do_plexd_matrix
       if (FQBEM(1:4)=='diag') then
         !> If diagonalize, diagonalize Perturbed Plexciton matrix    
         if(mdl(fmax(:,1)).gt.0.) call do_Hqm_int(fmax(:,1))
!         Hqm_evt(k,j)=Hqm+Hqm_int
         Hqm_evt=Hqm+Hqm_int
         if (myrank.eq.0)write(6,*) "Diagonalizing plexciton matrix"
         call diag_mat(Hqm_evt,Hqm_evl,Hqm_dim)
         if (myrank.eq.0)write(6,*) &
                "Perturbed Plexcitons matrix diagonalized"
         !> Print Energies and eigenstates.            
         if (myrank.eq.0) call out_QM_coupling
       else
         !> If propagate, transform Plexcitons integrals in plexciton basis
         call transform_plexd
         ! call do_Hqm_int(f(:))
         if (myrank.eq.0) write(6,*) "No QM propagation implemented"
       endif
       !> Deallocate matrices                                   
       call fin_QM_coupling 
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>    @brief Allocates and initializes matrices.
!>    @date Created   : S.Pipolo 02 May 2017
!>    @param[in] Hqm_dim
!>    @param[in,out] Hqm,Hqm_evt,Hqm_evl,occ
!----------------------------------------------------------------------------
      subroutine init_QM_coupling                         
       call do_BEM_quant
       FQBEM="diag-dip"
       if(FQBEM(6:8)=='all') then !< couple with all, but singly-occupied, modes.
         nmodes=nts_act
       else !< couple with the first "qmodes" singly-occupied modes. At present qmodes=1
!        nmodes=qmodes
         nmodes=4     
       endif
       allocate(g(nmodes,n_ci,n_ci))
       allocate(we(nmodes))
       allocate(omega_p(nmodes))
       allocate(qg(nmodes,nts_act))
       Hqm_dim=n_ci*(nmodes+1)
       allocate(occ(nmodes))
       occ=1.d0 !< all singly-occupied modes
       allocate(Hqm(Hqm_dim,Hqm_dim))
       allocate(Hqm_int(Hqm_dim,Hqm_dim))
       allocate(plexd(3,Hqm_dim,Hqm_dim))
       Hqm(:,:)=0.d0
       Hqm_int(:,:)=0.d0
       plexd(:,:,:)=0.d0
       allocate(Hqm_evt(Hqm_dim,Hqm_dim))
       allocate(Hqm_evl(Hqm_dim))
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>     @brief Finalize routine of QM_coupling  
!>     @date Created   : S.Pipolo 02 May 2017
!>     @param Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
      subroutine fin_QM_coupling                         
       call deallocate_BEM_public
       deallocate(we,omega_p,g,qg)
       deallocate(Hqm,Hqm_evt,Hqm_evl)
       if(allocated(Hqm_int)) deallocate(Hqm_int)
       deallocate(occ)
      return
      end subroutine
!     
!
!------------------------------------------------------------------------
!>     @brief Build Plexciton Hamiltonian: \f$ \mathcal{H}_{\text{M}}+\mathcal{H}_{\text{P}}+\mathcal{H}_{\text{MP}} \f$
!>     @date Created : 02 May 2017
!>     @author S.Pipolo 
!>     @note One mode coupled at a time
!>     @param 
!----------------------------------------------------------------------------
      subroutine do_Hqm_matrix  
       integer(4)::i,j,k,p,s !< indices    
       !> G0-G0 \f$ \mathbf{H}_{\text{M}} \f$ block: state energies, diagonal 
       do j=1,n_ci
         Hqm(j,j)=Hqm(j,j)+e_ci(j)
       enddo
       ! E1-E1,E1-G0
       do i=2,nmodes   
         do j=1,n_ci
           p=(i-1)*n_ci+j
           do k=j,n_ci
             s=(i-1)*n_ci+k
             !> E1-E1 \f$ \mathbf{H}_{\text{M}} \f$ off-diagonal subblocks
             Hqm(s,p)=Hqm(k,j)
             Hqm(p,s)=Hqm(s,p)
             !> G0-E1 \f$ \mathbf{H}_{\text{MP}} \f$ off-diagonal subblocks
             Hqm(k,p)=g(i,k,j)
             Hqm(j,s)=Hqm(k,p)
             Hqm(s,j)=Hqm(j,s)
             Hqm(p,k)=Hqm(s,j)
           enddo
           !> E1-E1 \f$ \mathcal{H}_{\text{P}} energies, diagonal 
           Hqm(p,p)=Hqm(p,p)+omega_p(i)*occ(i)    
         enddo
       enddo
      return
      end subroutine
!     
!
!------------------------------------------------------------------------
!>     @brief Build field-perturbation to Plexciton Hamiltonian: \f$ \mathcal{H}_{\text{MF}}+\mathcal{H}_{\text{PF}} \f$
!>     @date Created: 02 May 2017
!>     @author S.Pipolo 
!>     @note One mode coupled at a time
!>     @param 
!----------------------------------------------------------------------------
      subroutine do_plexd_matrix  
       integer(4)::i,j,k,p,s !< indices    
       real(dbl), allocatable:: gF(:) !< semiclassical particle-field couplings
       allocate(gF(3)) 
       !> Building \f$ \mathcal{H}_{\text{MF}} \f$ block
       plexd=mut
       do i=2,nmodes   
         !> Building \f$ \mathcal{H}_{\text{PF}} \f$ couplings
         gF(1)=-dot_product(BEM_Modes(i,:),cts_act(:)%x)*we(i)
         gF(2)=-dot_product(BEM_Modes(i,:),cts_act(:)%y)*we(i)
         gF(3)=-dot_product(BEM_Modes(i,:),cts_act(:)%z)*we(i)
         do j=1,n_ci
           p=(i-1)*n_ci+j
           do k=j,n_ci
             s=(i-1)*n_ci+k
             !> Adding \f$ \mathcal{H}_{\text{MF}} \f$ subblocks
             plexd(:,s,p)=plexd(:,k,j)
             plexd(:,p,s)=plexd(:,s,p)
           enddo
           !> Adding \f$ \mathcal{H}_{\text{PF}} \f$ subblocks
           plexd(:,p,j)=plexd(:,p,j)+gF(:)
           plexd(:,j,p)=plexd(:,j,p)+gF(:)
         enddo
       enddo
       deallocate(gF) 
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>    @brief transform plexd in Plexciton basis 
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!----------------------------------------------------------------------------
      subroutine transform_plexd 
       real(dbl), allocatable:: scr(:,:) 
       integer(4)::i
       allocate(scr(Hqm_dim,Hqm_dim))
       do i=1,3 
         scr(:,:)=matmul(transpose(Hqm_evt),plexd(i,:,:))
         plexd(i,:,:)=matmul(scr,Hqm_evt)                   
       enddo
       deallocate(scr)
      return
      end subroutine
!
!
!------------------------------------------------------------------------
!>    @brief build semiclassical interaction matrix
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!----------------------------------------------------------------------------
      subroutine do_Hqm_int(f)
       real(dbl), dimension(3), intent(in)  :: f  
       integer(4)::j,k
       do j=1,n_ci
         do k=1,n_ci
           Hqm_int(k,j)=-dot_product(plexd(:,k,j),f(:))
         enddo
       enddo
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
       open(7,file="Hqm.mat",status="unknown")
       open(8,file="Hqm.ene",status="unknown")
       !open(9,file="gCharges.mat",status="unknown")
       write(8,*) "Energies: "
       write(my_fmt,'(a,i0,a)') "(",Hqm_dim,"F10.6)"
       write(7,*) "Quantum-matrix: ", my_fmt
       do i=1,Hqm_dim
         write(7,my_fmt) (Hqm(i,j), j=1,Hqm_dim)
         if(i.le.n_ci) then
           write(8,"(i0,3F10.6)") i,Hqm_evl(i),e_ci(i),sqrt(BEM_W2(i))
         else
           write(8,"(i0,F10.6)") i, Hqm_evl(i)
         endif
       enddo
       !write(my_fmt,'(a,i0,a)') "(",nmodes,"E10.6)"
       !write(9,*) "# Nmodes = ",nmodes,"   Size = ", nts_act
       !do j=1,nts_act
       !  write(9,my_fmt) (qg(i,j), i=1,nmodes)
       !enddo
       close(7)
       close(8) 
       !close(9) 
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
       omega_p(1)=0.
       we(1)=0.
       do i=2,nmodes   
         omega_p(i)=sqrt(BEM_W2(i)) 
         we(i)=sqrt((omega_p(i)**2-eps_w0**2)/(two*omega_p(i)))
         qg(i,:)=BEM_Modes(i,:)*we(i)
         do j=1,n_ci
           do k=j,n_ci
             g(i,k,j)=dot_product(qg(i,:),vts(:,k,j))
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
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine test_QM_coupling                         
       character(len=32) :: my_fmt
       real(dbl):: r,d,mud,wl                
       integer(i4b) :: i,j,k   
       real(dbl), allocatable :: sp(:)               
       real(dbl), allocatable :: tot(:,:),ref(:,:)               
#ifndef MPI
       myrank=0
#endif
       allocate(sp(3),tot(n_ci,n_ci),ref(n_ci,n_ci))
       open(7,file="g.mat",status="unknown")
       write(7,*) "# Test for dipolar-mode couplings" 
       d=sqrt(sfe_act(1)%x**2+sfe_act(1)%y**2+sfe_act(1)%z**2)
       r=cts_act(1)%rsfe
       wl=sqrt(eps_A/3)
       sp(1)=sfe_act(1)%x 
       sp(2)=sfe_act(1)%y 
       sp(3)=sfe_act(1)%z 
       write(7,*) "# Sphere radius distance (bohr) and position"
       write(7,"(5F10.4)") cts_act(1)%rsfe,d,sp(1),sp(2),sp(3)
       write(7,*) "# g=dot_product(BEM_Modes(p,:),vts(:,i,j))*we" 
       write(7,*) "# g_ref=mu*sqrt(2*omega_p*r^3)/(d^3)"
       write(7,*) "#" 
       write(7,*) "#  p    i    j            g                 g_ref" 
       tot=zero
       ref=zero
       do i=2,4        
         do j=1,n_ci
           do k=j+1,n_ci
             mud=dot_product(mut(:,k,j),sp(:))/d
             tot(k,j)=tot(k,j)+g(i,k,j)*g(i,k,j)
             !ref: Garcia-Vidal PRL 112, 253601 (2014)
             !ref(k,j)=mud*sqrt(2*wl*r**3)/(d**3)
             ref(k,j)=2*mud*mud*wl*r**3/(d**6)
             !ref(k,j)=2*mud*mud*wl*r**3/(d+r)**6
           enddo
         enddo
       enddo
       do j=1,n_ci   
         do k=j+1,n_ci
           !write(7,"(3i5,3E20.12)") 2,j-1,k-1,sqrt(tot(k,j)),ref(k,j)
           write(7,"(3i5,3E20.12)") 2,j-1,k-1,sqrt(tot(k,j)),sqrt(ref(k,j))
         enddo
       enddo
       write(7,*) ""
       write(7,*) "# Dipolar resonance frequancy (a.u.)"
       write(7,*) "#  p          omega_p            sqrt(A/3)" 
       do i=2,4        
         write(7,"(i5, 2E20.12)")i, sqrt(BEM_W2(i)), sqrt(eps_A/3)
       enddo
       close(7)
       if (myrank.eq.0) then 
          write(6,*) "Test for dipolar-mode couplings...DONE" 
          write(6,*) "  Results in the g.mat file. " 
       endif
       open(9,file="gCharges.mat",status="unknown")
       write(my_fmt,'(a,i0,a)') "(",nmodes,"E20.6)"
       write(9,*) "# Nmodes = ",nmodes,"   Size = ", nts_act, my_fmt
       do j=1,nts_act
         write(9,my_fmt) (qg(i,j), i=1,nmodes)
       enddo
       close(9) 
       deallocate(sp,tot,ref)
       stop
      return
      end subroutine
!
!
      end module
