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
      use interface_classic
      use readio       
      use, intrinsic :: iso_c_binding
#ifdef OMP
      use omp_lib
#endif
#ifdef MPI
      use mpi
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
      integer(i4b) :: Hqm_dim  !<Dimension of \mathcal{H}_{\text{QM}} matrix 

      save
      private
      public do_QM_coupling
!
!
      contains
!
!------------------------------------------------------------------------
!>    @brief Driver routine of QM_coupling. 
!>    @date Created: 02 May 2017 
!>    @author S.Pipolo
!>    @note this routine should change to accomodate quantum external fields
!----------------------------------------------------------------------------
      subroutine do_QM_coupling(nstates_qm,energies_qm,trans_dipoles_qm,f0)
       implicit none 
       real(dbl), intent(in)    :: f0(3)     !< field at time0  
       integer(i4b), intent(in) :: nstates_qm  
       real(dbl), dimension(nstates_qm), intent(out) :: energies_qm
       real(dbl), dimension(3,nstates_qm,nstates_qm), intent(out)  :: &
                                                       trans_dipoles_qm 
#ifndef MPI
       myrank=0
#endif
       Hqm_dim=nstates_qm
       call initialize_qc_interface
       !> Allocate and initialize matrices here
       call init
       if (myrank.eq.0) write(6,*) "QM_coupling correcty initialized"
       !> compute charges associated to each mode
       call do_plasmon_charges(nmodes,omega_p,we)
       !> compute plexcitons couplings g
       call do_qm_couplings(nmodes,n_ci,g)
       if(this_Ftest.eq."qmt") then
         !> Test: performs the dipolar test on the spherical couplings and exit
         if (myrank.eq.0) call test_QM_coupling
       else 
         !> Build Plexcitons matrix: do_Hqm_matrix 
         call do_Hqm_matrix
         if (myrank.eq.0) write(6,*) "Plexcitons matrix built"
         !> Prepare perturbation integrals 
         call do_plexd_matrix(nmodes,n_ci,plexd)
         if (myrank.eq.0) write(6,*) &
                "Plexciton dipole integrals built"
         !> If a field is present at time zero, diagonalize Perturbed Plexciton matrix in presence of a static field.
         if(mdl(f0).gt.0.) call do_H_int(Hqm,plexd,f0,Hqm_dim)
         Hqm_evt=Hqm
         call diag_mat_in_wavet(Hqm_evt,Hqm_evl,Hqm_dim)
         if (myrank.eq.0.and.mdl(f0).le.0.) write(6,*) &
                "Plexcitons matrix diagonalized"
         if (myrank.eq.0.and.mdl(f0).gt.0.) write(6,*) &
                "Perturbed Plexcitons matrix diagonalized"
         !> If propagate, transform Plexcitons integrals in plexciton basis (pexciton dipole integrals)
         call transform_plexd
         if (myrank.eq.0) write(6,*) &
          "Plexciton dipole integrals transformed in plexciton basis"
         !> Print Energies and eigenstates.            
         if (myrank.eq.0) call out_QM_coupling
         !> Print out transition dipole integrals in plexciton basis
         if (myrank.eq.0) call out_plexd
       endif
#ifdef MPI
       call mpi_finalize(ierr_mpi)
#endif
       energies_qm=Hqm_evl
       trans_dipoles_qm=plexd
       !> Deallocate matrices                                   
       call fin_QM_coupling 
      return
      end subroutine do_QM_coupling
!
!
!------------------------------------------------------------------------
!     @brief Init routine of QM_coupling  
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm_dim,Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
      subroutine init
      implicit none
       !FQBEM='prop' !enforces use of correct QM_coupling flag
       allocate(g(nmodes,n_ci,n_ci))
       allocate(we(nmodes))
       allocate(omega_p(nmodes))
       allocate(occ(nmodes))
       occ=1.d0
       allocate(Hqm(Hqm_dim,Hqm_dim))
       allocate(Hqm_int(Hqm_dim,Hqm_dim))
       allocate(plexd(3,Hqm_dim,Hqm_dim))
       Hqm(:,:)=0.d0
       Hqm_int(:,:)=0.d0
       plexd(:,:,:)=0.d0
       allocate(Hqm_evt(Hqm_dim,Hqm_dim))
       allocate(Hqm_evl(Hqm_dim))
       write(*,*) "QM initialization done"
      return
      end subroutine init
!
!
!------------------------------------------------------------------------
!     @brief Finalize routine of QM_coupling  
!     @date Created   : S.Pipolo 02 May 2017
!     Modified  :
!     @param Hqm,Hqm_evt,Hqm_evl
!----------------------------------------------------------------------------
      subroutine fin_QM_coupling

       implicit none
       call deallocate_BEM_public_in_wavet
       deallocate(we,omega_p,g)
       deallocate(Hqm,Hqm_evt,Hqm_evl)
       if(allocated(Hqm_int)) deallocate(Hqm_int)
       deallocate(occ)
      return
      end subroutine fin_QM_coupling
!     
!------------------------------------------------------------------------
!>     @brief Build Plexciton Hamiltonian: \f$ \mathcal{H}_{\text{M}}+\mathcal{H}_{\text{P}}+\mathcal{H}_{\text{MP}} \f$
!>     @date Created : 02 May 2017
!>     @author S.Pipolo 
!>     @note One mode coupled at a time
!>     @param 
!----------------------------------------------------------------------------
      subroutine do_Hqm_matrix  
      implicit none
       integer(4)::i,j,k,p,s !< indices    
       !> G0-G0 \f$ \mathbf{H}_{\text{M}} \f$ block: state energies, diagonal 
       do j=1,n_ci
         Hqm(j,j)=Hqm(j,j)+e_ci(j)
       enddo
       ! E1-E1,E1-G0
       do i=1,nmodes   
         do j=1,n_ci
           p=i*n_ci+j
           do k=j,n_ci
             s=i*n_ci+k
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
!
!------------------------------------------------------------------------
!>    @brief transform plexd in Plexciton basis 
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!----------------------------------------------------------------------------
      subroutine transform_plexd 
      implicit none
       real(dbl), allocatable:: scr(:,:) 
       real(dbl), allocatable:: scr1(:,:) 
       integer(4)::i
       allocate(scr(Hqm_dim,Hqm_dim))
       allocate(scr1(Hqm_dim,Hqm_dim))
       do i=1,3 
         scr1(:,:)=plexd(i,:,:)
         scr=matmul(transpose(Hqm_evt),scr1)
         scr1=matmul(scr,Hqm_evt)                   
         plexd(i,:,:)=scr1(:,:) 
       enddo
       deallocate(scr)
       deallocate(scr1)
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

       implicit none
       integer(i4b) :: i,j   
       character(len=32) :: my_fmt
       open(7,file="Hqm_matrix.dat",status="unknown")
       open(8,file="Hqm_energies.dat",status="unknown")
       open(9,file="Hqm_vectors.dat",status="unknown")
       open(10,file="ci_mut_qc.inp",status="unknown",form="formatted")
       open(3,file="ci_energy_qc.inp",status="unknown",form="formatted")
       open(4,file="ci_ini_qc.inp",status="unknown",form="formatted")
       write(8,*) "Plexciton , molecule  and plasmon energies "
       write(my_fmt,'(a,i0,a)') "(",Hqm_dim,"F10.6)"
       write(7,*) "Super-matrix: ", my_fmt
       do i=1,Hqm_dim   
         write(7,my_fmt) (Hqm(i,j), j=1,Hqm_dim)
         if(i.le.n_ci.and.i.le.nmodes) then
           write(8,"(i0,3F10.6)") i,Hqm_evl(i),e_ci(i),omega_p(i)
         else if(i.le.n_ci.and.i.gt.nmodes) then
           write(8,"(i0,2F10.6,a)") i,Hqm_evl(i),e_ci(i), " - "
         else if(i.gt.n_ci.and.i.le.nmodes) then
           write(8,"(i0,F10.6,a,F10.6)") i,Hqm_evl(i)," - ",omega_p(i)
         else
           write(8,"(i0,F10.6,2a)") i, Hqm_evl(i), " - ", " - "
         endif
         ! ci_energy_qc.inp
         Hqm_evl(i)=Hqm_evl(i)-Hqm_evl(1)
         if(i.gt.1) write(7,'(A,I6,X,A,f15.8)') 'Root',i-1,':',Hqm_evl(i)/ev_to_au
         ! ci_ini_qc.inp
         if(i.eq.1) write(7,'(f15.8)') one
         if(i.gt.1) write(7,'(f15.8)') zero
         ! ci_mut_qc.inp 1st part
         write(7,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0, 'and',i-1, &
                & plexd(1,1,i),plexd(2,1,i),plexd(3,1,i)
         ! matrix
         write(9,my_fmt) (Hqm_evt(i,j), j=1,Hqm_dim)
       enddo
       ! ci_mut_qc.inp 2nd part
       do i=2,Hqm_dim
         do j=2,i   
           write(7,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1, & 
              & 'and',i-1,plexd(1,i,j),plexd(2,i,j),plexd(3,i,j)
         enddo
       enddo
       Hqm_evl(1)=zero                 
       close(7)
       close(8) 
       close(9) 
       close(10) 
       close(3) 
       close(4) 
       !call out_gcharges_in_wavet
      return
      end subroutine out_QM_coupling
!
!
!------------------------------------------------------------------------
!>    @brief Writes out dipole integrals in plexciton basis
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!>    @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine out_plexd

       implicit none
       integer(i4b) :: i,j   
       character(len=32) :: my_fmt
       open(7,file="ci_mut_plex.dat",status="unknown")
       do i=1,Hqm_dim
         write(7,'(a18,i4,3f15.6)') "States     0  and  ",i-1,plexd(1,1,i),plexd(2,1,i),plexd(3,1,i)
       enddo
       do i=2,Hqm_dim
         do j=2,i   
           write(7,'(a8,i4,a6,i4,3f15.6)') "States  ",i-1,"  and ",j-1,plexd(1,i,j),plexd(2,i,j),plexd(3,i,j)
         enddo
       enddo
       close(7)
      return
      end subroutine out_plexd
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
      implicit none
       real(dbl):: r,d,mud,wl,tmp                
       character(54):: my_fmt
       integer(i4b) :: i,j,k,pmax,p
       real(dbl), allocatable :: sp(:)               
       real(dbl), allocatable :: tot(:,:),ref(:,:)               
#ifndef MPI
       myrank=0
#endif
       allocate(sp(3),tot(n_ci,n_ci),ref(n_ci,n_ci))
       open(7,file="g.mat",status="unknown")
       write(7,*) "# Test for dipolar-mode couplings" 
       call grep_sphere_parameters(r,d,sp,wl)
       write(7,*) "# Sphere radius distance (bohr) and position"
       write(7,"(5F15.4)") r,d,sp(1),sp(2),sp(3)
       write(7,*) "# g=dot_product(BEM_Modes(p,:),vts(:,i,j))*we" 
       write(7,*) "# g_ref=mu*sqrt(2*omega_p*r^3)/(d^3)"
       write(7,*) "#" 
       write(7,*) "#  p    i    j            g                 g_ref" 
       tot=zero
       ref=zero
       do i=1,nmodes   
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
           write(7,"(3i5,3E20.12)") 2,j-1,k-1,sqrt(tot(k,j)),sqrt(ref(k,j))
         enddo
       enddo
       write(7,*) ""
       write(7,*) "# Dipolar resonance frequency (a.u.)"
       write(7,*) "#  p          omega_p            sqrt(A/3)" 
       do i=2,4        
         write(7,"(i5, 2E20.12)")i, omega_p2(i), wl 
       enddo
       close(7)
       if (myrank.eq.0) then 
          write(6,*) "Test for dipolar-mode couplings...DONE" 
          write(6,*) "  Results in the g.mat file. " 
       endif
       deallocate(sp,tot,ref)
       stop
      return
      end subroutine

      end module
