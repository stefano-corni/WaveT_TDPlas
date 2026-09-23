module interface_classic
      use constants
      use WTMathTools
      use readio
#ifdef TDPLAS
      use tdplas, only: set_qorf,set_qorf_pot,global_prop_Fmdm_relax,&
! used by dissipation
                        get_mdm_dip,get_gneq,init_mdm,prop_mdm,finalize_mdm,&
                        init_after_scf,set_charges,& 
! used by propagate
                        readio_and_init_tdplas_for_wt,&
! used by main and main_spectra
                        mpibcast_readio_mdm,quantum_init,global_qmodes_Fmop,global_qmodes_nprint,&
! used by main
                        global_sys_Fwrite,fr_0,BEM_Q0,mat_f0,global_prop_max_cycles,&
                        global_prop_threshold,global_prop_mix_coef,&
                        do_field_from_charges,out_gcharges,get_qr_fr,&
! used in scf
                        BEM_W2,global_sys_Ftest,drudel_eps_w0,BEM_Modes,pedra_surf_spheres,drudel_eps_A,&
                        do_BEM_quant,deallocate_BEM_public,global_qmodes_nmodes,&
                        global_qmodes_qmmodes,&
! used by QM_coupling 
                        q0,pedra_surf_n_tessere,global_prop_Fprop,global_prop_Fint,&
                        pedra_surf_tessere,global_prop_Finit_int,pedra_surf_n_spheres,global_medium_Fbem,&
                        global_sys_Fdeb,do_charges_from_pot,global_prop_n_q,init_mdm_prop,get_qorf,&
                        get_qorf0,do_Rfield_from_dip
! used only here in interface_tdplas
                        
#endif
#ifdef MPI
      use mpi
#endif
#ifdef OMP
      use omp_lib 
#endif

      implicit none

      type tess_pcm_in_wavet
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: area
       real(dbl) :: n(3)
       real(dbl) :: rsfe
      end type
!
      type sfera_in_wavet
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: r
      end type

      character(flg) :: this_Fmdm_relax, this_Fprop, this_Fint, this_Fwrite, &
                        this_Ftest, this_Finit_int, this_Fbem, this_Fmop

      real(dbl), allocatable :: quantum_vts(:,:,:) !< medium contribution to the hamiltonian
      real(dbl), allocatable :: quantum_vtsn(:) !< medium contribution to the hamiltonian
      real(dbl), allocatable :: h_mdm(:,:) !< medium contribution to the hamiltonian
      real(dbl), allocatable :: h_mdm_0(:,:) !< medium contribution to the hamiltonian

      integer(i4b) :: nts, this_nesf_act,this_nprint,this_max_mod_todiag

      type(tess_pcm_in_wavet), target, allocatable :: cts(:)
      integer(i4b) :: this_ncycmax !< maximum number of SCF cycles
      integer(i4b), allocatable :: this_imod(:) !<modes to print 
      real(dbl) :: this_thrshld    !< SCF threshold on (i) eigenvalues 10^-global_prop_thrshld (ii) eigenvectors 10^-(global_prop_thrshld+2)
      real(dbl), allocatable :: this_BEM_Q0(:,:)
      real(dbl), allocatable :: this_q0(:)             !< reaction BEM charges
      real(dbl), allocatable :: this_qx0(:)            !< local BEM charges 
      real(dbl), allocatable :: q0_fqfm(:),qx0_fqfm(:) !< reaction and local fwfm charges  
      real(dbl), allocatable :: m0_fqfm(:,:),mx0_fqfm(:,:) !< reaction and local fwfm dipoles  
      real(dbl) :: this_fr0(3)                         !< Reaction field at time 0 defined with Finit_mdm, here because used in scf
      real(dbl) :: this_fx0(3)                         !< Reaction field at time 0 defined with Finit_mdm, here because used in scf
      real(dbl), allocatable :: this_mat_f0(:,:) !< Onsager's total matrices needed for scf, free_energy and propagation
      real(dbl) :: this_mix_coef   !< SCF mixing ratio of old (1-global_prop_mix_coef) and new (global_prop_mix_coef) charges/field       
      real(dbl) :: this_eps_w0
      integer(i4b), allocatable :: qmmodes(:)
      integer(i4b) :: nmodes
! Atomistic medium
      integer(i4b) :: n_atoms
      real(dbl), allocatable :: r_atoms(:,:)
      real(dbl), allocatable :: vint_atoms(:,:,:)
      real(dbl), allocatable :: fint_atoms(:,:,:,:)
! Quantum plasmons
      real(dbl), allocatable :: qg(:,:)      !<Charges associated to each mode    
      real(dbl), allocatable :: wwe(:)        !<Plasmon coupling energy    

      public this_Fmdm_relax, &
! used by dissipation
             set_q0charges,get_medium_dip,get_energies,init_environment,prop_medium,finalize_medium,this_Finit_int,this_Fprop,&
             init_after_scf_in_wavet,init_env_prop,& ! used bypropagate (also this_mix_coef)
             read_medium_input,&
! used by main and main_spectra
             mpibcast_read_medium,set_global_tdplas_in_wavet,&
! used by main
             this_Fwrite,this_fr0,this_BEM_Q0,this_mat_f0,this_ncycmax,this_thrshld,quantum_vtsn,&
             this_mix_coef,&
             do_field_from_charges_in_wavet, nts, &
! used in scf
             update_environment_scf,out_environment_scf,&
             this_Ftest,this_eps_w0,&
             initialize_qc_interface,&
             this_Fmop,this_imod,this_nprint,this_max_mod_todiag,deallocate_BEM_public_in_wavet,&
             qmmodes,nmodes,get_m_or_v, do_plasmon_charges,do_qm_couplings
! used by QM_coupling 
! module variables this_q0 and this_fr0 contain the updated values of reaction cherges and filed during an SCF run
      contains
  
      ! begin - wrapper subroutines
!------------------------------------------------------------------------
! @brief Bridge subroutine to set charges qr_t to q0 during propagation 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine set_q0charges
        implicit none
#ifdef TDPLAS
        call set_charges(this_q0)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        return
      end subroutine set_q0charges

!------------------------------------------------------------------------
! @brief Set the dipole(t) in Sdip for spectra 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine get_medium_dip(mdm_dip)
        implicit none
        real(dbl), intent(inout) :: mdm_dip(3)
#ifdef TDPLAS
        call get_mdm_dip(mdm_dip)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        return
      end subroutine get_medium_dip
     
 
!------------------------------------------------------------------------
! @brief Read medium input 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine read_medium_input
        implicit none
        integer :: ii,shap(3), nthr
#ifdef TDPLAS
#ifdef OMP
        nthr=omp_get_max_threads()
#endif
        call readio_and_init_tdplas_for_wt(nthr)
        ! SP 18/05/20: added the following line to compute potentials from dipoles for the reaction field test
        !              better using global_sys_Fdeb
        this_Fprop=global_prop_Fprop
        this_Fwrite=global_sys_Fwrite
        this_Finit_int=global_prop_Finit_int
        this_Fmdm_relax = global_prop_Fmdm_relax
        this_Fint=global_prop_Fint
        this_Ftest=global_sys_Ftest
        nts=pedra_surf_n_tessere
        call read_gau_out_medium
        if(global_sys_Fdeb.eq."vmu") then
           call get_vts_from_dip
        elseif(this_Fprop.eq."chr-ief".or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") then 
           ! SP 17/05/20 shape is not a standard f90 function
           shap=shape(quantum_vts)
           if(shap(1).ne.nts) then
              write(6,*) "Error: the number of tesserae for the potential is different than those in the cavity/NP"
              write(6,*) shap(1)," vs ",nts
              write(6,*) "This is usually due to incoerent ci_pot.inp and cavity.inp files. I stop here" 
              stop
           endif
        end if
        this_ncycmax=global_prop_max_cycles
        this_thrshld=global_prop_threshold
        this_mix_coef=global_prop_mix_coef 
        this_eps_w0=drudel_eps_w0
        
        this_nesf_act=pedra_surf_n_spheres
        allocate(cts(nts))
        !cts=pedra_surf_tessere
        do ii=1, nts
         cts(ii)%x=pedra_surf_tessere(ii)%x
         cts(ii)%y=pedra_surf_tessere(ii)%y
         cts(ii)%z=pedra_surf_tessere(ii)%z
         cts(ii)%rsfe=pedra_surf_tessere(ii)%rsfe
         cts(ii)%n(:)=pedra_surf_tessere(ii)%n(:)
        end do
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        return
      end subroutine read_medium_input
      
      
!------------------------------------------------------------------------
! @brief Get energies 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)

        implicit none
        real(dbl), intent(inout) :: e_vac,g_neq_t,g_neq2_t,g_eq_t
#ifdef TDPLAS
        call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        return
      end subroutine get_energies
      
      
!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine init_environment(c,f)
        implicit none
        complex(cmp), intent(in)    :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in)    :: f(:)    !< (1:3)              - external field
        real(dbl), allocatable      :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl), allocatable      :: r(:)    !< auxiliary 3d vector array
        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external potential
        real(dbl), allocatable      :: h0(:,:) !< (1:n_ci,1:n_ci)     initial hamiltonian (deallocated)
        !> Allocating the interaction matrix h_mdm 
        allocate(h_mdm(n_ci,n_ci),h_mdm_0(n_ci,n_ci))
        allocate(h0(n_ci,n_ci))
        allocate(mu(3))
        h0=zero
        h_mdm=zero
        h_mdm_0=zero
        call do_dip_from_coeff(c,n_ci,mu,mut)
#ifdef TDPLAS
        !> Prepare for interaction with continuum medium
        if(this_Fprop.eq."dip") then
          !> build in principle the field coming from all dipoles
          !> initializing medium with molecular dipole and external field
          call init_mdm(mu_t = mu, f_tp = f, morv=mut(:,1,1))
          allocate(this_mat_f0(nts,nts))
          this_mat_f0=mat_f0
          call get_qorf(this_fr0)
        else
          !> build the potential on the surface from all charges and dipoles
          allocate(pot(nts),potf(nts))
          allocate(this_q0(nts),this_qx0(nts))
          pot(:)=zero
          potf(:)=zero
          call do_pot_tdplas(c,mu,f,pot,potf)
          call init_mdm(pot_t = pot, potf_t = potf, morv=quantum_vts(:,1,1))
          !SP07/05/26 initial charges set in tdplas
          call get_qorf(this_q0)
          deallocate(pot)
          deallocate(potf)
        end if
#endif
#ifdef fqfm  
        !> Prepare for interaction with discrete medium
        allocate(pot(n_atoms),potf(n_atoms))
        pot(:)=zero
        potf(:)=zero
        call do_pot_fld_fqfm(c,mu,f,pot,fld,potf)
        call init_fqfm(pot,fld,potf)
#endif
        !> update the OUT hamiltonian                 
        !> build the h_mdm interaction hamiltonian
        call do_interaction(h0)
        h_mdm_0=h_mdm
        deallocate(h0)
        return
      end subroutine init_environment


!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine init_env_prop(c,mu,f,h)
        implicit none
        complex(cmp), intent(in) :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(:)    !< (1:3)              - external field
        real(dbl)   , intent(inout) :: h(:,:)  !< (1:n_ci,1:n_ci) - interaction hamiltonian
        real(dbl), allocatable      :: r(:)  !< auxiliary 3d vector array
        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external potential
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
          !> build in principle the field coming from all dipoles
          !> initializing medium with molecular dipole and external field
          call init_mdm_prop(mu_t = mu, f_tp = f)
        else
          !> build the potential on the surface from all charges and dipoles
          allocate(pot(nts),potf(nts))
          pot(:)=zero
          potf(:)=zero
          call do_pot_tdplas(c,mu,f,pot,potf)
          call init_mdm_prop(pot_t = pot, potf_t = potf)
          deallocate(pot)
          deallocate(potf)
        end if
        !> build the h_mdm interaction hamiltonian
        call do_int_tdplas
#endif
#ifdef fqfm  
        allocate(pot(n_atoms),potf(n_atoms))
        pot(:)=zero
        potf(:)=zero
        call do_pot_fld_fqfm(c,mu,f,pot,fld,potf)
        call init_fqfm_prop(pot,fld,potf)
        ! WARNING: should we compute again the interaction???
        !> build the h_mdm interaction hamiltonian
        call do_int_fqfm
#endif
        !> update the OUT hamiltonian                 
        h=h+h_mdm
        !WARNING SP 15/05/24: CHECK THIS for restart, it was in propagate, is it really needed?
        !if (Fres.eq.'Nonr') then
        !   i=1
        !   call prop_medium(i,c,mu,f,h)
        !endif

        return
      end subroutine init_env_prop




!------------------------------------------------------------------------
! @brief Compute the potential at the BEM points from different sources:
!         - an external field f 
!         - the molecular charge density (coefficients c or dipome mu)
!         - wfqfm charges and dipoles
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_pot_tdplas(c,mu,f,pot,potf)
        implicit none
        complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)        - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(3)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(3)    !< (1:3)              - external field
        real(dbl)   , intent(OUT):: pot(nts)  !< (1:nts)   - potential on continuum surface
        real(dbl)   , intent(OUT):: potf(nts)  !< (1:nts)   - potential on continuum surface
        real(dbl), allocatable   :: r(:,:)  !< auxiliary 3d vector array
        allocate(r(3,nts))
        r(1,:)=cts(:)%x
        r(2,:)=cts(:)%y
        r(3,:)=cts(:)%z
        !> prepare the potetial acting on the medium for initialisation
        if(this_Fint.eq."ons") then
         !> computing molecular potential from its dipole and not the coefficients                     
         call do_pot_from_dip(1,mol_cc,mu,nts,r,pot)
        else
         !> computing molecular potential from the coefficients
         call do_pot_from_coeff(c,n_ci,nts,quantum_vts,pot)
        end if
        !> computing external potential in the long-wavelength limit
        call do_pot_from_field(f,nts,r,potf)
#ifdef fqfw 
        ! SP 6/5/2024: in the following routines pot is updated so that it contains the total pontential
        ! of the molecule charge density and the fwfm system that determines the apparent
        ! surface charges
        call do_pot_from_charges(n_atoms,r_atoms,q_disc,nts,r,pot)
        call do_pot_from_dip(n_atoms,r_atoms,m_disc,nts,r,pot)
#endif  
        deallocate(r)
      end subroutine do_pot_tdplas


!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_pot_fld_fqfm(c,mu,f,pot,fld,potf)
        implicit none
        complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(3)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(3)    !< (1:3)              - external field
        real(dbl)   , intent(OUT):: pot(n_atoms)  !< (1:nts)   - potential on continuum surface
        real(dbl)   , intent(OUT):: fld(3,n_atoms)  !< (1:nts)   - potential on continuum surface
        real(dbl)   , intent(OUT):: potf(3,n_atoms)  !< (1:nts)   - potential on continuum surface
        real(dbl), allocatable      :: r(:,:)  !< auxiliary 3d vector array
        real(dbl), allocatable      :: qorf(:) !< (1:pedra_surf_n_tessere)     - charges or field  
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
          allocate(qorf(3))
          stop "Coupling with fqfw not implemented for dipole propagation"
        else
          allocate(qorf(nts))
          allocate(r(3,nts))
          r(1,:)=cts(:)%x
          r(2,:)=cts(:)%y
          r(3,:)=cts(:)%z
          !> get charges from td_contmed               
          call get_qorf(qorf)
          !> compute potential and field from          
          call do_pot_from_charges(nts,r,qorf,n_atoms,r_atoms,pot)
          call do_fld_from_charges(nts,r,qorf,n_atoms,r_atoms,fld)
        end if
        deallocate(r)
#endif  
        !> prepare the potetial acting on the medium for initialisation
        ! WARNING CHECK THIS the following Flag should be defined for fqfm
        if(this_Fint.eq."ons") then
          !> computing molecular potential nad field from its dipole and not the coefficients                     
          call do_pot_from_dip(1,mol_cc,mu,n_atoms,r_atoms,pot)
          call do_fld_from_dip(1,mol_cc,mu,n_atoms,r_atoms,fld)
        else
          !> computing molecular potential from the coefficients
          call do_pot_from_coeff(c,n_ci,n_atoms,vint_atoms,pot)
          call do_fld_from_coeff(c,n_ci,n_atoms,fint_atoms,fld)
        end if
        !> computing external potential in the long-wavelength limit
        call do_pot_from_field(f,n_atoms,r_atoms,potf)

      end subroutine do_pot_fld_fqfm






!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_interaction(h)
        implicit none
        real(dbl), intent(INOUT) :: h(n_ci,n_ci) !<   

        h_mdm=zero
        !> compute that interaction with a continuum medium
#ifdef TDPLAS
        call do_int_tdplas
#endif
        !> compute that interaction with a discrete medium
#ifdef fqfw 
        call do_int_fqfm
#endif
        !> add the interaction term to the input hamiltonian
        h=h+h_mdm

        return
      end subroutine do_interaction


!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_int_tdplas
        implicit none
        real(dbl), allocatable      :: qorf(:) !< (1:pedra_surf_n_tessere)     - charges or field  
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
          !> allocate and initialise arrays and prepare for calls
          allocate(qorf(3))
        else
          !> allocate and initialise arrays and prepare for calls
          allocate(qorf(nts))
        end if
        !> get charges from td_contmed               
        call get_qorf(qorf)
        !> construct the interaction hamiltonian h_mdm
        call do_interaction_cont(qorf,h_mdm)
        deallocate(qorf)
        return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine do_int_tdplas


!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_int_fqfm  
        implicit none
        real(dbl), allocatable      :: q(:)   !< (1:N_atoms)     - fq  
        real(dbl), allocatable      :: m(:,:) !< (1:N_atoms)     - fm   
        allocate(q(n_atoms))
        allocate(m(3,n_atoms))
        !> get charges from fqfm               
        call get_q_and_m(q,m)
        !> construct the interaction hamiltonian h_mdm
        call do_interaction_discr(q,m,h_mdm)
        deallocate(q,m)
        return
      end subroutine do_int_fqfm 

      

!------------------------------------------------------------------------
! @brief Prepare medium for scf           
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
     subroutine prepare_mdm_for_scf
        implicit none
          allocate(this_BEM_Q0(nts,nts))
#ifdef TDPlas
          this_BEM_Q0=BEM_Q0
        return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
      end subroutine prepare_mdm_for_scf



!------------------------------------------------------------------------
! @brief update medium dof (q,f,etc) during a scf cycle
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
      subroutine update_environment_scf(c,f)
        implicit none
        complex(cmp), intent(in) :: c(n_ci) !> (1:n_ci)   - molecular wavefunction coefficients
        real(dbl), intent(in)    :: f(3)    !> (1:3)      - external field                     
        real(dbl), allocatable   :: mu(:)   !>            - molecular dipole 
        integer(i4b)::i    
        allocate(mu(3))
        mu=zero
#ifdef TDPLAS
        if((this_Fprop.eq."dip").or.(this_Fint.eq."ons")) then 
          call do_dip_from_coeff(c,n_ci,mu,mut)
        endif
        if(this_Fprop.eq."dip") then 
          call update_BEM_field(c,mu,f)
        else 
          call update_BEM_charges(c,mu,f)
        endif
#endif
#ifdef fqfm
        call update_fqfm_char_and_dip(c,mu,f)
#endif
        deallocate(mu)
        return
      end subroutine update_environment_scf

!------------------------------------------------------------------------
! @brief transform potentials            
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
      subroutine transform_environment_scf(c)
        implicit none
        real(dbl), intent(in) :: c(n_ci,n_ci)    !< (1:n_ci)           - molecular wavefunction coefficients
        integer(i4b)              :: its  

        if(this_Fprop.ne."dip") then 
        ! transform vts to the SCF state basis
         do its=1,nts
          quantum_vts(its,:,:)=matmul(quantum_vts(its,:,:),c)
          quantum_vts(its,:,:)=matmul(transpose(c),quantum_vts(its,:,:))
         enddo
        endif
        return
      end subroutine transform_environment_scf


     
!------------------------------------------------------------------------
! @brief write out potentials            
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
      subroutine out_environment_scf
        implicit none

        if(this_Fprop.ne."dip") then 
         call out_charges(this_q0)
         call out_vts
       endif
        return
      end subroutine out_environment_scf


     
!------------------------------------------------------------------------
! @brief Compute field from dipole 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_BEM_field(c,mu,f)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)   - molecular wavefunction coefficients
       real(dbl), intent(in) :: mu(3)         !< (1:3)      - molecular dipole
       real(dbl), intent(in)  :: f(:)         !< (1:3)      - field                               
       real(dbl), allocatable :: mf(:)        !< (1:3)      - field                               
#ifdef TDPLAS 
      ! SP 180226 the following should stay in tdplas     
       allocate(mf(3))
       call set_qorf_pot(mu,f)
       call get_qorf(mf)
       this_fr0=(1.-this_mix_coef)*this_fr0+this_mix_coef*mf
       deallocate(mf)
       call set_qorf(this_fr0)
       return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine update_BEM_field

!------------------------------------------------------------------------
! @brief Updates the value of this_q0 and this_qx0 by solving the BEM 
!         equations
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_BEM_charges(c,mu,f)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci)  !< (1:n_ci)        - molecular wavefunction coefficients
       real(dbl), intent(in)    :: mu(3)    !< (1:3)           - molecular dipole                   
       real(dbl), intent(in)    :: f(3)     !< (1:3)           - external field                    
       real(dbl), allocatable   :: pot(:)   !< (1:pedra_surf_n_tessere)     - molecular potential
       real(dbl), allocatable   :: potf(:)  !< (1:pedra_surf_n_tessere)     - field potential
       real(dbl), allocatable   :: q(:)     !< (1:pedra_surf_n_tessere)     - charges
       integer(i4b)::i  
#ifdef TDPLAS  
       allocate(pot(nts))
       allocate(potf(nts))
       allocate(q(nts))
       ! With the following call all potential (mol+fqfm) is included
       call do_pot_tdplas(c,mu,f,pot,potf)
       ! Set charges in TDPlas using the potential in pot                
       call set_qorf_pot(pot,potf)
       ! Get charges from TDPlas 
       call get_qorf(q)
       this_q0=(1.-this_mix_coef)*this_q0+this_mix_coef*q
       deallocate(pot,potf,q) 
       ! SC 12/8/2016: apparently for NP, charge compensation is needed
       !if (Fmdm.eq.'cnan'.or.Fmdm.eq.'qnan') then
       !  this_q0=this_q0-sum(this_q0)/nts
       !endif
       call set_qorf(this_q0)
       return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine update_BEM_charges


!------------------------------------------------------------------------
! @brief Compute fqfm charges and dipoles from potential and field 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_fqfm_char_and_dip(c,mu,f)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci)  !< (1:n_ci)       - molecular wavefunction coefficients
       real(dbl)   , intent(in) :: mu(3)    !< (1:3)          - molecular dipole
       real(dbl)   , intent(in) :: f(3)     !< (1:3)          - external field
       real(dbl)   , allocatable:: pot(:)   !< (1:n_atoms)    - potential on fqfw atoms       
       real(dbl)   , allocatable:: potf(:)  !< (1:n_atoms)    - external potential fqfw atoms        
       real(dbl)   , allocatable:: fld(:,:) !< (3,1:n_atoms)  - field on fqfw atoms
       real(dbl)   , allocatable:: q(:)     !< (1:n_atoms)    - charges on fqfw atoms
       real(dbl)   , allocatable:: m(:,:)   !< (3,1:n_atoms)  - dipoles on fqfw atoms
#ifdef TDPLAS
       allocate(pot(n_atoms))
       allocate(potf(n_atoms))
       allocate(fld(3,n_atoms))
       allocate(q(n_atoms))
       allocate(m(3,n_atoms))
       call do_pot_fld_fqfm(c,mu,f,pot,fld,potf)
       call do_qandm_fqfm(pot,fld,q,m)
       q0_fqfm=(1.-this_mix_coef)*q0_fqfm+this_mix_coef*q
       m0_fqfm=(1.-this_mix_coef)*m0_fqfm+this_mix_coef*m
       call do_qandm_fqfm(potf,f,q,m)
       call do_charges_from_pot(potf,q)
       qx0_fqfm=(1.-this_mix_coef)*qx0_fqfm+this_mix_coef*q
       mx0_fqfm=(1.-this_mix_coef)*mx0_fqfm+this_mix_coef*m
       deallocate(pot,potf,q,fld,m) 
       return
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine update_fqfm_char_and_dip


!------------------------------------------------------------------------
! @brief Write out the charges in the charges0_scf.dat file 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_charges(q)
       implicit none
       real(dbl), intent(IN):: q(nts)     
       integer(i4b) its
#ifndef MPI
       myrank=0
#endif
       open(unit=7,file="charges0_scf.inp",status="unknown", &
            form="formatted")
         write (7,*) nts
         do its=1,nts
          write (7,'(E22.8,F22.10)') q(its)
         enddo
       close(unit=7)
       if (myrank.eq.0) write(6,*) "Written out the SCF charges"
       return 
      end subroutine out_charges     


!------------------------------------------------------------------------
! @brief Transform to the new basis and write out potential integrals on
! tesserae (vts) 
!
! @date Created: S. Pipolo
! Modified: L. Biancorosso
!------------------------------------------------------------------------
      subroutine out_vts
       implicit none
       integer(i4b) :: its,i,j
#ifndef MPI
       myrank=0
#endif
       open(unit=7,file="ci_pot_scf.inp",status="unknown", &
          form="formatted")
       write(7,*) nts
       i=0
       j=0
       ! V00
       write(7,*) i,j 
       do its=1,nts
        write(7,*) quantum_vts(its,1,1)-quantum_vtsn(its),0.d0,quantum_vtsn(its)
       enddo

       do i=2,n_ci
          write(7,*) 0, i-1
          do its=1,nts
             write(7,*) quantum_vts(its,1,i) 
          enddo
       enddo

       do i=2,n_ci
           do j=2,i
              write(7,*)  i-1, j-1
              do its=1,nts
                 if (i.eq.j) then
                     write(7,*) quantum_vts(its,i,j)-quantum_vtsn(its)
                 else
                     write(7,*) quantum_vts(its,i,j) 
                 endif
              enddo
           enddo
       enddo
       if (myrank.eq.0) write(6,*) "Written out the SCF potentials"

       return 

      end subroutine out_vts      




!------------------------------------------------------------------------
! @brief Propagate medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine prop_medium(i,c,mu,f,h)
#ifdef MPI
      use mpi
#endif
        implicit none
        complex(cmp), intent(in) :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(:)    !< (1:3)              - external field
        real(dbl)   , intent(inout) :: h(:,:)  !< (1:n_ci,1:n_ci) - interaction hamiltonian
        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external  potential
        real(dbl), allocatable      :: r(:,:)  !< (3,1:pedra_surf_n_tessere)     - external  potential
        integer(i4b), intent(in) :: i
         ! To be more efficient this if should go in the propagate of waveT
         ! Propagate medium only every global_prop_n_q timesteps
#ifdef TDPlas
          if(mod(i,global_prop_n_q).ne.0) then
            ! Build the interaction Hamiltonian Reaction/Local with previous charges
            ! Update the interaction Hamiltonian
            call do_interaction(h)
            return
          endif
          if(this_Fprop.eq."dip") then
           ! propagating medium with molecular dipole and external field
           ! Get charges from external codes 
           call prop_mdm(i, mu_t = mu, f_tp = f)
          else
           allocate(pot(nts))
           allocate(potf(nts))
           pot=zero
           potf=zero
           allocate(r(3,nts))
           r(1,:)=cts(:)%x
           r(2,:)=cts(:)%y
           r(3,:)=cts(:)%z
           if(this_Fint.eq."ons") then
            ! computing molecular potential corresponding to a point-like dipole
            !call do_pot_from_dip(mu,pot)
            call do_pot_from_dip(1,mol_cc,mu,nts,r,pot)
           else
            ! computing molecular potential
            call do_pot_from_coeff(c,n_ci,nts,quantum_vts,pot)
           end if
           ! computing external potential in the long-wavelength limit
           call do_pot_from_field(f,nts,r,potf)
           ! propagating medium with molecular and external potentials
           ! SP 15/05/20 changed this_Ftest with global_sys_Ftest
           if(global_sys_Ftest.eq."n-r".or.global_sys_Ftest.eq."s-r") then
            call prop_mdm(i, mu_t = mu, pot_t = pot, potf_t = potf)
           else
            call prop_mdm(i, pot_t = pot, potf_t = potf)
#ifdef MPI
            call mpi_finalize(ierr_mpi)
            !WARNING Why stopping here?
            stop
#endif      
           end if
            deallocate(r)
            deallocate(pot)
            deallocate(potf)
          end if
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        ! Build the interaction Hamiltonian 
        ! WARNING Shall we grep charges/field here??
        ! Update the interaction Hamiltonian
        call do_interaction(h)
        return
      end subroutine prop_medium
      
     
!------------------------------------------------------------------------
! @brief Finalize medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine finalize_medium

        implicit none

#ifdef TDPLAS
        call finalize_mdm
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine finalize_medium

!------------------------------------------------------------------------
! @brief Broadcast input medium if parallel 
!
! @date Created   : E. Coccia 9/5/18 
! Modified  :  
!------------------------------------------------------------------------
      subroutine mpibcast_read_medium

        implicit none

#ifdef TDPLAS
        call mpibcast_readio_mdm 
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine mpibcast_read_medium 

!------------------------------------------------------------------------
! @brief Interface between TDPlas and QM code 
!
! @date Created   : 
! Modified  :  
!------------------------------------------------------------------------
      subroutine set_global_tdplas_in_wavet(this_dt,this_mdm,this_mol_cc,this_n_ci,this_n_ci_read,this_c_i,this_e_ci,&
                                            this_fmax,this_omega,this_Ffld,this_n_out,this_n_f,this_tdelay,this_pshift,&
                                            this_Fbin,this_Fopt,this_res,this_n_res)

        implicit none

        real(dbl)     , intent(in) :: this_dt                         ! time step
        character(3)  , intent(in) :: this_mdm                        ! kind of medium
        integer(i4b)  , intent(in) :: this_n_ci,this_n_ci_read        ! number of CIS states
        real(dbl)     , intent(in) :: this_e_ci(:)                    ! CIS energies
        real(dbl)     , intent(in) :: this_mol_cc(3)                  ! molecule center
        real(dbl)     , intent(in) :: this_fmax(3,10),this_omega(10)  ! field amplitude and frequency
        real(dbl)     , intent(in) :: this_tdelay(10),this_pshift(10) ! time delay and phase shift
        complex(cmp)  , intent(in) :: this_c_i(:)                     ! CIS coefficients
        character(3)  , intent(in) :: this_Ffld                       ! shape of impulse
        character(3)  , intent(in) :: this_Fbin                       ! binary output
        character(3)  , intent(in) :: this_Fopt                       ! matrix/vector multiplication 
        integer(i4b)  , intent(in) :: this_n_out,this_n_f             ! auxiliaries for output
        character(1)  , intent(in) :: this_res                        ! restart for medium 
        integer(i4b)  , intent(in) :: this_n_res                      ! frequency for restart

#ifdef TDPLAS
        call quantum_init(this_dt,this_mol_cc,this_fmax,this_omega,this_n_out,this_n_f,&
                          this_Fbin,this_Fopt, this_n_res)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine set_global_tdplas_in_wavet

      subroutine diag_mat_in_wavet(M,E,Md)

       implicit none

       integer(i4b), intent(in) :: Md
       real(dbl), intent(inout) :: M(Md,Md)
       real(dbl), intent(out) :: E(Md)

#ifdef TDPLAS
       call diag_mat(M,E,Md)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

       return

      end subroutine diag_mat_in_wavet




      subroutine do_field_from_charges_in_wavet(q,f)

       implicit none

       real(dbl),intent(out):: f(3)  
       real(dbl),intent(in):: q(nts)  

#ifdef TDPLAS
       call do_field_from_charges(q,f)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine do_field_from_charges_in_wavet









      subroutine deallocate_BEM_public_in_wavet
       implicit none
#ifdef TDPLAS
       call deallocate_BEM_public
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine deallocate_BEM_public_in_wavet







      subroutine get_vts_from_dip
       implicit none
       real(dbl),allocatable :: pos(:,:)
#ifdef TDPLAS
       allocate(pos(nts,3))
       pos(:,1)=pedra_surf_tessere(:)%x
       pos(:,2)=pedra_surf_tessere(:)%y
       pos(:,3)=pedra_surf_tessere(:)%z
       call do_vts_from_dip(quantum_vts,pos,mut,mol_cc,nts,n_ci)
       deallocate(pos)
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
      end subroutine get_vts_from_dip

!
!
!------------------------------------------------------------------------
!>     @brief initializes and allocates the variables for qm coupling 
!>     @date Created: 25 Sep 2019
!>     @author S.Pipolo
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine initialize_qc_interface
       implicit none
#ifdef TDPLAS
       call do_BEM_quant
       write(6,*) "Done BEM quantum"
       nmodes=global_qmodes_nmodes
       if(global_sys_Ftest.eq."qmt") nmodes=3
       allocate(qmmodes(nmodes))
       allocate(qg(nmodes,nts))
       allocate(wwe(nmodes))
       if(global_sys_Ftest.eq."qmt") then 
         nmodes=3
         qmmodes(1)=2 
         qmmodes(2)=3 
         qmmodes(3)=4 
       else
         qmmodes(:)=global_qmodes_qmmodes(:)       
       endif
       this_nprint=global_qmodes_nprint
       !> Test: use potentials from dipoles    
       if(this_Ftest.eq."qmt") then
         if (allocated(quantum_vts)) deallocate(quantum_vts)
         allocate (quantum_vts(nts,n_ci,n_ci))
         call get_vts_from_dip
         if (myrank.eq.0) write(6,*) "Integrals from dipoles computed"
       endif
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
      end subroutine initialize_qc_interface
!
!
!------------------------------------------------------------------------
!>     @brief computes the molecule-environment quantum couplig elements "g"
!>     @date Created: 25 Sep 2019
!>     @author S.Pipolo
!>     @param Hqm_evl  
!----------------------------------------------------------------------------
      subroutine do_qm_couplings(nmodes,n,g)
      implicit none
       integer(i4b),intent(in) :: nmodes,n
       real(dbl),intent(out) :: g(nmodes,n,n)
       integer(i4b) :: i,j,k   
#ifndef MPI
       myrank=0
#endif
       do i=1,nmodes   
         do j=1,n
           do k=j,n
             g(i,k,j)=dot_product(qg(i,:),quantum_vts(:,k,j))
           enddo
         enddo
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
      subroutine do_plasmon_charges(nmodes,omega_p,we)  
      implicit none
       integer(i4b),intent(in) :: nmodes
       real(dbl),intent(out) :: omega_p(nmodes)
       real(dbl),intent(out) :: we(nmodes)
       integer(i4b) :: i   
#ifndef MPI
       myrank=0
#endif
#ifdef TDPlas
       do i=1,nmodes  
         omega_p(i)=sqrt(BEM_W2(qmmodes(i))) 
         we(i)=sqrt((omega_p(i)**2-this_eps_w0**2)/(two*omega_p(i)))
         qg(i,:)=BEM_Modes(qmmodes(i),:)*we(i)
       enddo
       call out_gcharges
       wwe=we
      return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
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
      subroutine do_plexd_matrix(nmodes,n,plexd)  
      implicit none
       integer(i4b),intent(in) :: n,nmodes
       real(dbl),intent(out) :: plexd(3,n*nmodes,n*nmodes)
       integer(4)::i,j,k,p,s !< indices    
       real(dbl), allocatable:: gF(:) !< semiclassical particle-field couplings

#ifdef TDPlas
       allocate(gF(3)) 
       !> Building \f$ \mathcal{H}_{\text{MF}} \f$ block
       do j=1,n
         do k=j,n
           plexd(:,k,j)=mut(:,k,j)
           plexd(:,j,k)=plexd(:,k,j)
         enddo
       enddo
       do i=1,nmodes   
         !> Building \f$ \mathcal{H}_{\text{PF}} \f$ couplings
         gF(1)=-dot_product(BEM_Modes(qmmodes(i),:),cts(:)%x)*wwe(i)
         gF(2)=-dot_product(BEM_Modes(qmmodes(i),:),cts(:)%y)*wwe(i)
         gF(3)=-dot_product(BEM_Modes(qmmodes(i),:),cts(:)%z)*wwe(i)
         do j=1,n
           p=i*n+j
           do k=j,n
             s=i*n+k
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
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 

      end subroutine





      subroutine init_after_scf_in_wavet(pot_or_mut)

       implicit none
       real(dbl), intent(in) :: pot_or_mut(:)

#ifdef TDPLAS
       call init_after_scf(pot_or_mut)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine init_after_scf_in_wavet



! begin - subroutines to calculate dipoles, field and potentials from coefficients, dipoles and fields



!------------------------------------------------------------------------
! @brief Computes the modulus of a vector
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      function mdl(v) result(m)

        real(dbl), dimension(:), intent(in) :: v
        real(dbl) :: m
        integer(i4b) :: i

        m=zero

        do i=1,size(v)
          m=m+v(i)*v(i)
        enddo

        m=sqrt(m)

      end function mdl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!   INTERACTION  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
!> Build the interaction matrix in the molecular unperturbed state basis
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Build the interaction matrix 
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_interaction_discr(q,m,h)
       real(dbl),intent(IN)   :: q(n_atoms)
       real(dbl),intent(IN)   :: m(3,n_atoms)
       real(dbl),intent(INOUT):: h(n_ci,n_ci)
       integer(i4b):: i,j,k

#ifndef MPI
       myrank=0
#endif
       do j=1,n_ci
         do i=1,j
           h(i,j)=h(i,j)+dot_product(q(:),vint_atoms(:,i,j))
           do k=1,n_atoms
             h(i,j)=h(i,j)+dot_product(m(:,k),fint_atoms(:,k,i,j))
           enddo
           h(j,i)=h(i,j)
         enddo
       enddo
       return
      end subroutine do_interaction_discr

!------------------------------------------------------------------------
! @brief Build the interaction matrix 
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_interaction_cont(qorf,h)
       real(dbl),intent(IN)   :: qorf(:)
       real(dbl),intent(INOUT):: h(n_ci,n_ci)
       integer(i4b):: i,j

#ifndef MPI
       myrank=0
#endif
#ifdef TDPlas
       if (global_prop_Fint.eq.'ons') then
         h(:,:)=h(:,:)-mut(1,:,:)*qorf(1)-mut(2,:,:)*qorf(2)-mut(3,:,:)*qorf(3)
       elseif(global_prop_Fint.eq.'pcm') then
         do j=1,n_ci
           do i=1,j
             h(i,j)=h(i,j)+dot_product(qorf(:),quantum_vts(:,i,j))
             h(j,i)=h(i,j)
           enddo
         enddo
       else
         if (myrank.eq.0) write(*,*) "wrong interaction type "
#ifdef MPI
         call mpi_finalize(tp_ierr_mpi)
#endif
         stop
       endif
       return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
      end subroutine do_interaction_cont

!------------------------------------------------------------------------
! @brief send molecular dipole to other codes
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine get_m_or_v(m_or_v)     
       implicit none
       real(dbl), intent(out) :: m_or_v(:)
       ! CHECK THIS if one starts from a different quantum state.
#ifdef TDPLAS
       if(global_prop_Fprop.eq."dip") then
         m_or_v(:)=mut(:,1,1)
       else
         m_or_v(:)=quantum_vts(:,1,1)
       endif
       return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
      end subroutine get_m_or_v 


!------------------------------------------------------------------------
! Routines to be defined in fqfm code
!------------------------------------------------------------------------
      subroutine get_q_and_m(q,m)  
        implicit none
        real(dbl), intent(out)      :: q(n_atoms)   !< (1:N_atoms)     - fq  
        real(dbl), intent(out)      :: m(n_atoms) !< (1:N_atoms)     - fm   

        return
      end subroutine get_q_and_m 

      subroutine do_qandm_fqfm(pot,fld,q,m)
        implicit none
        real(dbl), intent(out) :: q(n_atoms)     !< (1:N_atoms)     - fq  
        real(dbl), intent(out) :: m(n_atoms)     !< (1:N_atoms)     - fm   
        real(dbl), intent(in) :: pot(n_atoms)   !< (1:N_atoms)     - potential   
        real(dbl), intent(in) :: fld(3,n_atoms) !< (1:N_atoms)     - field   

        return
      end subroutine do_qandm_fqfm 





!------------------------------------------------------------------------
! @brief Read transition potentials on tesserae
!
! @date Created: S. Pipolo
! Modified: E. Coccia
! Modified: S.Corni (27/06/2020): now the state pair is read from ci_pot,
!           we do not assume upper or lower triangular. Should work
!           for current gamess version as well
!------------------------------------------------------------------------
      subroutine read_gau_out_medium
       integer(i4b) :: i,j,its,nts
       real(dbl)  :: scr

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

      end subroutine read_gau_out_medium


!------------------------------------------------------------------------
! @brief Provide sphere parameters for test QM_coupling
!
! @date Created: S. Pipolo
!------------------------------------------------------------------------
      subroutine grep_sphere_parameters(r,d,sp,wl)
       real(dbl), intent(out)  :: r,d,wl
       real(dbl), intent(out)  :: sp(3)
#ifdef TDPLAS
       wl=sqrt(drudel_eps_A/3)
       sp(1)=pedra_surf_spheres(1)%x 
       sp(2)=pedra_surf_spheres(1)%y 
       sp(3)=pedra_surf_spheres(1)%z 
       r=cts(1)%rsfe
       d=sqrt(dot_product(sp,sp))
       return
#else   
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif 
      end subroutine grep_sphere_parameters



end module interface_classic
