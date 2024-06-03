module interface_tdplas
      use constants
      use readio
#ifdef TDPLAS
      use tdplas, only: set_charges,global_prop_Fmdm_relax,&
! used by dissipation
                        get_mdm_dip,get_gneq,init_mdm,prop_mdm,finalize_mdm,&
                        preparing_for_scf,init_after_scf,& 
! used by propagate
                        readio_and_init_tdplas_for_wt,&
! used by main and main_spectra
                        mpibcast_readio_mdm,quantum_init,global_qmodes_Fmop,global_qmodes_nprint,&
! used by main
                        global_sys_Fwrite,fr_0,BEM_Q0,mat_f0,global_prop_max_cycles,&
                        global_prop_threshold,quantum_vtsn,global_prop_mix_coef,diag_mat,&
                        do_field_from_charges,out_gcharges,get_qr_fr,&
! used in scf
                        BEM_W2,global_sys_Ftest,drudel_eps_w0,drudel_eps_A,BEM_Modes,pedra_surf_spheres,&
                        do_BEM_quant,do_vts_from_dip,deallocate_BEM_public,global_qmodes_nmodes,&
                        global_qmodes_qmmodes,&
! used by QM_coupling 
                        q0,quantum_vts,pedra_surf_n_tessere,global_prop_Fprop,global_prop_Fint,&
                        pedra_surf_tessere,global_prop_Finit_int,pedra_surf_n_spheres,global_medium_Fbem,&
                        global_sys_Fdeb 
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

      real(dbl), allocatable :: this_vts(:,:,:), this_vtsn(:) !<transition potentials on tesserae from cis
      real(dbl), allocatable :: h_mdm(:,:) !< medium contribution to the hamiltonian

      integer(i4b) :: this_nts_act, this_nesf_act,this_nprint,this_max_mod_todiag

      type(tess_pcm_in_wavet), target, allocatable :: this_cts_act(:)
      type(sfera_in_wavet), allocatable :: this_sfe_act(:)
      integer(i4b) :: this_ncycmax !< maximum number of SCF cycles
      integer(i4b), allocatable :: this_imod(:) !<modes to print 
      real(dbl) :: this_thrshld    !< SCF threshold on (i) eigenvalues 10^-global_prop_thrshld (ii) eigenvectors 10^-(global_prop_thrshld+2)
      real(dbl), allocatable :: this_BEM_Q0(:,:)
      real(dbl), allocatable :: this_BEM_W2(:)
      real(dbl), allocatable :: this_BEM_Modes(:,:)
      real(dbl), allocatable :: this_q0(:)
      real(dbl), allocatable :: qx0(:)                 !< local BEM charges 
      real(dbl), allocatable :: q0_fqfm(:),qx0_fqfm(:) !< reaction and local fwfm charges  
      real(dbl), allocatable :: m0_fqfm(:),mx0_fqfm(:) !< reaction and local fwfm dipoles  
      real(dbl) :: this_fr_0(3)                  !< Reaction field at time 0 defined with Finit_mdm, here because used in scf
      real(dbl), allocatable :: this_mat_f0(:,:) !< Onsager's total matrices needed for scf, free_energy and propagation
      real(dbl) :: this_mix_coef   !< SCF mixing ratio of old (1-global_prop_mix_coef) and new (global_prop_mix_coef) charges/field       
      real(dbl) :: this_eps_A,this_eps_w0
      integer(i4b), allocatable :: this_qmmodes(:)
      integer(i4b) :: this_nmodes
      public set_q0charges,this_Fmdm_relax,export_mdm_qmcoup, &
! used by dissipation
             get_medium_dip,get_energies,init_medium_prop,prop_medium,finalize_medium,this_Finit_int,this_Fprop,&
             preparing_for_scf_in_wavet,init_after_scf_in_wavet,& ! used bypropagate (also this_mix_coef)
             read_medium_input,&
! used by main and main_spectra
             mpibcast_read_medium,set_global_tdplas_in_wavet,&
! used by main
             this_Fwrite,this_fr_0,this_BEM_Q0,this_mat_f0,this_ncycmax,this_thrshld,this_vtsn,&
             this_mix_coef,diag_mat_in_wavet,&
             do_field_from_charges_in_wavet, this_nts_act, &
! used in scf
             this_BEM_W2,this_Ftest,this_eps_w0,this_eps_A,this_BEM_Modes,this_sfe_act,&
             do_BEM_quant_in_wavet,do_vts_from_dip_in_wavet,&
             this_Fmop,this_imod,this_nprint,this_max_mod_todiag,deallocate_BEM_public_in_wavet,&
             this_qmmodes,this_nmodes
! used by QM_coupling 
      contains
  
      ! begin - wrapper subroutines
      subroutine set_q0charges
!------------------------------------------------------------------------
! @brief Bridge subroutine to set charges qr_t to q0 during propagation 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
        implicit none
#ifdef TDPLAS
        call set_charges(q0)
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
        this_nts_act=pedra_surf_n_tessere
        if(global_sys_Fdeb.eq."vmu") then
                call do_vts_from_dip_in_wavet
                write(6,*) "Done replacing ci_pot with potential from dipole"
        elseif(this_Fprop.eq."chr-ief".or.this_Fprop.eq."chr-ied".or.this_Fprop.eq."chr-ons") then 
            ! SP 17/05/20 shape is not a standard f90 function
            shap=shape(quantum_vts)
            if(shap(1).eq.this_nts_act) then
               allocate(this_vts(this_nts_act,n_ci,n_ci))
               this_vts=quantum_vts
               allocate(this_vtsn(this_nts_act))
               this_vtsn=quantum_vtsn
            else
               write(6,*) "Error: the number of tesserae for the potential is different than those in the cavity/NP"
               write(6,*) shap(1)," vs ",this_nts_act
               write(6,*) "This is usually due to incoerent ci_pot.inp and cavity.inp files. I stop here" 
               stop
            endif
        end if
        this_ncycmax=global_prop_max_cycles
        this_thrshld=global_prop_threshold
        this_mix_coef=global_prop_mix_coef 
        this_eps_w0=drudel_eps_w0
        this_eps_A=drudel_eps_A
        this_nesf_act=pedra_surf_n_spheres
        allocate(this_sfe_act(this_nesf_act))
        !this_sfe_act=pedra_surf_spheres
        do ii=1, this_nesf_act
         this_sfe_act(ii)%x=pedra_surf_spheres(ii)%x
         this_sfe_act(ii)%y=pedra_surf_spheres(ii)%y
         this_sfe_act(ii)%z=pedra_surf_spheres(ii)%z
         this_sfe_act(ii)%r=pedra_surf_spheres(ii)%r
        end do
        allocate(this_cts_act(this_nts_act))
        !this_cts_act=pedra_surf_tessere
        do ii=1, this_nts_act
         this_cts_act(ii)%x=pedra_surf_tessere(ii)%x
         this_cts_act(ii)%y=pedra_surf_tessere(ii)%y
         this_cts_act(ii)%z=pedra_surf_tessere(ii)%z
         this_cts_act(ii)%rsfe=pedra_surf_tessere(ii)%rsfe
         this_cts_act(ii)%n(:)=pedra_surf_tessere(ii)%n(:)
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
      subroutine init_environment(c,mu,f,h,Fres)
        implicit none
        complex(cmp), intent(in) :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(:)    !< (1:3)              - external field
        real(dbl)   , intent(inout) :: h(:,:)  !< (1:n_ci,1:n_ci) - interaction hamiltonian
        real(dbl), allocatable      :: r(:)  !< auxiliary 3d vector array
        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external potential
        !> Allocating the interaction matrix h_mdm 
        allocate(h_mdm(quantum_n_ci,quantum_n_ci),h_mdm_0(quantum_n_ci,quantum_n_ci))
        h_mdm=zero

#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
          !> build in principle the field coming from all dipoles
          !> initializing medium with molecular dipole and external field
          call init_mdm(mu_t = mu, f_tp = f)
          allocate(this_mat_f0(this_nts_act,this_nts_act))
          this_mat_f0=mat_f0
          this_fr_0=fr_0
        else
          !> build the potential on the surface from all charges and dipoles
          allocate(pot(this_nts_act),potf(this_nts_act))
          pot(:)=zero
          potf(:)=zero
          call do_pot_tdplas(c,mu,f,pot,potf)
          call init_mdm(pot_t = pot, potf_t = potf)
          deallocate(pot)
          deallocate(potf)
          ! CHECK THIS: these two calls must go in scf and QM_coupling
          call prepare_mdm_for_scf
          call prepare_mdm_for_quantum
        end if
        !> build the h_mdm interaction hamiltonian
        call do_int_tdplas
#endif
#ifdef fqfm  
        allocate(pot(n_atoms),potf(n_atoms))
        pot(:)=zero
        potf(:)=zero
        call do_pot_fld_fqfm(c,mu,f,pot,fld,potf)
        call init_fqfm(pot,fld,potf)
        !> build the h_mdm interaction hamiltonian
        call do_int_fqfm
#endif
        !> update the OUT hamiltonian                 
        h=h+h_mdm
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
          allocate(pot(this_nts_act),potf(this_nts_act))
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
        !> build the h_mdm interaction hamiltonian
        call do_int_fqfm
#endif
        !> update the OUT hamiltonian                 
        h=h+h_mdm
        ! SP 15/05/24: CHECK THIS for restart, it was in propagate, is it really needed?
        if (Fres.eq.'Nonr') then
           i=1
           call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
        endif

        return
      end subroutine init_env_prop




!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_pot_tdplas(c,mu,f,pot,potf)
        implicit none
        complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)        - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(3)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(3)    !< (1:3)              - external field
        real(dbl)   , intent(OUT):: pot(this_nts_act)  !< (1:this_nts_act)   - potential on continuum surface
        real(dbl)   , intent(OUT):: potf(this_nts_act)  !< (1:this_nts_act)   - potential on continuum surface
        real(dbl), allocatable      :: r(:)  !< auxiliary 3d vector array
        allocate(r(3,this_nts_act))
        r(1,:)=this_cts_act(:)%x
        r(2,:)=this_cts_act(:)%y
        r(3,:)=this_cts_act(:)%z
        !> prepare the potetial acting on the medium for initialisation
        if(this_Fint.eq."ons") then
         !> computing molecular potential from its dipole and not the coefficients                     
         call do_pot_from_dip(1,this_mol_cc,mu,this_nts_act,r,pot)
        else
         !> computing molecular potential from the coefficients
         call do_pot_from_coeff(c,this_nts_act,this_vts,pot)
        end if
        !> computing external potential in the long-wavelength limit
        call do_pot_from_field(f,this_nts_act,r,potf)
#ifdef fqfw 
        ! SP 6/5/2024: in the following routines pot is updated so that it contains the total pontential
        ! of the molecule charge density and the fwfm system that determines the apparent
        ! surface charges
        call do_pot_from_charges(n_atoms,r_atoms,q_disc,this_nts_act,r,pot)
        call do_pot_from_dip(n_atoms,r_atoms,m_disc,this_nts_act,r,pot)
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
        real(dbl)   , intent(OUT):: pot(n_atoms)  !< (1:this_nts_act)   - potential on continuum surface
        real(dbl)   , intent(OUT):: fld(3,n_atoms)  !< (1:this_nts_act)   - potential on continuum surface
        real(dbl)   , intent(OUT):: potf(3,n_atoms)  !< (1:this_nts_act)   - potential on continuum surface
        real(dbl), allocatable      :: r(:)  !< auxiliary 3d vector array
        real(dbl), allocatable      :: qorf(:) !< (1:pedra_surf_n_tessere)     - charges or field  
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
          allocate(qorf(3)
          stop "Coupling with fqfw not implemented for dipole propagation"
        else
          allocate(qorf(this_nts_act))
          allocate(r(3,this_nts_act))
          r(1,:)=this_cts_act(:)%x
          r(2,:)=this_cts_act(:)%y
          r(3,:)=this_cts_act(:)%z
          !> get charges from td_contmed               
          call get_full_qorf(qorf)
          call do_pot_from_charges(this_nts_act,r,qorf,n_atoms,r_atoms,pot)
          call do_fld_from_charges(this_nts_act,r,qorf,n_atoms,r_atoms,fld)
          deallocate(r)
        end if
#endif  
        !> prepare the potetial acting on the medium for initialisation
        ! CHECK THIS the following Flag should be defined for fqfm
        if(this_Fint_fqfm.eq."ons") then
          !> computing molecular potential nad field from its dipole and not the coefficients                     
          call do_pot_from_dip(1,this_mol_cc,mu,n_atoms,r_atoms,pot)
          call do_fld_from_dip(1,this_mol_cc,mu,n_atoms,r_atoms,fld)
        else
          !> computing molecular potential from the coefficients
          call do_pot_from_coeff(c,n_atoms,vint_atoms,pot)
          call do_fld_from_coeff(c,n_atoms,fint_atoms,fld)
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
      subroutine do_int_tdplas
        implicit none
        real(dbl), allocatable      :: qorf(:) !< (1:pedra_surf_n_tessere)     - charges or field  
        if(this_Fprop.eq."dip") then
          !> allocate and initialise arrays and prepare for calls
          allocate(qorf(3))
        else
          !> allocate and initialise arrays and prepare for calls
          allocate(qorf(this_nts_act))
        end if
        !> get charges from td_contmed               
        call get_full_qorf(qorf)
        !> construct the interaction hamiltonian h_mdm
        call do_interaction_cont(qorf,h_mdm)
        deallocate(qorf)
        return
      end subroutine do_int_tdplas


!------------------------------------------------------------------------
! @brief Initialize medium for quantum states initialisation e.g. scf or 
!        quantum coupling.
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  
!------------------------------------------------------------------------
      subroutine do_int_fqfm  `
        implicit none
        real(dbl), allocatable      :: q(:) !< (1:pedra_surf_n_tessere)     - charges or field  
        real(dbl), allocatable      :: m(:,:) !< (1:pedra_surf_n_tessere)     - charges or field  
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
          allocate(this_BEM_Q0(this_nts_act,this_nts_act))
          this_BEM_Q0=BEM_Q0
        return
      end subroutine prepare_mdm_for_scf

!------------------------------------------------------------------------
! @brief Prepare medium for quantum coupling
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
     ! Now called by propagate, afterward should be called by initialize
     subroutine prepare_mdm_for_quantum
          if(this_Fbem.eq.'diag') then
            allocate(this_BEM_W2(this_nts_act))
            this_BEM_W2=BEM_W2
          end if
          if(Fmdm.eq."qnan") then
           allocate(this_BEM_Modes(this_nts_act,this_nts_act))
           this_BEM_Modes=BEM_Modes
          end if
        return
      end subroutine prepare_mdm_for_quantum


!------------------------------------------------------------------------
! @brief update medium dof (q,f,etc) during a scf cycle
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
      subroutine update_environment_scf(c,f)
        implicit none
        complex(cmp), intent(in) :: c(n_ci) !> (1:n_ci)           - molecular wavefunction coefficients
        complex(cmp), intent(in) :: f(3)    !> (1:3)           - external field                     
        real(dbl), allocatable   :: mu(:)   !> molecular dipole 
        integer(i4b)::i    
        allocate(mu(3))
        mu=zero
#ifdef TDPLAS
        if((this_Fprop.eq."dip").or.(this_Fint.eq."ons")) then 
          call do_dip_from_coeff(c,mu)
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
! @brief transform and write out potentials            
!
! @date Created   : S. Pipolo  
! Modified  :  
!------------------------------------------------------------------------
      subroutine out_environment_scf(c)
        implicit none
        complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)           - molecular wavefunction coefficients
        integer(i4b)              :: its  

        if(this_Fprop.ne."dip") then 
        ! transform vts to the SCF state basis
         do its=1,this_nts_act
          this_vts(its,:,:)=matmul(this_vts(its,:,:),c)
          this_vts(its,:,:)=matmul(transpose(c),this_vts(its,:,:))
         enddo
         if (myrank.eq.0) then
            call out_charges(q0)
            call out_vts
         endif
       endif
        return
      end subroutine out_environment_scf


     
!------------------------------------------------------------------------
! @brief Compute field from dipole 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_BEM_field(c,mu)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)   - molecular wavefunction coefficients
       real(dbl), intent(in) :: mu(3)    !< (1:n_ci)        - molecular wavefunction coefficients
       real(dbl), allocatable :: f(:)    !< (1:3)           - molecular wavefunction coefficients
       ! SP 04/0717 matmul for general spheroid orientation
       allocate(pot(this_nts_act))
       call do_rfield_from_dip(mu,f)
       !fr_0=(1.-this_mix_coef)*fr_0+this_mix_coef*matmul(this_mat_f0,mu)
       fr_0=(1.-this_mix_coef)*fr_0+this_mix_coef*f
       return
      end subroutine update_BEM_field

!------------------------------------------------------------------------
! @brief Compute charges from potential 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_BEM_charges(c,mu)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci)    !< (1:n_ci)           - molecular wavefunction coefficients
       complex(cmp), intent(in) :: mu(n_ci)    !< (1:n_ci)           - molecular dipole                   
       real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
       real(dbl), allocatable      :: potf(:)  !< (1:pedra_surf_n_tessere)     - field potential
       real(dbl), allocatable      :: q(:)  !< (1:pedra_surf_n_tessere)     - field potential
       integer(i4b)::i    
       allocate(pot(this_nts_act))
       allocate(potf(this_nts_act))
       allocate(q(this_nts_act))
       call do_pot_tdplas(c,mu,f,pot,potf)
       ! SP 11/05/24: shall we do the matmul operation in td_contmed?
       call do_charges_from_pot(pot,q)
       q0=(1.-this_mix_coef)*q0+this_mix_coef*q
       call do_charges_from_pot(potf,q)
       qx0=(1.-this_mix_coef)*qx0+this_mix_coef*q
       deallocate(pot,potf,q) 
       ! SC 12/8/2016: apparently for NP, charge compensation is needed
       if (Fmdm.eq.'cnan'.or.Fmdm.eq.'qnan') then
         q0=q0-sum(q0)/this_nts_act
         qx0=qx0-sum(qx0)/this_nts_act
       endif
       return
      end subroutine update_BEM_charges


!------------------------------------------------------------------------
! @brief Compute fqfm charges and dipoles from potential and field 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine update_fqfm_char_and_dip(c,mu,f)
       implicit none 
       complex(cmp), intent(in) :: c(n_ci) !< (1:n_ci)           - molecular wavefunction coefficients
       real(dbl)   , intent(in) :: mu(3)   !< (1:3)              - molecular dipole
       real(dbl)   , intent(in) :: f(3)    !< (1:3)              - external field
       real(dbl)   , intent(OUT):: pot(n_atoms)    !< (1:n_atoms)   - potential on fqfw atoms       
       real(dbl)   , intent(OUT):: fld(3,n_atoms)  !< (1:n_atoms)   - field on fqfw atoms
       real(dbl)   , intent(OUT):: potf(3,n_atoms) !< (1:n_atoms)   - external potential fqfw atoms        
       real(dbl)   , allocatable:: q(n_atoms)      !< (1:n_atoms)   - charges on fqfw atoms
       real(dbl)   , allocatable:: m(3,n_atoms)    !< (1:n_atoms)   - dipoles on fqfw atoms
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
      end subroutine update_fqfm_char_and_dip


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
       write(7,*) this_nts_act
       i=0
       j=0
       ! V00
       write(7,*) i,j 
       do its=1,this_nts_act
        write(7,*) this_vts(its,1,1)-this_vtsn(its),0.d0,this_vtsn(its)
       enddo

       do i=2,n_ci
          write(7,*) 0, i-1
          do its=1,this_nts_act
             write(7,*) this_vts(its,1,i) 
          enddo
       enddo

       do i=2,n_ci
           do j=2,i
              write(7,*)  i-1, j-1
              do its=1,this_nts_act
                 if (i.eq.j) then
                     write(7,*) this_vts(its,i,j)-this_vtsn(its)
                 else
                     write(7,*) this_vts(its,i,j) 
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
      subroutine prop_medium(i,c,mu,f)

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

        integer(i4b), intent(in) :: i

         ! To be more efficient this if should go in the propagate of waveT
         ! Propagate medium only every global_prop_n_q timesteps
          if(mod(i,global_prop_n_q).eq.0) then
            ! Build the interaction Hamiltonian Reaction/Local with previous charges
            call do_interaction
            ! Update the interaction Hamiltonian
            h(:,:)=h(:,:)+h_mdm(:,:)
            return
          endif
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
         ! propagating medium with molecular dipole and external field
         ! Get charges from external codes 
         call prop_mdm(i, mu_t = mu, f_tp = f)
        else
         allocate(pot(this_nts_act))
         allocate(potf(this_nts_act))
         if(this_Fint.eq."ons") then
          ! computing molecular potential corresponding to a point-like dipole
          call do_pot_from_dip(mu,pot)
         else
          ! computing molecular potential
          call do_pot_from_coeff(c,this_nts_act,this_vts,pot)
         end if
         ! computing external potential in the long-wavelength limit
         call do_pot_from_field(f,potf)

         ! propagating medium with molecular and external potentials
         ! SP 15/05/20 changed this_Ftest with global_sys_Ftest
         if(global_sys_Ftest.eq."n-r") then
          call prop_mdm(i, mu_t = mu, pot_t = pot, potf_t = potf, h_int = h)
         else
          !write(*,*) 'uella 1', myrank
          call prop_mdm(i, pot_t = pot, potf_t = potf, h_int = h)
          !write(*,*) 'uella 2', myrank
#ifdef MPI
         call mpi_finalize(ierr_mpi)
         stop
#endif
         end if
         deallocate(pot)
         deallocate(potf)


        end if
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

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
      subroutine set_global_tdplas_in_wavet(this_dt,this_mdm,this_mol_cc,this_n_ci,this_n_ci_read,this_c_i,this_e_ci,this_mut,&
				                                    this_fmax,this_omega,this_Ffld,this_n_out,this_n_f,this_tdelay,this_pshift,&
                                            this_Fbin,this_Fopt,this_res,this_n_res)

        implicit none

        real(dbl)     , intent(in) :: this_dt				         ! time step
        character(3)  , intent(in) :: this_mdm				         ! kind of medium
        integer(i4b)  , intent(in) :: this_n_ci,this_n_ci_read		 ! number of CIS states
        real(dbl)     , intent(in) :: this_e_ci(:)	        	     ! CIS energies
        real(dbl)     , intent(in) :: this_mut(:,:,:)			     ! CIS transition dipoles
        real(dbl)     , intent(in) :: this_mol_cc(3)			     ! molecule center
        real(dbl)     , intent(in) :: this_fmax(3,10),this_omega(10) ! field amplitude and frequency
        real(dbl)     , intent(in) :: this_tdelay(10),this_pshift(10)! time delay and phase shift
        complex(cmp)  , intent(in) :: this_c_i(:)                    ! CIS coefficients
        character(3)  , intent(in) :: this_Ffld			           	 ! shape of impulse
        character(3)  , intent(in) :: this_Fbin                      ! binary output
        character(3)  , intent(in) :: this_Fopt                      ! matrix/vector multiplication 
        integer(i4b)  , intent(in) :: this_n_out,this_n_f	         ! auxiliaries for output
        character(1)  , intent(in) :: this_res                      ! restart for medium 
        integer(i4b)  , intent(in) :: this_n_res                    ! frequency for restart

#ifdef TDPLAS
        call quantum_init(this_dt,this_mol_cc,this_n_ci,this_n_ci_read,this_mut,this_e_ci,this_c_i,&
			       this_fmax,this_omega,this_Ffld,this_n_out,this_n_f,&
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
       real(dbl),intent(in):: q(this_nts_act)  

#ifdef TDPLAS
       call do_field_from_charges(q,f)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine do_field_from_charges_in_wavet

      subroutine do_BEM_quant_in_wavet

       implicit none
       integer(i4b) :: i

#ifdef TDPLAS
       call do_BEM_quant
       this_nmodes=global_qmodes_nmodes
       allocate(this_qmmodes(this_nmodes))
       do i=1,this_nmodes
         this_qmmodes(i)=global_qmodes_qmmodes(i)
       enddo
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine do_BEM_quant_in_wavet

      subroutine deallocate_BEM_public_in_wavet

       implicit none

#ifdef TDPLAS
       call deallocate_BEM_public
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine deallocate_BEM_public_in_wavet

      subroutine do_vts_from_dip_in_wavet

       implicit none
#ifdef TDPLAS
       call do_vts_from_dip
       this_vts=quantum_vts
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine do_vts_from_dip_in_wavet

      subroutine init_environment_scf(c_c,f)
        implicit none
        complex(cmp), intent(in) :: c(n_ci) !> (1:n_ci)           - molecular wavefunction coefficients
        complex(cmp), intent(in) :: f(3)    !> (1:3)           - external field                     
        ! CHECK THIS:  are we forgetting about something for the initialisation? If not this is not needed and 
        ! update_environment_scf can be called at the beginning of the scf cycle
        call update_environment_scf(c_c,f)
        ! SP 15/05/24 before in propagate, now done in initialize
        !trans_dipoles = mut
        !energies=e_ci 
        !if(Fmag.eq.'mag') then 
        !    trans_mag=lt
        !endif
        ! compute the molecular dipole
        call do_dip_from_coeff(c_prev,mu_prev,nstates)
        if(this_Fprop.eq."dip") then
          call init_after_scf_in_wavet(mu_prev)
        else
          call do_pot_from_coeff(c_prev,pot_prev)
          call init_after_scf_in_wavet(pot_prev)
        endif

#ifdef TDPLAS
        call preparing_for_scf(mix, pot_or_mu)
        call get_qr_fr(q_or_f)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine init_environment_scf

      subroutine init_after_scf_in_wavet(pot_or_mut)

       implicit none
       real(dbl), intent(in) :: pot_or_mut(:)

#ifdef TDPLAS
! SC 16/10/2020: vts in tdplas must be updated
       quantum_vts=this_vts
       call init_after_scf(pot_or_mut)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine init_after_scf_in_wavet

      subroutine out_gcharges_in_wavet 

       implicit none

#ifdef TDPLAS
       call out_gcharges
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif                                                                  
      end subroutine out_gcharges_in_wavet                           


! begin - subroutines to calculate dipoles, field and potentials from coefficients, dipoles and fields

!------------------------------------------------------------------------
! @brief Compute dipole from CIS coefficients 
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_dip_from_coeff(c,dip,nc)

       implicit none

       integer(i4b), intent(IN)  :: nc  
       complex(cmp), intent(IN)  :: c(nc)
       real(dbl),    intent(OUT) :: dip(3)
       integer(i4b)              :: its,j,k  
       complex(cmp)              :: ctmp(nc) 

#ifndef OMP
       dip(1)=dot_product(c,matmul(mut(1,:,:),c))
       dip(2)=dot_product(c,matmul(mut(2,:,:),c))
       dip(3)=dot_product(c,matmul(mut(3,:,:),c))
#endif
#ifdef OMP
      if (Fopt.eq.'omp') then
         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(1,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(1)=dot_product(c,ctmp)

         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(2,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(2)=dot_product(c,ctmp)

         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(3,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(3)=dot_product(c,ctmp)
      else
         dip(1)=dot_product(c,matmul(mut(1,:,:),c))
         dip(2)=dot_product(c,matmul(mut(2,:,:),c))
         dip(3)=dot_product(c,matmul(mut(3,:,:),c))
      endif 
#endif

      end subroutine do_dip_from_coeff

!------------------------------------------------------------------------
! @brief Compute potential (pot) on n points from CIS coefficientes (c) 
!        and potential integrals (v)
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_pot_from_coeff(c,n,v,pot)

       implicit none

       complex(cmp), intent(IN)        :: c(n_ci)
       real(dbl),    intent(IN)        :: n             
       real(dbl),    intent(IN)        :: v(n,n_ci,n_ci)
       real(dbl),    intent(INOUT)       :: pot(n)

       integer(i4b)                       :: i,k,j  
       complex(cmp), save, allocatable    :: ctmp(:)
       complex(cmp), save                 :: cc

#ifndef OMP
       do i=1,this_nts_act
          pot(i)=pot(i)+dot_product(c,matmul(v(i,:,:),c))
       enddo
#endif

#ifdef OMP
       if (Fopt.eq.'omp') then
          allocate(ctmp(n*n_ci))
!$OMP PARALLEL REDUCTION (+:cc)
!$OMP DO 
          do i=1,n
             do k=1,n_ci
                cc=0.d0
                do j=1,n_ci
                   cc = cc + v(i,k,j)*c(j)
                enddo
                ctmp(k+(i-1)*n_ci) = cc
             enddo
          enddo
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,ctmp((i-1)*n_ci+1:i*n_ci))
          enddo
!$OMP END PARALLEL
          deallocate(ctmp)
       else
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,matmul(v(i,:,:),c))
          enddo 
!$OMP END PARALLEL
       endif
#endif

      end subroutine do_pot_from_coeff

!------------------------------------------------------------------------
! @brief Compute the fiel (fld) on n points from CIS coefficientes (c) 
!        and field integrals (f)
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_fld_from_coeff(c,n,f,fld)

       implicit none

       complex(cmp), intent(IN)        :: c(n_ci)
       real(dbl),    intent(IN)        :: n             
       real(dbl),    intent(IN)        :: v(3,n,n_ci,n_ci)
       real(dbl),    intent(INOUT)       :: fld(n)

       integer(i4b)                       :: i,k,j,l  
       complex(cmp), save, allocatable    :: ctmp(:)
       complex(cmp), save                 :: cc

#ifndef OMP
       do i=1,this_nts_act
         do j=1,3
           fld(i)=fld(i)+dot_product(c,matmul(f(j,i,:,:),c))
         enddo
       enddo
#endif

#ifdef OMP
       if (Fopt.eq.'omp') then
          allocate(ctmp(n*n_ci))
!$OMP PARALLEL REDUCTION (+:cc)
!$OMP DO 
          do i=1,n
             do k=1,n_ci
                cc=0.d0
                do j=1,n_ci
                   cc = cc + v(i,k,j)*c(j)
                enddo
                ctmp(k+(i-1)*n_ci) = cc
             enddo
          enddo
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,ctmp((i-1)*n_ci+1:i*n_ci))
          enddo
!$OMP END PARALLEL
          deallocate(ctmp)
       else
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,matmul(v(i,:,:),c))
          enddo 
!$OMP END PARALLEL
       endif
#endif

      end subroutine do_fld_from_coeff

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on n points od coordinates r 
!        generated by an external electric field (fld)
! (fld) 
!
! @date Created: S. Pipolo
! Modified: 
!------------------------------------------------------------------------
      subroutine do_pot_from_field(fld,n,r,pot)

       implicit none

       real(dbl), intent(IN):: fld(3) 
       integer(i4b), intent(IN):: n 
       real(dbl), intent(IN):: r(n) 
       real(dbl), intent(INOUT):: pot(n) 
       integer(i4b) :: i  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO 
#endif
        do i=1,n
          pot(i)=pot(i)-dot_product(fld,r(:,i))           
        enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_field

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nd dipoles (dip) at positions rd
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_dip(nd,rd,dip,n,r,pot)

       integer(i4b), intent(IN) :: nd
       real(dbl), intent(IN) :: rd(3,nd)
       real(dbl), intent(IN) :: dip(3,nd)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: pot(n)
       real(dbl):: diff(3)  
       real(dbl):: distm1 
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nd
            diff(:)=-(rd(:,j)-r(:,i))
            distm1=1/sqrt(dot_product(diff,diff))
            pot(i)=pot(i)+dot_product(diff,dip(:,j))*distm1*distm1*distm1
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_dip

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nd dipoles (dip) at positions rd
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_fld_from_dip(nd,rd,dip,n,r,fld)

       integer(i4b), intent(IN) :: nd
       real(dbl), intent(IN) :: rd(3,nd)
       real(dbl), intent(IN) :: dip(3,nd)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: fld(3,n)
       real(dbl):: diff(3),f(3)  
       real(dbl):: distm1 
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nd
            diff(:)=-(rd(:,j)-r(:,i))
            distm1=1/sqrt(dot_product(diff,diff))
            diff(:)=diff(:)*distm1
            f(:)=3*dot_product(diff,dip(:,j))*diff(:)-dip(:,j)
            fld(:,i)=fld(:,i)+f(:)*distm1*distm1*distm1
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_fld_from_dip

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nq charges q at positions rq
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_charges(nq,rq,q,n,r,pot)

       integer(i4b), intent(IN) :: nq
       real(dbl), intent(IN) :: rq(3,nq)
       real(dbl), intent(IN) :: q(nq)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: pot(n)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nq
            diff(:)=(rq(:,j)-r(:,i))
            dist=sqrt(dot_product(diff,diff))
            pot(i)=pot(i)+q(j)/dist
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_charges

!------------------------------------------------------------------------
! @brief Compute the field (fld) on a number (n) of points of 
!        coordinates r generated by nq charges q at positions rq
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_fld_from_charges(nq,rq,q,n,r,fld)

       integer(i4b), intent(IN) :: nq
       real(dbl), intent(IN) :: rq(3,nq)
       real(dbl), intent(IN) :: q(nq)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: fld(3:n)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nq
            diff(:)=-(rq(:,j)-r(:,i))
            dist=sqrt(dot_product(diff,diff))
            fld(i)=fld(i)+q(j)*diff/dist/dist/dist
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_fld_from_charges

subroutine export_mdm_qmcoup
   implicit none
   integer(i4b) :: i
         this_nprint=global_qmodes_nprint
         allocate(this_BEM_W2(this_nts_act))
         this_BEM_W2=BEM_W2
         allocate(this_BEM_Modes(this_nts_act,this_nts_act))
         this_BEM_Modes=BEM_Modes
end subroutine

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
       real(dbl),intent(IN):: q(n_atoms)
       real(dbl),intent(IN):: m(3,n_atoms)
       real(dbl),intent(IN):: h(n_ci,n_ci)
       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif
       do j=1,quantum_n_ci
         do i=1,j
           h(i,j)=h(i,j)+dot_product(q(:),vint_atoms(:,i,j))
           do k=1,_n_atoms
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
       real(dbl),intent(IN):: qorf(:)
       real(dbl),intent(IN):: h(n_ci,n_ci)
       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif
       if (global_prop_Fint.eq.'ons') then
         h(:,:)=h(:,:)-mut(1,:,:)*qorf(1)-mut(2,:,:)*qorf(2)-mut(3,:,:)*qorf(3)
       elseif(global_prop_Fint.eq.'pcm') then
         do j=1,quantum_n_ci
           do i=1,j
             h(i,j)=h(i,j)+dot_product(qorf(:),quantum_vts(:,i,j))
             h(j,i)=h(i,j)
           enddo
         enddo
       else
         if (tp_myrank.eq.0) write(*,*) "wrong interaction type "
#ifdef MPI
         call mpi_finalize(tp_ierr_mpi)
#endif
         stop
       endif
       return
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
       if(global_prop_Fprop.eq."dip") then
         m_or_v(:)=mut(:,1,1)
       else
         m_or_v(:)=quantum_vts(:,1,1)
       endif
       return
      end subroutine get_mol_dipole     
!
!------------------------------------------------------------------------
! @brief compute charges from potential
! This shoud stay in TDPLAS td_contmed????
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_charges_from_pot(pot,q)
       real(dbl),intent(IN):: pot(this_nts_act)
       real(dbl),intent(OUT):: q(this_nts_act)
       q=matmul(this_BEM_Q0,pot)
       return
      end subroutine do_interaction_cont


!
!   Deal with restart
     if(global_prop_Fmdm_res.eq.'yesr') then
        h_int=0
        if(global_sys_Fdeb.ne."off") call do_interaction_h
        h_int(:,:)=h_int(:,:)+h_mdm(:,:)

! Check where to put this
#ifndef MPI
       if (i.eq.1.and.global_sys_Fwrite.eq."high") then
        write(6,*) "h_mdm at the first propagation step"
        do j=1,quantum_n_ci
         do k=1,quantum_n_ci
          write (6,*) j,k,h_mdm(j,k)
         enddo
        enddo
       endif
#endif
! this is after do_c_oldbasis in 


end module interface_tdplas
