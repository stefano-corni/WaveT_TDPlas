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
      real(dbl) :: this_fr_0(3)                       !< Reaction field at time 0 defined with Finit_mdm, here because used in scf
      real(dbl), allocatable :: this_mat_f0(:,:) !< Onsager's total matrices needed for scf, free_energy and propagation
      real(dbl) :: this_mix_coef   !< SCF mixing ratio of old (1-global_prop_mix_coef) and new (global_prop_mix_coef) charges/field       
      real(dbl) :: this_eps_A,this_eps_w0
      integer(i4b), allocatable :: this_qmmodes(:)
      integer(i4b) :: this_nmodes
      public set_q0charges,this_Fmdm_relax,export_mdm_qmcoup, &
! used by dissipation
             get_medium_dip,get_energies,init_medium,prop_medium,finalize_medium,this_Finit_int,this_Fprop,&
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
! @brief Initialize medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
     subroutine init_medium(c,mu,f,h)

        implicit none
        complex(cmp), intent(in) :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(:)    !< (1:3)              - external field
        real(dbl)   , intent(inout) :: h(:,:)  !< (1:n_ci,1:n_ci) - interaction hamiltonian


        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external potential
#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
         ! initializing medium with molecular dipole and external field
         call init_mdm(mu_t = mu, f_tp = f, h_int = h)
         allocate(this_mat_f0(this_nts_act,this_nts_act))
         this_mat_f0=mat_f0
         this_fr_0=fr_0
        else
         allocate(pot(this_nts_act))
         allocate(potf(this_nts_act))
         if(this_Fint.eq."ons") then
          ! computing molecular potential corresponding to a point-like dipole
          call do_pot_from_dip(mu,pot)
         else
          ! computing molecular potential
          call do_pot_from_coeff(c,pot)
         end if
         ! computing external potential in the long-wavelength limit
         call do_pot_from_field(f,potf)
         ! initializing medium with molecular and external potentials
         call init_mdm(pot_t = pot, potf_t = potf, h_int = h)
         deallocate(pot)
         deallocate(potf)
         allocate(this_q0(this_nts_act))
         this_q0=q0
         allocate(this_BEM_Q0(this_nts_act,this_nts_act))
         this_BEM_Q0=BEM_Q0
         if(this_Fbem.eq.'diag') then
           allocate(this_BEM_W2(this_nts_act))
           this_BEM_W2=BEM_W2
         end if
         if(Fmdm.eq."qnan") then
          allocate(this_BEM_Modes(this_nts_act,this_nts_act))
          this_BEM_Modes=BEM_Modes
         end if
        end if
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine init_medium
      
      
!------------------------------------------------------------------------
! @brief Propagate medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------
      subroutine prop_medium(i,c,mu,f,h)

        implicit none

        complex(cmp), intent(in) :: c(:)    !< (1:n_ci)           - molecular wavefunction coefficients
        real(dbl)   , intent(in) :: mu(:)   !< (1:3)              - molecular dipole
        real(dbl)   , intent(in) :: f(:)    !< (1:3)              - external field
        real(dbl)   , intent(inout) :: h(:,:)  !< (1:n_ci,1:n_ci) - interaction hamiltonian

        real(dbl), allocatable      :: pot(:)  !< (1:pedra_surf_n_tessere)     - molecular potential
        real(dbl), allocatable      :: potf(:) !< (1:pedra_surf_n_tessere)     - external  potential

        integer(i4b), intent(in) :: i

#ifdef TDPLAS
        if(this_Fprop.eq."dip") then
         ! propagating medium with molecular dipole and external field
         call prop_mdm(i, mu_t = mu, f_tp = f, h_int = h)
        else
         allocate(pot(this_nts_act))
         allocate(potf(this_nts_act))
         if(this_Fint.eq."ons") then
          ! computing molecular potential corresponding to a point-like dipole
          call do_pot_from_dip(mu,pot)
         else
          ! computing molecular potential
          call do_pot_from_coeff(c,pot)
         end if
         ! computing external potential in the long-wavelength limit
         call do_pot_from_field(f,potf)
         ! propagating medium with molecular and external potentials
         ! SP 15/05/20 changed this_Ftest with global_sys_Ftest
         if(global_sys_Ftest.eq."n-r") then
          call prop_mdm(i, mu_t = mu, pot_t = pot, potf_t = potf, h_int = h)
         else
          call prop_mdm(i, pot_t = pot, potf_t = potf, h_int = h)
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

      subroutine preparing_for_scf_in_wavet(mix,pot_or_mu,q_or_f)

       implicit none

       real(dbl), intent(in) :: mix
       real(dbl), intent(in) :: pot_or_mu(:)
       real(dbl), intent(out) :: q_or_f(:)

#ifdef TDPLAS
        call preparing_for_scf(mix, pot_or_mu)
        call get_qr_fr(q_or_f)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

      end subroutine preparing_for_scf_in_wavet

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
! @brief Compute potential on BEM surface from CIS coefficientes 
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_pot_from_coeff(c,pot)

       implicit none

       complex(cmp), intent(IN)        :: c(n_ci)
       real(dbl),    intent(OUT)       :: pot(this_nts_act)

       integer(i4b)                       :: its,k,j  
       complex(cmp), save, allocatable    :: ctmp(:)
       complex(cmp), save                 :: cc

#ifndef OMP
       do its=1,this_nts_act
          pot(its)=dot_product(c,matmul(this_vts(its,:,:),c))
       enddo
#endif

#ifdef OMP
       if (Fopt.eq.'omp') then
          allocate(ctmp(this_nts_act*n_ci))
!$OMP PARALLEL REDUCTION (+:cc)
!$OMP DO 
          do its=1,this_nts_act
             do k=1,n_ci
                cc=0.d0
                do j=1,n_ci
                   cc = cc + this_vts(its,k,j)*c(j)
                enddo
                ctmp(k+(its-1)*n_ci) = cc
             enddo
          enddo
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO
          do its=1,this_nts_act
             pot(its)=dot_product(c,ctmp((its-1)*n_ci+1:its*n_ci))
          enddo
!$OMP END PARALLEL
          deallocate(ctmp)
       else
!$OMP PARALLEL
!$OMP DO
          do its=1,this_nts_act
             pot(its)=dot_product(c,matmul(this_vts(its,:,:),c))
          enddo 
!$OMP END PARALLEL
       endif
#endif

      end subroutine do_pot_from_coeff

!------------------------------------------------------------------------
! @brief Compute the potential on the BEM surface generated by a field
! (fld) 
!
! @date Created: S. Pipolo
! Modified: 
!------------------------------------------------------------------------
      subroutine do_pot_from_field(fld,pot)

       implicit none

       real(dbl), intent(in):: fld(3) 
       real(dbl), intent(out):: pot(this_nts_act) 
       integer(i4b) :: its  

       ! Field
       pot(:)=zero
#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO 
#endif
do its=1,this_nts_act
          pot(its)=pot(its)-fld(1)*this_cts_act(its)%x           
          pot(its)=pot(its)-fld(2)*this_cts_act(its)%y          
          pot(its)=pot(its)-fld(3)*this_cts_act(its)%z         
        enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_field

!------------------------------------------------------------------------
! @brief Compute potential on BEM surface from CIS coefficients (ons) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_dip(dip,pot)


       real(dbl), intent(IN) :: dip(3)
       real(dbl), intent(OUT) :: pot(this_nts_act)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: its  

       pot(:)=zero
#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do its=1,this_nts_act
          diff(1)=-(mol_cc(1)-this_cts_act(its)%x)
          diff(2)=-(mol_cc(2)-this_cts_act(its)%y)
          diff(3)=-(mol_cc(3)-this_cts_act(its)%z)
          dist=sqrt(dot_product(diff,diff))
          pot(its)=pot(its)+dot_product(diff,dip)/(dist**3)
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_dip

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


end module interface_tdplas
