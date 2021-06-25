      Module BEM_medium
      use tdplas_constants 
      use cavity_types
      use global_tdplas
      use global_quantum
      use readfile_freq

      use pedra_friends
      use sphe_surface
      use MathTools
      use readfile_epsilon
      use drudel_epsilon
      use debye_epsilon
      use dielectric_function

#ifdef OMP
      use omp_lib
#endif

#ifdef MPI
      use mpi
#endif

!      use, intrinsic :: iso_c_binding

      implicit none

! SP 25/06/17: variables starting with "BEM_" are public.
      real(dbl), allocatable :: BEM_S(:,:),BEM_D(:,:)    !< Calderon S and D matrices
      real(dbl), allocatable :: BEM_L(:),BEM_T(:,:)      !< $\Lambda$ and T eigenMatrices
      complex(cmp), allocatable :: BEM_VLc(:,:), BEM_VRc(:,:), BEM_Lc(:)
      real(dbl), allocatable :: BEM_W2(:),BEM_Modes(:,:) !< BEM squared frequencies and Modes ($T*S^{1/2}$)
      real(dbl), allocatable :: BEM_ADtm1(:,:)           !< BEM inverse matrix for general eps-prop 3
! SP 25/06/17: K0 and Kd are still common to 'deb' and 'drl' cases
      real(dbl), allocatable :: K0(:),Kd(:)              !< Diagonal $K_0$ and $K_d$ matrices
      real(dbl), allocatable :: K0x(:),Kdx(:)
      real(dbl), allocatable :: fact1(:),fact2(:)        !< Diagonal vectors for propagation matrices
      real(dbl), allocatable :: fact2x(:)
      real(dbl), allocatable :: Sm12T(:,:),TSm12(:,:)    !< $S^{-1/2}T$ and $T*S^{-1/2}$ matrices
      real(dbl), allocatable :: TSp12(:,:)               !< $T*S^{1/2}$ matrices
      real(dbl), allocatable :: BEM_Sm12(:,:),Sp12(:,:)  !< $S^{-1/2}$ and $S^{1/2}$ matrices
      real(dbl), allocatable :: BEM_Q0(:,:),BEM_Qd(:,:)  !< Static and Dyanamic BEM matrices $Q_0$ and $Q_d$
      real(dbl), allocatable :: BEM_Qt(:,:),BEM_R(:,:)   !< Debye propagation matrices $\tilde{Q}$ and $R$
      real(dbl), allocatable :: BEM_Qw(:,:),BEM_Qf(:,:)  !< Drude-Lorents propagation matrices $Q_\omega$ and $Q_f$
      real(dbl), allocatable :: BEM_Qg(:,:)              !< General propagation matrix $Q_\gamma$
      real(dbl), allocatable :: BEM_2G(:)
      real(dbl), allocatable :: BEM_Q0x(:,:),BEM_Qdx(:,:)  !< Static and Dyanamic BEM matrices $Q_0$ and $Q_d$ for local (x) field
      real(dbl), allocatable :: BEM_Qtx(:,:)
      real(dbl), allocatable :: BEM_Qfx(:,:)
      real(dbl), allocatable :: MPL_Ff(:,:,:,:),MPL_Fw(:,:)!< Onsager's Matrices with factors for reaction and local (x) field
      real(dbl), allocatable :: MPL_F0(:,:,:,:),MPL_Fx0(:,:,:) !< Onsager's Matrices with factors for reaction and local (x) field
      real(dbl), allocatable :: MPL_Ft0(:,:,:),MPL_Ftx0(:,:,:) !< Onsager's Matrices with factors/tau for reaction and local (x) field
      real(dbl), allocatable :: MPL_Fd(:,:,:,:),MPL_Fxd(:,:,:) !< Onsager's Matrices with factors for local(x) field
      real(dbl), allocatable :: MPL_Taum1(:),MPL_Tauxm1(:)  !< Onsager's Matrices with factors for local(x) field
      real(dbl), allocatable :: mat_f0(:,:),mat_fd(:,:)     !< Onsager's total matrices needed for scf, free_energy and propagation
      real(dbl), allocatable :: ONS_ff(:)                   !< Onsager spherical propagation factors corresponding to $Q_f$
      real(dbl):: ONS_fw                                 !< Onsager spherical propagation factors corresponding to $Q_\omega$
      real(dbl):: ONS_f0,ONS_fd                          !< Onsager spherical propagation factors corresponding to $Q_0$ and $Q_d$
      real(dbl):: ONS_fx0,ONS_fxd                        !< Onsager spherical propagation factors for local field
      real(dbl):: ONS_taum1,ONS_tauxm1                   !< Onsager spherical time constants
      real(dbl) :: sgn                                   !< discriminates between BEM equations for solvent and nanoparticle
      !complex(cmp) :: eps,eps_f                          !< drl complex eps(\omega) and (eps(\omega)-1)/(eps(\omega)+2)
      real(dbl), allocatable :: lambda(:,:)              !< depolarizing factors (3,sph_surf_nsph)
      complex(cmp), allocatable :: q_omega(:) !< Medium charges in frequency domain
      complex(cmp), allocatable :: Kdiag_omega(:)         !< Diagonal K matrix in frequency domain
      real(dbl), allocatable :: scrd3(:) ! Scratch vector dim=3

      real(dbl), allocatable :: BEM_2ppDA(:,:),BEM_2ppDAx(:,:)
      real(dbl), allocatable :: BEM_Sm1(:,:)
      real(dbl), allocatable :: BEM_ADt(:,:)

      real(dbl), allocatable :: gg(:), w2(:), kf(:)

      real(dbl), allocatable :: sin_delta(:), cos_delta(:), kf_prime(:), kf0(:)

      real(dbl), allocatable :: fact3(:),fact3x(:)
      integer(i4b)           :: npoles    !< number of poles when general dielectric function is used
      real(dbl), allocatable :: BEM_Qdf(:,:), BEM_Qdfx(:,:)
      real(dbl), allocatable :: BEM_Qdf_2g(:,:), BEM_Qdfx_2g(:,:)
      type(poles_t) :: poles_eps
      type(poles_t), allocatable :: poles(:)

      save
      private
      public BEM_L,BEM_T,ONS_ff,ONS_fw,                                &
             BEM_Sm12,MPL_F0,MPL_Ft0,MPL_Fd,MPL_Fx0,MPL_Ftx0,MPL_Fxd,  &
             MPL_Tauxm1,MPL_Taum1,mat_f0,mat_fd,MPL_Ff,MPL_Fw,         &
             ONS_f0,ONS_fd,ONS_taum1,ONS_fx0,ONS_fxd,ONS_tauxm1,       &
             BEM_Qt,BEM_R,BEM_Qw,BEM_Qf,BEM_Qd,BEM_Q0,BEM_W2,BEM_Modes,&
             BEM_Qtx,BEM_Qfx,BEM_Qdx,BEM_Q0x,mpibcast_pole_file,       &
             do_BEM_prop,do_BEM_freq,do_BEM_quant,do_MPL_prop,npoles,  &
             do_charge_freq,out_gcharges,read_pole_file,poles_eps,     &
             deallocate_BEM_public,deallocate_MPL_public,BEM_Qg,BEM_2G,&
             BEM_ADt,kf,w2,gg,kf_prime,BEM_Qdf,BEM_Qdfx,BEM_Qdf_2g,    &
             BEM_Qdfx_2g,kf0,deallocate_BEM_end_propagation,BEM_ADtm1, &
             clean_all_ocpy_BEM

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  DRIVER  ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------
! @brief BEM driver routine for propagation
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine do_BEM_prop

       real(dbl), allocatable :: Sm1(:,:)  !< $S^{-1}$ Onsager matrix

#ifndef MPI
       tp_myrank=0
#endif
       !Cavity read/write and S D matrices
       call init_BEM
       if(global_eps_Feps.eq."gen".and.global_medium_Fbem.eq."stan") &
          call mpibcast_pole_file()
       if(global_medium_Fbem.eq.'diag') then
         call init_BEM_diagonal
       endif
       if(global_out_Fgamess.eq.'yes') then
         allocate(BEM_Qd(pedra_surf_n_tessere,pedra_surf_n_tessere))
         allocate(BEM_Q0(pedra_surf_n_tessere,pedra_surf_n_tessere))
         if(global_medium_Floc=='loc'.and.global_medium_Fmdm.eq.'csol') then
             allocate(BEM_Qdx(pedra_surf_n_tessere,pedra_surf_n_tessere))
             allocate(BEM_Q0x(pedra_surf_n_tessere,pedra_surf_n_tessere))
         endif
         !Standard or Diagonal BEM
         if(global_medium_Fbem.eq.'stan') then
           call init_BEM_standard
           call do_BEM_standard
           if (tp_myrank.eq.0)write(6,*) "Standard BEM is experimental"
#ifdef MPI
       call mpi_finalize(tp_ierr_mpi)
#endif
         elseif(global_medium_Fbem.eq.'diag') then
             call do_BEM_diagonal
         endif
         !Write out matrices for gamess
         call out_BEM_gamess
         if(global_medium_Floc=='loc') then
             call out_BEM_lf
         end if
         call finalize_BEM
       endif
       if((global_prop_Fprop.eq."chr-ief").or.&
          (global_prop_Fprop.eq."chr-ied").or.&
          (global_prop_Fprop.eq."chr-ons")) then
         if(.not.allocated(BEM_Qd)) allocate(BEM_Qd(pedra_surf_n_tessere,pedra_surf_n_tessere))
         if(.not.allocated(BEM_Q0)) allocate(BEM_Q0(pedra_surf_n_tessere,pedra_surf_n_tessere))
         if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
           allocate(BEM_Qdx(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Q0x(pedra_surf_n_tessere,pedra_surf_n_tessere))
         endif
       endif
       if((global_prop_Fprop.eq."chr-ief").or.&
          (global_prop_Fprop.eq."chr-ied")) then
         !Standard or Diagonal BEM
         if(global_medium_Fbem.eq.'stan') then
           call init_BEM_standard
           call do_BEM_standard
           if (tp_myrank.eq.0)write(6,*) "Standard BEM is experimental"
#ifdef MPI
       call mpi_finalize(tp_ierr_mpi)
#endif
         elseif(global_medium_Fbem.eq.'diag') then
           call do_BEM_diagonal
           !Save Modes for quantum BEM
           if(global_medium_Fmdm.eq."qnan") then
             allocate(BEM_Modes(pedra_surf_n_tessere,pedra_surf_n_tessere))
             BEM_Modes=TSm12
           endif
         endif
       endif
         !Write out matrices for gamess
         if(global_sys_Fwrite.eq."high") call out_BEM_gamess
         if(global_sys_Fwrite.eq."high") call out_BEM_mat
         !Build propagation Matrices
       if(global_prop_Fprop.ne."non") then
         if(global_eps_Feps.eq."deb") then
           allocate(BEM_R(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qt(pedra_surf_n_tessere,pedra_surf_n_tessere))
           if((global_medium_Floc.eq.'loc').and.&
              (global_medium_Fmdm.eq.'csol')) then
                 allocate(BEM_Qtx(pedra_surf_n_tessere,pedra_surf_n_tessere))
           endif
           if(global_medium_Fbem.eq.'stan') call do_propBEM_std_deb
           if(global_medium_Fbem.eq.'diag') call do_propBEM_dia_deb
         elseif(global_eps_Feps.eq."drl") then
           allocate(BEM_Qw(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qf(pedra_surf_n_tessere,pedra_surf_n_tessere))
           if((global_medium_Floc.eq.'loc').and.&
              (global_medium_Fmdm.eq.'csol')) then
                 allocate(BEM_Qfx(pedra_surf_n_tessere,pedra_surf_n_tessere))
           endif
           if(global_medium_Fbem.eq.'stan') call do_propBEM_std_drl
           if(global_medium_Fbem.eq.'diag') call do_propBEM_dia_drl
         elseif(global_eps_Feps.eq."gen") then
           allocate(BEM_Qg(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qw(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qf(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qdf(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qdf_2g(pedra_surf_n_tessere,pedra_surf_n_tessere))
           if((global_medium_Floc.eq.'loc').and.&
              (global_medium_Fmdm.eq.'csol')) then
            allocate(BEM_Qfx(pedra_surf_n_tessere,pedra_surf_n_tessere))
            allocate(BEM_Qdfx(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Qdfx_2g(pedra_surf_n_tessere,pedra_surf_n_tessere))
           endif
           if(global_medium_Fbem.eq.'stan') call do_propBEM_std_gen
           if(global_medium_Fbem.eq.'diag') call do_propBEM_dia_gen
         endif
       endif
       !Write out propagation matrices
         if(global_sys_Fwrite.eq."high") call out_BEM_propmat
         if(global_prop_Fprop.eq."chr-ons") then
          allocate(Sm1(pedra_surf_n_tessere,pedra_surf_n_tessere))
          ! Form $S^{-1}$ matrix
          Sm1=inv(BEM_S)
          pedra_surf_n_spheres=1
          if(global_eps_Feps.eq."deb") then
            call do_propfact_ons_deb
          elseif(global_eps_Feps.eq."drl") then
            allocate(ONS_ff(pedra_surf_n_spheres))
            call do_propfact_ons_drl
          endif
          ! SP: Computing matrices needed for initialization and propgation eq.2 JPCA 2015
          BEM_Q0=-ONS_f0*Sm1
          BEM_Qd=-ONS_fd*Sm1
          if(global_eps_Feps.eq."drl") then
           !allocate(BEM_Qf(pedra_surf_n_tessere,pedra_surf_n_tessere))
           BEM_Qf=-ONS_ff(1)*Sm1
          endif
          if(global_medium_Floc.eq.'loc') then
            BEM_Q0x=ONS_fx0*Sm1
            BEM_Qdx=ONS_fxd*Sm1
          endif
         endif
       if(global_medium_Fmdm.eq."qnan") then
          if (global_medium_Fbem.eq.'diag') then
             allocate(BEM_Qd(pedra_surf_n_tessere,pedra_surf_n_tessere))
             allocate(BEM_Q0(pedra_surf_n_tessere,pedra_surf_n_tessere))
             call do_BEM_diagonal
             allocate(BEM_Modes(pedra_surf_n_tessere,pedra_surf_n_tessere))
             BEM_Modes=TSm12
             call out_gcharges
             write(*,*) "Printing quantum plasmons charges"
          else
             write(*,*) "Quantum nanoparticle requires BEM diagonal"
             write(*,*) "Please specify bem_type=""diag"""
             if (allocated(BEM_Qd)) deallocate(BEM_Qd)
             if (allocated(BEM_Q0)) deallocate(BEM_Q0)
             stop
          endif
       endif

       if (allocated(Sm1)) deallocate(Sm1)
       !Deallocate private arrays
       call finalize_BEM
       return

      end subroutine


!------------------------------------------------------------------------
! @brief BEM driver routine for frequency calculation (old do_freq_mat)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_BEM_freq

       real(dbl), allocatable :: pot(:),pot2(:)
       complex(cmp) :: mu_omega(3)
       integer(i4b):: i,its

       ! Cavity read/write and S D matrices
       call init_BEM
       allocate(BEM_Qd(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(BEM_Q0(pedra_surf_n_tessere,pedra_surf_n_tessere))
       if(global_medium_Floc=='loc'.and.global_medium_Fmdm.eq.'csol') then
         allocate(BEM_Qdx(pedra_surf_n_tessere,pedra_surf_n_tessere))
         allocate(BEM_Q0x(pedra_surf_n_tessere,pedra_surf_n_tessere))
       end if
       ! Using Diagonal BEM
       call init_BEM_diagonal
       call do_BEM_diagonal
       ! Write out matrices
       if(global_out_Fgamess.eq.'yes') then
           call out_BEM_gamess
           ! Write out local-field matrices
           if(global_medium_Floc=='loc'.and.global_medium_Fmdm.eq.'csol') then
               call out_BEM_lf
           end if
       endif
       ! Calculate potential on tesserae
       allocate(pot(pedra_surf_n_tessere))
       !call do_pot_from_field(fmax(:,1),pot)

       pot(:)=zero
       !SC 7/12/2020: modified to either use n homogeneous field...
       if (global_ext_pert_Ftyp.eq."field") then
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
        do its=1,pedra_surf_n_tessere
          pot(its)=pot(its)-global_ext_pert_direction(1)*pedra_surf_tessere(its)%x
          pot(its)=pot(its)-global_ext_pert_direction(2)*pedra_surf_tessere(its)%y
          pot(its)=pot(its)-global_ext_pert_direction(3)*pedra_surf_tessere(its)%z
        enddo
!$OMP ENDDO
!$OMP END PARALLEL
       !SC 7/12/2020: or read in the transition potential for a molecule
       elseif (global_ext_pert_Ftyp.eq."molecule")  then
        open(888,file="ci_trans.inp",status="old")
         do its=1,pedra_surf_n_tessere
          read(888,*) pot(its)
         enddo
        close(888)
! SC: if it is a eet="yes" calculation read in the trans potential
!  of a second molecule, and place in pot2
         if (global_ext_pert_Feet.eq."yes") then
           open(888,file="ci_trans_2.inp",status="old")
           allocate(pot2(pedra_surf_n_tessere))
           do its=1,pedra_surf_n_tessere
             read(888,*) pot2(its)
           enddo
           close(888)
         endif
       elseif (global_ext_pert_Ftyp.eq."from_cipot") then
                call read_molecule_file
                call read_gau_out_medium(global_ext_pert_n_ci+1)
                pot=quantum_vts(:,1,global_ext_pert_nstate+1)
       elseif (global_ext_pert_Ftyp.eq."dipole") then
            call read_molecule_file
            call do_pot_from_dip(mu_trans,pot)

       endif
! 
       allocate(Kdiag_omega(pedra_surf_n_tessere))
       allocate(q_omega(pedra_surf_n_tessere))
       call do_charge_freq(pot,pot2,mu_omega)
       deallocate(pot,q_omega,Kdiag_omega)
       if (global_ext_pert_Feet.eq."yes") deallocate(pot2)
       !Deallocate private arrays
       call finalize_BEM

       return

      end subroutine do_BEM_freq


!------------------------------------------------------------------------
! @brief BEM driver routine for quantum BEM calculation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_BEM_quant

       ! Cavity read/write and S D matrices
       call init_BEM
       allocate(BEM_Qd(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(BEM_Q0(pedra_surf_n_tessere,pedra_surf_n_tessere))
       ! Diagonal BEM
       call init_BEM_diagonal
       call do_BEM_diagonal
       !Save Modes for quantum BEM

       allocate(BEM_Modes(pedra_surf_n_tessere,pedra_surf_n_tessere))
       BEM_Modes=TSm12
       call finalize_BEM
       call out_gcharges
       return

      end subroutine do_BEM_quant

!------------------------------------------------------------------------
! @brief BEM initialization routine: cavity read/write and S and D
! matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_BEM

       integer(i4b)              :: its

#ifndef MPI
       tp_myrank=0
#endif
       allocate(scrd3(3))
       sgn=one
       if(global_medium_Fmdm.eq."cnan".or.global_medium_Fmdm.eq."qnan") sgn=-one
       if (global_medium_read_write.eq.'wri') then
       ! Write out geometric info and stop
         ! Build the cavity/nanoparticle surface
         !if(pedra_surf_Fcav.eq.'fil') then
         !  call read_cavity_full_file
         !elseif(pedra_surf_Fcav.eq.'gms') then
         !  call read_gmsh_file(pedra_surf_Finv)
         !else
         !  if(global_medium_Fmdm(2:4).eq.'sol') call pedra_int('act')
         !  if(global_medium_Fmdm(2:4).eq.'nan') call pedra_int('met')
         !endif
         ! write out the cavity/nanoparticle surface
         !call output_surf
         if (tp_myrank.eq.0) then
            call output_surf
            write(6,*) "Created output file with surface points"
         endif
         ! Build and write out Calderon SD matrices
         allocate(BEM_S(pedra_surf_n_tessere,pedra_surf_n_tessere))
         if (global_prop_Fprop.ne.'chr-ons') allocate(BEM_D(pedra_surf_n_tessere,pedra_surf_n_tessere))
         call do_BEM_SD
         if (tp_myrank.eq.0) call write_BEM_SD
         if (tp_myrank.eq.0) then
           write(6,*) "Matrixes S D have been written out"
         endif
       elseif (global_medium_read_write.eq.'rea') then
       !Read in geometric info and proceed
         !call read_cavity_file
         allocate(BEM_S(pedra_surf_n_tessere,pedra_surf_n_tessere))
         if (global_prop_Fprop.ne.'chr-ons') allocate(BEM_D(pedra_surf_n_tessere,pedra_surf_n_tessere))
         call read_BEM_SD
         if (tp_myrank.eq.0) write(6,*) &
         "BEM surface and Matrixes S D have been read in"
       endif
       if(global_eps_Feps.eq."gen".and.global_medium_Fbem.eq.'stan') then
                    allocate(BEM_ADtm1(pedra_surf_n_tessere,pedra_surf_n_tessere))
                    call read_pole_file
       endif
       if (tp_myrank.eq.0) write(6,*) "BEM correctly initialized"

       return

      end subroutine init_BEM

!------------------------------------------------------------------------
! @brief BEM finalized and deallocation routine
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine finalize_BEM

       if(allocated(scrd3))deallocate(scrd3)
       if(allocated(BEM_S))deallocate(BEM_S)
       if(allocated(BEM_D)) deallocate(BEM_D)
       if(allocated(fact1)) deallocate(fact1)
       if(allocated(fact2)) deallocate(fact2)
       if(allocated(K0)) deallocate(K0)
       if(allocated(Kd)) deallocate(Kd)
       if(allocated(Sp12)) deallocate(Sp12)
       if(allocated(Sm12T)) deallocate(Sm12T)
       if(allocated(TSm12)) deallocate(TSm12)
       if(allocated(TSp12)) deallocate(TSp12)
       if(allocated(fact2x)) deallocate(fact2x)
       if(allocated(K0x)) deallocate(K0x)
       if(allocated(Kdx)) deallocate(Kdx)
       if(allocated(poles)) deallocate(poles)
       if(allocated(BEM_2ppDA)) deallocate(BEM_2ppDA)
       if(allocated(BEM_2ppDAx)) deallocate(BEM_2ppDAx)
       if(allocated(BEM_Sm1)) deallocate(BEM_Sm1)

       return

      end subroutine finalize_BEM


!------------------------------------------------------------------------
! @brief BEM finalized and deallocation routine
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine deallocate_BEM_public

       if((global_prop_Fprop.eq.'chr-ief').or.&
          (global_prop_Fprop.eq.'chr-ied').or.&
          (global_prop_Fprop.eq.'chr-ons')) then
         if(allocated(BEM_Qd)) deallocate(BEM_Qd)
         if(allocated(BEM_Q0)) deallocate(BEM_Q0)
         if(allocated(BEM_Qt)) deallocate(BEM_Qt)
         if(allocated(BEM_R)) deallocate(BEM_R)
         if(allocated(BEM_Qw)) deallocate(BEM_Qw)
         if(allocated(BEM_Qf)) deallocate(BEM_Qf)
         if(allocated(BEM_L)) deallocate(BEM_L)
         if(allocated(BEM_W2)) deallocate(BEM_W2)
         if(allocated(BEM_T)) deallocate(BEM_T)
         if(allocated(BEM_Sm12)) deallocate(BEM_Sm12)
         if(allocated(BEM_Qdx)) deallocate(BEM_Qdx)
         if(allocated(BEM_Q0x)) deallocate(BEM_Q0x)
         if(allocated(BEM_Qtx)) deallocate(BEM_Qtx)
         if(allocated(BEM_Qf)) deallocate(BEM_Qfx)
         if(allocated(BEM_Qg)) deallocate(BEM_Qg)
         if(allocated(BEM_2G)) deallocate(BEM_2G)
         if(allocated(BEM_VLc)) deallocate(BEM_VLc)
         if(allocated(BEM_Lc)) deallocate(BEM_Lc)
         if(allocated(BEM_VRc)) deallocate(BEM_VRc)
       if(allocated(BEM_Qdf)) deallocate(BEM_Qdf)
       if(allocated(BEM_Qdfx)) deallocate(BEM_Qdfx)
       if(allocated(BEM_Qdf)) deallocate(BEM_Qdf_2g)
       if(allocated(BEM_Qdfx)) deallocate(BEM_Qdfx_2g)

       if(allocated(BEM_ADt)) deallocate(BEM_ADt)
       if(allocated(BEM_ADtm1)) deallocate(BEM_ADtm1)

       if(allocated(kf)) deallocate(kf)
       if(allocated(w2)) deallocate(w2)
       if(allocated(gg)) deallocate(gg)
       if(allocated(fact3)) deallocate(fact3)
       if(allocated(fact1)) deallocate(fact1)
       if(allocated(fact2)) deallocate(fact2)


       endif

       return

      end subroutine deallocate_BEM_public



subroutine clean_all_ocpy_BEM
 if(allocated(BEM_S)) deallocate(BEM_S)
 if(allocated(BEM_D)) deallocate(BEM_D)
 if(allocated(BEM_L)) deallocate(BEM_L)
 if(allocated(BEM_T)) deallocate(BEM_T)
 if(allocated(BEM_W2)) deallocate(BEM_W2)
 if(allocated(BEM_Modes)) deallocate(BEM_Modes)
 if(allocated(K0)) deallocate(K0)
 if(allocated(Kd)) deallocate(Kd)
 if(allocated(K0x)) deallocate(K0x)
 if(allocated(Kdx)) deallocate(Kdx)
 if(allocated(fact1)) deallocate(fact1)
 if(allocated(fact2)) deallocate(fact2)
 if(allocated(fact2x)) deallocate(fact2x)
 if(allocated(Sm12T)) deallocate(Sm12T)
 if(allocated(TSm12)) deallocate(TSm12)
 if(allocated(TSp12)) deallocate(TSp12)
 if(allocated(BEM_Sm12)) deallocate(BEM_Sm12)
 if(allocated(Sp12)) deallocate(Sp12)
 if(allocated(BEM_Q0)) deallocate(BEM_Q0)
 if(allocated(BEM_Qd)) deallocate(BEM_Qd)
 if(allocated(BEM_Qt)) deallocate(BEM_Qt)
 if(allocated(BEM_R)) deallocate(BEM_R)
 if(allocated(BEM_Qw)) deallocate(BEM_Qw)
 if(allocated(BEM_Qf)) deallocate(BEM_Qf)
 if(allocated(BEM_Qg)) deallocate(BEM_Qg)
 if(allocated(BEM_2G)) deallocate(BEM_2G)
 if(allocated(BEM_Qtx)) deallocate(BEM_Qtx)
 if(allocated(BEM_Qfx)) deallocate(BEM_Qfx)
 if(allocated(MPL_Ff)) deallocate(MPL_Ff)
 if(allocated(MPL_Fw)) deallocate(MPL_Fw)
 if(allocated(MPL_F0)) deallocate(MPL_F0)
 if(allocated(MPL_Fx0)) deallocate(MPL_Fx0)
 if(allocated(MPL_Ft0)) deallocate(MPL_Ft0)
 if(allocated(MPL_Ftx0)) deallocate(MPL_Ftx0)
 if(allocated(MPL_Fd)) deallocate(MPL_Fd)
 if(allocated(MPL_Fxd)) deallocate(MPL_Fxd)
 if(allocated(MPL_Taum1)) deallocate(MPL_Taum1)
 if(allocated(MPL_Tauxm1)) deallocate(MPL_Tauxm1)
 if(allocated(mat_f0)) deallocate(mat_f0)
 if(allocated(mat_fd)) deallocate(mat_fd)
 if(allocated(ONS_ff)) deallocate(ONS_ff)
 if(allocated(lambda)) deallocate(lambda)
 if(allocated(q_omega)) deallocate(q_omega)
 if(allocated(Kdiag_omega)) deallocate(Kdiag_omega)
 if(allocated(scrd3)) deallocate(scrd3)
 if(allocated(BEM_2ppDA)) deallocate(BEM_2ppDA)
 if(allocated(BEM_2ppDAx)) deallocate(BEM_2ppDAx)
 if(allocated(BEM_Sm1)) deallocate(BEM_Sm1)
 if(allocated(BEM_ADt)) deallocate(BEM_ADt)
 if(allocated(gg)) deallocate(gg)
 if(allocated(w2)) deallocate(w2)
 if(allocated(kf)) deallocate(kf)
 if(allocated(sin_delta)) deallocate(sin_delta)
 if(allocated(cos_delta)) deallocate(cos_delta)
 if(allocated(kf_prime)) deallocate(kf_prime)
 if(allocated(kf0)) deallocate(kf0)
 if(allocated(fact3)) deallocate(fact3)
 if(allocated(fact3x)) deallocate(fact3x)
 if(allocated(BEM_Qdf)) deallocate(BEM_Qdf)
 if(allocated(BEM_Qdfx)) deallocate(BEM_Qdfx)
 if(allocated(BEM_Qdf_2g)) deallocate(BEM_Qdf_2g)
 if(allocated(BEM_Qdfx_2g)) deallocate(BEM_Qdfx_2g)
 if(allocated(poles)) deallocate(poles)
 if(allocated(poles_eps%omega_p)) deallocate(poles_eps%omega_p)
 if(allocated(poles_eps%gamma_p)) deallocate(poles_eps%gamma_p)
 if(allocated(poles_eps%eps_omega_p))deallocate(poles_eps%eps_omega_p)
 if(allocated(poles_eps%re_deps_domega_p))deallocate(poles_eps%re_deps_domega_p)
 if(allocated(poles_eps%im_deps_domega_p))deallocate(poles_eps%im_deps_domega_p)
 if(allocated(poles_eps%A_coeff_p))deallocate(poles_eps%A_coeff_p)
 ONS_fw = 0
 ONS_f0 = 0
 ONS_fd = 0
 ONS_fx0 = 0
 ONS_fxd = 0
 ONS_taum1 = 0
 ONS_tauxm1 = 0
 sgn = 0
 npoles = 0
                                                                                 
end subroutine

!------------------------------------------------------------------------
! @brief Calculate propagation Onsager matrices from factors including
! depolarization
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_MPL_prop

       real(dbl):: tmp(3),m
       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif

       ! SP 05/07/17 Only one cavity!!!
       if(global_medium_Fmdm.eq."csol") sph_surf_nsph=1
       call init_MPL
       if(global_MPL_ord.eq.1) then
       !SPHEROID
         if(sph_surf_Fshape.eq."spho") then
           do i=1,sph_surf_nsph
           ! Determine minor axis versors for spheroids
             ! build a vector not parallel to major unit vector
             m=sph_surf_maj(i)/sph_surf_min(i)
             tmp=zero
             if(sph_surf_vrs(1,1,i).ne.one) tmp(1)=one
             if(sph_surf_vrs(1,1,i).eq.one) tmp(2)=one
             ! build the two minor unit vectors
             sph_surf_vrs(:,2,i)=vprod(sph_surf_vrs(:,1,i),tmp)
             sph_surf_vrs(:,2,i)=sph_surf_vrs(:,2,i)/mdl(sph_surf_vrs(:,2,i))
             sph_surf_vrs(:,3,i)=vprod(sph_surf_vrs(:,1,i),sph_surf_vrs(:,2,i))
             sph_surf_vrs(:,3,i)=sph_surf_vrs(:,3,i)/mdl(sph_surf_vrs(:,3,i))
            ! Determine depolarization factors: Osborn Phys. Rev. 67 (1945), 351.
             if(int(m*100).gt.100) then   ! Prolate spheroid
               lambda(1,i)=-4.d0*pi/(m*m-one)*(m/two/sqrt(m*m-one)* &
                     log((m+sqrt(m*m-one))/(m-sqrt(m*m-one)))-one)
             elseif(int(m*100).lt.100) then    ! Oblate spheroid
               m=1/m
               lambda(1,i)=-4.d0*pi*m*m/(m*m-one)*&
                       (one-one/sqrt(m*m-one)*asin(sqrt(m*m-one)/m))
             else
               if (tp_myrank.eq.0) write(6,*) "This is a Sphere"
               lambda(1,i)=one/three
             endif
             lambda(2,i)=pt5*(one-lambda(1,i))
             lambda(3,i)=pt5*(one-lambda(1,i))
           enddo
           if(global_eps_Feps.eq."deb") call do_propMPL_deb
           if(global_eps_Feps.eq."drl") call do_propMPL_drl
           ! mat_f0 and mat_fd for scf, g_ and prop
           mat_f0=zero
           mat_fd=zero
           do j=1,sph_surf_nsph
             do i=1,3
               mat_f0(:,:)=mat_f0(:,:)+MPL_F0(:,:,i,j)
               mat_fd(:,:)=mat_fd(:,:)+MPL_Fd(:,:,i,j)
             enddo
           enddo
         else
       !SPHERE
         !Build propagation Matrices with factors
           ! mat_f0 and mat_fd for scf, g_ and prop
           if(global_eps_Feps.eq."drl")then
             call do_propfact_ons_drl
             do i=1,sph_surf_nsph
               ONS_ff=ONS_ff*sph_surf_maj(i)**3
               ONS_f0=ONS_f0*sph_surf_maj(i)**3
               ONS_fd=ONS_fd*sph_surf_maj(i)**3
             enddo
            if (tp_myrank.eq.0) then
              write (6,*) "Onsager"
              write (6,*) "eps_0,eps_d",drudel_eps_0,drudel_eps_d
              write (6,*) "lambda",1/three
              write (6,*) "f0",ONS_f0
              write (6,*) "fd",ONS_fd
            endif
           elseif(global_eps_Feps.eq."deb") then
             call do_propfact_ons_deb
             do i=1,sph_surf_nsph
               ONS_f0=ONS_f0/sph_surf_maj(i)**3
               ONS_fd=ONS_fd/sph_surf_maj(i)**3
             enddo
            if (tp_myrank.eq.0) then
              write (6,*) "Onsager"
              write (6,*) "eps_0,eps_d",debye_eps_0,debye_eps_d
              write (6,*) "lambda",1/three
              write (6,*) "f0",ONS_f0
              write (6,*) "fd",ONS_fd
              write (6,*) "tau",1./ONS_taum1
            endif
           endif
           mat_fd=zero
           mat_f0=zero
           do i=1,3
             mat_fd(i,i)=ONS_fd
             mat_f0(i,i)=ONS_f0
           enddo
         endif
       else
         if (tp_myrank.eq.0)write(6,*) "Higher multipoles not implemented "
#ifdef MPI
         call mpi_finalize(tp_ierr_mpi)
#endif
         stop
       endif
       call finalize_MPL

       return

      end subroutine do_MPL_prop


!------------------------------------------------------------------------
! @brief Allocate arrays for multipolar routines
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_MPL

       allocate(mat_fd(3,3))
       allocate(mat_f0(3,3))
       if(sph_surf_Fshape.eq."spho") then
         allocate(lambda(3,sph_surf_nsph))
         allocate(MPL_F0(3,3,3,sph_surf_nsph))
         allocate(MPL_Fd(3,3,3,sph_surf_nsph))
         if(global_eps_Feps.eq."drl") then
           allocate(MPL_Fw(3,sph_surf_nsph))
           allocate(MPL_Ff(3,3,3,sph_surf_nsph))
         elseif(global_eps_Feps.eq."deb") then
           allocate(MPL_Ft0(3,3,3))
           allocate(MPL_Taum1(3))
           if(global_medium_Floc.eq.'loc') then
             allocate(MPL_Fx0(3,3,3))
             allocate(MPL_Fxd(3,3,3))
             allocate(MPL_Ftx0(3,3,3))
             allocate(MPL_Tauxm1(3))
           endif
         endif
       else
         allocate(ONS_ff(sph_surf_nsph))
       endif

       return

      end subroutine init_MPL

!------------------------------------------------------------------------
! @brief Deallocate lambda
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine finalize_MPL

       if(allocated(lambda)) deallocate(lambda)

       return

      end subroutine finalize_MPL


!------------------------------------------------------------------------
! @brief Deallocate MPL arrays
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine deallocate_MPL_public

       if(allocated(ONS_ff)) deallocate(ONS_ff)
       if(allocated(MPL_Fw)) deallocate(MPL_Fw)
       if(allocated(MPL_Ff)) deallocate(MPL_Ff)
       if(allocated(MPL_F0)) deallocate(MPL_F0)
       if(allocated(MPL_Fd)) deallocate(MPL_Fd)
       if(allocated(MPL_Fx0)) deallocate(MPL_Fx0)
       if(allocated(MPL_Fxd)) deallocate(MPL_Fxd)
       if(allocated(MPL_Ft0)) deallocate(MPL_Ft0)
       if(allocated(MPL_Ftx0)) deallocate(MPL_Ftx0)
       if(allocated(MPL_Taum1)) deallocate(MPL_Taum1)
       if(allocated(MPL_Tauxm1)) deallocate(MPL_Tauxm1)

       return

      end subroutine deallocate_MPL_public
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  CORE ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Compute Calderon's S and D matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_BEM_SD

       real(dbl) :: temp
       integer(i4b) :: i,j

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
        do j=1,pedra_surf_n_tessere
          call green_s(i,j,temp)
          BEM_S(i,j)=temp
          if (global_prop_Fprop.ne.'chr-ons') then
            call green_d(i,j,temp)
            BEM_D(i,j)=temp
          endif
        enddo
       enddo
!$OMP enddo
!$OMP END PARALLEL

       return

      end subroutine do_BEM_SD


!------------------------------------------------------------------------
! @brief Calderon D matrix with Purisima Dii elements
! SC: changed to a diagonal value of D_ii that should be more general
! than those for the sphere
!
! @date Created: S. Pipolo
! Modified: S. Corni 30/5/17
!------------------------------------------------------------------------
      subroutine green_d (i,j,value)

       integer(i4b), intent(in):: i,j
       real(dbl), intent(out) :: value
       real(dbl):: dist,sum_d
       integer(i4b) :: k

       if (i.ne.j) then
          scrd3(1)=(pedra_surf_tessere(i)%x-pedra_surf_tessere(j)%x)
          scrd3(2)=(pedra_surf_tessere(i)%y-pedra_surf_tessere(j)%y)
          scrd3(3)=(pedra_surf_tessere(i)%z-pedra_surf_tessere(j)%z)
          dist=sqrt(dot_product(scrd3,scrd3))
          value=dot_product(pedra_surf_tessere(j)%n,scrd3)/dist**3
       else
          sum_d=0.d0

          do k=1,i-1
             scrd3(1)=(pedra_surf_tessere(i)%x-pedra_surf_tessere(k)%x)
             scrd3(2)=(pedra_surf_tessere(i)%y-pedra_surf_tessere(k)%y)
             scrd3(3)=(pedra_surf_tessere(i)%z-pedra_surf_tessere(k)%z)
             dist=sqrt(dot_product(scrd3,scrd3))
             sum_d=sum_d+dot_product(pedra_surf_tessere(k)%n,scrd3)/dist**3*pedra_surf_tessere(k)%area
          enddo

          do k=i+1,pedra_surf_n_tessere
             scrd3(1)=(pedra_surf_tessere(i)%x-pedra_surf_tessere(k)%x)
             scrd3(2)=(pedra_surf_tessere(i)%y-pedra_surf_tessere(k)%y)
             scrd3(3)=(pedra_surf_tessere(i)%z-pedra_surf_tessere(k)%z)
             dist=sqrt(dot_product(scrd3,scrd3))
             sum_d=sum_d+dot_product(pedra_surf_tessere(k)%n,scrd3)/dist**3*pedra_surf_tessere(k)%area
          enddo

          sum_d=-(2.0*pi+sum_d)/pedra_surf_tessere(i)%area
          value=sum_d
          !value=-1.0694*sqrt(4.d0*pi*pedra_surf_tessere(i)%area)/(2.d0* &
          !       pedra_surf_tessere(i)%rsfe)/pedra_surf_tessere(i)%area
       endif

       return

      end subroutine green_d


!------------------------------------------------------------------------
! @brief Calderon S matrix
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine green_s(i,j,value)

       integer(i4b), intent(in):: i,j
       real(dbl), intent(out) :: value
       real(dbl):: dist

       if (i.ne.j) then
         scrd3(1)=(pedra_surf_tessere(i)%x-pedra_surf_tessere(j)%x)
         scrd3(2)=(pedra_surf_tessere(i)%y-pedra_surf_tessere(j)%y)
         scrd3(3)=(pedra_surf_tessere(i)%z-pedra_surf_tessere(j)%z)
         dist=sqrt(dot_product(scrd3,scrd3))
         value=one/dist
       else
         value=1.0694*sqrt(4.d0*pi/pedra_surf_tessere(i)%area)
       endif

       return

      end subroutine green_s

!------------------------------------------------------------------------
! @brief Compute BEM matrices within diagonal approach
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine do_BEM_diagonal

       integer(i4b) :: i,j,k
       real(8), allocatable :: scr1(:,:),scr2(:,:),scr3(:,:)
       real(8), allocatable :: eigt(:,:),eigt_t(:,:)
       real(8), allocatable :: eigv(:),BEM_WI(:)
       complex(cmp), allocatable :: scrc(:,:)
       complex(cmp) :: vr_tot, vl_tot
       integer(i4b) :: init
       real(dbl) :: fac_eps0,fac_epsd,re,im


#ifndef MPI
       tp_myrank=0
#endif
       if (global_eps_Feps.eq."non") then
           call mpi_error("ERROR: in input epsilon_omega = non but you are trying to", &
                          "calculate bem matrices. Something is wrong in your input", &
                          " and there should be a check. Apologies! ")
       endif
       allocate(scr1(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(scr2(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(scr3(pedra_surf_n_tessere,pedra_surf_n_tessere))
   if(global_medium_bem_sym.eq."yes") then
       allocate(eigv(pedra_surf_n_tessere))
       allocate(eigt(pedra_surf_n_tessere,pedra_surf_n_tessere),eigt_t(pedra_surf_n_tessere,pedra_surf_n_tessere))

       ! Form S^1/2 and S^-1/2
       ! Copy the matrix in the eigenvector matrix

       eigt = BEM_S
       call diag_mat(eigt,eigv,pedra_surf_n_tessere)
       if(global_sys_Fwrite.eq."high") then
          if (tp_myrank.eq.0) write(6,*) "S matrix diagonalized "
       endif

       do i=1,pedra_surf_n_tessere
          if(eigv(i).le.0.d0) then
            write(6,*) "WARNING:",i," eig of S is negative or zero!"
            write(6,*) "    It is set to 1e-8"
            write(6,*) "Please use a more symmetric mesh or set bem_symmetric = 'non' "
            eigv(i)=1.d-8
          endif
          scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
       enddo
       eigt_t=transpose(eigt)

       Sp12=matmul(scr1,eigt_t)

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
       enddo
!$OMP ENDDO
!$OMP END PARALLEL
       deallocate(eigv)

       BEM_Sm12=matmul(scr1,eigt_t)

!      Form the S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 , and diagonalize it
       !S^-1/2 D A S^1/2
       deallocate(eigt)
       deallocate(eigt_t)

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=BEM_D(:,i)*pedra_surf_tessere(i)%area
       enddo
!$OMP ENDDO
!$OMP END PARALLEL

       scr3=matmul(BEM_Sm12,scr1)
       scr2=matmul(scr3,Sp12)

       !S^-1/2 D A S^1/2+S^1/2 A D* S^-1/2 and diagonalize

!$OMP PARALLEL
!$OMP DO
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
           BEM_T(i,j)=0.5*(scr2(i,j)+scr2(j,i))
        enddo
       enddo
!$OMP ENDDO
!$OMP END PARALLEL
       deallocate(scr2,scr3)
       call diag_mat(BEM_T,BEM_L,pedra_surf_n_tessere)

       if(global_sys_Fwrite.eq."high") then
          if (tp_myrank.eq.0) then
          write(6,*) "S^-1/2DAS^1/2+S^1/2AD*S^-1/2 matrix diagonalized"
          endif
       endif
       if (global_eps_Feps.eq."deb") then
!       debye dielectric function
         if(debye_eps_0.ne.one) then
           fac_eps0=(debye_eps_0+one)/(debye_eps_0-one)
           K0(:)=(twp-sgn*BEM_L(:))/(twp*fac_eps0-sgn*BEM_L(:))
           ! GG: analogous to K_0 matrix in the case of local-field for solvent external medium
           if((global_medium_Floc.eq.'loc').and.&
              (global_medium_Fmdm.eq.'csol')) then
                  K0x(:)=-(twp+BEM_L(:))/(twp*fac_eps0-BEM_L(:))
           endif
         else
           K0(:)=zero
         endif
         if(debye_eps_d.ne.one) then
           fac_epsd=(debye_eps_d+one)/(debye_eps_d-one)
           Kd(:)=(twp-sgn*BEM_L(:))/(twp*fac_epsd-sgn*BEM_L(:))
           ! GG: analogous to K_d matrix in the case of local-field for solvent external medium
           if((global_medium_Floc.eq.'loc').and.&
              (global_medium_Fmdm.eq.'csol')) then
                  Kdx(:)=-(twp+BEM_L(:))/(twp*fac_epsd-BEM_L(:))
           endif
         else
           Kd=zero
         endif
         ! SP: Need to check the signs of the second part for a debye medium localized in space
         fact1(:)=((twp-sgn*BEM_L(:))*debye_eps_0+twp+BEM_L(:))/ &
                 ((twp-sgn*BEM_L(:))*debye_eps_d+twp+BEM_L(:))/debye_eps_tau
         fact2(:)=K0(:)*fact1(:)
         ! GG: analogous to \tau K_0 matrix in the case of local-field for solvent external medium
         if((global_medium_Floc.eq.'loc').and.&
            (global_medium_Fmdm.eq.'csol')) then
                 fact2x(:)=K0x(:)*fact1(:)
         endif
       elseif (global_eps_Feps.eq."drl") then
!        Drude-Lorentz dielectric function
         Kd=zero
         fact2(:)=(twp-sgn*BEM_L(:))*drudel_eps_A/(two*twp)
! SC: the first eigenvector should be 0 for the NP
         if ((global_medium_Fmdm.eq.'cnan').or.&
             (global_medium_Fmdm.eq.'qnan')) then
                fact2(1)=0.d0
         endif
         ! SC: no spurious negative square frequencies

         do i=1,pedra_surf_n_tessere
           if(fact2(i).lt.0.d0) then
             write(6,*) "WARNING: BEM_W2(",i,") is ", fact2(i)+drudel_eps_w0*drudel_eps_w0
             write(6,*) "   I put it to 1e-8"
             fact2(i)=1.d-8
             BEM_L(i)=-twp
           endif
         enddo

         if (drudel_eps_w0.eq.zero) drudel_eps_w0=1.d-8
         BEM_W2(:)=fact2(:)+drudel_eps_w0*drudel_eps_w0
         K0(:)=fact2(:)/BEM_W2(:)
         ! GG: analogous to K_f and K_0 matrices in the case of local-field for solvent external medium
         if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
           fact2x(:)=-(twp+BEM_L(:))*drudel_eps_A/(two*twp)
           K0x(:)=fact2x(:)/BEM_W2(:)
         endif
       endif

       if(global_sys_Fwrite.eq."high") then
         if (tp_myrank.eq.0) &
             write(6,*) "Done BEM eigenmodes"
       endif
       Sm12T=matmul(BEM_Sm12,BEM_T)

       TSm12=transpose(Sm12T)

       TSp12=matmul(transpose(BEM_T),Sp12)

      ! SC 05/11/2016 write out the transition charges in pqr format
       if(global_sys_Fwrite.eq."high".and.tp_myrank.eq.0) call output_charge_pqr

      ! Do BEM_Q0 and and BEM_Qd


!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*K0(i)

       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Q0=-matmul(scr1,TSm12)
!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*Kd(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qd=-matmul(scr1,TSm12)

       ! GG: analogous to Q_0 and Q_d matrices in the case of
       ! local-field for solvent external medium
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*K0x(i)
        enddo
        BEM_Q0x=-matmul(scr1,TSm12)
        !BEM_Q0x=-mat_mat_mult(scr1,TSm12)
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*Kdx(i)
        enddo
        BEM_Qdx=-matmul(scr1,TSm12)
       endif
       !Print matrices in output
       if(global_sys_Fwrite.eq."high".and.tp_myrank.eq.0) call out_BEM_diagmat
    elseif(global_medium_bem_sym.eq."non") then

        allocate(BEM_WI(pedra_surf_n_tessere))   
        do i=1,pedra_surf_n_tessere
         scr1(:,i)=BEM_D(:,i)*pedra_surf_tessere(i)%area
        enddo

        call diag_mat_nosym(scr1,BEM_L,BEM_WI,scr2,scr3,pedra_surf_n_tessere)
        i=1
       do while (i.le.pedra_surf_n_tessere)
          if (abs(BEM_WI(i)).gt.1.0d-7) then
             k=i
             BEM_VLc(:,k)=cmplx(scr2(:,k),scr2(:,k+1))
             BEM_VRc(:,k)=cmplx(scr3(:,k),scr3(:,k+1))
             k=i+1
             BEM_VLc(:,k)=cmplx(scr2(:,k-1),-scr2(:,k))
             BEM_VRc(:,k)=cmplx(scr3(:,k-1),-scr3(:,k))
             i=i+2
          else
             k=i
             BEM_VLc(:,k)=cmplx(scr2(:,k),0)
             BEM_VRc(:,k)=cmplx(scr3(:,k),0)
             i=i+1
          endif
       enddo

       BEM_VLc=transpose(BEM_VLc)

       do i=1,pedra_surf_n_tessere
           if (BEM_L(i).lt.(-twp)) BEM_L(i)=-twp
           BEM_Lc(i)=cmplx(BEM_L(i),BEM_WI(i))
       enddo

       BEM_Sm1=inv(BEM_S)
       scrc=inv_cmp(BEM_VRc)
       BEM_VLc=scrc
       BEM_VRc=matmul(BEM_Sm1,BEM_VRc)
       if(pedra_surf_n_particles.eq.1.or.global_medium_Fnorm.eq."tot") then
           do i=2,pedra_surf_n_tessere
              vl_tot=0
              vr_tot=0
              do j=1,pedra_surf_n_tessere
                  vl_tot=vl_tot+BEM_VLc(i,j)
                  vr_tot=vr_tot+BEM_VRc(j,i)
              enddo
              do j=1,pedra_surf_n_tessere
                  re=real(BEM_VLc(i,j))-real(vl_tot)/pedra_surf_n_tessere
                  im=aimag(BEM_VLc(i,j))-aimag(vl_tot)/pedra_surf_n_tessere
                  BEM_VLc(i,j)=cmplx(re,im)
                  re=real(BEM_VRc(j,i))-real(vr_tot)/pedra_surf_n_tessere
                  im=aimag(BEM_VRc(j,i))-aimag(vr_tot)/pedra_surf_n_tessere
                  BEM_VRc(j,i)=cmplx(re,im)
              enddo
           enddo
       else
           do k=1,pedra_surf_n_particles
              if (k.eq.1) then 
                      init=pedra_surf_n_particles+1
              else
                      init=pedra_surf_comp(k,2)
              endif
              do i=init,pedra_surf_comp(k,3)
                   vl_tot=0
                   vr_tot=0
                   do j=pedra_surf_comp(k,2),pedra_surf_comp(k,3)
                         vl_tot=vl_tot+BEM_VLc(i,j)
                         vr_tot=vr_tot+BEM_VRc(j,i)
                   enddo
                   do j=1,pedra_surf_n_tessere
                        re=real(BEM_VLc(i,j))-real(vl_tot)/pedra_surf_comp(k,1)
                        im=aimag(BEM_VLc(i,j))-aimag(vl_tot)/pedra_surf_comp(k,1)
                        BEM_VLc(i,j)=cmplx(re,im)
                        re=real(BEM_VRc(j,i))-real(vr_tot)/pedra_surf_comp(k,1)
                        im=aimag(BEM_VRc(j,i))-aimag(vr_tot)/pedra_surf_comp(k,1)
                        BEM_VRc(j,i)=cmplx(re,im)
                   enddo
              enddo
           enddo
       endif

       deallocate(scrc,scr2,scr3,BEM_WI)
    endif

       deallocate(scr1)
       if (tp_myrank.eq.0) write(6,*) "Done BEM diagonal"


       return

      end subroutine

      subroutine deallocate_BEM_end_propagation

          if(allocated(sin_delta)) deallocate (sin_delta)
          if(allocated(cos_delta)) deallocate (cos_delta)
          if(allocated(kf)) deallocate (kf)
          if(allocated(w2)) deallocate (w2)
          if(allocated(gg)) deallocate (gg)
          if(allocated(kf0)) deallocate (kf0)
          if(allocated(kf_prime)) deallocate (kf_prime)

      end subroutine deallocate_BEM_end_propagation




      subroutine do_BEM_standard
!------------------------------------------------------------------------
! @brief Compute BEM matrices within diagonal approach
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------

       integer(i4b) :: i,j
       real(dbl), allocatable :: scr1(:,:),scr2(:,:),scr3(:,:),scr4(:,:)
       real(dbl) :: scrd3(3), dist


#ifndef MPI
       tp_myrank=0
#endif
       if (global_eps_Feps.eq."non") then
           call mpi_error("ERROR: in input epsilon_omega = non but you are trying to", &
                          "calculate bem matrices. Something is wrong in your input", &
                          " and there should be a check. Apologies! ")
       endif
       allocate(scr1(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                scr2(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                scr3(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                scr4(pedra_surf_n_tessere,pedra_surf_n_tessere))
       if (global_eps_Feps.eq."gen") then
               allocate(kf(npoles),w2(npoles),gg(npoles),kf0(npoles),kf_prime(npoles))
               allocate(sin_delta(npoles),cos_delta(npoles))
               sin_delta(:) = zero
               cos_delta(:) = one
               kf_prime(:) = zero
               kf(:) = poles_eps%A_coeff_p(:)/(two*pi)

               w2(:) = poles_eps%omega_p(:)**2+poles_eps%gamma_p(:)**2
               gg(:) = two*poles_eps%gamma_p(:)
               kf0(:) = kf(:) / w2(:)
       endif
       ! Form S^-1 matrix

       BEM_Sm1=inv(BEM_S)

       ! Form -DA
       BEM_ADt = zero
       scr1=zero
       do i=1,pedra_surf_n_tessere
!         scr1(i,i)= -sgn * BEM_D(i,i)*pedra_surf_tessere(i)%area
         scr1(:,i)= -sgn * BEM_D(:,i)*pedra_surf_tessere(i)%area
       enddo

       ! Form transpose DA
       BEM_ADt= -transpose(scr1)
       
       if (global_eps_Feps.eq."gen") then
            scr2=-BEM_ADt
            scr4=kf0(npoles)*scr2
            do i=1,pedra_surf_n_tessere
               scr4(i,i)=1-scr4(i,i)
            enddo
            BEM_ADtm1= inv(scr4)
       else
            scr2 = scr1
       endif

       ! Form 2 pi - DA

       BEM_2ppDA = scr1
       do i=1,pedra_surf_n_tessere
         BEM_2ppDA(i,i)= BEM_2ppDA(i,i) + twp
       enddo

       ! Form eps0 dependent matrix term
       
       if(global_eps_Feps.eq."deb") then
         do i=1,pedra_surf_n_tessere
           scr2(i,i)= scr2(i,i) + twp * (debye_eps_0+one) / (debye_eps_0-one)
         enddo
       elseif(global_eps_Feps.eq."drl") then
         do i=1,pedra_surf_n_tessere
           scr2(i,i)= scr2(i,i) + twp * (drudel_eps_0+one) / (drudel_eps_0-one)
         enddo
       elseif(global_eps_Feps.eq."gen") then
         do i=1,pedra_surf_n_tessere
           scr2(i,i)= scr2(i,i) + one/sum(kf0)
         enddo
       endif
       ! inverse

       scr2 = inv(scr2)

       ! Form Q0
       if ( global_eps_Feps.eq."gen" ) then
               BEM_Q0=-matmul(scr2,matmul(BEM_Sm1,BEM_2ppDA))
       else
               BEM_Q0=-matmul(BEM_Sm1,matmul(scr2,BEM_2ppDA))
       endif

       ! Form epsd dependent matrix term
       if (global_eps_Feps.eq."gen") then
            scr3=-BEM_ADt
       else
            scr3 = scr1
       endif

       if(global_eps_Feps.eq."deb") then
         do i=1,pedra_surf_n_tessere
           scr3(i,i)= scr3(i,i) + twp * (debye_eps_d+one) / (debye_eps_d-one)
         enddo
       elseif(global_eps_Feps.eq."drl") then
         do i=1,pedra_surf_n_tessere
           scr3(i,i)= scr3(i,i) + twp * (drudel_eps_d+one) / (drudel_eps_d-one)
         enddo
       elseif (global_eps_Feps.eq."gen") then
          do i=1,pedra_surf_n_tessere
           scr3(i,i)= scr3(i,i) + twp * (readf_eps_d+one) / (readf_eps_d-one)
         enddo
       endif

         

       ! inverse

       scr3 = inv(scr3)

       ! Form Qd
       if(global_eps_Feps.eq."gen") then
               BEM_Qd=-matmul(scr3,matmul(BEM_Sm1,BEM_2ppDA))
       else
               BEM_Qd=-matmul(BEM_Sm1,matmul(scr3,BEM_2ppDA))
       endif
       ! GG: analogous to Q_0 and Q_d matrices in the case of
       ! local-field for solvent external medium
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        BEM_2ppDAx = scr1
        do i=1,pedra_surf_n_tessere
          BEM_2ppDAx(i,i)= -BEM_2ppDAx(i,i) + twp
        enddo
        if (global_eps_Feps.eq."gen") then
                BEM_Q0x=matmul(scr2,matmul(BEM_Sm1,BEM_2ppDAx))
                BEM_Qdx=matmul(scr3,matmul(BEM_Sm1,BEM_2ppDAx))
        else
                BEM_Q0x=matmul(BEM_Sm1,matmul(scr2,BEM_2ppDAx))
                BEM_Qdx=matmul(BEM_Sm1,matmul(scr3,BEM_2ppDAx))
        endif
       endif

       deallocate(scr1,scr2,scr3,scr4)

       if (tp_myrank.eq.0) write(6,*) "Done BEM general"


       return

      end subroutine

!------------------------------------------------------------------------
! @brief Initialize diagonal BEM
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine init_BEM_diagonal

       allocate(fact1(pedra_surf_n_tessere),fact2(pedra_surf_n_tessere))
       allocate(Kd(pedra_surf_n_tessere),K0(pedra_surf_n_tessere))
        if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        allocate(fact2x(pedra_surf_n_tessere),fact3x(pedra_surf_n_tessere))
        allocate(Kdx(pedra_surf_n_tessere),K0x(pedra_surf_n_tessere))
       endif
       allocate(BEM_L(pedra_surf_n_tessere))
       allocate(BEM_W2(pedra_surf_n_tessere))
       allocate(BEM_2G(pedra_surf_n_tessere))
       allocate(BEM_T(pedra_surf_n_tessere,pedra_surf_n_tessere))
       if(global_medium_bem_sym.eq."yes") then
           allocate(BEM_Sm12(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(Sp12(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(Sm12T(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(TSm12(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(TSp12(pedra_surf_n_tessere,pedra_surf_n_tessere))
       elseif(global_medium_bem_sym.eq."non") then
           allocate(BEM_VRc(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_VLc(pedra_surf_n_tessere,pedra_surf_n_tessere))
           allocate(BEM_Lc(pedra_surf_n_tessere))
       endif

       allocate(fact3(pedra_surf_n_tessere))

       if(global_eps_Feps.eq."gen") then
              allocate(poles(pedra_surf_n_tessere))
       endif
       return

      end subroutine init_BEM_diagonal

      subroutine init_BEM_standard
!------------------------------------------------------------------------
! @brief Initialize standard BEM
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------

       allocate(BEM_Sm1(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                BEM_2ppDA(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                BEM_ADt(pedra_surf_n_tessere,pedra_surf_n_tessere))
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        allocate(BEM_2ppDAx(pedra_surf_n_tessere,pedra_surf_n_tessere))
       endif

       return

      end subroutine


!------------------------------------------------------------------------
! @brief Compute charges in the frequency domain
!
! @date Created: S. Pipolo
! Modified: E. Coccia 4/12/18
!------------------------------------------------------------------------
      subroutine do_charge_freq(pot,pot2,mu_omega)

       real(dbl),       intent(in)  :: pot(:)
       real(dbl),       intent(in)  :: pot2(:)
       complex(cmp),    intent(out) :: mu_omega(3)
       complex(cmp) :: eps, v_eet, mu_ind_abs(3), mu_ind_emi(3)
       real(dbl)  :: gamma_met, shift_met, re_q, im_q, abs_val
       real(dbl) :: a,b,pl_omega_abs,pl_omega_emi, gamma_met_emi, E_tot, E_0
       complex(dbl)  :: field_point(3),q_tot
       integer(4) :: its, i, j, k,grid_dim, idum
       real, allocatable :: grid(:,:)
       character(30) :: namef
       Kdiag_omega(1)=0
       sgn=-one
! SP 16/07/17: added drudel_eps_w0 to Kdiag_omega
       if(drudel_eps_w0.eq.zero) drudel_eps_w0=1.d-8
       if (global_ext_pert_print_field.eq."yes") then
           open(8,file="grid.inp", status="old")
           read(8,*) grid_dim
           allocate(grid(grid_dim,3))
           do i=1,grid_dim
              read(8,*) idum, grid(i,1), grid(i,2), grid(i,3)
              write(namef,'(a12,i0,a4)')"field_point_",i,".dat"
              open(20+i,file=namef,status="unknown")
              write (20+i,*)"freq re(Ex) re(Ey) re(Ez) im(Ex) im(Ey) &
                            im(Ez)  E_tot  E_tot/E_0"
           enddo
           close(8)
       endif
       if (global_ext_pert_print_surface_charges.eq."yes") then
               open(9, file="surface_charges_real.pqr", status="unknown")
               open(10, file="surface_charges_imag.pqr", status="unknown")
               open(11, file="surface_charges_abs.pqr", status="unknown")
       endif

       open(7,file="dipole_freq.dat",status="unknown")
       if(global_ext_pert_Ftyp.eq."field") then
        write (7,*)"freq re(mux) re(muy) re(muz) im(mux) im(muy) im(muz)"
       elseif (global_ext_pert_Ftyp.eq."molecule") then
        if (global_ext_pert_Feet.eq."non") then
           write (7,*)"freq re(mux) re(muy) re(muz) im(mux) im(muy) im(muz) Gamma_met Shift_met "
        elseif (global_ext_pert_Feet.eq."yes") then
           write (7,*) "freq re(mux) re(muy) re(muz) im(mux) im(muy) im(muz) Gamma_met Shift_met Re(V)  Im(V)"
        endif
       endif
!
!
       do i=1,dielectric_func_n_omega

!$OMP PARALLEL
           select case( global_eps_Feps )
           case('deb')
            ! debye eps
            eps = dielectric_func_eps_deb(dielectric_func_omegas(i))
           case('drl')
            ! drude-lorentz eps
            eps = dielectric_func_eps_drl(dielectric_func_omegas(i))
           case('gen')
            eps = dielectric_func_eps_fromfile(dielectric_func_omegas(i))
           case('gold')
            eps = dielectric_func_eps_gold(dielectric_func_omegas(i))
           end select
          do its=2,pedra_surf_n_tessere 
              Kdiag_omega(its)=(twp-sgn*BEM_L(its))/( ((eps+onec)/(eps-onec))*twp -sgn*BEM_L(its))
          enddo
!$OMP END PARALLEL
          
          if(global_medium_bem_sym.eq."yes") then
              q_omega=matmul(BEM_Sm12,pot)
              q_omega=matmul(transpose(BEM_T),q_omega)
              q_omega=Kdiag_omega*q_omega
              q_omega=matmul(BEM_T,q_omega)
              q_omega=-matmul(BEM_Sm12,q_omega)
           elseif(global_medium_bem_sym.eq."non") then
              q_omega=matmul(BEM_VLc,pot)
              q_omega=Kdiag_omega*q_omega
              q_omega=matmul(BEM_VRc,q_omega)
              q_omega=-matmul(BEM_Sm1,q_omega)
           endif
           if(global_medium_Fnorm.eq."tot") then
              q_tot=0
              do its=1,pedra_surf_n_tessere
                 q_tot=q_tot+q_omega(its)
              enddo
              do its=1,pedra_surf_n_tessere
                 re_q=real(q_omega(its))-real(q_tot)/pedra_surf_n_tessere
                 im_q=aimag(q_omega(its))-aimag(q_tot)/pedra_surf_n_tessere
                 q_omega(its)=cmplx(re_q,im_q)
              enddo
           elseif(global_medium_Fnorm.eq."sep") then
              do k=1,pedra_surf_n_particles
                 q_tot=0
                 do its=pedra_surf_comp(k,2),pedra_surf_comp(k,3)
                       q_tot=q_tot+q_omega(its)
                 enddo
                 do its=pedra_surf_comp(k,2),pedra_surf_comp(k,3)
                       re_q=real(q_omega(its))-real(q_tot)/pedra_surf_comp(k,1)
                       im_q=aimag(q_omega(its))-aimag(q_tot)/pedra_surf_comp(k,1)
                       q_omega(its)=cmplx(re_q,im_q)
                 enddo
              enddo
           endif

          if(global_ext_pert_print_surface_charges.eq."yes") then
              write (9,'("Frequency: ",1e15.6)') dielectric_func_omegas(i)
              write (10,'("Frequency: ",1e15.6)') dielectric_func_omegas(i)
              write (11,'("Frequency: ",1e15.6)') dielectric_func_omegas(i)
              do its=1,pedra_surf_n_tessere
                  write (9,'("ATOM ",I6,"  H  ",I6,"  H  ",3F11.3,F15.5,"  1.5")') &
                  its,its,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y, &
                  pedra_surf_tessere(its)%z, real(q_omega(its))
                  write (10,'("ATOM ",I6,"  H  ",I6,"  H  ",3F11.3,F15.5,"  1.5")') &
                  its,its,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y, &
                  pedra_surf_tessere(its)%z, aimag(q_omega(its))
                  abs_val=sqrt(real(q_omega(its))**2+aimag(q_omega(its))**2)
                  write (11,'("ATOM ",I6,"  H  ",I6,"  H  ",3F11.3,F15.5,"  1.5")') &
                  its,its,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y, &
                  pedra_surf_tessere(its)%z, abs_val
              enddo
          endif

          mu_omega=0.d0

!$OMP PARALLEL REDUCTION(+:mu_omega)
!$OMP DO
          do its=1,pedra_surf_n_tessere
             mu_omega(1)=mu_omega(1)+q_omega(its)*(pedra_surf_tessere(its)%x)
             mu_omega(2)=mu_omega(2)+q_omega(its)*(pedra_surf_tessere(its)%y)
             mu_omega(3)=mu_omega(3)+q_omega(its)*(pedra_surf_tessere(its)%z)
          enddo
!$OMP enddo
!$OMP END PARALLEL
          if(global_ext_pert_print_field.eq."yes") then
                do j=1,grid_dim
                    quantum_mol_cc=grid(j,:)
                    call do_field_from_charges_cmp(q_omega, field_point)
                    E_tot = (real(field_point(1))+global_ext_pert_direction(1))**2+(aimag(field_point(1)))**2
                    E_tot = E_tot+(real(field_point(2))+global_ext_pert_direction(2))**2+(aimag(field_point(2)))**2
                    E_tot = E_tot+(real(field_point(3))+global_ext_pert_direction(3))**2+(aimag(field_point(3)))**2
                    E_tot = sqrt(E_tot)
                    E_0 = sqrt(global_ext_pert_direction(1)**2+global_ext_pert_direction(2)**2+global_ext_pert_direction(3)**2)
                    write(20+j,'(9e15.6)') dielectric_func_omegas(i),real(field_point(:)),aimag(field_point(:)), E_tot, E_tot/E_0
                enddo
          endif


          if(global_ext_pert_Ftyp.eq."field") then
              write (7,'(7e15.6)') dielectric_func_omegas(i),real(mu_omega(:)),aimag(mu_omega(:))
          else
              gamma_met=-2.d0*dot_product(aimag(q_omega),pot)
              shift_met=dot_product(real(q_omega),pot)
              if (global_ext_pert_Feet.eq."non") then
               write (7,'(9e15.6)') dielectric_func_omegas(i),real(mu_omega(:)), & 
                                aimag(mu_omega(:)), gamma_met,shift_met
              elseif (global_ext_pert_Feet.eq."yes") then
               v_eet=dot_product(q_omega,pot2)   
               write (7,'(9e15.6)') dielectric_func_omegas(i),real(mu_omega(:)), & 
                                aimag(mu_omega(:)),gamma_met,shift_met, &
                                real(v_eet),aimag(v_eet)
              endif      
              if (global_ext_pert_Ftyp.eq."from_cipot".or.global_ext_pert_Ftyp.eq."dipole") then
                       if (global_ext_pert_pl_type.eq."from_cienergy") then
                           if (dielectric_func_omegas(i).gt.(tomega-0.0002).and.dielectric_func_omegas(i).lt.(tomega+0.0002)) then
                               mu_ind_abs= mu_omega
                               mu_ind_emi= mu_omega
                               call do_pl_intensity(gamma_met,mu_ind_abs,mu_ind_emi)
                           endif
                       elseif (global_ext_pert_pl_type.eq."from_omega") then
                            pl_omega_emi=global_ext_pert_pl_omega_emi
                            pl_omega_abs=global_ext_pert_pl_omega_abs
                            tomega=pl_omega_emi
                 if (dielectric_func_omegas(i).gt.(pl_omega_emi-0.0002).and.dielectric_func_omegas(i).lt.(pl_omega_emi+0.0002)) then
                               gamma_met_emi=gamma_met
                               mu_ind_emi = mu_omega
                             endif
                 if (dielectric_func_omegas(i).gt.(pl_omega_abs-0.0002).and.dielectric_func_omegas(i).lt.(pl_omega_abs+0.0002)) then
                               mu_ind_abs = mu_omega
                               call do_pl_intensity(gamma_met_emi,mu_ind_abs,mu_ind_emi)
                           endif
                       endif
               endif    
          endif

       enddo

       do i=1,grid_dim
            close(20+i)
       enddo
       close(7)
       if(allocated(grid)) deallocate(grid)

       return

      end subroutine do_charge_freq


!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (debye)
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine do_propBEM_dia_deb

       integer(i4b) :: i
       real(8), allocatable :: scr1(:,:)

       allocate(scr1(pedra_surf_n_tessere,pedra_surf_n_tessere))
!      Form the \tilde{Q} and R for debye propagation

!$OMP PARALLEL
!$OMP DO
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact1(i)
        enddo
!$OMP enddo
!$OMP END PARALLEL

        BEM_R=matmul(scr1,TSp12)

!$OMP PARALLEL
!$OMP DO
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact2(i)
        enddo
!$OMP enddo
!$OMP END PARALLEL

        BEM_Qt=-matmul(scr1,TSm12)
        ! GG: analogous to \tilde{Q} matrix in the case of local-field
        ! for solvent external medium
        if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
         do i=1,pedra_surf_n_tessere
           scr1(:,i)=Sm12T(:,i)*fact2x(i)
         enddo
         BEM_Qtx=-matmul(scr1,TSm12)
        endif

        deallocate(scr1)

        return

      end subroutine do_propBEM_dia_deb

      subroutine do_propBEM_std_deb
!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (debye)
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------

       real(dbl) :: factor

!      Form the \tilde{Q} and R for debye propagation

       factor = (debye_eps_0-one)/(debye_eps_d-one)/debye_eps_tau

       BEM_R=  factor * matmul(BEM_Qd,inv(BEM_Q0))
       BEM_Qt= factor * matmul(BEM_Q0,matmul(inv(BEM_Qd),BEM_Q0))

       ! GG: analogous to \tilde{Q} matrix in the case of local-field
       ! for solvent external medium
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
          BEM_Qtx= factor * matmul(BEM_Q0x,matmul(inv(BEM_Qdx),BEM_Q0x))
       endif
       return

      end subroutine


!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (drude-lorentz)
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine do_propBEM_dia_drl


       integer(i4b) :: i
       real(8), allocatable :: scr1(:,:)

       allocate(scr1(pedra_surf_n_tessere,pedra_surf_n_tessere))
!      Form the Q_w and Q_f for drude-lorentz propagation

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_W2(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qw=matmul(scr1,TSp12)

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*fact2(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qf=-matmul(scr1,TSm12)
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact2x(i)
        enddo
        BEM_Qfx=-matmul(scr1,TSm12)
       endif

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_2G(i)*fact3(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qdf_2g=-matmul(scr1,TSm12)
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact3x(i)
        enddo
        BEM_Qdfx_2g=-matmul(scr1,TSm12)
       endif

      ! addition with respect to do_propBEM_dia_drl
!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_2G(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qg=matmul(scr1,TSp12)

       deallocate(scr1)

       return

      end subroutine do_propBEM_dia_drl

      subroutine do_propBEM_std_drl
!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (drude-lorentz)
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------

      real(dbl) :: factor
      integer(i4b) :: i

!      Form the Q_w and Q_f for drude-lorentz propagation

      factor = -drudel_eps_A/(twp*two)

       BEM_Qw= factor * BEM_2ppDA
       do i=1,pedra_surf_n_tessere
        BEM_Qw(i,i)=BEM_Qw(i,i) + drudel_eps_w0*drudel_eps_w0
       enddo

       BEM_Qf= factor * matmul(BEM_Sm1,BEM_2ppDA)

       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') BEM_Qfx= -factor * matmul(BEM_Sm1,BEM_2ppDAx)

       return

      end subroutine


!------------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (general dielectric function)
!
! @date Created: G. Gil
! Modified:
! Notes: Taken from do_propBEM_dia_drl and building up also BEM_Qg
!------------------------------------------------------------------------------
      subroutine do_propBEM_dia_gen

       integer(i4b) :: i
       real(8), allocatable :: scr1(:,:)

       allocate(scr1(pedra_surf_n_tessere,pedra_surf_n_tessere))
!      Form the Q_w and Q_f for general dielectric function propagation

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_W2(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qw=matmul(scr1,TSp12)

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*fact2(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qf=-matmul(scr1,TSm12)
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact2x(i)
        enddo
        BEM_Qfx=-matmul(scr1,TSm12)
       endif

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*fact3(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qdf=-matmul(scr1,TSm12)
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact3x(i)
        enddo
        BEM_Qdfx=-matmul(scr1,TSm12)
       endif

!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_2G(i)*fact3(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qdf_2g=-matmul(scr1,TSm12)
       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') then
        do i=1,pedra_surf_n_tessere
          scr1(:,i)=Sm12T(:,i)*fact3x(i)
        enddo
        BEM_Qdfx_2g=-matmul(scr1,TSm12)
       endif

      ! addition with respect to do_propBEM_dia_drl
!$OMP PARALLEL
!$OMP DO
       do i=1,pedra_surf_n_tessere
         scr1(:,i)=Sm12T(:,i)*BEM_2G(i)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       BEM_Qg=matmul(scr1,TSp12)

       deallocate(scr1)

       return

      end subroutine do_propBEM_dia_gen

      subroutine do_propBEM_std_gen
!------------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (general dielectric function)
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------------

!      Form the Q_f for general dielectric function propagation

       BEM_Qf= -matmul(BEM_Sm1,BEM_2ppDA)

       if(global_medium_Floc.eq.'loc'.and.global_medium_Fmdm.eq.'csol') BEM_Qfx= matmul(BEM_Sm1,BEM_2ppDAx)

       return

      end subroutine

!------------------------------------------------------------------------
! @brief Initialize factors for debye dipole propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_propMPL_deb

       real(dbl):: k1,f1,g(3,3),m
       integer(i4b):: i,j,k,l

       ! build geometric factors due to orientations
       MPL_Fd=zero
       MPL_F0=zero
       MPL_Taum1=zero
       MPL_Ft0=zero
       if(global_medium_Floc.eq."loc") then
         MPL_Fx0=zero
         MPL_Fxd=zero
         MPL_Tauxm1=zero
         MPL_Ftx0=zero
       endif
       do i=1,sph_surf_nsph
         m=sph_surf_maj(i)/sph_surf_min(i)
         do j=1,3 !Spheroids maj/min
           g=zero
           do l=1,3 !xyz
             do k=1,3 !xyz
               g(k,l)=sph_surf_vrs(k,j,i)*sph_surf_vrs(l,j,i)
             enddo
           enddo
           ! Add depolarized reaction factors 3l/ab^2*(eps-1)/(eps+l/(1-l))
           ! Onsager JACS 58 (1936), Buckingham Trans.Far.Soc. 49 (1953), Abbott Trans.Far.Soc. 48 (1952),
           k1=lambda(j,i)/(one-lambda(j,i))
           f1=three*lambda(j,i)/sph_surf_maj(i)/sph_surf_min(i)**2
           MPL_F0(:,:,j,i)=g(:,:)*f1*(debye_eps_0-one)/(debye_eps_0+k1)
           MPL_Fd(:,:,j,i)=g(:,:)*f1*(debye_eps_d-one)/(debye_eps_d+k1)
           ! f0/tau for propagation equations
           MPL_Ft0(:,:,j)=g(:,:)*f1*(debye_eps_0-one)/(debye_eps_d+k1)/debye_eps_tau
           ! SP 29/06/17: changed from debye_eps_tau to 1/debye_eps_tau
           MPL_Taum1(j)=(debye_eps_0+k1)/(debye_eps_d+k1)/debye_eps_tau
           if(global_medium_Floc.eq."loc") then
             ! Add depolarized local factors eps/(eps+(1-eps)*l), 3eps/(2eps+1) for spheres
             ! Sihvola Journal of Nanomaterials 2007
             MPL_Fx0(:,:,j)=g(:,:)*debye_eps_0/(debye_eps_0+(one-debye_eps_0)*lambda(j,i))
             MPL_Fxd(:,:,j)=g(:,:)*debye_eps_d/(debye_eps_d+(one-debye_eps_d)*lambda(j,i))
             MPL_Ftx0(:,:,j)=g(:,:)*debye_eps_0/(debye_eps_d+(1-debye_eps_d)*lambda(j,i))&
                                         /debye_eps_tau
             MPL_Tauxm1(j)=(debye_eps_0+(1-debye_eps_0)*lambda(j,i))/ &
                           (debye_eps_d+(1-debye_eps_d)*lambda(j,i))/debye_eps_tau
           endif
         enddo
       enddo

       return

      end subroutine do_propMPL_deb

!------------------------------------------------------------------------
! @brief Onsager propagation matrix for debye
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_propfact_ons_deb

       ONS_f0=(debye_eps_0-one)/(debye_eps_0+pt5)
       ONS_fd=(debye_eps_d-one)/(debye_eps_d+pt5)
       ONS_fx0=three*debye_eps_0/(two*debye_eps_0+one)
       ONS_fxd=three*debye_eps_d/(two*debye_eps_d+one)
       ONS_taum1=(two*debye_eps_0+one)/(two*debye_eps_d+one)/debye_eps_tau
       !ONS_taum1=(debye_eps_0+pt5)/(debye_eps_d+pt5)/debye_eps_tau
       !ONS_tauxm1=(two*debye_eps_0+one)/(two*debye_eps_d+one)/debye_eps_tau

       return

      end subroutine do_propfact_ons_deb


!------------------------------------------------------------------------
! @brief Initialize factors for drude-lorentz dipole propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_propMPL_drl

       real(dbl):: f1,g(3,3),m
       integer(i4b):: i,j,k,l

       ! build geometric factors due to orientations
       MPL_Fw=zero
       MPL_Ff=zero
       MPL_F0=zero
       MPL_Fd=zero
       ! Add depolarized Onsager factors:
       ! POLARIZABILITY ANALYSIS OF CANONICAL DIELECTRIC AND BIANISOTROPIC SCATTERERS
       ! Juha Avelin Phd thesis
       do i=1,sph_surf_nsph
         m=sph_surf_maj(i)/sph_surf_min(i)
         ! SP 4\pi\epsilon_0=1
         !f1=two*two*pi*sph_surf_maj(i)*sph_surf_min(i)**2
         f1=sph_surf_maj(i)*sph_surf_min(i)**2
         do j=1,3 !Spheroids maj/min
           MPL_Fw(j,i)=drudel_eps_A*lambda(j,i)+drudel_eps_w0*drudel_eps_w0
           g=zero
           do l=1,3 !xyz
             do k=1,3 !xyz
               g(k,l)=sph_surf_vrs(k,j,i)*sph_surf_vrs(l,j,i)
             enddo
           enddo
           MPL_Ff(:,:,j,i)=g(:,:)*f1*drudel_eps_A/three
           MPL_F0(:,:,j,i)=g(:,:)*f1/three*(drudel_eps_0-one)/(one+&
                                     (drudel_eps_0-one)*lambda(j,i))
           MPL_Fd(:,:,j,i)=g(:,:)*f1/three*(drudel_eps_d-one)/(one+&
                                     (drudel_eps_d-one)*lambda(j,i))
         enddo
       enddo

       return

      end subroutine do_propMPL_drl


!------------------------------------------------------------------------
! @brief Onsager propagation matrix for drude-lorentz
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_propfact_ons_drl

       integer(i4b)::i

       ! Form the ONS_fw=-1/tau+w_0^2 and ONS_ff=A/3
       ! SP 4\pi\epsilon_0=1
       ! May be extended to spheres with different \epsilon(\omega)
       ONS_fw=drudel_eps_A/three+drudel_eps_w0*drudel_eps_w0
       ONS_f0=(drudel_eps_0-one)/(drudel_eps_0+two)
       ONS_fd=(drudel_eps_d-one)/(drudel_eps_d+two)
       ONS_ff(:)=drudel_eps_A/three

       return

      end subroutine do_propfact_ons_drl



      subroutine do_poles(sol,npoles,const,j)
!------------------------------------------------------------------------------
! @brief Compute the real part of the poles of the PCM response kernel
!   * real part of the poles - frequencies
!   * imaginary part of the poles - damping parameters
!   * dielectric function at the poles frequencies
!   * derivative of the dielectric function at the poles frequencies
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------------

       implicit none

       type(poles_t), intent(out) :: sol
       integer(i4b),  intent(out) :: npoles
       real(dbl),     intent(in)  :: const
       integer(i4b),  intent(in)  :: j

       real(dbl) :: val_omega

       integer(i4b), parameter :: niter = 1000
       integer(i4b) :: i, k, kp, end_idx(readf_eps_n_omega), iter, ipole

        ! General notes:
        ! 1) We assume first-order Taylor expansion for eps around frequency of the pole.
        !    This means that the function whose roots we are interested in is not Re{eps} but Re{eps} + Im{d eps/d omega} x gamma
        !    Also this means that gamma can be written as gamma =  Im{eps} / Re{d eps/d omega}.
        !    Since gamma is positive by definition, Re{d eps/d omega} > 0.
        ! 2) This is consistent with other assumptions: simple poles and only one gamma per omega.

        ipole = 0
        write(3,*) const
        end_idx = 0
        do i=1,readf_eps_n_omega-1
         ! first approximation to the root / function crossing zero
         if(.not.(((readf_eps_func(i+1)+const.gt.zero) .and. (readf_eps_func(i)+const.lt.zero)) .or. &
                  ((readf_eps_func(i+1)+const.lt.zero) .and. (readf_eps_func(i)+const.gt.zero))     )) cycle
          k = i
          ! Newton method to find roots of discrete functions
          do iter = 1, niter
            !write(*,*) "pt=", i, "iter=", iter, "final pt=", k, abs(readf_eps_dfunc(k)), readf_eps_re_deps(k)
            val_omega = readf_eps_omegas(k) -(readf_eps_func(k)+const)/readf_eps_dfunc(k)
            kp = k
            k = minloc(abs(readf_eps_omegas(:)-val_omega),1)
            if( k == kp ) exit
          enddo
          ! constraint to estimate gamma from first-order Taylor
          if( nint(sign(one,readf_eps_re_deps(k))) .eq. -1 ) cycle
          !write(*,*) "pt=", i, "iter=", iter, "final pt=", k
          ipole = ipole + 1
          end_idx(i) = k
        end do
        npoles=ipole
        !write(55,*) "How many poles per PCM matrix kernel component", j,"?", npoles
        if( npoles .ge. 1 ) then
         ! Since Newton method from different initial points can arrive to the same root, we remove roots considered multiply.
         if( npoles .ge. 2) then
           ipole = 0
           do i=1, readf_eps_n_omega-1
             if( end_idx(i) .eq. 0 ) cycle
             ipole = ipole + 1
             if( ipole .eq. 1 ) then
               k = end_idx(i)
               cycle
             endif
             if( end_idx(i) .eq. k ) then
               ipole = ipole - 1
               end_idx(i) = 0
             else
               k = end_idx(i)
             endif
           enddo
           npoles=ipole
         endif
         !write(55,*) "How many poles per PCM matrix kernel component", j,"?", npoles
         allocate(sol%omega_p(npoles),sol%gamma_p(npoles),sol%A_coeff_p(npoles))
         allocate(sol%eps_omega_p(npoles),sol%re_deps_domega_p(npoles),sol%im_deps_domega_p(npoles))
         ipole = 0
         do i=1,readf_eps_n_omega-1
           if( end_idx(i) .eq. 0 ) cycle
           ipole = ipole + 1
           k = end_idx(i)
           ! computing omega, gamma, eps and derivative of eps
           sol%omega_p(ipole) = readf_eps_omegas(k) -(readf_eps_func(k)+const)/readf_eps_dfunc(k)
           sol%eps_omega_p(ipole) = cmplx(readf_eps_re_deps(k)*(val_omega-readf_eps_omegas(k)) + real(readf_eps_epsilons(k),dbl),&
                                          readf_eps_im_deps(k)*(val_omega-readf_eps_omegas(k)) + dimag(readf_eps_epsilons(k))    )
           sol%re_deps_domega_p(ipole) = readf_eps_re_deps(k)
           sol%im_deps_domega_p(ipole) = readf_eps_im_deps(k)
           sol%gamma_p(ipole) = dimag(readf_eps_epsilons(k))/readf_eps_re_deps(k)
           sol%A_coeff_p(ipole)=one
           write(2,*) j, const, k, sol%omega_p(ipole), sol%gamma_p(ipole), sol%eps_omega_p(ipole),&
                                   sol%re_deps_domega_p(ipole), sol%im_deps_domega_p(ipole), readf_eps_func(k), i
         end do
        else
     !write(55,*) "Warning! No poles for the PCM matrix kernel component", j
        endif

      end subroutine



!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   INPUT/OUTPUT ROUTINES   !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Write out Calderon's D and S matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine write_BEM_SD

       integer(i4b) :: i,j

       open(7,file="mat_SD.inp",status="unknown")
       write(7,*) pedra_surf_n_tessere
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
          if (global_prop_Fprop.eq.'chr-ons') then
            write(7,'(2E26.16)')BEM_S(i,j)
          else
            write(7,'(2E26.16)')BEM_S(i,j),BEM_D(i,j)
          endif
        enddo
       enddo

       close(7)

       return

      end subroutine write_BEM_SD

!------------------------------------------------------------------------
! @brief Read Calderon's D and S matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine read_BEM_SD

       integer(i4b) :: i,j

#ifndef MPI
       tp_myrank=0
#endif

       if (tp_myrank.eq.0) then
          open(7,file="mat_SD.inp",status="old")
          read(7,*) pedra_surf_n_tessere
          do j=1,pedra_surf_n_tessere
             do i=1,pedra_surf_n_tessere
                if (global_prop_Fprop.eq.'chr-ons') then
                   read(7,*) BEM_S(i,j)
                else
                   read(7,*) BEM_S(i,j), BEM_D(i,j)
               endif
             enddo
          enddo
       endif

#ifdef MPI
      call mpi_bcast(pedra_surf_n_tessere,  1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
      call mpi_bcast(BEM_S,    pedra_surf_n_tessere*pedra_surf_n_tessere,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
      if (global_prop_Fprop.ne.'chr-ons') then
         call mpi_bcast(BEM_D, pedra_surf_n_tessere*pedra_surf_n_tessere,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
      endif
#endif

       close(7)

       return

      end subroutine read_BEM_SD

!------------------------------------------------------------------------
! @brief Output BEM matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_BEM_mat

       integer(i4b):: i,j

       open(7,file="BEM_matrices.mat",status="unknown")
       write(7,*) "# Q_0 , Q_d "
       write(7,*) pedra_surf_n_tessere
       do j=1,pedra_surf_n_tessere
        do i=j,pedra_surf_n_tessere
         write(7,'(2E26.16)')BEM_Q0(i,j),BEM_Qd(i,j)
        enddo
       enddo
       close(7)
       write(6,*) "Written out the propagation BEM matrixes"

       return

      end subroutine out_BEM_mat

!------------------------------------------------------------------------
! @brief Output propagation matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_BEM_propmat

       integer(i4b):: i,j

       open(7,file="BEM_propmat.mat",status="unknown")
       if(global_eps_Feps.eq."deb")write(7,*) "# \tilde{Q} , R "
       if(global_eps_Feps.eq."drl")write(7,*) "# Q_w , Q_t "
       write(7,*) pedra_surf_n_tessere
       do j=1,pedra_surf_n_tessere
        do i=j,pedra_surf_n_tessere
         if(global_eps_Feps.eq."deb")write(7,'(2E26.16)')BEM_Qt(i,j),BEM_R(i,j)
         if(global_eps_Feps.eq."drl")write(7,'(2E26.16)')BEM_Qw(i,j),BEM_Qt(i,j)
        enddo
       enddo
       close(7)
       write(6,*) "Written out the propagation BEM matrixes"

       return

      end subroutine out_BEM_propmat


!------------------------------------------------------------------------
! @brief Output propagation matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_BEM_gamess

       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif
       open(7,file="np_bem.mat",status="unknown")
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
         write(7,'(D20.12)') BEM_Q0(i,j)/pedra_surf_tessere(i)%area
        enddo
       enddo
       close(7)
       if (tp_myrank.eq.0)write(6,*) "Written the static matrix for gamess"
       open(7,file="np_bem.mdy",status="unknown")
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
         write(7,'(D20.12)') BEM_Qd(i,j)/pedra_surf_tessere(i)%area
        enddo
       enddo
       close(7)
       if (tp_myrank.eq.0)write(6,*)"Written the dynamic matrix for gamess"

       return

      end subroutine out_BEM_gamess


!------------------------------------------------------------------------
! @brief Output diagonal matrices and frequencies
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_BEM_diagmat

       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif

       open(7,file="BEM_TSpm12.mat",status="unknown")
       write(7,*) "# T , S^1/2, S^-1/2 "
       write(7,*) pedra_surf_n_tessere
       do j=1,pedra_surf_n_tessere
        do i=j,pedra_surf_n_tessere
         write(7,'(3D26.16)')BEM_T(i,j),Sp12(i,j),BEM_Sm12(i,j)
        enddo
       enddo
       close(7)
       open(7,file="BEM_W2L.mat",status="unknown")
       write(7,*) "# \omega^2 (drude-lorentz) , \Lambda, K_0, K_d "
       write(7,*) pedra_surf_n_tessere
       do j=1,pedra_surf_n_tessere
         write(7,'(i10, 4D20.12)') j,BEM_W2(j),BEM_L(j),K0(j),Kd(j)
       enddo
       close(7)
       if (tp_myrank.eq.0) then
       write(6,*) "Written out BEM diagonal matrices in BEM_TSpm12.mat"
       write(6,*) "Written out BEM squared frequencies in BEM_W2L.mat"
       endif

       return

      end subroutine out_BEM_diagmat

!------------------------------------------------------------------------
! @brief Output propagation matrices
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_BEM_lf

       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif

       open(7,file="np_bem.mlf",status="unknown")
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
         write(7,'(D20.12)') BEM_Q0x(i,j)/pedra_surf_tessere(i)%area
        enddo
       enddo
       close(7)
       if (tp_myrank.eq.0)write(6,*) 'Static loc-field mat in gamess form'
       open(7,file="np_bem.mld",status="unknown")
       do j=1,pedra_surf_n_tessere
        do i=1,pedra_surf_n_tessere
         write(7,'(D20.12)') BEM_Qdx(i,j)/pedra_surf_tessere(i)%area
        enddo
       enddo
       close(7)
       if (tp_myrank.eq.0)write(6,*)'Dynamic loc-field mat in gamess form'

       return

      end subroutine out_BEM_lf


!------------------------------------------------------------------------
! @brief Output cavity files
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
     subroutine output_surf

       integer(i4b) :: i
! This routine creates the same files also created by Gaussian

       open(unit=7,file="cavity.inp",status="unknown",form="formatted")
       write (7,*) pedra_surf_n_tessere,pedra_surf_n_spheres
       do i=1,pedra_surf_n_spheres
         write (7,'(3F22.10)') pedra_surf_spheres(i)%x,pedra_surf_spheres(i)%y, &
                               pedra_surf_spheres(i)%z
       enddo
       do i=1,pedra_surf_n_tessere
         write (7,'(4F22.10,D14.5)') pedra_surf_tessere(i)%x,pedra_surf_tessere(i)%y,pedra_surf_tessere(i)%z, &
                               pedra_surf_tessere(i)%area,pedra_surf_tessere(i)%rsfe
       enddo
       close(unit=7)
! SC 28/9/2016: write out cavity in GAMESS ineq format
! WHICH UNITS ARE ASSUMED IN GAMESS? CHECK!!
       open(unit=7,file="np_bem.cav",status="unknown",form="formatted")
       write (7,*) pedra_surf_n_tessere
       do i=1,pedra_surf_n_tessere
         write (7,'(F22.10)') pedra_surf_tessere(i)%x
       enddo
       do i=1,pedra_surf_n_tessere
         write (7,'(F22.10)') pedra_surf_tessere(i)%y
       enddo
       do i=1,pedra_surf_n_tessere
         write (7,'(F22.10)') pedra_surf_tessere(i)%z
       enddo
       do i=1,pedra_surf_n_tessere
         write (7,'(F22.10)') pedra_surf_tessere(i)%area
       enddo
       close(unit=7)

       return

      end subroutine output_surf

      subroutine output_charge_pqr
!------------------------------------------------------------------------
! @brief Output charges in pqr files
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
       integer :: its,i
       character(30) :: fname
       real(dbl) :: area, maxt
       area=sum(pedra_surf_tessere(:)%area)

#ifndef MPI
       tp_myrank=0
#endif

       do i=1,10
        write (fname,'("charge_freq_",I0,".pqr")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) pedra_surf_n_tessere
        do its=1,pedra_surf_n_tessere
          write (7,'("ATOM ",I6," H    H  ",I6,3F11.3,F15.5,"  1.5")') &
               its,its,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y,pedra_surf_tessere(its)%z, &
               Sm12T(its,i+1)/pedra_surf_tessere(its)%area*area
        enddo
        close(unit=7)
        ! SP : eigenvectors
        write (fname,'("eigenvector_",I0,".pdb")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) "MODEL        1"
        maxt=zero
        do its=1,pedra_surf_n_tessere
          if(sqrt(TSm12(i+1,its)*TSm12(i+1,its)).gt.maxt)&
             maxt=sqrt(TSm12(i+1,its)*TSm12(i+1,its))
        enddo
        do its=1,pedra_surf_n_tessere
          write (7,'("ATOM ",I6,"  H   HHH H ",I3,"    ",3F8.3,2F6.2)') &
               its,i,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y,pedra_surf_tessere(its)%z, &
               TSm12(i+1,its)/maxt*10.,TSm12(i+1,its)/maxt*10.
        enddo
        close(unit=7)
       enddo
       do i=pedra_surf_n_tessere-10,pedra_surf_n_tessere-1
        write (fname,'("eigenvector_",I0,".pdb")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) "MODEL        1"
        maxt=zero
        do its=1,pedra_surf_n_tessere
          if(sqrt(TSm12(i+1,its)*TSm12(i+1,its)).gt.maxt)&
             maxt=sqrt(TSm12(i+1,its)*TSm12(i+1,its))
        enddo
        do its=1,pedra_surf_n_tessere
          write (7,'("ATOM ",I6,"  H   HHH H ",I3,"    ",3F8.3,2F6.2)') &
               its,i,pedra_surf_tessere(its)%x,pedra_surf_tessere(its)%y,pedra_surf_tessere(its)%z, &
               TSm12(i+1,its)/maxt*10.,TSm12(i+1,its)/maxt*10.
        enddo
        close(unit=7)
       enddo
       ! SC 31/10/2016: print out total charges associated with eigenvectors
       open(unit=7,file="charge_eigv.dat",status="unknown", &
            form="formatted")
       if (tp_myrank.eq.0) then
       write (6,*) "Total charge associated to each eigenvector written"
       write (6,*) "   in file charge_eigv.dat"
       endif
       do i=1,pedra_surf_n_tessere
         write(7,'(i20, 2D20.12)') i,sum(Sm12T(:,i))
       enddo
       close(unit=7)

       return

      end subroutine output_charge_pqr

      subroutine out_gcharges
!------------------------------------------------------------------------
! @brief Output charges for quantum plasmons
!
! @date Created: J. Fregoni
! Modified:
!------------------------------------------------------------------------
       integer(i4b) :: i,j,p
       character(len=52) :: my_fmt, my_fmt1
       real(dbl),allocatable,dimension(:) :: omega_p,we
       real(dbl), allocatable :: qg(:,:)        !<Charges associated to each mode

       allocate(we(pedra_surf_n_tessere))
       allocate(qg(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(omega_p(pedra_surf_n_tessere))
       do i = 1, global_qmodes_nprint
           omega_p(i)=sqrt(BEM_W2(i))
           we(i)=sqrt((omega_p(i)**2-drudel_eps_w0**2)/(two*omega_p(i)))
           qg(i,:)=BEM_Modes(i,:)*we(i)
       enddo
       we(1)=zero
       qg(1,:)=zero
       open(9,file="qnp_charges.dat",status="unknown")
       write(my_fmt,'(a,i0,a)') "(",global_qmodes_nprint+3,"E15.6)"
       write(my_fmt1,'(a,i0,a)') "(A22,",global_qmodes_nprint+3,"E15.6)"
       write(9,my_fmt1) "# Plasmon_Frequencies",(omega_p(i),i=2,global_qmodes_nprint)
       write(9,*) "# Modes: x y z q_m1 q_m2 .... q_mN   with   N = ", global_qmodes_nprint
       do j=1,pedra_surf_n_tessere
         write(9,my_fmt) pedra_surf_tessere(j)%x,pedra_surf_tessere(j)%y,pedra_surf_tessere(j)%z,(qg(i,j),i=1,global_qmodes_nprint)
       enddo
         close(9)

       !JF 13/11/2019 Output charges for external QM coupling in .pqr,
       !trajectory like
       open(23,file="qnp_charges.pqr",status="unknown")
       write(my_fmt,'(a,i0,a)') "(",global_qmodes_nprint,"E20.6)"
       do p=2,global_qmodes_nprint
         write(23,*) "mode number = ",p,"   Size = ", pedra_surf_n_tessere, my_fmt
         do j=1,pedra_surf_n_tessere
         write(23,'("ATOM ",I6," H    H  ",I6,3F11.3,3X,E13.6,2X,"1.5")')&
             j,j,pedra_surf_tessere(j)%x,pedra_surf_tessere(j)%y,pedra_surf_tessere(j)%z,qg(p,j)
         enddo
       enddo
       write(*,*) "Charges printed ok"
       close(23)
       !JF Includes mopac print format for charges in gmop.mat file
       if(global_qmodes_Fmop.eq."yes") then
        open(20,file="qnp_mop.dat",status="unknown")
        write(20,*) "#sphere center ?"
        write(20,*) "xcoord   ycoord   zcoord  area   ",(p,p=2,global_qmodes_nprint)
          do j=1,pedra_surf_n_tessere
            write(20,'(4F11.3,3X,100000(ES16.6E3,3X))')&
            pedra_surf_tessere(j)%x,pedra_surf_tessere(j)%y,pedra_surf_tessere(j)%z,&
            pedra_surf_tessere(j)%area,(qg(p,j),p=2,global_qmodes_nprint)
          enddo
          write(*,*) "Mopac Charges printed"
       endif
       close(20)
       deallocate(qg,we,omega_p)
      return
      end subroutine

      subroutine read_pole_file
!------------------------------------------------------------------------
! @brief Read pole file when general dielectric function is used
!
! @date Created: 11/09/2020 G. Dall'Osto
! Modified:
!------------------------------------------------------------------------
      integer(i4b) :: i
#ifndef MPI
      tp_myrank=0
#endif
       if (tp_myrank.eq.0) then
        open(4,file="poles.inp")
        read(4,*) npoles

        allocate(poles_eps%omega_p(npoles),poles_eps%gamma_p(npoles),&
                 poles_eps%re_deps_domega_p(npoles),poles_eps%im_deps_domega_p(npoles),poles_eps%A_coeff_p(npoles))

        do i=1,npoles
           read(4,*) poles_eps%omega_p(i),poles_eps%gamma_p(i),poles_eps%A_coeff_p(i)
           !        poles_eps%re_deps_domega_p(i),poles_eps%im_deps_domega_p(i)
        enddo
        close(4)
       endif
       write(*,*) 'Pole file has been read'

       end subroutine

       subroutine mpibcast_pole_file
!------------------------------------------------------------------------
! @brief MPI BCAST pole information when general dielectric function is used
!
! @date Created: 11/09/2020 G. Dall'Osto
! Modified:
!------------------------------------------------------------------------

#ifdef MPI
           call mpi_bcast(npoles,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
           if(tp_myrank.ne.0) then
               allocate(poles_eps%omega_p(npoles),poles_eps%gamma_p(npoles),&
                      poles_eps%re_deps_domega_p(npoles),poles_eps%im_deps_domega_p(npoles),poles_eps%A_coeff_p(npoles))
           endif
           call mpi_bcast(poles_eps%omega_p,    npoles,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
           call mpi_bcast(poles_eps%gamma_p,    npoles,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
           call mpi_bcast(poles_eps%A_coeff_p,    npoles,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
#endif

      end subroutine


!------------------------------------------------------------------------
! @brief Computation of photoluminescence quantities
!
! @date Created: G. Dall'Osto
! Modified:
!------------------------------------------------------------------------

       subroutine do_pl_intensity(gamma_met,mu_ind_abs, mu_ind_emi)

         real(dbl), intent(in) :: gamma_met
         complex(cmp), intent(in) :: mu_ind_abs(3),mu_ind_emi(3)
         complex(cmp) :: mu_tot(3)
         real(dbl) ::  sp_gam, quantum_yield, tmom2, tmom0, sp_gam_0, quantum_yield_0
         real(dbl) ::  intensity, absorbance, intensity_0, absorbance_0, ratio
         integer(i4b) :: n_ci, nstate
         real(dbl) :: mu_vac(3)

         n_ci=global_ext_pert_n_ci+1
         nstate=global_ext_pert_nstate+1
         mu_vac=global_ext_pert_mu_vac

         mu_tot(1) = mu_trans(1)+mu_ind_emi(1)
         mu_tot(2) = mu_trans(2)+mu_ind_emi(2)
         mu_tot(3) = mu_trans(3)+mu_ind_emi(3)
         tmom2 = mu_tot(1)*conjg(mu_tot(1))+mu_tot(2)*conjg(mu_tot(2))+mu_tot(3)*conjg(mu_tot(3))
         sp_gam = 4.d0/(3.d0*clight**3)*tomega**3*tmom2
         quantum_yield = sp_gam/(sp_gam+gamma_met+global_ext_pert_gamma_vac)
  
         mu_tot(1) = mu_trans(1)+mu_ind_abs(1)
         mu_tot(2) = mu_trans(2)+mu_ind_abs(2)
         mu_tot(3) = mu_trans(3)+mu_ind_abs(3)
         tmom2 = mu_tot(1)*conjg(mu_tot(1))+mu_tot(2)*conjg(mu_tot(2))+mu_tot(3)*conjg(mu_tot(3))
         absorbance = twp/(3.d0*clight)*tmom2
         intensity = quantum_yield*absorbance
         tmom0 = mu_vac(1)**2 + mu_vac(2)**2 + mu_vac(3)**2
         sp_gam_0 = 4.d0/(3.d0*clight**3)*tomega**3*tmom0
         quantum_yield_0 = sp_gam_0/(sp_gam_0+global_ext_pert_gamma_vac)
         absorbance_0 = twp/(3.d0*clight)*tmom0
         intensity_0 = quantum_yield_0*absorbance_0
         ratio = intensity / intensity_0
         write(*,*) "Quantum yield: ", quantum_yield
         write(*,*) "Quantum yield in vacuum: ", quantum_yield_0
         write(*,*) "Molecular absorbance: ", absorbance
         write(*,*) "Molecular absorbance in vacuum: ", absorbance_0
         write(*,*) "Photoluminescence intensity: ",intensity
         write(*,*) "Photoluminescence intensity in vacuum: ",intensity_0
         write(*,*) "Radiative decay rate: ", sp_gam
         write(*,*) "Radiative decay rate in vacuum : ", sp_gam_0
         write(*,*) "Ratio intensity/intensity in vacuum: ", ratio

      end subroutine do_pl_intensity

      end module
