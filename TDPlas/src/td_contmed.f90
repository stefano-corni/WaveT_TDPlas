      Module td_ContMed
      use tdplas_constants
      use global_tdplas

      use pedra_friends
      use sphe_surface
      use readfile_epsilon
      use drudel_epsilon
      use debye_epsilon
      use MathTools
      use BEM_medium
      use readfile_freq
      !xxx just for a test use scf
      use global_quantum
#ifdef OMP
      use omp_lib
#endif

#ifdef MPI
      use mpi
#endif

      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), allocatable :: q0(:)            !< moved from readio, used only here. Charges at time 0 defined with Finit_mdm, here because used in scf
      real(dbl) :: fr_0(3)                       !< moved from readio, used only here. Reaction field at time 0 defined with Finit_mdm, here because used in scf

!
      real(dbl) :: t !< time variable. Subscripts _t _tp _tp2 ... _0 refer respectively to dynamic variables at times t, t-dt, t-2dt, ... 0.
! Molecular observables and Maxwell potential/field
      ! c_tp: input from propagate, coefficients of states at time tp
      ! f_tp: input from propagate, Maxwell field at time tp

      ! DOWN - THIS BLOCK IS NEEDS TO ELMINATED
      real(dbl) :: f_tp2(3)                                   !< old Maxwell field stored in td_contmed
      real(dbl) :: mu_tp(3),mu_tp2(3),mu_0(3)                 !< molecular dipole
      real(dbl), allocatable :: pot_tp(:),pot_tp2(:),pot_0(:) !< molecular potential on BEM points
      real(dbl), allocatable :: potf_tp(:),potf_tp2(:)        !< Maxwell potential on BEM points
      ! UP   - THIS BLOCK IS NEEDS TO ELMINATED

      !SP mu_0 contains quantum_mut(:,1,1)
      !SP pot_0 contains quantum_vts(:,1,1) for global_prop_Fint="ief" and potential of mu_0 if global_prop_Fint='ons'
! Interaction and medium description
      real(dbl), allocatable :: h_mdm(:,:),h_mdm_0(:,:) !< medium contribution to the hamiltonian
      real(dbl), allocatable :: q_mdm(:)                !< medium total charges at current time
      real(dbl), allocatable :: mu_mdm(:,:)             !< medium total dipole for each spheroid/cavity at current time
      real(dbl) :: mu_mdm_p(3,1)                        !< medium total dipole for each pole at current time
      real(dbl) ::  f_mdm(3)                            !< medium total field on the molecule center of charge
! Medium Propagation varibles: charges and dipoles
      ! Dipoles
      real(dbl), allocatable :: mr_0(:,:)               !< reaction dipole (mr) of the (BEM) medium, one vector for each spheroid
      real(dbl), allocatable :: mr_t(:,:),mr_tp(:,:)    !< reaction dipole (mr) of the (BEM) medium, one vector for each spheroid
      real(dbl), allocatable :: dmr_t(:,:),dmr_tp(:,:)  !< reaction dipole difference mr_t-mr_tp
      real(dbl), allocatable :: fmr_t(:,:),fmr_tp(:,:)  !< force on the reaction dipole (vv propagator)
      real(dbl), allocatable :: mx_t(:,:),mx_tp(:,:)    !< dipole induced by the Maxwell field ("external" dipole - mx)
      real(dbl), allocatable :: dmx_t(:,:),dmx_tp(:,:)  !< external dipole difference mx_t-mx_tp
      real(dbl), allocatable :: fmx_t(:,:),fmx_tp(:,:)  !< force on the external dipole (vv propagator)
      ! Charges
      real(dbl), allocatable :: qr_t(:),qr_tp(:)        !< reaction BEM charges (qr)
      real(dbl), allocatable :: dqr_t(:),dqr_tp(:)      !< reaction charge difference qr_t-qr_tp
      real(dbl), allocatable :: fqr_t(:),fqr_tp(:)      !< force on the reaction chares (vv propagator)
      real(dbl), allocatable :: qx_t(:),qx_tp(:)        !< charges induced by the Maxwell field ("external" charges - qx)
      real(dbl), allocatable :: dqx_t(:),dqx_tp(:)      !< external charge difference qx_t-qx_tp
      real(dbl), allocatable :: fqx_t(:),fqx_tp(:)      !< force on the external medium dipole (vv propagator)
      ! charges per pole
      real(dbl), allocatable :: qr_t_p(:,:),qr_tp_p(:,:)        !< reaction BEM charges (qr)
      real(dbl), allocatable :: dqr_t_p(:,:),dqr_tp_p(:,:)      !< reaction charge difference qr_t-qr_tp
      real(dbl), allocatable :: fqr_t_p(:,:),fqr_tp_p(:,:)      !< force on the reaction chares (vv propagator)
      real(dbl), allocatable :: fqr_tp2_p(:),fqr_tp3_p(:)       !< force on the reaction chares (vv propagator type 2)
      real(dbl), allocatable :: dfqr_t_p(:),dfqr_tp_p(:)        !< derivative of force on reaction chares (vv propagator type 2)
      real(dbl), allocatable :: qx_t_p(:,:),qx_tp_p(:,:)        !< charges induced by the Maxwell field ("external" charges - qx)
      real(dbl), allocatable :: dqx_t_p(:,:),dqx_tp_p(:,:)      !< external charge difference qx_t-qx_tp
      real(dbl), allocatable :: fqx_t_p(:,:),fqx_tp_p(:,:)      !< force on the external medium dipole (vv propagator)
      real(dbl), allocatable :: fqx_tp2_p(:),fqx_tp3_p(:)       !< force on the external medium dipole (vv propagator type 2)
      real(dbl), allocatable :: dfqx_t_p(:),dfqx_tp_p(:)        !< derivative of force on external medium dipole (vv propagator type 2)
      real(dbl), allocatable :: sum_r(:),sum_x(:)               !< sum of external and reaction charges on poles and tesserae
      ! Fields and potentials
      real(dbl) :: fr_t(3),fr_tp(3)                     !< reaction field on the molecule centre of charge
      real(dbl) :: dfr_t(3)                             !< reaction field difference (fr_t-fr_tp)
      real(dbl) :: fx_t(3),fx_tp(3)                     !< external field (fx) on the molecule centre of charge
!SC 06/02/16: eq and neq free energies;
      real(dbl) :: e_vac,g_eq,g_eq_gs                   !< energy/equilibrium Free energies
      real(dbl) :: g_neq2,g_neq2_0,g_neq_0,g_neq1,g_neq1_part !< Non Equilibrium Free energies
! Working variables
      real(dbl) :: f1,f2,f3,f4,f5 !< constant factors (calculated once and for all) for velocity-verlet propagator (Drude-Lorentz)
      real(dbl), allocatable :: BEM_f1(:,:),BEM_f3(:,:),BEM_f5(:,:) !< matrices (calculated once and for all) for velocity-verlet propagator (gral dielec func)
      real(dbl), allocatable :: std_f1(:),std_f3(:),std_f5(:)
      real(dbl) :: dip(3)                               !< used in a test
      real(dbl) :: ref                                  !< reference value in different debug tests
      real(dbl) :: qtot                                 !< total BEM charge
      real(dbl) :: taum1                                !< onsager tau-1 for charges single-tau propagation
      real(dbl) :: qtot0                                !< total charge at time 0
      integer(i4b) :: file_med=11                       !< medium output file number

      real(dbl), allocatable :: qr_tp2(:),qx_tp2(:)     !< reaction and external BEM charges in the iteration before the last

      save
      private
!SC 07/02/16: added output_gneq
      public init_mdm,prop_mdm,finalize_mdm,qtot,ref,get_gneq, &
             get_ons,get_mdm_dip,set_charges,preparing_for_scf,&
             init_after_scf,mpibcast_readio_mdm,fr_0,q0, &
             set_potential, init_potential, prop_chr, init_charges,&
             get_propagated_charges, get_corrected_propagated_charges,&
             init_vv_propagator,get_qr_fr,deallocate_potential,       &
             finalize_prop, clean_all_ocpy_tdcont

      contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   INTERFACE ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Medium initialization called by WaveT or other programs
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_mdm(mu_t, f_tp, pot_t, potf_t, h_int)

      implicit none

      real(dbl), optional, intent(in)    ::  mu_t(:)         !< (1:3)           - molecular dipole
      real(dbl), optional, intent(in)    ::  f_tp(:)         !< (1:3)           - external field
      real(dbl), optional, intent(in)    :: pot_t(:)         !< (1:pedra_surf_n_tessere)     - molecular potential
      real(dbl), optional, intent(in)    :: potf_t(:)        !< (1:pedra_surf_n_tessere)     - external  potential
      real(dbl)          , intent(inout) :: h_int(quantum_n_ci,quantum_n_ci) !< (1:quantum_n_ci,1:quantum_n_ci) - interaction hamiltonian
      integer(i4b) :: its,i,j
      character(20) :: name_f


! OPEN FILES
      write(name_f,'(a9,i0,a4)') "medium_t_",quantum_n_f,".dat"
      if (global_prop_Fmdm_res.eq.'nonr') then
      if (quantum_Fbin.ne.'bin') then
         open (file_med,file=name_f,status="unknown")
         write(file_med,*) "# step  time  dipole(x) dipole(y) dipole(z)"
      else
         open (file_med,file=name_f,status="unknown",form="unformatted")
      endif
      elseif (global_prop_Fmdm_res.eq.'yesr') then
         open (file_med,file=name_f,status="unknown",position='append')
      !   if (global_prop_Fint.eq.'pcm') then
      !      allocate(qr_i(pedra_surf_n_tessere))
      !      if (global_medium_Floc.eq.'loc') allocate(qx_i(pedra_surf_n_tessere))
      !   endif
      !   call read_medium_restart()
      endif
      allocate(h_mdm(quantum_n_ci,quantum_n_ci),h_mdm_0(quantum_n_ci,quantum_n_ci))
      h_mdm=zero
      h_mdm_0=zero
      if(global_prop_Fprop.eq."dip") then
        f_tp2=f_tp
! SC: First dipole propagation...
        call do_MPL_prop  !in BEM_medium
        call init_dip_and_field(mu_t)
      else
        if(.not.allocated(mu_mdm))allocate(mu_mdm(3,1))
! SC: ...then charges propagation
        call do_BEM_prop !in BEM_medium
! SC 03/05/2016: create a new BEM_Q0=BEM_Qw^-1*BEM_Qf that should avoid
!                spurious charge dynamics for stationary states
!        call init_BEM_Q0
        call init_potential(pot_t,potf_t)
        call init_charges
!EC: restart values
        !if (global_prop_Fmdm_res.eq.'Yesr') then
         !if (global_prop_Fint.eq.'ons') then
        !    fr_t=fr_i
        !    if (global_medium_Floc.eq.'loc') fx_t=fx_i
         !elseif (global_prop_Fint.eq.'pcm') then
        !    qr_t=qr_i
        !    qr_tp=qr_i
        !    if (global_medium_Floc.eq.'loc') then
        !       qx_t=qx_i
        !       qx_tp=qx_i
        !    endif
         !endif
        !endif
! SC: predifine the factors used in the VV propagator, used for
! Drude-Lorentz
      endif
      call init_vv_propagator
      if (global_medium_Fmdm.ne.'vac') call correct_hamiltonian
! SC set the initial values of the solvent component of the 
! neq free energies
     g_neq1=zero
     g_neq2=zero
      if(global_prop_Fmdm_res.eq.'yesr') then
                qr_t=qr_tp2
                qx_t=qx_tp2
                h_int=0
                if(global_sys_Fdeb.ne."off") call do_interaction_h
                h_int(:,:)=h_int(:,:)+h_mdm(:,:)
                qr_t=qr_tp
                qx_t=qx_tp
        endif
      return


      end subroutine init_mdm

!------------------------------------------------------------------------
! @brief Medium propagation called by WaveT or other programs
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine prop_mdm(i, mu_t, f_tp, pot_t, potf_t, h_int)

      implicit none

      real(dbl), optional, intent(in)    ::  mu_t(:)         !< (1:3)           - molecular dipole
      real(dbl), optional, intent(in)    ::  f_tp(:)         !< (1:3)           - external field
      real(dbl), optional, intent(in)    :: pot_t(:)         !< (1:pedra_surf_n_tessere)     - molecular potential
      real(dbl), optional, intent(in)    :: potf_t(:)        !< (1:pedra_surf_n_tessere)     - external  potential
      real(dbl)          , intent(inout) :: h_int(quantum_n_ci,quantum_n_ci) !< (1:quantum_n_ci,1:quantum_n_ci) - interaction hamiltonian
       integer(i4b), intent(IN) :: i
       integer(i4b) :: its,k,j

       ! Propagate medium only every global_prop_n_q timesteps
      !  otherwise, just resum the interaction hamiltonian and exit
       if(mod(i,global_prop_n_q).ne.0) then
         ! SP 230916: added to perform tests on the local field
         h_int(:,:)=h_int(:,:)+h_mdm(:,:)
         return
       endif
       t=(i-1)*quantum_dt

       ! Build the interaction Hamiltonian Reaction/Local with previous charges
       if(global_sys_Fdeb.ne."off") call do_interaction_h
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
       ! SP 18/05/20 test purposes
       if(global_sys_Ftest.eq."n-r") h_mdm=0
       ! Update the interaction Hamiltonian
       h_int(:,:)=h_int(:,:)+h_mdm(:,:)

       if (global_prop_Fprop.eq."dip") then
       ! Dipole propagation:
         mu_tp = mu_t
         !call prop_dip(f_tp)
         call do_gneq(mu_t,quantum_mut,dfr_t,fr_t,fr_0,mat_fd,3,-1)
         if(global_medium_Fmdm.eq."cnan".or.global_medium_Fmdm.eq."qnan") then
           mu_mdm=mr_t
           if((global_sys_Ftest.eq."n-l".or.global_sys_Ftest.eq."s-l").or.(global_sys_Fdeb.eq."off")) mu_mdm=zero
           if(global_medium_Floc.eq."loc") mu_mdm=mu_mdm+mx_t
         endif
       else
       ! Charges propagation:
         call set_potential(pot_t, potf_t)
         qtot=zero
         mu_mdm=zero
         ! SP 26/06/17: MathUtils, do_dip_from_charges updates the value in mu_mdm
         call do_dip_from_charges(qr_t,mu_mdm(:,1),qtot)
         if((global_sys_Ftest.eq."n-l".or.global_sys_Ftest.eq."s-l").or.(global_sys_Fdeb.eq."off")) mu_mdm=zero
         if (global_medium_Floc.eq."loc") then
           call do_dip_from_charges(qx_t,mu_mdm(:,1),qtot)
         endif
         ! Calculate Reaction Field from charges
         if (global_prop_Fint.eq.'ons') dfr_t=-fr_t
         call do_field_from_charges(qr_t,fr_t)
         if (global_prop_Fint.eq.'ons') dfr_t=dfr_t+fr_t
         ! Calculate Local Field from external charges
         if(global_medium_Floc.eq."loc") call do_field_from_charges(qx_t,fx_t)
       ! SC calculate free energy:
         if (global_prop_Fint.eq.'ons') then
! SP 10/07/17: commented the following, with global_prop_Fint=ons do_gneq should be changed
           !call do_gneq(mu_t,quantum_mut,dfr_t,fr_t,fr_0,mat_fd,3,-1)
         else
           call do_gneq(pot_t,quantum_vts,dqr_t,qr_t,q0,BEM_Qd,pedra_surf_n_tessere,1)
         endif
       endif
       ! SP 230916: added to perform tests on the local/reaction field
       if(global_sys_Ftest.eq."s-r".or.global_sys_Ftest.eq."n-r") mu_tp=quantum_mut(:,1,1) !Only for debug purposes
       if(global_sys_Ftest.eq."n-r") then
         call do_ref(mu_tp)
       elseif(global_sys_Ftest.eq."n-l".or.global_sys_Ftest.eq."s-r".or.global_sys_Ftest.eq."s-l") then
         call do_ref
       end if
       ! SP 24/02/16  Write output
       if (mod(i,quantum_n_out).eq.0.or.i.eq.1) then
          if (quantum_Fbin.eq.'bin') then
             call out_mdm_bin(i)
          else
             call out_mdm(i)
          endif
       endif
       ! Calculate medium's dipole from charges
       if (global_prop_Fprop.eq."dip") then
         call prop_dip(f_tp)
       else
         call prop_chr
       endif      
       ! EC 28/11/17 Write restart
       if (mod(i,quantum_n_res).eq.0) call wrt_restart_mdm(i)
       return

      end subroutine prop_mdm


!------------------------------------------------------------------------
! @brief Medium finalization called by WaveT or other programs
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------
      subroutine finalize_mdm

       if (allocated(h_mdm))deallocate(h_mdm)
       if (allocated(h_mdm_0))deallocate(h_mdm_0)
       call finalize_prop
       if((global_prop_Fprop.eq."chr-ief").or.&
          (global_prop_Fprop.eq."chr-ied").or.&
          (global_prop_Fprop.eq."chr-ons")) then
         call deallocate_BEM_public
       elseif (global_prop_Fprop.eq."dip".or.global_prop_Fint.eq."ons") then
         call deallocate_MPL_public
       endif
       !if (global_prop_Fmdm_res.eq.'Yesr') then
       !   if (global_prop_Fint.eq.'pcm') then
       !      deallocate(qr_i)
       !      if (global_medium_Floc.eq.'loc') then
       !         deallocate(qx_i)
       !      endif
       !   endif
       !endif

       close (file_med)

       return

      end subroutine finalize_mdm

!------------------------------------------------------------------------
! @brief Solvent contribution to free energies
!
! @date Created: S. Corni
! Modified:
!------------------------------------------------------------------------
      subroutine get_gneq(e_vac_t,g_eqr_t,g_neqr_t,g_neq2_t)

       real(dbl),intent(inout):: e_vac_t,g_eqr_t,g_neqr_t,g_neq2_t

       e_vac_t=e_vac_t+e_vac
       g_eqr_t=g_eq+g_eqr_t
       g_neqr_t=g_neq1+g_neqr_t
       g_neq2_t=g_neq2+g_neq2_t

       return

      end subroutine get_gneq


!------------------------------------------------------------------------
! @brief Medium dipole for spectra
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine get_mdm_dip(mdm_dip)

       real(dbl),intent(out):: mdm_dip(3)

       mdm_dip(:)=mu_mdm(:,1)

       return

      end subroutine get_mdm_dip

!------------------------------------------------------------------------
! @brief Set charges
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine set_charges(q)

       real(dbl),intent(in):: q(pedra_surf_n_tessere)

       qr_t=q

       return

      end subroutine set_charges

!------------------------------------------------------------------------
! @brief  Get the current reaction field charges or Onsager reaction fields
!
! @date Created: S. Corni
! Modified:
!------------------------------------------------------------------------
      subroutine get_qr_fr(q)

       real(dbl),intent(out):: q(:)
  
         if (global_prop_Fint.eq."ons") then
          q=fr_t
         else 
          q=qr_t
         endif
         return
        end subroutine get_qr_fr

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!   INTERACTION  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
!> Build the interaction matrix in the molecular unperturbed state basis
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Build the interaction matrix in the molecular unperturbed state
! basis
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_interaction_h

       integer(i4b):: i,j

#ifndef MPI
       tp_myrank=0
#endif

! SC 19/03/2016: added a temporarty array to avoid summing the local field
! to ther reaction field for onsager
! SC 19/03/2016: moved from prop_mdm (so that h_mdm has values assined only here

       h_mdm(:,:)=zero

       if (global_prop_Fint.eq.'ons') then
         f_mdm(:)=fr_t(:)
         if(global_medium_Floc.eq."loc") f_mdm(:)=f_mdm(:)+fx_t(:)
           h_mdm(:,:)=-h_mdm_0(:,:)-quantum_mut(1,:,:)*f_mdm(1) &
                     -quantum_mut(2,:,:)*f_mdm(2) -quantum_mut(3,:,:)*f_mdm(3)
       elseif(global_prop_Fint.eq.'pcm') then
         q_mdm=qr_t
         if(global_medium_Floc.eq."loc") q_mdm(:)=qr_t(:)+qx_t(:)
! SC 31/10/2016: avoid including interaction with an unwanted net charge
         q_mdm=q_mdm+(qtot0-sum(q_mdm))/pedra_surf_n_tessere
         do j=1,quantum_n_ci
            do i=1,j
               h_mdm(i,j)=-h_mdm_0(i,j)+dot_product(q_mdm(:),quantum_vts(:,i,j))
               h_mdm(j,i)=h_mdm(i,j)
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

      end subroutine do_interaction_h

!------------------------------------------------------------------------
! @brief Build interaction matrix at time 0
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine correct_hamiltonian

       implicit none

       integer(4) :: its,i,j

       if (global_prop_Fint.eq.'ons') then
         h_mdm_0(:,:)=-quantum_mut(1,:,:)*fr_0(1)-quantum_mut(2,:,:)*fr_0(2) &
                                         -quantum_mut(3,:,:)*fr_0(3)
       elseif (global_prop_Fint.eq.'pcm') then
!SC 04/05/2016: updated to be coerent with both GS and SCF initialization
!SC 27/09/2016: corrected bug introduced previously
         q_mdm=q0+(qtot0-sum(q0))/pedra_surf_n_tessere

!$OMP PARALLEL REDUCTION(+:h_mdm_0)
!$OMP DO
         do its=1,pedra_surf_n_tessere
           h_mdm_0(:,:)=h_mdm_0(:,:)+q_mdm(its)*quantum_vts(its,:,:)
         enddo
!$OMP END DO
!$OMP END PARALLEL

       endif

#ifndef MPI
       if(global_sys_Fwrite.eq."high") then
         do i=1,quantum_n_ci
          do j=1,quantum_n_ci
           write(6,*) i,j,h_mdm_0(i,j)
          enddo
         enddo
       endif
#endif
       return

      end subroutine correct_hamiltonian
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialization/deallocation !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Initialize potentials for propaagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_potential(pot_t,potf_t)

       real(dbl), intent(IN) :: pot_t(:)
       real(dbl), intent(IN) :: potf_t(:)
       complex(cmp) :: c_gs(quantum_n_ci)
       integer(i4b) :: its

       allocate(pot_tp(pedra_surf_n_tessere))
       pot_tp = pot_t

       allocate(pot_0(pedra_surf_n_tessere))
! SC 08/04/2016: a routine to test by calculating the potentials from the dipoles
! SP 18/05/20: now done in interface_tdplas.f90 (WaveT module)
       !if (global_sys_Fdeb.eq.'vmu') call do_vts_from_dip
       c_gs(:)=zeroc
       c_gs(1)=onec
! SP 26/06/17: changed to use general MathTools
       if(global_prop_Fint.eq.'ons') then
         call do_dip_from_coeff(c_gs,dip,quantum_n_ci)
         call do_pot_from_dip(dip,pot_0)
       else
         call do_pot_from_coeff(c_gs,pot_0)
       endif
       ! Test for a sudden switch-on of an unpolarizable molecular dipole
       if(global_sys_Ftest.eq."s-r") pot_0=zero
       if(global_sys_Ftest.eq."s-r") pot_tp=zero
       allocate(pot_tp2(pedra_surf_n_tessere))
       pot_tp2=pot_tp
       if(global_medium_Floc.eq."loc") then
       allocate(potf_tp(pedra_surf_n_tessere))
       allocate(potf_tp2(pedra_surf_n_tessere))
         potf_tp=potf_t
         potf_tp2=potf_tp

! SP 09/07/16 commented the following
         !fx_t(:)=zero
       endif

       return

      end subroutine init_potential

      subroutine deallocate_potential
        if(allocated(pot_tp)) deallocate(pot_tp)
        if(allocated(pot_tp2)) deallocate(pot_tp2)
        if(allocated(potf_tp)) deallocate(potf_tp)
        if(allocated(potf_tp2)) deallocate(potf_tp2)
        if(allocated(pot_0)) deallocate(pot_0)
        if(allocated(qr_t)) deallocate(qr_t)
        if(allocated(q_mdm)) deallocate(q_mdm)
        if(allocated(qr_tp)) deallocate(qr_tp)
        if(allocated(dqr_tp)) deallocate(dqr_tp)
        if(allocated(dqr_t)) deallocate(dqr_t)
        if(allocated(qr_tp2)) deallocate(qr_tp2)
        if(allocated(qx_t_p))deallocate(qx_t_p)   
        if(allocated(qx_tp_p))deallocate(qx_tp_p)
        if(allocated(dqx_t_p))deallocate(dqx_t_p)
        if(allocated(dqx_tp_p))deallocate(dqx_tp_p)
        if(allocated(fqx_tp_p))deallocate(fqx_tp_p)
        if(allocated(fqx_t_p))deallocate(fqx_t_p)
        if(allocated(qr_t_p))deallocate(qr_t_p)
        if(allocated(qr_tp_p))deallocate(qr_tp_p)
        if(allocated(dqr_t_p))deallocate(dqr_t_p)
        if(allocated(dqr_tp_p))deallocate(dqr_tp_p)
        if(allocated(fqr_tp_p))deallocate(fqr_tp_p)
        if(allocated(fqr_t_p))deallocate(fqr_t_p)
        if(allocated(sum_r))deallocate(sum_r)
        if(allocated(sum_x))deallocate(sum_x)
        if(allocated(dqr_tp))deallocate(dqr_tp)
        if(allocated(fqr_tp))deallocate(fqr_tp)
        if(allocated(fqr_t))deallocate(fqr_t)
        if(allocated(dqx_tp))deallocate(dqx_tp)
        if(allocated(fqx_tp))deallocate(fqx_tp)
        if(allocated(fqx_t))deallocate(fqx_t)
        if(allocated(qr_tp2))deallocate(qr_tp2)
        if(allocated(qx_tp2))deallocate(qx_tp2)
        if(allocated(std_f1))deallocate(std_f1)
        if(allocated(std_f3))deallocate(std_f3)
        if(allocated(std_f5))deallocate(std_f5)
      end subroutine deallocate_potential

      subroutine clean_all_ocpy_tdcont
       if(allocated(q0))deallocate(q0)   
       if(allocated(pot_tp))deallocate(pot_tp)
       if(allocated(pot_tp2))deallocate(pot_tp2)
       if(allocated(pot_0))deallocate(pot_0)
       if(allocated(potf_tp))deallocate(potf_tp)
       if(allocated(potf_tp2))deallocate(potf_tp2)
       if(allocated(h_mdm))deallocate(h_mdm)
       if(allocated(h_mdm_0))deallocate(h_mdm_0)
       if(allocated(q_mdm))deallocate(q_mdm)
       if(allocated(mu_mdm))deallocate(mu_mdm)
       if(allocated(mr_0))deallocate(mr_0)
       if(allocated(mr_t))deallocate(mr_t)
       if(allocated(mr_tp))deallocate(mr_tp)
       if(allocated(dmr_t))deallocate(dmr_t)
       if(allocated(dmr_tp))deallocate(dmr_tp)
       if(allocated(fmr_t))deallocate(fmr_t)
       if(allocated(fmr_tp))deallocate(fmr_tp)
       if(allocated(mx_t))deallocate(mx_t)
       if(allocated(mx_tp))deallocate(mx_tp)
       if(allocated(dmx_t))deallocate(dmx_t)
       if(allocated(dmx_tp))deallocate(dmx_tp)
       if(allocated(fmx_t))deallocate(fmx_t)
       if(allocated(fmx_tp))deallocate(fmx_tp)
       if(allocated(qr_t))deallocate(qr_t)
       if(allocated(qr_tp))deallocate(qr_tp)
       if(allocated(dqr_t))deallocate(dqr_t)
       if(allocated(dqr_tp))deallocate(dqr_tp)
       if(allocated(fqr_t))deallocate(fqr_t)
       if(allocated(fqr_tp))deallocate(fqr_tp)
       if(allocated(qx_t))deallocate(qx_t)
       if(allocated(qx_tp))deallocate(qx_tp)
       if(allocated(dqx_t))deallocate(dqx_t)
       if(allocated(dqx_tp))deallocate(dqx_tp)
       if(allocated(fqx_t))deallocate(fqx_t)
       if(allocated(fqx_tp))deallocate(fqx_tp)
       if(allocated(qr_t_p))deallocate(qr_t_p)
       if(allocated(qr_tp_p))deallocate(qr_tp_p)
       if(allocated(dqr_t_p))deallocate(dqr_t_p)
       if(allocated(dqr_tp_p))deallocate(dqr_tp_p)
       if(allocated(fqr_t_p))deallocate(fqr_t_p)
       if(allocated(fqr_tp_p))deallocate(fqr_tp_p)
       if(allocated(qx_t_p))deallocate(qx_t_p)
       if(allocated(qx_tp_p))deallocate(qx_tp_p)
       if(allocated(dqx_t_p))deallocate(dqx_t_p)
       if(allocated(dqx_tp_p))deallocate(dqx_tp_p)
       if(allocated(fqx_t_p))deallocate(fqx_t_p)
       if(allocated(fqx_tp_p))deallocate(fqx_tp_p)
       if(allocated(sum_r))deallocate(sum_r)
       if(allocated(sum_x))deallocate(sum_x)
       if(allocated(BEM_f1))deallocate(BEM_f1)
       if(allocated(BEM_f3))deallocate(BEM_f3)
       if(allocated(BEM_f5))deallocate(BEM_f5)
       if(allocated(std_f1))deallocate(std_f1)
       if(allocated(std_f3))deallocate(std_f3)
       if(allocated(std_f5))deallocate(std_f5)
       if(allocated(qr_tp2))deallocate(qr_tp2)
       if(allocated(qx_tp2))deallocate(qx_tp2)
       fr_t = 0
       fr_tp = 0
       dfr_t = 0         
       fx_t= 0
       fx_tp= 0
       e_vac= 0
       g_eq= 0
       g_eq_gs= 0
       g_neq2= 0
       g_neq2_0= 0
       g_neq_0= 0
       g_neq1= 0
       g_neq1_part= 0
       f1= 0
       f2= 0
       f3= 0
       f4= 0
       f5 = 0
       fr_0= 0
       t = 0
       f_tp2= 0
       mu_tp= 0
       mu_tp2= 0
       mu_0= 0
       mu_mdm_p= 0
       f_mdm= 0
       dip= 0
       ref= 0
       qtot = 0
       taum1= 0
       qtot0= 0
       file_med=11 
      end subroutine

      

!------------------------------------------------------------------------
! @brief Initialize charges for propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_charges

       implicit none

       integer(i4b):: its
       real(dbl), allocatable :: qd(:)

       integer(i4b) :: ipoles

#ifndef MPI
       tp_myrank=0
#endif

       allocate(qd(pedra_surf_n_tessere))
       allocate(qr_t(pedra_surf_n_tessere))
       allocate(q_mdm(pedra_surf_n_tessere))
       allocate(qr_tp(pedra_surf_n_tessere))
       allocate(dqr_t(pedra_surf_n_tessere))
       if (.not.allocated(q0)) allocate (q0(pedra_surf_n_tessere))
       ! init the state and the RF before propagation
!SP 29/05/16: pot_0 replaces quantum_vts(:,1,1) to allow treating global_prop_Fprop=ief and global_prop_Fint=ons
       select case(global_medium_Finit)
        case ('vac')
          q0(:)=zero
          qtot0=zero
        case ('fro')
          q0(:)=matmul(BEM_Q0,pot_0)
          qtot0=sum(q0)
        case ('rea')
          call read_charges_gau
       end select
! SC 31/10/2016: in case of nanoparticle, normalize initial charges to zero
       if(global_medium_Fmdm.eq.'cnan'.or.global_medium_Fmdm.eq.'qnan') then
         !q0=q0-qtot0/pedra_surf_n_tessere
         qtot0=zero
       endif
       g_eq_gs=0.5d0*dot_product(q0,pot_0)
       if (tp_myrank.eq.0) then
           write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       endif
       ! see readio_medium for global_prop_Finit_int
       if(global_prop_Finit_int.eq.'nsc') then

!         g_neq_0=0.5*MPL_fd*dot_product(mu_tp,mu_tp) &
!                -MPL_fd*dot_product(mu_tp,mu_0) &
!                +0.5*MPL_fd*dot_product(mu_0,mu_0)
!         g_neq_0=-g_neq_0

! SC in principle a non equilibrium self consistency if debye_eps_d=1 is needed
!    here we use the non self-consitent density but use the correct RF
         qd=matmul(BEM_Qd,pot_0)
         g_neq_0=0.5d0*dot_product(qd,pot_0)
         g_neq2_0=g_neq_0
         g_neq_0=g_neq_0-dot_product(qd,pot_tp)
         qd=matmul(BEM_Qd,pot_tp)
         g_neq_0=g_neq_0+0.5d0*dot_product(qd,pot_tp)
         qr_tp=q0+matmul(BEM_Qd,(pot_tp-pot_0))
       end if
       if (tp_myrank.eq.0) write(6,*) 'Medium contribution to G_neq at &
                                           & t=0:',g_neq_0
       qr_t(:)=qr_tp(:)
       allocate(qr_tp2(pedra_surf_n_tessere))
       qr_tp2(:)=qr_tp(:)
       dqr_t(:)=zero
       if(global_prop_Fint.eq."ons") call do_field_from_charges(qr_t,fr_0)
       if(global_medium_Floc.eq."loc") then
         allocate(qx_t(pedra_surf_n_tessere))
         allocate(qx_tp(pedra_surf_n_tessere))
         allocate(dqx_t(pedra_surf_n_tessere))
         ! SP 09/07/17 changed the following for coherence
         qx_t(:)=zero
         qx_tp(:)=zero
         allocate(qx_tp2(pedra_surf_n_tessere))
         qx_tp2(:)=zero
         dqx_t(:)=zero
       endif
       if(global_eps_Feps.eq."drl") then
         allocate(dqr_tp(pedra_surf_n_tessere))
         allocate(fqr_tp(pedra_surf_n_tessere))
         allocate(fqr_t(pedra_surf_n_tessere))
         dqr_tp(:)=zero
         fqr_tp=-mat_mult(BEM_Qw,qr_t)+mat_mult(BEM_Qf,pot_tp)
         if(global_medium_Floc.eq."loc") then
           allocate(dqx_tp(pedra_surf_n_tessere))
           allocate(fqx_tp(pedra_surf_n_tessere))
           allocate(fqx_t(pedra_surf_n_tessere))
           dqx_tp(:)=zero
           fqx_tp=zero
         endif
       endif
       if(global_eps_Feps.eq."gen".and.global_medium_Fbem.eq."stan" ) then
          allocate(qr_t_p(pedra_surf_n_tessere,npoles))
          allocate(qr_tp_p(pedra_surf_n_tessere,npoles))
          allocate(dqr_t_p(pedra_surf_n_tessere,npoles))
          allocate(dqr_tp_p(pedra_surf_n_tessere,npoles))
          allocate(fqr_tp_p(pedra_surf_n_tessere,npoles))
          allocate(fqr_t_p(pedra_surf_n_tessere,npoles))
          allocate(fqr_tp2_p(pedra_surf_n_tessere))
          allocate(fqr_tp3_p(pedra_surf_n_tessere))
          allocate(dfqr_t_p(pedra_surf_n_tessere))
          allocate(dfqr_tp_p(pedra_surf_n_tessere))
          allocate(sum_r(npoles))
          allocate(sum_x(npoles))
          qr_tp(:)=zero
          do ipoles=1,npoles
             qr_tp_p(:,ipoles)=kf0(ipoles)*(matmul(BEM_Qf,pot_0)+matmul(BEM_ADt,q0) )
             qr_tp(:)=qr_tp(:)+qr_tp_p(:,ipoles)
          enddo
          if (npoles.ne.1) then
                  do ipoles=1,npoles-1
                      fqr_tp_p(:,ipoles)=-w2(ipoles)*qr_tp_p(:,ipoles)+kf(ipoles)*(matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp))
                      fqr_tp_p(:,ipoles)=fqr_tp_p(:,ipoles)-sum(fqr_tp_p(:,ipoles))/pedra_surf_n_tessere
                  enddo
                  if (typ_prop.eq."0") then
                          fqr_tp_p(:,npoles)=-w2(npoles)*qr_tp_p(:,npoles)+kf(npoles)*(matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp))
                          fqr_tp_p(:,npoles)=fqr_tp_p(:,npoles)-sum(fqr_tp_p(:,npoles))/pedra_surf_n_tessere
                  else
                          fqr_tp_p(:,npoles)=matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp)
                  endif
          else
                  ipoles=1
                  fqr_tp_p(:,ipoles)=-w2(ipoles)*qr_tp_p(:,ipoles)+kf(ipoles)*(matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp))
                  fqr_tp_p(:,ipoles)=fqr_tp_p(:,ipoles)-sum(fqr_tp_p(:,ipoles))/pedra_surf_n_tessere
          endif
          fqr_tp2_p(:)=zero
          fqr_tp3_p(:)=zero
          dqr_tp_p(:,:)=zero
          if(global_medium_Floc.eq."loc") then
            allocate(qx_t_p(pedra_surf_n_tessere,npoles))
            allocate(qx_tp_p(pedra_surf_n_tessere,npoles))
            allocate(dqx_t_p(pedra_surf_n_tessere,npoles))
            allocate(dqx_tp_p(pedra_surf_n_tessere,npoles))
            allocate(fqx_tp_p(pedra_surf_n_tessere,npoles))
            allocate(fqx_t_p(pedra_surf_n_tessere,npoles))
            allocate(fqx_tp2_p(pedra_surf_n_tessere))
            allocate(fqx_tp3_p(pedra_surf_n_tessere))
            allocate(dfqx_t_p(pedra_surf_n_tessere))
            allocate(dfqx_tp_p(pedra_surf_n_tessere))
            qx_tp_p(:,:)=zero
            dqx_tp_p(:,:)=zero
            fqx_tp_p(:,:)=zero
            fqx_t_p(:,:)=zero
            fqx_tp2_p(:)=zero
            fqx_tp3_p(:)=zero
            dfqx_t_p(:)=zero
            dfqx_tp_p(:)=zero
          endif
       endif
       deallocate(qd)

       if (global_prop_Fmdm_res .eq.'yesr') then
          call read_medium_restart()
       endif

       return

      end subroutine init_charges

!------------------------------------------------------------------------
! @brief Initialize dipoles and field for propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_dip_and_field(mu_t)

       implicit none

       real(dbl), intent(in) :: mu_t(:)

       mu_tp=mu_t

       allocate(mu_mdm(3,sph_surf_nsph))
       allocate(mr_0(3,sph_surf_nsph))
       allocate(mr_t(3,sph_surf_nsph))
       allocate(mr_tp(3,sph_surf_nsph))
       allocate(dmr_t(3,sph_surf_nsph))
       if(global_medium_Floc.eq."loc") then
         allocate(mx_t(3,sph_surf_nsph))
         allocate(mx_tp(3,sph_surf_nsph))
         allocate(dmx_t(3,sph_surf_nsph))
       endif
       if(global_eps_Feps.eq."drl") then
         allocate(dmr_tp(3,sph_surf_nsph))
         allocate(fmr_t(3,sph_surf_nsph))
         allocate(fmr_tp(3,sph_surf_nsph))
         dmr_tp=zero
         fmr_tp=zero
         if(global_medium_Floc.eq."loc") then
           allocate(dmx_tp(3,sph_surf_nsph))
           allocate(fmx_t(3,sph_surf_nsph))
           allocate(fmx_tp(3,sph_surf_nsph))
           dmx_tp=zero
           fmx_tp=zero
         endif
       endif
       ! The following subroutines should be merged!!
       if(sph_surf_Fshape.eq."sphe") call init_dip_sphe
       if(sph_surf_Fshape.eq."spho") call init_dip_spho

       return

      end subroutine init_dip_and_field
!
!------------------------------------------------------------------------
!> Initialize Dipole for propagation with spherical object
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Initialize dipole for propagation with spherical object
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_dip_sphe

       implicit none

       real(dbl) :: fld(3),fld0(3)
       integer(i4b) :: i

       mu_0(:)=quantum_mut(:,1,1)
       mr_0=zero
       fr_0=zero
       if(global_sys_Ftest.eq."s-r") mu_0=zero
       if(global_sys_Ftest.eq."s-r") mu_tp=zero
!SC 7/5/2018
       mu_tp2=mu_tp
       if (global_medium_Fmdm.eq."cnan".or.global_medium_Fmdm.eq."qnan") then
       ! SP 06/07/17: need to implement a proper scf cycle for more then a nanoparticle
       !              only non-self consistent initialization is done
         ! SP compute the field acting on spheres/oids given m_0=zero
         fr_t=zero
         do i=1,sph_surf_nsph
           call do_field_from_dip(mu_0,quantum_mol_cc,fld0,sph_surf_centre(:,i))
           if(global_medium_Finit.eq."fro") mr_0(:,i)=ONS_f0*fld0
           call do_field_from_dip(mu_tp,quantum_mol_cc,fld,sph_surf_centre(:,i))
! SP 11/07/17: to have the same initial RF for charges and dipoles
           mr_t(:,i)=mr_0(:,i)+ONS_fd*(fld-fld0)
           call do_field_from_dip(mr_t(:,i),sph_surf_centre(:,i),fld,quantum_mol_cc)
           fr_t=fr_t+fld
         enddo
         mr_tp=mr_t
         if(global_medium_Floc.eq."loc") then
           mx_t=zero
           mx_tp=mx_t
         endif
       elseif (global_medium_Fmdm.eq."csol") then
         if(global_medium_Finit.eq."fro") fr_0(:)=ONS_f0*mu_0(:)
         call init_dip_sphe_sol
       endif
       fr_tp=fr_t
       if(global_medium_Floc.eq."loc") then
         fx_t=zero
         fx_tp=fx_t
       endif

       return

      end subroutine init_dip_sphe
!
!------------------------------------------------------------------------
!> Initialize Dipole for propagation with spheroidal object
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Initialize dipole for propagation with spheroidal object
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_dip_spho

      ! inizialize onsager
       implicit none

       real(dbl) :: fld(3),fld0(3)
       integer(i4b) :: i

       ! Init molecular dipole
       mu_0(:)=quantum_mut(:,1,1)
       mr_0=zero
       fr_0=zero
       if(global_sys_Ftest.eq."s-r") mu_0=zero
       if(global_sys_Ftest.eq."s-r") mu_tp=zero
       if (global_medium_Fmdm.eq."cnan".or.global_medium_Fmdm.eq.'qnan') then
       ! SP 06/07/17: need to implement a proper scf cycle for more then a nanoparticle
       !              only non-self consistent initialization is done
         ! SP compute the field acting on spheres/oids given m_0=zero
         fr_t=zero
         do i=1,sph_surf_nsph
           call do_field_from_dip(mu_0,quantum_mol_cc,fld0,sph_surf_centre(:,i))
           if(global_medium_Finit.eq."fro") mr_0(:,i)=matmul(mat_f0,fld0)
           call do_field_from_dip(mu_tp,quantum_mol_cc,fld,sph_surf_centre(:,i))
           mr_t(:,i)=mr_0(:,i)+matmul(mat_fd,fld-fld0)
           call do_field_from_dip(mr_t,sph_surf_centre(:,i),fld,quantum_mol_cc)
           fr_t=fr_t+fld
         enddo
         mr_tp=mr_t
         if(global_medium_Floc.eq."loc") then
           mx_t=zero
           mx_tp=mx_t
         endif
       elseif (global_medium_Fmdm.eq."csol") then
         if(global_medium_Finit.eq."fro") fr_0=matmul(mat_f0,mu_0)
         call init_dip_spho_sol
         fr_tp=fr_t
       endif
       if(global_medium_Floc.eq."loc") then
         fx_t=zero
         fx_tp=fx_t
       endif

       return

      end subroutine init_dip_spho


!------------------------------------------------------------------------
! @brief Initialize free energy (spheroidal)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_dip_spho_sol

       implicit none

#ifndef MPI
       tp_myrank=0
#endif

       g_eq_gs=-0.5d0*dot_product(fr_0,mu_0)
       if (tp_myrank.eq.0) then
          write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       endif
       if (global_prop_Finit_int.eq.'nsc') then
! SC in principle a non equilibrium self consistency if debye_eps_d=1 is needed
!    here we use the non self-consitent dipole but use the correct RF
         fr_t=fr_0+matmul(mat_fd,mu_tp-mu_0)
         g_neq_0=0.5*dot_product(mu_tp,matmul(mat_fd,mu_tp)) &
! SP 04/0717 WARNING: NOT SURE OF THE FOLLOWING LINE
                -dot_product(matmul(mat_fd,mu_tp),mu_0) &
                !-f_d*dot_product(mu_tp,mu_0) &
                +0.5*dot_product(mu_0,matmul(mat_fd,mu_0))
         g_neq_0=-g_neq_0
         g_neq2_0=-0.5*dot_product(mu_0,matmul(mat_fd,mu_0))
       end if
       if (tp_myrank.eq.0) write(6,*) 'Medium contribution to G_neq at t=0:',g_neq_0

       return

      end subroutine init_dip_spho_sol

!------------------------------------------------------------------------
! @brief Initialize free energy (spherical)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_dip_sphe_sol

       implicit none

#ifndef MPI
       tp_myrank=0
#endif

       g_eq_gs=-0.5d0*dot_product(fr_0,mu_0)
       if (tp_myrank.eq.0) then
           write(6,*) 'Medium contribution to ground state free energy:', &
                    g_eq_gs
       endif
       if(global_prop_Finit_int.eq.'nsc') then
! SC in principle a non equilibrium self consistency if debye_eps_d=1 is needed
!    here we use the non self-consitent dipole but use the correct RF
         fr_t=fr_0+ONS_fd*(mu_tp-mu_0)
         g_neq_0=0.5*ONS_fd*dot_product(mu_tp,mu_tp) &
                -ONS_fd*dot_product(mu_tp,mu_0) &
                +0.5*ONS_fd*dot_product(mu_0,mu_0)
         g_neq_0=-g_neq_0
         g_neq2_0=-0.5*ONS_fd*dot_product(mu_0,mu_0)
       end if
       if (tp_myrank.eq.0) write(6,*) 'G_neq at t=0:',g_neq_0

       return

      end subroutine init_dip_sphe_sol

      subroutine preparing_for_scf(mix,pot_or_mut)

       implicit none

       real(dbl) :: mix

       real(dbl), intent(in) :: pot_or_mut(:) 

       if (global_prop_Fprop.eq."dip") then
        if(sph_surf_Fshape.eq."sphe") then
         fr_t=mix*ONS_f0*pot_or_mut+(1.-mix)*fr_0
        else if(sph_surf_Fshape.eq."spho") then
         fr_t=mix*matmul(mat_f0,pot_or_mut)+(1.-mix)*fr_0
        end if
       else
        qr_t=mix*matmul(BEM_Q0,pot_or_mut)+(1.-mix)*q0
       endif

      end subroutine preparing_for_scf

      subroutine init_after_scf(pot_or_mut)

       implicit none
       integer(4) :: ipoles
       real(dbl), intent(in) :: pot_or_mut(:)

       if (global_prop_Fprop.eq."dip") then
        if(sph_surf_Fshape.eq."sphe") then
         g_neq_0=-0.5*ONS_f0*dot_product(pot_or_mut,pot_or_mut)
         g_neq2_0=-0.5*ONS_fd*dot_product(pot_or_mut,pot_or_mut)
        else if(sph_surf_Fshape.eq."spho") then
         g_neq_0=-0.5*dot_product(pot_or_mut,matmul(mat_f0,pot_or_mut))
         g_neq2_0=-0.5*dot_product(pot_or_mut,matmul(mat_fd,pot_or_mut))
        end if
       else
        qr_tp=matmul(BEM_Q0,pot_or_mut)
        g_neq_0=0.5*dot_product(qr_tp,pot_or_mut)
        qr_tp=matmul(BEM_Qd,pot_or_mut)
        g_neq2_0=0.5*dot_product(qr_tp,pot_or_mut)
        qr_tp=matmul(BEM_Q0,pot_or_mut)
        qr_t(:)=qr_tp(:)
        dqr_t(:)=zero
        if (global_eps_Feps.eq."drl") then
            fqr_tp=-mat_mult(BEM_Qw,qr_t)+mat_mult(BEM_Qf,pot_or_mut)
            dqr_tp(:)=zero
        endif
        if(global_eps_Feps.eq."gen".and.global_medium_Fbem.eq."stan") then 
            qr_tp(:)=zero
          do ipoles=1,npoles
           qr_tp_p(:,ipoles)=kf0(ipoles)*(matmul(BEM_Qf,pot_or_mut)+matmul(BEM_ADt,qr_t) )
           qr_tp(:)=qr_tp(:)+qr_tp_p(:,ipoles)
          enddo
          dqr_tp_p(:,:)=zero
          if (npoles.ne.1) then
           do ipoles=1,npoles-1
             fqr_tp_p(:,ipoles)=-w2(ipoles)*qr_tp_p(:,ipoles)+kf(ipoles)*(matmul(BEM_Qf,pot_or_mut)+matmul(BEM_ADt,qr_tp))
             fqr_tp_p(:,ipoles)=fqr_tp_p(:,ipoles)-sum(fqr_tp_p(:,ipoles))/pedra_surf_n_tessere
           enddo
             fqr_tp_p(:,npoles)=matmul(BEM_Qf,pot_or_mut)+matmul(BEM_ADt,qr_tp)
          else
             ipoles=1
             fqr_tp_p(:,ipoles)=-w2(ipoles)*qr_tp_p(:,ipoles)+kf(ipoles)*(matmul(BEM_Qf,pot_or_mut)+matmul(BEM_ADt,qr_tp))
             fqr_tp_p(:,ipoles)=fqr_tp_p(:,ipoles)-sum(fqr_tp_p(:,ipoles))/pedra_surf_n_tessere
          endif
           qr_t=qr_tp
        endif
        if(global_prop_Fint.eq."ons") call do_field_from_charges(qr_t,fr_0)
       endif
       ! SC 27/06/2020: set q0 to the charge after the SCF inside tdplas as well
       q0=qr_t
       ! SC 27/06/2020: recalculate the hamiltonian from q0
       call correct_hamiltonian

       if (tp_myrank.eq.0) write(6,*) 'Medium contribution to G_neq at t=0:',g_neq_0

      end subroutine init_after_scf

!------------------------------------------------------------------------
! @brief Deallocate propagation variables
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine finalize_prop

       integer(i4b) :: its

       ! Molecule
       if(allocated(pot_0)) deallocate (pot_0)
       if(allocated(pot_tp)) deallocate (pot_tp)
       if(allocated(pot_tp2)) deallocate (pot_tp2)
       if(allocated(potf_tp)) deallocate (potf_tp)
       if(allocated(potf_tp2)) deallocate (potf_tp2)
       ! Interaction
       if(allocated(q_mdm)) deallocate(q_mdm)
       if(allocated(mu_mdm)) deallocate(mu_mdm)
       ! Dipole
       if(allocated(mr_0)) deallocate(mr_0)
       if(allocated(mr_t)) deallocate(mr_t)
       if(allocated(mr_tp)) deallocate(mr_tp)
       if(allocated(dmr_t)) deallocate(dmr_t)
       if(allocated(dmr_tp)) deallocate(dmr_tp)
       if(allocated(fmr_t)) deallocate(fmr_t)
       if(allocated(fmr_tp)) deallocate(fmr_tp)
       if(allocated(mx_t)) deallocate(mx_t)
       if(allocated(mx_tp)) deallocate(mx_tp)
       if(allocated(dmx_t)) deallocate(dmx_t)
       if(allocated(dmx_tp)) deallocate(dmx_tp)
       if(allocated(fmx_t)) deallocate(fmx_t)
       if(allocated(fmx_tp)) deallocate(fmx_tp)
       ! Charges
       if(allocated(qr_t)) deallocate(qr_t)
       if(allocated(qr_tp)) deallocate(qr_tp)
       if(allocated(dqr_t)) deallocate(dqr_t)
       if(allocated(dqr_tp)) deallocate(dqr_tp)
       if(allocated(fqr_t)) deallocate(fqr_t)
       if(allocated(fqr_tp)) deallocate(fqr_tp)
       if(allocated(qx_t)) deallocate(qx_t)
       if(allocated(qx_tp)) deallocate(qx_tp)
       if(allocated(dqx_t)) deallocate(dqx_t)
       if(allocated(dqx_tp)) deallocate(dqx_tp)
       if(allocated(fqx_t)) deallocate(fqx_t)
       if(allocated(fqx_tp)) deallocate(fqx_tp)
       return

      end subroutine finalize_prop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! Propagation routines !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SP 22/02/16: A general remark: matmul and DGEMV used randomly,
!             DGEMV shoud be efficient for big matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! @brief Potential propagation                                                                         
!                                                                                                      
! @date Created: M. Rosa                                                                               
! Modified:                                                                                            
!------------------------------------------------------------------------                              
      subroutine set_potential(pot_t,potf_t)                                                           
          implicit none                                                                                
          real(dbl), intent(IN) :: pot_t(:)                                                            
          real(dbl), optional, intent(IN) :: potf_t(:)                                                 
                                                                                                       
          pot_tp  = pot_t                                                                              
          if(global_medium_Floc.eq."loc") then                                                         
              potf_tp = potf_t                                                                         
          endif                 
      end subroutine
!------------------------------------------------------------------------
! @brief Potential and charges propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_chr

       implicit none

#ifndef MPI
       tp_myrank=0
#endif

       ! Propagate
       if(global_sys_Ftest.eq."s-r") pot_tp=quantum_vts(:,1,1) !Only for debug purposes
       if(global_eps_Feps.eq."deb") then
         if(global_prop_Fprop.eq."chr-ief") then
           call prop_ief_deb
         elseif(global_prop_Fprop.eq."chr-ied") then
           taum1=(debye_eps_0+one)/(debye_eps_d+one)/debye_eps_tau
           call prop_ied_deb
         elseif(global_prop_Fprop.eq."chr-ons") then
           call prop_ons_deb
         else
           if (tp_myrank.eq.0) write(*,*) "not implemented yet"
#ifdef MPI
           call mpi_finalize(tp_ierr_mpi)
#endif
           stop
         endif
       ! SP 22/02/16 dqr_t calculated for g_neq
         dqr_t=qr_t-qr_tp
       elseif(global_eps_Feps.eq."drl") then
         if(global_prop_Fprop.eq."chr-ief") then
           call prop_vv_ief_drl
         elseif(global_prop_Fprop.eq."chr-ons") then
           call prop_ons_drl ! uses vv
         else
           if (tp_myrank.eq.0) write(*,*) "Wrong propagation method"
#ifdef MPI
           call mpi_finalize(tp_ierr_mpi)
#endif
           stop
         endif
       elseif(global_eps_Feps.eq."gen") then
         if(global_prop_Fprop.eq."chr-ief") then
           if(global_medium_Fbem.eq."stan") then
            call prop_vv_ief_gen_std
           else
            write(*,*) "Wrong propagation method"
           stop
           endif
         elseif(global_prop_Fprop.eq."chr-ons") then
           ! FIXME: add C-PCM propagation type
           !call prop_ons_gen ! uses vv - not yet
           write(*,*) "Wrong propagation method"
           stop
         else
           write(*,*) "Wrong propagation method"
           stop
         endif
       endif
       if(global_sys_Fdeb.eq."equ") then
!EC:  mat_mult optimizes quantum_n_ci**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of quantum_n_ci
         qr_t=mat_mult(BEM_Q0,pot_tp)
       endif
       qr_tp2=qr_tp
       qr_tp=qr_t
       pot_tp2=pot_tp
       if(global_medium_Floc.eq."loc") then
         qx_tp2=qx_tp
         qx_tp=qx_t
         potf_tp2=potf_tp
       endif
       return

      end subroutine prop_chr


!------------------------------------------------------------------------
! @brief Dipole and field propagation
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_dip(f_tp)

      ! evolve onsager field/dipole
       implicit none

       real(dbl), intent(IN):: f_tp(3)

       ! calculate molecular dipole from CI coefficients
       if(global_sys_Ftest.eq."s-r") mu_tp=quantum_mut(:,1,1) !Only for debug purposes
       ! propagate onsager reaction and local fields
       if(global_eps_Feps.eq."deb") then
         if(sph_surf_Fshape.eq."sphe") call prop_dip_sphe_deb(f_tp)
         if(sph_surf_Fshape.eq."spho") call prop_dip_spho_deb(f_tp)
       elseif(global_eps_Feps.eq."drl") then
         if(sph_surf_Fshape.eq."sphe") call prop_dip_sphe_drl(f_tp)
         if(sph_surf_Fshape.eq."spho") call prop_dip_spho_drl(f_tp)
       endif
       dfr_t=fr_t-fr_tp
       fr_tp=fr_t
       mu_tp2=mu_tp
       if(global_medium_Floc.eq."loc") then
         fx_tp=fx_t
         f_tp2=f_tp
       endif

       return

      end subroutine prop_dip

      subroutine init_vv_propagator

        if(global_eps_Feps.eq."drl") then            
          call init_vv_propagator_drl                
        elseif (global_eps_Feps.eq."gen") then       
           if( global_medium_Fbem.eq."stan" ) then   
            call init_vv_propagator_gen_std          
           else                                      
            call init_vv_propagator_gen              
           endif                                     
        endif 


        end subroutine init_vv_propagator
!------------------------------------------------------------------------
! @brief Initialization for velocity Verlet propagation (vv)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_vv_propagator_drl

#ifndef MPI
       tp_myrank=0
#endif

       f1=quantum_dt*(1.d0-quantum_dt*0.5d0*drudel_eps_gm)
       f2=quantum_dt*quantum_dt*0.5d0
       f3=1.d0-quantum_dt*drudel_eps_gm*(1.d0-quantum_dt*0.5*drudel_eps_gm)
       f4=0.5d0*quantum_dt
       f5=drudel_eps_gm*f2
       if (tp_myrank.eq.0) write(6,*) "Initiated VV propagator"

       return

      end subroutine init_vv_propagator_drl

!------------------------------------------------------------------------
! @brief Initialization for velocity Verlet propagation (vv)
!
! @date Created: G. Gil
! Modified:
! Notes: Taken from init_vv_propagator and replacing drudel_eps_gm by BEM_Qg
!------------------------------------------------------------------------
      subroutine init_vv_propagator_gen

       real(dbl), allocatable :: BEM_I(:,:)
       integer :: j

       allocate(BEM_I(pedra_surf_n_tessere,pedra_surf_n_tessere))
       BEM_I = zero
       forall(j = 1:pedra_surf_n_tessere) BEM_I(j,j) = one

       ! FIXME: need to be deallocated somewhere
       allocate(BEM_f1(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                BEM_f3(pedra_surf_n_tessere,pedra_surf_n_tessere),&
                BEM_f5(pedra_surf_n_tessere,pedra_surf_n_tessere))

       BEM_f1=quantum_dt*(BEM_I-quantum_dt*0.5d0*BEM_Qg)
       f4=0.5d0*quantum_dt
       f2=quantum_dt*f4
       BEM_f3=BEM_I-matmul(BEM_Qg,BEM_f1)
       BEM_f5=BEM_Qg*f2
       write(6,*) "Initiated VV propagator"


       deallocate(BEM_I)

       return

      end subroutine

      subroutine init_vv_propagator_gen_std
!------------------------------------------------------------------------
! @brief Initialization for velocity Verlet propagation (vv)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------


#ifndef MPI
       tp_myrank=0
#endif

       allocate(std_f1(npoles),std_f3(npoles),std_f5(npoles))

       std_f1(:)=quantum_dt*(1.d0-quantum_dt*0.5d0*gg(:))
       f2=quantum_dt*quantum_dt*0.5d0
       std_f3(:)=1.d0-gg(:)*std_f1(:)
       f4=0.5d0*quantum_dt
       std_f5(:)=gg(:)*f2
       if (tp_myrank.eq.0) write(6,*) "Initiated VV propagator"

       return

      end subroutine


!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of charges with Onsager-BEM equations
! ! SP 07/07/17 vv propagaton
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_ons_drl

      ! charge propagation with drude/lorentz and cosmo/onsager equations

!EC:  mat_mult optimizes nts**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of nts

      ! Reaction Field
       qr_t=qr_tp+f1*dqr_tp+f2*fqr_tp
       ! SP: changed using $BEM_Qf=-ONS_ff(1)*S{-1}$
       fqr_t=-ONS_fw*qr_t+mat_mult(BEM_Qf,pot_tp)
       dqr_t=f3*dqr_tp+f4*(fqr_t+fqr_tp)-f5*fqr_tp
       fqr_tp=fqr_t
       dqr_tp=dqr_t
       ! Local Field
       if(global_medium_Floc.eq."loc") then
         qx_t=qx_tp+f1*dqx_tp+f2*fqx_tp
         fqx_t=-ONS_fw*qx_t+mat_mult(BEM_Qf,potf_tp)
         dqx_t=f3*dqx_tp+f4*(fqx_t+fqx_tp)-f5*fqx_tp
         fqx_tp=fqx_t
         dqx_tp=dqx_t
       endif

       return

      end subroutine prop_ons_drl

!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of charges with Onsager-BEM equations
! ! SP 07/07/17 updated with vv propagaton
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine prop_ons_deb

      ! charge propagation with drude/lorentz and onsager equations

!EC:  mat_mult optimizes nts**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of nts

      ! SP: Reaction Field eq.46 Corni et al. JPCA 2015
      qr_t=qr_tp-quantum_dt*ONS_taum1*qr_tp+quantum_dt*ONS_taum1*mat_mult(BEM_Q0,pot_tp2)&
                                 +mat_mult(BEM_Qd,pot_tp-pot_tp2)
      ! Local Field analogous equation, BEM_Q0x=ONS_fx0*Sm1 and BEM_Qdx=ONS_fxd*Sm1
       if(global_medium_Floc.eq."loc") then
         call do_field_from_charges(qx_tp,fx_tp)
         qx_t=qx_tp-quantum_dt*ONS_taum1*qx_tp+quantum_dt*ONS_taum1*mat_mult(BEM_Q0x,potf_tp2) &
                                 +mat_mult(BEM_Qdx,potf_tp-potf_tp2)
       endif

       return

      end subroutine prop_ons_deb
!
!------------------------------------------------------------------------
!> Drude-Lorentz propagation global_prop_Fprop=chr-ief standard algorithm
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation global_prop_Fprop=chr-ief standard algorithm
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine prop_ief_drl

      ! Charge propagation with drude/lorentz and IEF equations

      ! Reaction Field
       qr_t=qr_tp+quantum_dt*dqr_tp
       call DGEMV('N',pedra_surf_n_tessere,pedra_surf_n_tessere,quantum_dt,BEM_Qf,pedra_surf_n_tessere,pot_tp,one_i, &
                     zero,dqr_t,one_i)
       call DGEMV('N',pedra_surf_n_tessere,pedra_surf_n_tessere,-quantum_dt,BEM_Qw,pedra_surf_n_tessere,qr_tp,one_i,  &
                     one,dqr_t,one_i)
! SC 17/8/2016: changed the following,
!      dqr_t=dqr_t+(1-quantum_dt*drudel_eps_gm)*dqr_tp
       dqr_t=(dqr_t+dqr_tp)/(1.d0+drudel_eps_gm*quantum_dt)
! SC avoid developing a total charge
       dqr_t=dqr_t-sum(dqr_t)/pedra_surf_n_tessere
       dqr_tp=dqr_t
      ! Local Field
       if(global_medium_Floc.eq."loc") then
         qx_t=qx_tp+quantum_dt*dqx_tp
         if(global_medium_Fmdm.eq.'csol') then
          call DGEMV('N',pedra_surf_n_tessere,pedra_surf_n_tessere,quantum_dt,BEM_Qfx,pedra_surf_n_tessere,potf_tp,one_i,&
                        zero,dqx_t,one_i)
         else if((global_medium_Fmdm.eq.'cnan').or.&
                 (global_medium_Fmdm.eq.'qnan')) then
          call DGEMV('N',pedra_surf_n_tessere,pedra_surf_n_tessere,quantum_dt,BEM_Qf,pedra_surf_n_tessere,potf_tp,one_i,&
                        zero,dqx_t,one_i)
         endif
         call DGEMV('N',pedra_surf_n_tessere,pedra_surf_n_tessere,-quantum_dt,BEM_Qw,pedra_surf_n_tessere,qx_tp,one_i,&
                       one,dqx_t,one_i)
! SC 17/8/2016: changed the following,
!        dqx_t=dqx_t+(1.d0-quantum_dt*drudel_eps_gm)*dqx_tp
         dqx_t=(dqx_t+dqx_tp)/(1.d0+drudel_eps_gm*quantum_dt)
! SC avoid developing a total charge
         dqx_t=dqx_t-sum(dqx_t)/pedra_surf_n_tessere
         dqx_tp=dqx_t
       endif

       return

      end subroutine prop_ief_drl

!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation global_prop_Fprop=chr-ief velocity-verlet
! algorithm
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine prop_vv_ief_drl

      ! Charge propagation with drude/lorentz and IEF equations

!      qr_t=qr_tp+quantum_dt*(1.d0-quantum_dt*0.5d0*drudel_eps_gm)*dqr_tp+quantum_dt*quantum_dt*0.5d0*fqr_tp
!      fqr_t=-matmul(BEM_Qw,qr_t)+matmul(BEM_Qf,pot_tp)
!      dqr_t=(1.d0-quantum_dt*0.5d0*drudel_eps_gm)*dqr_tp+0.5d0*quantum_dt*(fqr_t+fqr_tp)
!      dqr_t=dqr_t/(1.d0+quantum_dt*0.5d0*drudel_eps_gm)
! SC integrator from E. Vanden-Eijnden, G. Ciccotti CPL 429 (2006) 310–316

!EC:  mat_mult optimizes nts**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of nts

       !qr_t=qr_tp+f1*dqr_tp+f2*fqr_tp 
       fqr_t=-mat_mult(BEM_Qw,qr_t)+mat_mult(BEM_Qf,pot_tp)
       dqr_t=f3*dqr_tp+f4*(fqr_t+fqr_tp)-f5*fqr_tp
       qr_t=qr_tp+f1*dqr_t+f2*fqr_t
       fqr_tp=fqr_t
       dqr_tp=dqr_t


      ! Local Field
       if(global_medium_Floc.eq."loc") then
         !qx_t=qx_tp+f1*dqx_tp+f2*fqx_tp
        if(global_medium_Fmdm.eq.'csol') then
         !fqx_t=-matmul(BEM_Qw,qx_t)+matmul(BEM_Qfx,potf_tp)
         fqx_t=-mat_mult(BEM_Qw,qx_t)+mat_mult(BEM_Qfx,potf_tp)
        else if((global_medium_Fmdm.eq.'cnan').or.&
                (global_medium_Fmdm.eq.'qnan')) then
         !fqx_t=-matmul(BEM_Qw,qx_t)+matmul(BEM_Qf,potf_tp)
         fqx_t=-mat_mult(BEM_Qw,qx_t)+mat_mult(BEM_Qf,potf_tp)
        endif
        dqx_t=f3*dqx_tp+f4*(fqx_t+fqx_tp)-f5*fqx_tp
        qx_t=qx_tp+f1*dqx_t+f2*fqx_t
        fqx_tp=fqx_t
        dqx_tp=dqx_t
       endif

       return

      end subroutine prop_vv_ief_drl

!------------------------------------------------------------------------
! @brief General dielectric function propagation global_prop_Fprop=chr-ief
! velocity-verlet
! algorithm
!
! @date Created: G. Gil
! Modified:
! Notes: Taken from prop_vv_ief_drl and changing to the vv with matrix
! factors
!        N.B. that matrix BEM_Qg is inside the BEM_f1, BEM_f3 and BEM_f5
!        matrices
!------------------------------------------------------------------------
      subroutine prop_vv_ief_gen

       integer(i4b) :: i
       real(dbl) :: threshold

       threshold = 1d-9


      ! Charge propagation with general dielectric function and IEF
      ! equations
       qr_t=qr_tp+matmul(BEM_f1,dqr_tp)+f2*fqr_tp
       qr_t=qr_t+quantum_dt*0.5d0*matmul(BEM_Qdf,pot_tp-pot_tp2)
       fqr_t=-matmul(BEM_Qw,qr_t)+matmul(BEM_Qf,pot_tp)
       fqr_t=fqr_t+matmul((BEM_Qdf-quantum_dt*0.5d0*BEM_Qdf_2g),pot_tp-pot_tp2)
       dqr_t=matmul(BEM_f3,dqr_tp)+f4*(fqr_t+fqr_tp)-matmul(BEM_f5,fqr_tp)
       fqr_tp=fqr_t
       dqr_tp=dqr_t

      ! Local Field
       if(global_medium_Floc.eq."loc") then
        qx_t=qx_tp+matmul(BEM_f1,dqx_tp)+f2*fqx_tp
        if(global_medium_Fmdm.eq.'csol') then
         qx_t=qx_t+quantum_dt*0.5d0*matmul(BEM_Qdfx,potf_tp-potf_tp2)
         fqx_t=-matmul(BEM_Qw,qx_t)+matmul(BEM_Qfx,potf_tp)
         fqx_t=fqx_t+matmul((BEM_Qdfx-quantum_dt*0.5d0*BEM_Qdfx_2g),potf_tp-potf_tp2)
        else if((global_medium_Fmdm.eq.'cnan').or.&
                (global_medium_Fmdm.eq.'qnan')) then
         qx_t=qx_t+quantum_dt*0.5d0*matmul(BEM_Qdf,potf_tp-potf_tp2)
         fqx_t=-matmul(BEM_Qw,qx_t)+matmul(BEM_Qf,potf_tp)
         fqx_t=fqx_t+matmul((BEM_Qdf-quantum_dt*0.5d0*BEM_Qdf_2g),potf_tp-potf_tp2)
        endif
        dqx_t=matmul(BEM_f3,dqx_tp)+f4*(fqx_t+fqx_tp)-matmul(BEM_f5,fqx_tp)
        fqx_tp=fqx_t
        dqx_tp=dqx_t
       endif

       return

      end subroutine

      subroutine prop_vv_ief_gen_std
!------------------------------------------------------------------------
! @brief General dielectric function propagation global_prop_Fprop=chr-ief
! velocity-verlet
! algorithm
!
! @date Created: G. Gil
! Modified:
!------------------------------------------------------------------------

      ! Charge propagation with general dielectric function and IEF
      ! equations
       integer(i4b) :: pidx,ncycle
       integer(i4b) :: i,j,its
       real(dbl) :: threshold, q_tot

       threshold = .5d-8
       if(global_medium_Fnorm.eq."non") global_medium_Fnorm = "tot"
       if (typ_prop.eq."0") then
          ncycle=npoles
       else
            if(npoles.eq.1) then
               ncycle=npoles
            else
               ncycle=npoles-1
            endif
       endif

       qr_tp(:)=qr_t(:)
       qr_t(:) = zero
       if(global_medium_Floc.eq."loc") qx_t(:) = zero
       do pidx = 1, ncycle
        fqr_t_p(:,pidx)=-w2(pidx)*qr_tp_p(:,pidx)+kf(pidx)*(matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp))
        fqr_t_p(:,pidx)=fqr_t_p(:,pidx)-sum(fqr_t_p(:,pidx))/pedra_surf_n_tessere
        dqr_t_p(:,pidx)=std_f3(pidx)*dqr_tp_p(:,pidx)+f4*(fqr_t_p(:,pidx)+fqr_tp_p(:,pidx))-std_f5(pidx)*fqr_tp_p(:,pidx)
        qr_t_p(:,pidx)=qr_tp_p(:,pidx)+std_f1(pidx)*dqr_t_p(:,pidx)+f2*fqr_t_p(:,pidx)
        fqr_tp_p(:,pidx)=fqr_t_p(:,pidx)
        dqr_tp_p(:,pidx)=dqr_t_p(:,pidx)
        qr_tp_p(:,pidx)=qr_t_p(:,pidx)
        qr_t(:) = qr_t(:) + qr_t_p(:,pidx)
      ! Local Field
       if(global_medium_Floc.eq."loc") then
        if(global_medium_Fmdm.eq.'csol') then
         fqx_t_p(:,pidx)=-w2(pidx)*qx_tp_p(:,pidx)+kf(pidx)*(matmul(BEM_Qfx,potf_tp)+matmul(BEM_ADt,qx_tp))
        else if((global_medium_Fmdm.eq.'cnan').or.&
                (global_medium_Fmdm.eq.'qnan')) then
         fqx_t_p(:,pidx)=-w2(pidx)*qx_tp_p(:,pidx)+kf(pidx)*matmul(BEM_Qf,potf_tp)+kf(pidx)*matmul(BEM_ADt,qx_tp)
        endif
        dqx_t_p(:,pidx)=std_f3(pidx)*dqx_tp_p(:,pidx)+f4*(fqx_t_p(:,pidx)+fqx_tp_p(:,pidx))-std_f5(pidx)*fqx_tp_p(:,pidx)
        qx_t_p(:,pidx)=qx_tp_p(:,pidx)+std_f1(pidx)*dqx_t_p(:,pidx)+f2*fqx_t_p(:,pidx)
        if(global_medium_Fnorm.eq."tot") then 
                qx_t_p(:,pidx)=qx_t_p(:,pidx)-sum(qx_t_p(:,pidx))/pedra_surf_n_tessere
        elseif(global_medium_Fnorm.eq."sep") then
             do i=1,pedra_surf_n_particles
                 q_tot=0
                 do its=pedra_surf_comp(i,2),pedra_surf_comp(i,3)
                       q_tot=q_tot+qx_t_p(its,pidx)
                 enddo
                 do its=pedra_surf_comp(i,2),pedra_surf_comp(i,3)
                       qx_t_p(its,pidx)=qx_t_p(its,pidx)-q_tot/pedra_surf_comp(i,1)
                 enddo
             enddo 
        endif
        fqx_tp_p(:,pidx)=fqx_t_p(:,pidx)
        dqx_tp_p(:,pidx)=dqx_t_p(:,pidx)
        qx_tp_p(:,pidx)=qx_t_p(:,pidx)
        qx_t(:) = qx_t(:) + qx_t_p(:,pidx)
       endif
       enddo
       if (typ_prop.ne."0".and.npoles.ne.1) then
           pidx=npoles
           if (typ_prop.eq."1") then
                !derivative method
                fqr_t_p(:,pidx)=matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp)
                qr_t_p(:,pidx)=qr_tp_p(:,pidx)+kf0(pidx)*(fqr_t_p(:,pidx)-fqr_tp_p(:,pidx))
                fqr_tp_p(:,pidx)=fqr_t_p(:,pidx)
           elseif (typ_prop.eq."2") then
                !second order method
                fqr_t_p(:,pidx)=matmul(BEM_Qf,pot_tp)+matmul(BEM_ADt,qr_tp)
                dfqr_t_p(:)=3*fqr_t_p(:,pidx)-8*fqr_tp_p(:,pidx)+7*fqr_tp2_p(:)-2*fqr_tp3_p(:)
                qr_t_p(:,pidx)=qr_tp_p(:,pidx)+quantum_dt*dqr_tp_p(:,pidx)+kf(pidx)/(w2(pidx))*dfqr_t_p(:)
                dqr_t_p(:,pidx)=dqr_tp_p(:,pidx)+kf(pidx)/(quantum_dt*w2(pidx))*(dfqr_t_p(:))
                fqr_tp3_p(:)=fqr_tp2_p(:)
                fqr_tp2_p(:)=fqr_tp_p(:,pidx)
                dqr_tp_p(:,pidx)=dqr_t_p(:,pidx)
                dfqr_tp_p(:)=dfqr_t_p(:)
                fqr_tp_p(:,pidx)=fqr_t_p(:,pidx)
           elseif (typ_prop.eq."3") then
                   fqr_t_p(:,pidx)=matmul(BEM_ADt,qr_t)+matmul(BEM_Qf,pot_tp)
                   qr_t_p(:,pidx)=kf0(pidx)*matmul(BEM_ADtm1,fqr_tp_p(:,pidx))
           endif
           qr_tp_p(:,pidx)=qr_t_p(:,pidx)
           qr_t(:) = qr_t(:) + qr_t_p(:,pidx)
           if(global_medium_Floc.eq."loc") then
               if (typ_prop.eq."1") then
                   !derivative method
                    if(global_medium_Fmdm.eq.'csol') then
                         fqx_t_p(:,pidx)=matmul(BEM_Qfx,potf_tp)+matmul(BEM_ADt,qx_tp)
                    else if(global_medium_Fmdm.eq.'cnan'.or.global_medium_Fmdm.eq.'qnan') then
                         fqx_t_p(:,pidx)=matmul(BEM_Qf,potf_tp)+matmul(BEM_ADt,qx_tp)
                    endif
                    qx_t_p(:,pidx)=qx_tp_p(:,pidx)+kf0(pidx)*(fqx_t_p(:,pidx)-fqx_tp_p(:,pidx))
               elseif (typ_prop.eq."2") then
                   !second order method
                    if(global_medium_Fmdm.eq.'csol') then
                         fqx_t_p(:,pidx)=matmul(BEM_Qfx,potf_tp)+matmul(BEM_ADt,qx_tp)
                    else if(global_medium_Fmdm.eq.'cnan'.or.global_medium_Fmdm.eq.'qnan') then
                         fqx_t_p(:,pidx)=matmul(BEM_Qf,potf_tp)+matmul(BEM_ADt,qx_tp)
                    endif
                   dfqx_t_p(:)=3*fqx_t_p(:,pidx)-8*fqx_tp_p(:,pidx)+7*fqx_tp2_p(:)-2*fqx_tp3_p(:)
                   qx_t_p(:,pidx)=qx_tp_p(:,pidx)+quantum_dt*dqx_tp_p(:,pidx)+kf(pidx)/(2*w2(pidx))*dfqx_t_p(:)
                   dqx_t_p(:,pidx)=dqx_tp_p(:,pidx)+kf(pidx)/(2*quantum_dt*w2(pidx))*(dfqx_t_p(:)+dfqx_tp_p(:))
                   fqx_tp3_p(:)=fqx_tp2_p(:)
                   fqx_tp2_p(:)=fqx_tp_p(:,pidx)
                   dqx_tp_p(:,pidx)=dqx_t_p(:,pidx)
                   dfqx_tp_p(:)=dfqx_t_p(:)
               elseif (typ_prop.eq."3") then
                   if(global_medium_Fmdm.eq.'csol') then
                         fqx_t_p(:,pidx)= matmul(BEM_Qfx,potf_tp)+matmul(BEM_ADt,qx_t)
                         qx_t_p(:,pidx)=kf0(pidx)*matmul(BEM_ADtm1,fqx_tp_p(:,pidx))
                    else if(global_medium_Fmdm.eq.'cnan'.or.global_medium_Fmdm.eq.'qnan') then
                         fqx_t_p(:,pidx)=matmul(BEM_ADt,qx_t)+matmul(BEM_Qf,potf_tp)
                         qx_t_p(:,pidx)=kf0(pidx)*matmul(BEM_ADtm1,fqx_tp_p(:,pidx))
                    endif
               endif
               qx_t_p(:,pidx)=qx_t_p(:,pidx)-sum(qx_t_p(:,pidx))/pedra_surf_n_tessere
               fqx_tp_p(:,pidx)=fqx_t_p(:,pidx)
               qx_tp_p(:,pidx)=qx_t_p(:,pidx)
               qx_t(:) = qx_t(:) + qx_t_p(:,pidx)
           endif
       endif
       return

      end subroutine

!------------------------------------------------------------------------
! @brief Debye propagation global_prop_Fprop=chr-ief
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine prop_ief_deb

      ! Charge propagation with debye and IEF equations

!EC:  mat_mult optimizes nts**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of nts

      ! Reaction Field
       qr_t=qr_tp-quantum_dt*mat_mult(BEM_R,qr_tp)+quantum_dt*mat_mult(BEM_Qt,pot_tp) &
                                    +mat_mult(BEM_Qd,pot_tp-pot_tp2)
      ! Local Field eq.47 JPCA 2015
       if(global_medium_Floc.eq."loc") then
        if(global_medium_Fmdm.eq.'csol') then
         ! GG: BEM matrices (except R) are different in the case of local-field for solvent external medium
         qx_t=qx_tp-quantum_dt*mat_mult(BEM_R,qx_tp)+quantum_dt*mat_mult(BEM_Qtx,potf_tp) &
                                     +mat_mult(BEM_Qdx,potf_tp-potf_tp2)
        else if((global_medium_Fmdm.eq.'cnan').or.&
                (global_medium_Fmdm.eq.'qnan')) then
         qx_t=qx_tp-quantum_dt*mat_mult(BEM_R,qx_tp)+quantum_dt*mat_mult(BEM_Qt,potf_tp) &
                                     +mat_mult(BEM_Qd,potf_tp-potf_tp2)
        endif
       endif

       return

      end subroutine prop_ief_deb

!------------------------------------------------------------------------
! @brief  Debye propagation global_prop_Fprop=chr-ied
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------
      subroutine prop_ied_deb

      ! Charge propagation with debye and IEF equations one taud

!EC:  mat_mult optimizes nts**2-based statements
!     mat_mult uses matmul or explicit loops (with OMP),
!     according to the value of nts

      ! Reaction Field
       qr_t=qr_tp-quantum_dt*taum1*qr_tp+quantum_dt*taum1*mat_mult(BEM_Q0,pot_tp) &
                                  +mat_mult(BEM_Qd,pot_tp-pot_tp2)
      ! Local Field eq.47 JPCA 2015
       if(global_medium_Floc.eq."loc") then
        if(global_medium_Fmdm.eq.'csol') then
         ! GG: BEM matrices (except taum1) are different in the case of local-field for solvent external medium
         qx_t=qx_tp-quantum_dt*taum1*qx_tp+quantum_dt*taum1*mat_mult(BEM_Q0x,potf_tp) &
                                   +mat_mult(BEM_Qdx,potf_tp-potf_tp2)
        else if((global_medium_Fmdm.eq.'cnan').or.&
                (global_medium_Fmdm.eq.'qnan')) then
         qx_t=qx_tp-quantum_dt*taum1*qx_tp+quantum_dt*taum1*mat_mult(BEM_Q0,potf_tp) &
                                   +mat_mult(BEM_Qd,potf_tp-potf_tp2)
        endif
       endif

       return

      end subroutine prop_ied_deb

!------------------------------------------------------------------------
! @brief  Debye propagation of the spheroid reaction and local fields
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_dip_spho_deb(f_tp)

       real(dbl), intent(IN):: f_tp(3)
       real(dbl) :: f(3),fp(3),mp(3),mp2(3)
       integer(i4b) :: i

       ! Equation (23) JCP 146, 064116 (2017)
       ! SP 05/07/17: this is easily genealizable to multipoles see JCP 146, 064116 (2017)
       fr_t=zero
       do i=1,3
         ! Determine previous field component along spheroid axis:
         fp=dot_product(fr_tp,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
         mp=dot_product(mu_tp,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
         mp2=dot_product(mu_tp2,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
         f=fp+matmul(MPL_Fd(:,:,i,1),mp-mp2)+quantum_dt*global_prop_n_q*  &
             (matmul(MPL_Ft0(:,:,i),mp)-MPL_Taum1(i)*fp)
         if(global_sys_Fdeb.eq."equ") f=matmul(MPL_F0(:,:,i,1),mu_tp) !Equilibrium RF for debug purposes
         fr_t=fr_t+f
       enddo
       ! Local Field
       if(global_medium_Floc.eq."loc") then
         ! Equation (23) JCP 146, 064116 (2017)
         fx_t=zero
         do i=1,3
           fp=dot_product(fx_tp,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
           mp=dot_product(mu_tp,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
           mp2=dot_product(mu_tp2,sph_surf_vrs(:,i,1))*sph_surf_vrs(:,i,1)
           f=fp+matmul(MPL_Fxd(:,:,i),mp-mp2)+quantum_dt*global_prop_n_q*  &
               (matmul(MPL_Ftx0(:,:,i),mp)-MPL_Taum1(i)*fp)
           fx_t=fx_t+f
         enddo
       endif

       return

      end subroutine prop_dip_spho_deb

!------------------------------------------------------------------------
! @brief  Debye propagation of the sphere reaction and local fields
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_dip_sphe_deb(f_tp)

       real(dbl), intent(IN):: f_tp(3)
       real(dbl) :: f(3),fp(3),mt(3),mp(3)
       integer(i4b) :: i

       fr_t=fr_tp+ONS_fd*(mu_tp-mu_tp2)+quantum_dt*global_prop_n_q*ONS_taum1* &
           (ONS_f0*mu_tp-fr_tp)
       if(global_sys_Fdeb.eq."equ") fr_t=ONS_f0*mu_tp
       ! Local Field
       if(global_medium_Floc.eq."loc") then
        fx_t=fx_tp+ONS_fxd*(f_tp-f_tp2)+quantum_dt*global_prop_n_q*ONS_tauxm1* &
           (ONS_fx0*f_tp-fx_tp)
       endif

       return

      end subroutine prop_dip_sphe_deb

!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of the spheroid reaction and local fields
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_dip_spho_drl(f_tp)

       real(dbl), intent(IN):: f_tp(3)
       real(dbl):: f(3),fld(3),mr(3)
       integer(i4b):: i,j,k

       ! Equation ....  paper ....
       fr_t=zero
       do i=1,sph_surf_nsph
         ! SP 07/07/17: Determine previous field component along spheroid axis:
         mr_t(:,i)=mr_tp(:,i)+f1*dmr_tp(:,i)+f2*fmr_tp(:,i)
         ! Determine total field acting on each dipole from:
         ! other particles
         fld=zero
         do k=1,sph_surf_nsph
           if (k.ne.i) then
             call do_field_from_dip(mr_tp(:,k),sph_surf_centre(:,k),f, &
                                              sph_surf_centre(:,i))
           fld=fld+f
           endif
         enddo
         ! and the molecule
         call do_field_from_dip(mu_tp,quantum_mol_cc,f,sph_surf_centre(:,i))
         fld=fld+f
         fmr_t=zero
         do j=1,3
           mr(:)=dot_product(mr_t(:,i),sph_surf_vrs(:,j,i))*sph_surf_vrs(:,j,i)
           fmr_t(:,i)=fmr_t(:,i)-MPL_Fw(j,i)*mr(:) &
                                +matmul(MPL_Ff(:,:,j,i),fld)
         enddo
         ! propagate dipole using vv
         dmr_t(:,i)=f3*dmr_tp(:,i)+f4*(fmr_t(:,i)+fmr_tp(:,i))- &
                                                f5*fmr_tp(:,i)
         call do_field_from_dip(mr_t(:,i),sph_surf_centre(:,i),f,quantum_mol_cc)
         fr_t=fr_t+f
      enddo
       fmr_tp=fmr_t
       dmr_tp=dmr_t
       mr_tp=mr_t
       ! Local Field
       if(global_medium_Floc.eq."loc") then
         ! Equation ....  paper ....
         fx_t=zero
         do i=1,sph_surf_nsph
           ! SP 07/07/17: Determine previous field component along spheroid axis:
           mx_t(:,i)=mx_tp(:,i)+f1*dmx_tp(:,i)+f2*fmx_tp(:,i)
           fld=zero
           do k=1,sph_surf_nsph
             if (k.ne.i) then
               call do_field_from_dip(mx_tp(:,k),sph_surf_centre(:,k),f, &
                                                sph_surf_centre(:,i))
             fld=fld+f
             endif
           enddo
           ! and the external field
           fld=fld+f_tp
           fmx_t=zero
           do j=1,3
             fmx_t(:,i)=fmx_t(:,i)-MPL_Fw(j,i)*mx_t(:,i) &
                                  +matmul(MPL_Ff(:,:,j,i),fld)
           enddo
           ! propagate dipole using vv
           dmx_t(:,i)=f3*dmx_tp(:,i)+f4*(fmx_t(:,i)+fmx_tp(:,i))-&
                                                 f5*fmx_tp(:,i)
           ! Update contribution to reaction field on the molecule
           call do_field_from_dip(mx_t(:,i),sph_surf_centre(:,i),f,quantum_mol_cc)
           fx_t=fx_t+f
         enddo
         fmx_tp=fmx_t
         dmx_tp=dmx_t
         mx_tp=mx_t
       endif

       return

      end subroutine prop_dip_spho_drl

!------------------------------------------------------------------------
! @brief  Drude-Lorentz propagation of the sphere reaction and local fields
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine prop_dip_sphe_drl(f_tp)

       real(dbl), intent(IN):: f_tp(3)
       real(dbl):: f(3),fld(3)
       integer(i4b):: i,j,k

       ! Equation ....  paper ....
       fr_t=zero
       do i=1,sph_surf_nsph
         mr_t(:,i)=mr_tp(:,i)+f1*dmr_tp(:,i)+f2*fmr_tp(:,i)
         fld=zero
         do k=1,sph_surf_nsph
           if (k.ne.i) then
             call do_field_from_dip(mr_tp(:,k),sph_surf_centre(:,k),f, &
                                              sph_surf_centre(:,i))
           fld=fld+f
           endif
         enddo
         call do_field_from_dip(mu_tp,quantum_mol_cc,f,sph_surf_centre(:,i))
         fld=fld+f
         fmr_t(:,i)=-ONS_fw*mr_t(:,i)+ONS_ff(i)*fld
         dmr_t(:,i)=f3*dmr_tp(:,i)+f4*(fmr_t(:,i)+fmr_tp(:,i))- &
                                                f5*fmr_tp(:,i)
         call do_field_from_dip(mr_t(:,i),sph_surf_centre(:,i),f,quantum_mol_cc)
         fr_t=fr_t+f
       enddo
       fmr_tp=fmr_t
       dmr_tp=dmr_t
       mr_tp=mr_t
       if(global_medium_Floc.eq."loc") then
         fx_t=zero
         do i=1,sph_surf_nsph
           mx_t(:,i)=mx_tp(:,i)+f1*dmx_tp(:,i)+f2*fmx_tp(:,i)
           fld=zero
           do k=1,sph_surf_nsph
             if (k.ne.i) then
               call do_field_from_dip(mx_tp(:,k),sph_surf_centre(:,k),f, &
                                                sph_surf_centre(:,i))
             fld=fld+f
             endif
           enddo
           ! and the external field
           fld=fld+f_tp
           fmx_t=zero
           fmx_t(:,i)=-ONS_fw*mx_t(:,i)+ONS_ff(i)*fld
           dmx_t(:,i)=f3*dmx_tp(:,i)+f4*(fmx_t(:,i)+fmx_tp(:,i))-&
                                                 f5*fmx_tp(:,i)
           call do_field_from_dip(mx_t(:,i),sph_surf_centre(:,i),f,quantum_mol_cc)
           fx_t=fx_t+f
         enddo
         fmx_tp=fmx_t
         dmx_tp=dmx_t
         mx_tp=mx_t
       endif

       return

      end subroutine prop_dip_sphe_drl

!------------------------------------------------------------------------
! @brief Fetch the current value of the onsager field from elsewhere
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine get_ons(f_med)

       real(dbl):: f_med(3)

       f_med=fr_t

       return

      end subroutine get_ons
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Free-energy   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------
! @brief Update non-equilibrium free energy
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_gneq(v_avg,mu_or_v,df_or_dq,f_or_q,f_or_q0,fact_d, &
                           n_coor_or_ts,sig)

        implicit none

        integer(i4b),  intent(in) :: n_coor_or_ts,sig
        real(dbl),     intent(in) :: mu_or_v(n_coor_or_ts,quantum_n_ci,quantum_n_ci)
        real(dbl),     intent(in) :: df_or_dq(n_coor_or_ts), &
                                     f_or_q(n_coor_or_ts),&
                                     f_or_q0(n_coor_or_ts),&
                                     fact_d(n_coor_or_ts,n_coor_or_ts)
        real(dbl),     intent(in) :: v_avg(:) !< molecular potential or dipole depending on the case (changing name!)
!       real(dbl) :: de_a
        integer(i4b)              :: its,k,j

        g_eq=0.d0
!       de_a=0.d0

        g_neq1_part=g_neq1_part+sig*dot_product(v_avg,df_or_dq)
        g_eq=sig*dot_product(v_avg,f_or_q)

        g_eq=0.5*g_eq
! SC 27/09/2016: corrected bug in expression of e_vac
!       e_vac=2.*g_eq_gs-de_a
        e_vac=-2.*g_eq+2.*g_eq_gs
!        g_neq1=g_neq_0+g_neq1_part-g_eq
        g_neq1=g_neq_0-g_neq1_part
        g_neq2=-sig*(dot_product((mu_or_v(:,1,1)-v_avg), &
               (f_or_q-matmul(fact_d,v_avg)))- &
               0.5*dot_product(v_avg,matmul(fact_d,v_avg)))- &
               g_neq2_0+e_vac
!       write(6,*) 'g_neq2 a, g_neq2 b', &
!             -sig*dot_product((mu_or_v(:,1,1)-v_avg), &
!               (f_or_q-matmul(fact_d,v_avg))), &
!             -sig*0.5*dot_product(v_avg,matmul(fact_d,v_avg))
!       write (6,*) 'g_eq,g_neq1_part',g_eq,g_neq1_part
        g_eq=-g_eq_gs+g_eq+e_vac
! SC to be completed with other means to calculate gneq
! SC: Caricato et al. JCP 2006, currently only for Onsager

        return

      end subroutine do_gneq
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Test/Debug    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Test for Solvent/Nanoparticle reaction/local fields
!
! @date Created: S. Pipolo
! Modified: G. Gil, M. Rosa
!------------------------------------------------------------------------
      subroutine do_ref(mu)

       real(dbl), optional, intent(IN) :: mu(3)
       complex(cmp) :: refc, E0
       complex(cmp) :: eps,eps_f

       integer(i4b) :: its
       real(dbl):: dist,pos(3),dp,f(3),rr,facd,fac0,tau

       select case (global_sys_Ftest)
         ! Spherical Nanoparticle reaction field
         case ("n-r")
           if(global_prop_Fprop.eq."dip") then
             pos(1)=sph_surf_centre(1,1)
             pos(2)=sph_surf_centre(2,1)
             pos(3)=sph_surf_centre(3,1)
             rr=sph_surf_maj(1)
           else
             pos=zero
             do its=1,pedra_surf_n_tessere
               pos(1)=pos(1)+pedra_surf_tessere(its)%x
               pos(2)=pos(2)+pedra_surf_tessere(its)%y
               pos(3)=pos(3)+pedra_surf_tessere(its)%z
             enddo
             pos=pos/pedra_surf_n_tessere
             rr=pedra_surf_tessere(1)%rsfe
           endif
           call do_field_from_dip(mu,quantum_mol_cc,f,pos)
           pos=pos-quantum_mol_cc
           dist=sqrt(dot_product(pos,pos))
           dp=dot_product(mu,pos)
           f(:)=(3*dp*pos(:)/dist**2-mu(:))/dist**3
           !write(6,*) "pos(1), mu(1), dp, f(1)", pos(1), mu(1), dp, f(1)
           ref=f(1)*rr*rr*rr
         ! Spherical Nanoparticle local field
         case ("n-l")
           if(global_prop_Fprop.eq."dip") then
             rr=sph_surf_maj(1)
           else
             rr=pedra_surf_tessere(1)%rsfe
           endif
           !15/04/2020 MR: do_eps_drl don't exist anymore, has been moved to do_epsilon module and takes omega as argument
           !things that were done inside do_eps_drl are now done here explicitly
           !17/05/20: SP: why can't we use do_epsilon here?
           !call do_eps_drl
           eps=dcmplx(drudel_eps_A,zero)/dcmplx(drudel_eps_w0**2-quantum_omega(1)**2,-quantum_omega(1)*drudel_eps_gm)
           eps=eps+onec
           eps_f=(eps-onec)/(eps+twoc)
           E0=dcmplx(zero,0.5d0*sqrt(dot_product(quantum_fmax(:,1),quantum_fmax(:,1))))
           refc=eps_f*E0*exp(-ui*quantum_omega(1)*t)
           ref=real(refc+conjg(refc))*rr**3
           !stop
         case ("s-r")
         !Corni JCPA 2015, fig.1
           if(global_prop_Fprop.eq."dip") then
             rr=sph_surf_maj(1)
           else
             rr=sqrt(pedra_surf_tessere(1)%x**2+pedra_surf_tessere(1)%y**2+pedra_surf_tessere(1)%z**2)
           endif
           facd=(two*debye_eps_d-two)/(two*debye_eps_d+one)
           fac0=two*three*(debye_eps_0-debye_eps_d)/(two*debye_eps_d-two)/(two*debye_eps_0+one)
           tau=(two*debye_eps_d+one)/(two*debye_eps_0+one)*debye_eps_tau
           ref=facd*(one-fac0*(exp(-t/tau)-one))/rr**3*mu_tp(3)
         case ("s-l")
           if(global_prop_Fprop.eq."dip") then
             rr=sph_surf_maj(1)
           else
             rr=sqrt(pedra_surf_tessere(1)%x**2+pedra_surf_tessere(1)%y**2+pedra_surf_tessere(1)%z**2)
           endif
           !15/04/2020 MR: do_eps_deb don't exist anymore, has been moved to do_epsilon module and takes omega as argument
           !things that were done inside do_eps_deb are now done here explicitly
           !call do_eps_deb
           eps=dcmplx(debye_eps_d,zero)+dcmplx(debye_eps_0-debye_eps_d,zero)/ &
               dcmplx(one,-quantum_omega(1)*debye_eps_tau)
           eps_f=(three*eps)/(two*eps+onec)
           if((global_medium_Fmdm.eq.'csol').and.&
             (global_prop_Fprop.ne.'dip').and.&
             (global_prop_Fprop.ne.'chr-ons')) then
                 eps_f=(eps-onec)/(two*eps+onec)
           endif
           E0=dcmplx(zero,0.5d0*sqrt(dot_product(quantum_fmax(:,1),quantum_fmax(:,1))))
           refc=eps_f*E0*exp(-ui*quantum_omega(1)*t)
           ref=real(refc+conjg(refc))
       end select

       return

      end subroutine do_ref
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Input/Output  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
! @brief Output medium reaction/local field/dipoles
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine out_mdm(i)

       integer(i4b),intent(IN) :: i
       real(dbl):: fm(3)

       select case(global_sys_Ftest)
       case ('n-r','n-l')
          if(global_eps_Feps.eq."gen") then 
             write (file_med,'(i8,f12.2,4e22.10e3)') i,t,mu_mdm(:,1)
          else
             write (file_med,'(i8,f12.2,4e22.10e3)') i,t,mu_mdm(:,1),ref
          endif
       case ('s-r')
       if((global_prop_Fprop.eq."chr-ief").or.&
          (global_prop_Fprop.eq."chr-ied").or.&
          (global_prop_Fprop.eq."chr-ons")) then
                call do_field_from_charges(qr_t,fr_t)
        endif
         write (file_med,'(i8,f12.2,4e22.10e3)') i,t,fr_t(:),ref
       case ('s-l')
         write (file_med,'(i8,f12.2,4e22.10e3)') i,t,fx_t(:),ref
       case default
         if(global_prop_Fprop.eq."dip") then
           write (file_med,'(i8,f12.2,3e22.10e3)') i,t,fr_t(:)
         else
           call do_field_from_charges(qr_t,fr_t)
           write (file_med,'(i8,f12.2,9e22.10e3)')i,t,mu_mdm(:,1),&
                                                    fr_t(:),qtot,qtot0
         endif
       end select

       return

      end subroutine out_mdm

!------------------------------------------------------------------------
! @brief Output medium reaction/local field/dipoles in binary format
!
! @date Created: E. Coccia 18/6/18
! Modified:
!------------------------------------------------------------------------
      subroutine out_mdm_bin(i)

       integer(i4b),intent(IN) :: i
       real(dbl):: fm(3)

       select case(global_sys_Ftest)
       case ('n-r','n-l')
         write (file_med) i,t,mu_mdm(:,1),ref
       case ('s-r')
       if((global_prop_Fprop.eq."chr-ief").or.&
          (global_prop_Fprop.eq."chr-ied").or.&
          (global_prop_Fprop.eq."chr-ons")) then
              call do_field_from_charges(qr_t,fr_t)
       endif
       write (file_med) i,t,fr_t(:),ref
       case ('s-l')
         write (file_med) i,t,fx_t(:),ref
       case default
         if(global_prop_Fprop.eq."dip") then
           write (file_med) i,t,fr_t(:)
         else
           call do_field_from_charges(qr_t,fr_t)
           write (file_med) i,t,mu_mdm(:,1),fr_t(:),qtot,qtot0
         endif
       end select

       return

      end subroutine out_mdm_bin

!------------------------------------------------------------------------
! @brief Read charges from gaussian "charges0.inp" output
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine read_charges_gau

        integer(4) :: i,nts

#ifndef MPI
        tp_myrank=0
#endif

        open(7,file="charges0.inp",status="old",err=10)
        if (tp_myrank.eq.0) then
           write (6,*) "Initial charges read from charges0.inp"
        endif
          read(7,*) nts
          if(pedra_surf_n_tessere.eq.0.or.nts.eq.pedra_surf_n_tessere) then
            pedra_surf_n_tessere=nts
          else
           if (tp_myrank.eq.0) write(*,*) "Tesserae number conflict"
#ifdef MPI
           call mpi_finalize(tp_ierr_mpi)
#endif
            stop
          endif
          qtot0=zero
          do i=1,pedra_surf_n_tessere
            read(7,*) q0(i)
            qtot0=qtot0+q0(i)
          enddo
        close(7)

        return

10      write (6,*) "No file charges0.inp found,", &
                    "charges calculated for the initial state"
        return

      end subroutine read_charges_gau


! @brief Return propagated charges
!
! @date Created: M. Rosa
! Modified:
!------------------------------------------------------------------------
      subroutine get_propagated_charges(q_out)

           real(dbl), intent(out)  ::  q_out(pedra_surf_n_tessere)
           q_mdm=qr_t
           if(global_medium_Floc.eq."loc") q_mdm(:)=qr_t(:)+qx_t(:)
           q_out=q_mdm

       end subroutine


      subroutine get_corrected_propagated_charges(q_out)

           real(dbl), intent(out)  ::  q_out(pedra_surf_n_tessere)
           q_mdm=qr_t
           if(global_medium_Floc.eq."loc") q_mdm(:)=qr_t(:)+qx_t(:)
! SC 31/10/2016: avoid including interaction with an unwanted net charge
           q_out=q_mdm+(qtot0-sum(q_mdm))/pedra_surf_n_tessere

       end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Miscellaneous !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------
!> create a new BEM_Q0=BEM_Qw^-1*BEM_Qf that should avoid
!!               spurious charge dynamics for stationary states
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! @brief Create a new BEM_Q0=BEM_Qw^-1*BEM_Qf that should avoid
!               spurious charge dynamics for stationary states
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine init_BEM_Q0

       implicit none

       real(dbl),allocatable :: mat_t(:,:),mat_t2(:,:)
       integer(i4b) :: info,its
       integer(i4b),allocatable :: ipiv(:)

       allocate(mat_t(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(mat_t2(pedra_surf_n_tessere,pedra_surf_n_tessere))
       allocate(ipiv(pedra_surf_n_tessere))
       mat_t=BEM_Qw
       mat_t2=BEM_Qf
       call dgetrf(pedra_surf_n_tessere,pedra_surf_n_tessere,mat_t,pedra_surf_n_tessere,ipiv,info)
       call dgetrs('N',pedra_surf_n_tessere,pedra_surf_n_tessere,mat_t,pedra_surf_n_tessere,ipiv,mat_t2,pedra_surf_n_tessere,info)
       BEM_Q0=mat_t2
       deallocate (mat_t)
       deallocate (mat_t2)
       deallocate (ipiv)

       return

      end subroutine init_BEM_Q0

!------------------------------------------------------------------------
! @brief write restart
!
! @date Created   : E. Coccia 28 Nov 2017
! Modified  : G. Gil 02 Jul 2018
!------------------------------------------------------------------------
      subroutine wrt_restart_mdm(istep)

       implicit none

       integer(i4b), intent(IN)     :: istep
       integer(i4b)     :: i
       open(778, file='restart_mdm')

!        if (global_prop_Fint.eq.'ons') then
!          write(778,*) 'Dipoles for Onsager'
!          do i=1,3
!             write(778,*) fr_t(1), fr_t(2), fr_t(3)
!          enddo
!          if (global_medium_Floc.eq.'loc') then
!             write(778,*) 'Dipoles for Onsager (local)'
!             do i=1,3
!                write(778,*) fx_t(1), fx_t(2), fx_t(3)
!             enddo
!          endif
!       elseif (global_prop_Fint.eq.'pcm') then
       if (global_prop_Fint.eq.'pcm') then
          write(778,*) 'Reaction-field polarization charges (PCM)'
          do i=1,pedra_surf_n_tessere
             write(778,*) qr_tp(i)
          enddo
          write(778,*) 'Reaction-field polarization charges prev (PCM)'
          do i=1,pedra_surf_n_tessere
             write(778,*) qr_tp2(i)
          enddo
          write(778,*) 'Molecular potential'
          do i=1,pedra_surf_n_tessere
             write(778,*) pot_tp(i), pot_tp2(i)
          enddo
          if (global_medium_Floc.eq.'loc') then
             write(778,*) 'Local-field polarization charges (PCM)'
             do i=1,pedra_surf_n_tessere
                write(778,*) qx_tp(i)
             enddo
             write(778,*) 'Local-field polarization charges prev (PCM)'
             do i=1,pedra_surf_n_tessere
                write(778,*) qx_tp2(i)
             enddo
             write(778,*) 'External-field potential'
             do i=1,pedra_surf_n_tessere
                write(778,*) potf_tp(i), potf_tp2(i)
             enddo
          endif
          if(global_eps_Feps.eq."gen") then
             if( global_medium_Fbem.eq."stan" ) then
                 write(778,*) 'Reaction-field pol-charges poles(PCM)'
                 write(778,*) npoles
                 do i=1,pedra_surf_n_tessere
                     write(778,*) qr_t_p(i,:)
                 enddo
                 write(778,*) 'Dereaction-field pol-charges poles(PCM)'
                 do i=1,pedra_surf_n_tessere
                     write(778,*) dqr_t_p(i,:)
                enddo
                write(778,*) 'Reaction-field force poles(PCM)'
                do i=1,pedra_surf_n_tessere
                     write(778,*) fqr_t_p(i,:)
                enddo
                if (global_medium_Floc.eq.'loc') then
                     write(778,*) 'Local-field pol-charges poles(PCM)'
                     do i=1,pedra_surf_n_tessere
                         write(778,*) qx_t_p(i,:)
                     enddo
                     write(778,*) 'Delocal-field pol-charges poles(PCM)'
                     do i=1,pedra_surf_n_tessere
                         write(778,*) dqx_t_p(i,:)
                     enddo
                     write(778,*) 'Local-field force poles (PCM)'
                     do i=1,pedra_surf_n_tessere
                         write(778,*) fqx_t_p(i,:)
                     enddo
                endif
             endif
          elseif(global_eps_Feps.eq.'drl') then
             write(778,*) 'De reaction-field polarization charges (PCM)'
                do i=1,pedra_surf_n_tessere
                   write(778,*) dqr_t(i)
                enddo
             write(778,*) 'Reaction-field force (PCM)'
                do i=1,pedra_surf_n_tessere
                   write(778,*) fqr_t(i)
                enddo
             if (global_medium_Floc.eq.'loc') then
                write(778,*) 'De local-field polarization charges (PCM)'
                do i=1,pedra_surf_n_tessere
                   write(778,*) dqx_t(i)
                enddo
                write(778,*) 'Local-field force (PCM)'
                do i=1,pedra_surf_n_tessere
                   write(778,*) fqx_t(i)
                enddo
             endif
          endif
       else
         if (istep.eq.quantum_n_res) then 
           write(*,*) 'Error: restart is not implemented yet for other ', &
                    'interaction types other than PCM. wrt'
         endif
       endif

       close(778)

       return

      end subroutine wrt_restart_mdm

!------------------------------------------------------------------------
! @brief Read restart
!
! @date Created   : E. Coccia 28 Nov 2017
! Modified  : G. Gil 02 Jul 2018
!------------------------------------------------------------------------
      subroutine read_medium_restart()

       implicit none

       integer(i4b)     :: i,j
       character(3)     :: cdum
       logical          :: exist

       inquire(file='restart_mdm', exist=exist)
       if (exist) then
          open(779, file='restart_mdm', status="old")
       else
          write(*,*) 'ERROR:  file restart_mdm is missing'
          stop
       endif

!       if (global_prop_Fint.eq.'ons') then
!          read(779,*) cdum
!          do i=1,3
!             read(779,*) fr_t(1), fr_t(2), fr_t(3)
!          enddo
!          if (global_medium_Floc.eq.'loc') then
!             read(779,*) cdum
!             do i=1,3
!                read(779,*) fx_t(1), fx_t(2), fx_t(3)
!             enddo
!          endif
!       elseif (global_prop_Fint.eq.'pcm') then
       if (global_prop_Fint.eq.'pcm') then
          read(779,*) cdum
          do i=1,pedra_surf_n_tessere
             read(779,*) qr_tp(i)
          enddo
          read(779,*) cdum
          do i=1,pedra_surf_n_tessere
             read(779,*) qr_tp2(i)
          enddo
          read(779,*) cdum
          do i=1,pedra_surf_n_tessere
             read(779,*) pot_tp(i), pot_tp2(i)
          enddo
          if (global_medium_Floc.eq.'loc') then
             read(779,*) cdum
             do i=1,pedra_surf_n_tessere
                read(779,*) qx_tp(i)
             enddo
             read(779,*) cdum
             do i=1,pedra_surf_n_tessere
                read(779,*) qx_tp2(i)
             enddo
             read(779,*) cdum
             do i=1,pedra_surf_n_tessere
                read(779,*) potf_tp(i), potf_tp2(i)
             enddo
          endif
          if(global_eps_Feps.eq."gen") then
             if (global_medium_Fbem.eq.'stan') then
                read(779,*) cdum
                read(779,*) npoles
                do i=1,pedra_surf_n_tessere
                   read(779,*) (qr_tp_p(i,j), j=1,npoles)
                enddo
                read(779,*) cdum
                do i=1,pedra_surf_n_tessere
                   read(779,*) (dqr_tp_p(i,j), j=1,npoles)
                enddo
                read(779,*) cdum
                do i=1,pedra_surf_n_tessere
                   read(779,*) (fqr_tp_p(i,j), j=1,npoles)
                enddo
                if (global_medium_Floc.eq.'loc') then
                   read(779,*) cdum
                   do i=1,pedra_surf_n_tessere
                      read(779,*) (qx_tp_p(i,j), j=1,npoles)
                   enddo
                   read(779,*) cdum
                   do i=1,pedra_surf_n_tessere
                      read(779,*) (dqx_tp_p(i,j), j=1,npoles)
                   enddo
                   read(779,*) cdum
                   do i=1,pedra_surf_n_tessere
                      read(779,*) (fqx_tp_p(i,j), j=1,npoles)
                   enddo
                endif
             endif
          elseif(global_eps_Feps.eq."drl") then
             read(779,*) cdum
             do i=1,pedra_surf_n_tessere
                read(779,*) dqr_tp(i)
             enddo
             read(779,*) cdum
             do i=1,pedra_surf_n_tessere
                read(779,*) fqr_tp(i)
             enddo
             if (global_medium_Floc.eq.'loc') then
                read(779,*) cdum
                do i=1,pedra_surf_n_tessere
                   read(779,*) dqx_tp(i)
                enddo
                read(779,*) cdum
                do i=1,pedra_surf_n_tessere
                   read(779,*) fqx_tp(i)
                enddo
             endif
          endif
       else
         write(*,*) 'Error: restart is not implemented yet for other ',&
                    'interaction types other than PCM. read'
       endif

       close(779)

       return

      end subroutine read_medium_restart

!------------------------------------------------------------------------
! @brief Broadcast input data
!
! @date Created: E. Coccia 24/4/18
! Modified:
!------------------------------------------------------------------------
      subroutine mpibcast_readio_mdm()

#ifdef MPI

       call mpi_bcast(global_sys_Fwrite,           flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_sys_Ftest,            flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_sys_Fdeb,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(quantum_Ffld,                flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_Fopt_chr,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_Fcalc,                flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_MPL_ord,              flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_out_Fgamess,          flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_out_Fprint_lf_matrix, flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_prop_Fprop,           flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(blobal_prop_Fmdm_relax,      flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_prop_Finit_int,       flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_prop_Fint,            flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_eps_Feps,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Floc,          flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Fpol,          flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Fmdm,          flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Fbem,          flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Fqbem,         flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_read_write,      flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_medium_Finit,         flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(sph_surf_Fshape,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(global_surf_Fsurf,           flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(pedra_surf_Finv,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(pedra_surf_Fcav,             flg,MPI_CHARACTER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       if((global_prop_Fprop.eq.'chr-ief').or.&
          (global_prop_Fprop.eq.'chr-ied').or.&
          (global_prop_Fprop.eq.'chr-ons')) then
              call mpi_bcast(pedra_surf_n_tessere,    1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
      endif
       if (global_prop_Fprop.eq.'dip') then
              call mpi_bcast(sph_surf_nsph,       1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       endif
       if (tp_myrank.ne.0) then
       if((global_prop_Fprop.eq.'chr-ief').or.&
          (global_prop_Fprop.eq.'chr-ied').or.&
          (global_prop_Fprop.eq.'chr-ons')) then
             allocate(quantum_vts(pedra_surf_n_tessere,quantum_n_ci,quantum_n_ci))
             allocate(quantum_vtsn(pedra_surf_n_tessere))
          elseif (global_prop_Fprop.eq.'dip') then
             allocate(sph_surf_maj(sph_surf_nsph))
             allocate(sph_surf_min(sph_surf_nsph))
             allocate(sph_surf_vrs(3,3,sph_surf_nsph))
             allocate(sph_surf_centre(3,sph_surf_nsph))
          endif
       endif

       !Send input variables to the slaves
       call mpi_bcast(global_prop_n_q,                  1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(ntst,                             1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)

       call mpi_bcast(max_cycles,                       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(debye_eps_0,                      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(drudel_eps_gm,                    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(debye_eps_d,                      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(drudel_eps_A,                     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(drudel_eps_w0,                    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(drudel_f_vel,                     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(debye_tau_deb,                    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(thrshld,                   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       call mpi_bcast(mix_coef,                  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       if (general_prop_Fprop.eq.'dip') then
          call mpi_bcast(sph_surf_min,           sph_surf_nsph,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
          call mpi_bcast(sph_surf_maj,           sph_surf_nsph,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
          call mpi_bcast(sph_surf_centre,        3*sph_surf_nsph,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
          call mpi_bcast(sph_surf_vrs,           3*3*sph_surf_nsph,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       endif
       if((global_prop_Fprop.eq.'chr-ief').or.&
          (global_prop_Fprop.eq.'chr-ied').or.&
          (global_prop_Fprop.eq.'chr-ons')) then
          call mpi_bcast(quantum_vtsn,   pedra_surf_n_tessere,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
          call mpi_bcast(quantum_vts,    pedra_surf_n_tessere*quantum_n_ci*quantum_n_ci,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
       endif

#endif

       return

      end subroutine mpibcast_readio_mdm



      end module
