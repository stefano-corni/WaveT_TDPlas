      Module td_ContMed       
      use constants    
      use global_tdplas
      use readio_medium
      use pedra_friends
      use MathTools
      use BEM_medium
      use scf         
      use, intrinsic :: iso_c_binding

      implicit none
! 
      real(dbl) :: t !< time variable. Subscripts _t _tp _tp2 ... _0 refer respectively to dynamic variables at times t, t-dt, t-2dt, ... 0.
! Molecular observables and Maxwell potential/field
      ! c_tp: input from propagate, coefficients of states at time tp
      ! f_tp: input from propagate, Maxwell field at time tp
      real(dbl) :: f_tp2(3)                                   !< old Maxwell field stored in td_contmed  
      real(dbl) :: mu_tp(3),mu_tp2(3),mu_0(3)                 !< molecular dipole 
      real(dbl), allocatable :: pot_tp(:),pot_tp2(:),pot_0(:) !< molecular potential on BEM points
      real(dbl), allocatable :: potf_tp(:),potf_tp2(:)        !< Maxwell potential on BEM points
      !SP mu_0 contains mut(:,1,1)  
      !SP pot_0 contains vts(:,1,1) for Fint="ief" and potential of mu_0 if Fint='ons'
! Interaction and medium description
      real(dbl), allocatable :: h_mdm(:,:),h_mdm_0(:,:) !< medium contribution to the hamiltonian 
      real(dbl), allocatable :: q_mdm(:)                !< medium total charges at current time            
      real(dbl), allocatable :: mu_mdm(:,:)             !< medium total dipole for each spheroid/cavity at current time
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
      ! Fields and potentials
      real(dbl) :: fr_t(3),fr_tp(3)                     !< reaction field on the molecule centre of charge
      real(dbl) :: dfr_t(3)                             !< reaction field difference (fr_t-fr_tp) 
      real(dbl) :: fx_t(3),fx_tp(3)                     !< external field (fx) on the molecule centre of charge
!SC 06/02/16: eq and neq free energies;  
      real(dbl) :: e_vac,g_eq,g_eq_gs                   !< energy/equilibrium Free energies
      real(dbl) :: g_neq2,g_neq2_0,g_neq_0,g_neq1,g_neq1_part !< Non Equilibrium Free energies
! Working variables
      real(dbl) :: f1,f2,f3,f4,f5 !< constant factors (calculated once and for all) for velocity-verlet propagator (Drude-Lorentz) 
      real(dbl) :: dip(3)                               !< used in a test 
      real(dbl) :: ref                                  !< reference value in different debug tests    
      real(dbl) :: qtot                                 !< total BEM charge    
      real(dbl) :: taum1                                !< onsager tau-1 for charges single-tau propagation 
      real(dbl) :: qtot0                                !< total charge at time 0                           
      integer(i4b) :: file_med=11                       !< medium output file number 
      save
      private
!SC 07/02/16: added output_gneq
      public init_mdm,prop_mdm,finalize_mdm,qtot,ref,get_gneq, &
             get_ons,get_mdm_dip,set_charges
      contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   INTERFACE ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine init_mdm(c_tp,f_tp,h_int)   
!------------------------------------------------------------------------
! @brief Medium initialization called by WaveT or other programs 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      real(dbl), intent(INOUT) :: f_tp(3)
      complex(cmp), intent(INOUT) :: c_tp(n_ci)
      real(dbl), intent(INOUT):: h_int(n_ci,n_ci)
      integer(i4b) :: its,i,j                   
      character(20) :: name_f

 
! OPEN FILES
      write(name_f,'(a9,i0,a4)') "medium_t_",n_f,".dat"
      !if (Fmdm_res.eq.'Nonr') then
      open (file_med,file=name_f,status="unknown")
      write(file_med,*) "# step  time  dipole  field  qtot  qtot0"
      !elseif (Fmdm_res.eq.'Yesr') then
      !   open (file_med,file=name_f,status="unknown",position='append')
      !   if (Fint.eq.'pcm') then
      !      allocate(qr_i(nts_act))
      !      if (Floc.eq.'loc') allocate(qx_i(nts_act))
      !   endif 
      !   call read_medium_restart()
      !endif
      allocate(h_mdm(n_ci,n_ci),h_mdm_0(n_ci,n_ci))
      h_mdm=zero
      h_mdm_0=zero
      f_tp2=f_tp
      if(Fprop(1:3).eq."dip") then
! SC: First dipole propagation...
        call do_MPL_prop  !in BEM_medium
        call init_dip_and_field(c_tp,f_tp)
      else 
        if(.not.allocated(mu_mdm))allocate(mu_mdm(3,1))
! SC: ...then charges propagation        
        call do_BEM_prop !in BEM_medium
! SC 03/05/2016: create a new BEM_Q0=BEM_Qw^-1*BEM_Qf that should avoid
!                spurious charge dynamics for stationary states
!        call init_BEM_Q0
        call init_potential(c_tp,f_tp)
        call init_charges(c_tp)
!EC: restart values
        !if (Fmdm_res.eq.'Yesr') then
         !if (Fint.eq.'ons') then
        !    fr_t=fr_i
        !    if (Floc.eq.'loc') fx_t=fx_i
         !elseif (Fint.eq.'pcm') then
        !    qr_t=qr_i
        !    qr_tp=qr_i
        !    if (Floc.eq.'loc') then
        !       qx_t=qx_i
        !       qx_tp=qx_i
        !    endif
         !endif 
        !endif
! SC: predifine the factors used in the VV propagator, used for
! Drude-Lorentz
      endif
      if(Feps.eq."drl") call init_vv_propagator
      if (Fmdm(1:3).ne.'vac') call correct_hamiltonian
! SC set the initial values of the solvent component of the 
! neq free energies
      g_neq1=zero
      g_neq2=zero

      return


      end subroutine init_mdm


      subroutine prop_mdm(i,c_tp,f_tp,h_int)   
!------------------------------------------------------------------------
! @brief Medium propagation called by WaveT or other programs 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(INOUT) :: f_tp(3)
       complex(cmp), intent(IN) :: c_tp(n_ci)
       real(dbl), intent(INOUT):: h_int(n_ci,n_ci)
       integer(i4b), intent(IN) :: i                    
       integer(i4b) :: its,k,j                    

 
      ! Propagate medium only every n_q timesteps  
      !  otherwise, just resum the interaction hamiltonian and exit
       if(mod(i,n_q).ne.0) then
         ! SP 230916: added to perform tests on the local field
         h_int(:,:)=h_int(:,:)+h_mdm(:,:)
         return
       endif
       t=(i-1)*dt
       call do_dip_from_coeff(c_tp,mu_tp,n_ci)

       if (Fprop(1:3).eq."dip") then
       ! Dipole propagation: 
         call prop_dip(c_tp,f_tp)
         call do_gneq(c_tp,mut,dfr_t,fr_t,fr_0,mat_fd,3,-1)
         if(Fmdm(2:4).eq."nan") then
           mu_mdm=mr_t
           if((Ftest(2:3).eq."-l").or.(Fdeb.eq."off")) mu_mdm=zero
           if(Floc.eq."loc") mu_mdm=mu_mdm+mx_t
         endif
       else
       ! Charges propagation: 
         ! Calculate external potential on tesserae for local field       
         if(Floc.eq."loc") call do_pot_from_field(f_tp,potf_tp)
         ! Calculate the molecule potential on tesserae
! SP 26/06/17: changed to use general MathTools  
         if(Fint.eq.'ons') then
           call do_dip_from_coeff(c_tp,dip,n_ci)
           call do_pot_from_dip(dip,pot_tp)
         else 
           call do_pot_from_coeff(c_tp,pot_tp)
         endif

         call prop_chr
         ! Calculate medium's dipole from charges 
         qtot=zero
         mu_mdm=zero
         ! SP 26/06/17: MathUtils, do_dip_from_charges updates the value in mu_mdm
         call do_dip_from_charges(qr_t,mu_mdm(:,1),qtot)        
         if((Ftest(2:3).eq."-l").or.(Fdeb.eq."off")) mu_mdm=zero
         if(Floc.eq."loc") call do_dip_from_charges(qx_t,mu_mdm(:,1),qtot)
         ! Calculate Reaction Field from charges
         if (Fint.eq.'ons') dfr_t=-fr_t
         call do_field_from_charges(qr_t,fr_t)
         if (Fint.eq.'ons') dfr_t=dfr_t+fr_t
         ! Calculate Local Field from external charges
         if(Floc.eq."loc") call do_field_from_charges(qx_t,fx_t)
       ! SC calculate free energy:
         if (Fint.eq.'ons') then 
! SP 10/07/17: commented the following, with Fint=ons do_gneq should be changed
           !call do_gneq(c_tp,mut,dfr_t,fr_t,fr_0,mat_fd,3,-1)
         else
           call do_gneq(c_tp,vts,dqr_t,qr_t,q0,BEM_Qd,nts_act,1)
         endif
       endif
       ! Build the interaction Hamiltonian Reaction/Local
       if(Fdeb.ne."off") call do_interaction_h
       if (i.eq.1.and.Fwrite.eq."high") then
        write(6,*) "h_mdm at the first propagation step"
        do j=1,n_ci
         do k=1,n_ci
          write (6,*) j,k,h_mdm(j,k)
         enddo
        enddo
       endif
       ! SP 230916: added to perform tests on the local/reaction field
       if(Ftest(2:2).eq."-") call do_ref(c_tp)
       ! Update the interaction Hamiltonian 
       h_int(:,:)=h_int(:,:)+h_mdm(:,:)

       ! SP 24/02/16  Write output
       if (mod(i,n_out).eq.0.or.i.eq.1) call out_mdm(i)

       ! EC 28/11/17 Write restart
       !if (mod(i,n_res).eq.0) call wrt_restart_mdm()


       return

      end subroutine


      subroutine finalize_mdm   
!------------------------------------------------------------------------
! @brief Medium finalization called by WaveT or other programs 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       deallocate(h_mdm,h_mdm_0)
       call finalize_prop
       if (Fprop(1:3).eq."chr") then
         call deallocate_BEM_public    
       elseif (Fprop(1:3).eq."dip".or.Fint.eq."ons") then
         call deallocate_MPL_public    
       endif
       !if (Fmdm_res.eq.'Yesr') then
       !   if (Fint.eq.'pcm') then
       !      deallocate(qr_i)
       !      if (Floc.eq.'loc') then
       !         deallocate(qx_i)
       !      endif
       !   endif 
       !endif
 
       close (file_med)

       return

      end subroutine


      subroutine get_gneq(e_vac_t,g_eqr_t,g_neqr_t,g_neq2_t)
!------------------------------------------------------------------------
! @brief Solvent contribution to free energies 
!
! @date Created: S. Corni 
! Modified:
!------------------------------------------------------------------------

       real(dbl),intent(inout):: e_vac_t,g_eqr_t,g_neqr_t,g_neq2_t

       e_vac_t=e_vac_t+e_vac
       g_eqr_t=g_eq+g_eqr_t
       g_neqr_t=g_neq1+g_neqr_t
       g_neq2_t=g_neq2+g_neq2_t

       return

      end subroutine


      subroutine get_mdm_dip(mdm_dip)
!------------------------------------------------------------------------
! @brief Medium dipole for spectra 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl),intent(out):: mdm_dip(3)

       mdm_dip(:)=mu_mdm(:,1)

       return

      end subroutine


      subroutine set_charges(q)
!------------------------------------------------------------------------
! @brief Set charges 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl),intent(in):: q(nts_act)

       qr_t=q

       return

      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!   INTERACTION  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!------------------------------------------------------------------------
!> Build the interaction matrix in the molecular unperturbed state basis 
!------------------------------------------------------------------------
      subroutine do_interaction_h
!------------------------------------------------------------------------
! @brief Build the interaction matrix in the molecular unperturbed state
! basis 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b):: i,j
! SC 19/03/2016: added a temporarty array to avoid summing the local field
! to ther reaction field for onsager
! SC 19/03/2016: moved from prop_mdm (so that h_mdm has values assined only here 
       h_mdm(:,:)=zero

       if (Fint.eq.'ons') then
         f_mdm(:)=fr_t(:)
         if(Floc.eq."loc") f_mdm(:)=f_mdm(:)+fx_t(:) 
         h_mdm(:,:)=-h_mdm_0(:,:)-mut(1,:,:)*f_mdm(1) &
                          -mut(2,:,:)*f_mdm(2) -mut(3,:,:)*f_mdm(3)
       elseif(Fint.eq.'pcm') then
         q_mdm=qr_t
         if(Floc.eq."loc") q_mdm(:)=qr_t(:)+qx_t(:)
! SC 31/10/2016: avoid including interaction with an unwanted net charge
         q_mdm=q_mdm+(qtot0-sum(q_mdm))/nts_act
         do j=1,n_ci   
           do i=1,j       
             h_mdm(i,j)=-h_mdm_0(i,j)+dot_product(q_mdm(:),vts(:,i,j))
             h_mdm(j,i)=h_mdm(i,j)
           enddo
         enddo
       else
         write(*,*) "wrong interaction type "
         stop
       endif

       return

      end subroutine


      subroutine correct_hamiltonian
!------------------------------------------------------------------------
! @brief Build interaction matrix at time 0 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       integer(4) :: its,i,j 

       if (Fint.eq.'ons') then
         h_mdm_0(:,:)=-mut(1,:,:)*fr_0(1)-mut(2,:,:)*fr_0(2) &
                                         -mut(3,:,:)*fr_0(3)
       elseif (Fint.eq.'pcm') then
!SC 04/05/2016: updated to be coerent with both GS and SCF initialization
!SC 27/09/2016: corrected bug introduced previously
         q_mdm=q0+(qtot0-sum(q0))/nts_act
         do its=1,nts_act     
           h_mdm_0(:,:)=h_mdm_0(:,:)+q_mdm(its)*vts(its,:,:)
         enddo
       endif
       if(Fwrite.eq."high") then
         do i=1,n_ci
          do j=1,n_ci
           write(6,*) i,j,h_mdm_0(i,j)
          enddo
         enddo
       endif

       return

      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialization/deallocation !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine init_potential(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Initialize potentials for propaagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       complex(cmp), intent(IN) :: c_tp(n_ci)
       real(dbl), intent(IN) :: f_tp(3)
       complex(cmp) :: c_gs(n_ci)
       integer(i4b) :: its  

       allocate(pot_tp(nts_act))
       allocate(pot_0(nts_act))
! SC 08/04/2016: a routine to test by calculating the potentials from the dipoles
       if (Fdeb.eq.'vmu') call do_vts_from_dip
! SP 26/06/17: changed to use general MathTools  
       if(Fint.eq.'ons') then
         call do_dip_from_coeff(c_tp,dip,n_ci)
         call do_pot_from_dip(dip,pot_tp)
! SP 10/07/17: commented the following, is it really needed here? Should stay in the init_charges
         !call do_field_from_charges(q0,fr_tp)
       else 
         call do_pot_from_coeff(c_tp,pot_tp)
       endif
       c_gs(:)=zeroc
       c_gs(1)=onec
! SP 26/06/17: changed to use general MathTools  
       if(Fint.eq.'ons') then
         call do_dip_from_coeff(c_gs,dip,n_ci)
         call do_pot_from_dip(dip,pot_0)
       else 
         call do_pot_from_coeff(c_gs,pot_0)
       endif
       ! Test for a sudden switch-on of an unpolarizable molecular dipole
       if(Ftest.eq."s-r") pot_0=zero
       if(Ftest.eq."s-r") pot_tp=zero
       allocate(pot_tp2(nts_act))
       pot_tp2=pot_tp
       if(Floc.eq."loc") then
         allocate(potf_tp(nts_act))
         allocate(potf_tp2(nts_act))      
         call do_pot_from_field(f_tp,potf_tp)
         potf_tp2=potf_tp    
! SP 09/07/16 commented the following 
         !fx_t(:)=zero
       endif

       return

      end subroutine


      subroutine init_charges(c_tp)
!------------------------------------------------------------------------
! @brief Initialize charges for propagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       complex(cmp), intent(INOUT) :: c_tp(n_ci)
       integer(i4b):: its
       real(dbl), allocatable :: qd(:)

       allocate(qd(nts_act))
       allocate(qr_t(nts_act))
       allocate(q_mdm(nts_act))
       allocate(qr_tp(nts_act))
       allocate(dqr_t(nts_act))
       if (.not.allocated(q0)) allocate (q0(nts_act))
       ! init the state and the RF before propagation
!SP 29/05/16: pot_0 replaces vts(:,1,1) to allow treating Fprop=ief and Fint=ons
       select case(Finit_mdm)
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
       if(Fmdm(2:4).eq.'nan') then
         !q0=q0-qtot0/nts_act
         qtot0=zero
       endif
       g_eq_gs=0.5d0*dot_product(q0,pot_0)
       write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       ! see readio_medium for Finit_int
       select case (Finit_int)
        case ('nsc')

!         g_neq_0=0.5*MPL_fd*dot_product(mu_tp,mu_tp) &
!                -MPL_fd*dot_product(mu_tp,mu_0) &
!                +0.5*MPL_fd*dot_product(mu_0,mu_0)
!         g_neq_0=-g_neq_0

! SC in principle a non equilibrium self consistency if eps_d=1 is needed
!    here we use the non self-consitent density but use the correct RF
         qd=matmul(BEM_Qd,pot_0)
         g_neq_0=0.5d0*dot_product(qd,pot_0)
         g_neq2_0=g_neq_0
         g_neq_0=g_neq_0-dot_product(qd,pot_tp)
         qd=matmul(BEM_Qd,pot_tp)
         g_neq_0=g_neq_0+0.5d0*dot_product(qd,pot_tp)
         qr_tp=q0+matmul(BEM_Qd,(pot_tp-pot_0))
        case ('sce')
         qr_tp=mix_coef*matmul(BEM_Q0,pot_tp)+(1.-mix_coef)*q0
         call do_scf(qr_tp,c_tp)
! update the potential
!         do its=1,nts_act
!          pot_tp(its)=dot_product(c_tp,matmul(vts(its,:,:),c_tp))
!         enddo
! SP 26/06/17: changed to use general MathTools  
         if(Fint.eq.'ons') then
           call do_dip_from_coeff(c_tp,dip,n_ci)
           call do_pot_from_dip(dip,pot_tp)
         else 
           call do_pot_from_coeff(c_tp,pot_tp)
         endif
         g_neq_0=0.5*dot_product(qr_tp,pot_tp)
         qr_tp=matmul(BEM_Qd,pot_tp)
         g_neq2_0=0.5*dot_product(qr_tp,pot_tp)
         qr_tp=matmul(BEM_Q0,pot_tp)
       end select
       write(6,*) 'G_neq at t=0:',g_neq_0
       qr_t(:)=qr_tp(:)
       dqr_t(:)=zero  
       if(Fint.eq."ons") call do_field_from_charges(qr_t,fr_0)
       if(Floc.eq."loc") then
         allocate(qx_t(nts_act))
         allocate(qx_tp(nts_act))
         allocate(dqx_t(nts_act))
         ! SP 09/07/17 changed the following for coherence  
         qx_t(:)=zero                          
         qx_tp(:)=zero    
         dqx_t(:)=zero 
       endif
       if(Feps.eq."drl") then
         allocate(dqr_tp(nts_act))
         allocate(fqr_tp(nts_act))
         allocate(fqr_t(nts_act))
         dqr_tp(:)=zero
         fqr_tp=zero
         if(Floc.eq."loc") then
           allocate(dqx_tp(nts_act))
           allocate(fqx_tp(nts_act))
           allocate(fqx_t(nts_act))
           dqx_tp(:)=zero
           fqx_tp=zero
         endif
       endif
       deallocate(qd)

       return

      end subroutine


      subroutine init_dip_and_field(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Initialize dipoles and field for propagation  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       complex(cmp), intent(INOUT) :: c_tp(n_ci)
       real(dbl), intent(INOUT) :: f_tp(3)

       allocate(mu_mdm(3,nsph))
       allocate(mr_0(3,nsph))
       allocate(mr_t(3,nsph))
       allocate(mr_tp(3,nsph))
       allocate(dmr_t(3,nsph))
       if(Floc.eq."loc") then
         allocate(mx_t(3,nsph))
         allocate(mx_tp(3,nsph))
         allocate(dmx_t(3,nsph))
       endif
       if(Feps.eq."drl") then
         allocate(dmr_tp(3,nsph))
         allocate(fmr_t(3,nsph))
         allocate(fmr_tp(3,nsph))
         dmr_tp=zero
         fmr_tp=zero
         if(Floc.eq."loc") then
           allocate(dmx_tp(3,nsph))
           allocate(fmx_t(3,nsph))
           allocate(fmx_tp(3,nsph))
           dmx_tp=zero
           fmx_tp=zero
         endif
       endif
       ! The following subroutines should be merged!!
       if(Fshape.eq."sphe") call init_dip_sphe(c_tp,f_tp)
       if(Fshape.eq."spho") call init_dip_spho(c_tp,f_tp)

       return

      end subroutine
!     
!------------------------------------------------------------------------
!> Initialize Dipole for propagation with spherical object                
!------------------------------------------------------------------------
      subroutine init_dip_sphe(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Initialize dipole for propagation with spherical object
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       complex(cmp), intent(INOUT) :: c_tp(:)
       real(dbl), intent(IN) :: f_tp(3)
       real(dbl) :: fld(3),fld0(3)
       integer(i4b) :: i

       mu_0(:)=mut(:,1,1)
       mr_0=zero
       fr_0=zero
       call do_dip_from_coeff(c_tp,mu_tp,n_ci)
       if(Ftest.eq."s-r") mu_0=zero
       if(Ftest.eq."s-r") mu_tp=zero
       if (Fmdm(2:4).eq."nan") then 
       ! SP 06/07/17: need to implement a proper scf cycle for more then a nanoparticle
       !              only non-self consistent initialization is done
         ! SP compute the field acting on spheres/oids given m_0=zero
         fr_t=zero                   
         do i=1,nsph
           call do_field_from_dip(mu_0,mol_cc,fld0,sph_centre(:,i)) 
           if(Finit_mdm.eq."fro") mr_0(:,i)=ONS_f0*fld0
           call do_field_from_dip(mu_tp,mol_cc,fld,sph_centre(:,i)) 
! SP 11/07/17: to have the same initial RF for charges and dipoles           
           mr_t(:,i)=mr_0(:,i)+ONS_fd*(fld-fld0)
           call do_field_from_dip(mr_t(:,i),sph_centre(:,i),fld,mol_cc) 
           fr_t=fr_t+fld
         enddo
         mr_tp=mr_t
         if(Floc.eq."loc") then
           mx_t=zero
           mx_tp=mx_t
         endif
       elseif (Fmdm(2:4).eq."sol") then 
         if(Finit_mdm.eq."fro") fr_0(:)=ONS_f0*mu_0(:)
         call init_dip_sphe_sol(c_tp)
       endif
       fr_tp=fr_t
       if(Floc.eq."loc") then
         fx_t=zero         
         fx_tp=fx_t
       endif

       return

      end subroutine init_dip_sphe
!     
!------------------------------------------------------------------------
!> Initialize Dipole for propagation with spheroidal object                
!------------------------------------------------------------------------
      subroutine init_dip_spho(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Initialize dipole for propagation with spheroidal object 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! inizialize onsager 
       implicit none

       complex(cmp), intent(INOUT) :: c_tp(n_ci)
       real(dbl), intent(IN) :: f_tp(3)
       real(dbl) :: fld(3),fld0(3)
       integer(i4b) :: i

       ! Init molecular dipole 
       mu_0(:)=mut(:,1,1)
       mr_0=zero
       fr_0=zero
       call do_dip_from_coeff(c_tp,mu_tp,n_ci)
       if(Ftest.eq."s-r") mu_0=zero
       if(Ftest.eq."s-r") mu_tp=zero
       if (Fmdm(2:4).eq."nan") then 
       ! SP 06/07/17: need to implement a proper scf cycle for more then a nanoparticle
       !              only non-self consistent initialization is done
         ! SP compute the field acting on spheres/oids given m_0=zero
         fr_t=zero                   
         do i=1,nsph
           call do_field_from_dip(mu_0,mol_cc,fld0,sph_centre(:,i)) 
           if(Finit_mdm.eq."fro") mr_0(:,i)=matmul(mat_f0,fld0)
           call do_field_from_dip(mu_tp,mol_cc,fld,sph_centre(:,i)) 
           mr_t(:,i)=mr_0(:,i)+matmul(mat_fd,fld-fld0)
           call do_field_from_dip(mr_t,sph_centre(:,i),fld,mol_cc) 
           fr_t=fr_t+fld
         enddo
         mr_tp=mr_t
         if(Floc.eq."loc") then
           mx_t=zero                
           mx_tp=mx_t
         endif
       elseif (Fmdm(2:4).eq."sol") then 
         if(Finit_mdm.eq."fro") fr_0=matmul(mat_f0,mu_0)
         call init_dip_spho_sol(c_tp,f_tp)
         fr_tp=fr_t
       endif
       if(Floc.eq."loc") then
         fx_t=zero                            
         fx_tp=fx_t
       endif

       return

      end subroutine


      subroutine init_dip_spho_sol(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Initialize free energy (spheroidal) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       complex(cmp), intent(INOUT) :: c_tp(n_ci)
       real(dbl), intent(IN) :: f_tp(3)

       g_eq_gs=-0.5d0*dot_product(fr_0,mu_0)
       write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       select case (Finit_int)
        case ('nsc')
! SC in principle a non equilibrium self consistency if eps_d=1 is needed
!    here we use the non self-consitent dipole but use the correct RF
         fr_t=fr_0+matmul(mat_fd,mu_tp-mu_0)
         g_neq_0=0.5*dot_product(mu_tp,matmul(mat_fd,mu_tp)) &
! SP 04/0717 WARNING: NOT SURE OF THE FOLLOWING LINE 
                -dot_product(matmul(mat_fd,mu_tp),mu_0) &
                !-f_d*dot_product(mu_tp,mu_0) &
                +0.5*dot_product(mu_0,matmul(mat_fd,mu_0))
         g_neq_0=-g_neq_0
         g_neq2_0=-0.5*dot_product(mu_0,matmul(mat_fd,mu_0))
        case ('sce')
         fr_t=mix_coef*matmul(mat_f0,mu_tp)+(1.-mix_coef)*fr_0
         call do_scf(fr_t,c_tp)
         call do_dip_from_coeff(c_tp,mu_tp,n_ci)
         g_neq_0=-0.5*dot_product(mu_tp,matmul(mat_f0,mu_tp))
         g_neq2_0=-0.5*dot_product(mu_tp,matmul(mat_fd,mu_tp))
       end select
       write(6,*) 'G_neq at t=0:',g_neq_0

       return

      end subroutine


      subroutine init_dip_sphe_sol(c_tp)
!------------------------------------------------------------------------
! @brief Initialize free energy (spherical) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       complex(cmp), intent(INOUT) :: c_tp(n_ci)

       g_eq_gs=-0.5d0*dot_product(fr_0,mu_0)
       write(6,*) 'Medium contribution to ground state free energy:', &
                    g_eq_gs
       select case (Finit_int)
       case ('nsc')
! SC in principle a non equilibrium self consistency if eps_d=1 is needed
!    here we use the non self-consitent dipole but use the correct RF
         fr_t=fr_0+ONS_fd*(mu_tp-mu_0)
         g_neq_0=0.5*ONS_fd*dot_product(mu_tp,mu_tp) &
                -ONS_fd*dot_product(mu_tp,mu_0) &
                +0.5*ONS_fd*dot_product(mu_0,mu_0)
         g_neq_0=-g_neq_0
         g_neq2_0=-0.5*ONS_fd*dot_product(mu_0,mu_0)
       case ('sce')
         fr_t=mix_coef*ONS_f0*mu_tp+(1.-mix_coef)*fr_0
         call do_scf(fr_t,c_tp)
         call do_dip_from_coeff(c_tp,mu_tp,n_ci)
         g_neq_0=-0.5*ONS_f0*dot_product(mu_tp,mu_tp)
         g_neq2_0=-0.5*ONS_fd*dot_product(mu_tp,mu_tp)
       end select
       write(6,*) 'G_neq at t=0:',g_neq_0

       return

      end subroutine


      subroutine finalize_prop
!------------------------------------------------------------------------
! @brief Deallocate propagation variables 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

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

      end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! Propagation routines !!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SP 22/02/16: A general remark: matmul and DGEMV used randomly, 
!             DGEMV shoud be efficient for big matrices              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine prop_chr
!------------------------------------------------------------------------
! @brief Potential and charges propagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------


       ! Propagate
       if(Ftest.eq."s-r") pot_tp=vts(:,1,1) !Only for debug purposes
       if(Feps.eq."deb") then
         if(Fprop.eq."chr-ief") then
           call prop_ief_deb
         elseif(Fprop.eq."chr-ied") then
           taum1=(eps_0+one)/(eps_d+one)/tau_deb 
           call prop_ied_deb
         elseif(Fprop.eq."chr-ons") then
           call prop_ons_deb 
         else
           write(*,*) "not implemented yet"
           stop
         endif
       ! SP 22/02/16 dqr_t calculated for g_neq
         dqr_t=qr_t-qr_tp
       elseif(Feps.eq."drl") then
         if(Fprop.eq."chr-ief") then
           call prop_vv_ief_drl 
         elseif(Fprop.eq."chr-ons") then
           call prop_ons_drl ! uses vv
         else
           write(*,*) "Wrong propagation method"
           stop
         endif
       endif
       if(Fdeb.eq."equ") qr_t=matmul(BEM_Q0,pot_tp)
       qr_tp=qr_t
       pot_tp2=pot_tp
       if(Floc.eq."loc") then
         qx_tp=qx_t
         potf_tp2=potf_tp
       endif

       return

      end subroutine


      subroutine prop_dip(c_tp,f_tp)
!------------------------------------------------------------------------
! @brief Dipole and field propagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! evolve onsager field/dipole  
       implicit none

       complex(cmp), intent(IN) :: c_tp(n_ci)
       real(dbl), intent(IN):: f_tp(3)  

       ! calculate molecular dipole from CI coefficients
       if(Ftest.eq."s-r") mu_tp=mut(:,1,1) !Only for debug purposes
       ! propagate onsager reaction and local fields 
       if(Feps.eq."deb") then
         if(Fshape.eq."sphe") call prop_dip_sphe_deb(f_tp) 
         if(Fshape.eq."spho") call prop_dip_spho_deb(f_tp) 
       elseif(Feps.eq."drl") then
         if(Fshape.eq."sphe") call prop_dip_sphe_drl(f_tp) 
         if(Fshape.eq."spho") call prop_dip_spho_drl(f_tp) 
       endif
       dfr_t=fr_t-fr_tp
       fr_tp=fr_t
       mu_tp2=mu_tp
       if(Floc.eq."loc") then
         fx_tp=fx_t
         f_tp2=f_tp
       endif

       return

      end subroutine


      subroutine init_vv_propagator
!------------------------------------------------------------------------
! @brief Initialization for velocity Verlet propagation (vv) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       f1=dt*(1.d0-dt*0.5d0*eps_gm)
       f2=dt*dt*0.5d0
       f3=1.d0-dt*eps_gm*(1.d0-dt*0.5*eps_gm)
       f4=0.5d0*dt
       f5=eps_gm*f2
       write(6,*) "Initiated VV propagator"

       return

      end subroutine


      subroutine prop_ons_drl
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of charges with Onsager-BEM equations
! ! SP 07/07/17 vv propagaton 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! charge propagation with drude/lorentz and cosmo/onsager equations
       integer(i4b) :: its  

      ! Reaction Field
       qr_t=qr_tp+f1*dqr_tp+f2*fqr_tp
       ! $Q_f=-ONS_ff(1)*S{-1}$
       fqr_t=-ONS_fw*qr_t-ONS_ff(1)*matmul(BEM_sm1,pot_tp)
       dqr_t=f3*dqr_tp+f4*(fqr_t+fqr_tp)-f5*fqr_tp
       fqr_tp=fqr_t
       dqr_tp=dqr_t
       ! Local Field
       if(Floc.eq."loc") then
         qx_t=qx_tp+f1*dqx_tp+f2*fqx_tp
         fqx_t=-ONS_fw*qx_t-ONS_ff(1)*matmul(BEM_sm1,potf_tp)
         dqx_t=f3*dqx_tp+f4*(fqx_t+fqx_tp)-f5*fqx_tp
         fqx_tp=fqx_t
         dqx_tp=dqx_t
       endif

       return

      end subroutine


      subroutine prop_ons_deb
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of charges with Onsager-BEM equations
! ! SP 07/07/17 updated with vv propagaton 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! charge propagation with drude/lorentz and onsager equations
       integer(i4b) :: its  

      ! Reaction Field eq.47 Corni et al. JPCA 2015
       qr_t=qr_tp-dt*taum1*qr_tp-dt*taum1*ONS_f0*matmul(BEM_Sm1,pot_tp2)&
                                 -ONS_fd*matmul(BEM_Sm1,pot_tp-pot_tp2)
      ! Local Field 
       if(Floc.eq."loc") then
         qx_t=qx_tp-dt*taum1*qx_tp-dt*taum1*ONS_fx0*matmul(BEM_Sm1 &
              ,potf_tp2)-ONS_fxd*matmul(BEM_Sm1,potf_tp-potf_tp2)
                                  
       endif

       return

      end subroutine
!     
!------------------------------------------------------------------------
!> Drude-Lorentz propagation Fprop=chr-ief standard algorithm 
!------------------------------------------------------------------------
      subroutine prop_ief_drl
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation Fprop=chr-ief standard algorithm 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! Charge propagation with drude/lorentz and IEF equations 
       integer(i4b) :: its  

      ! Reaction Field
       qr_t=qr_tp+dt*dqr_tp
       call DGEMV('N',nts_act,nts_act,dt,BEM_Qf,nts_act,pot_tp,one_i, &
                     zero,dqr_t,one_i)
       call DGEMV('N',nts_act,nts_act,-dt,BEM_Qw,nts_act,qr_tp,one_i,  &
                     one,dqr_t,one_i)
! SC 17/8/2016: changed the following, 
!      dqr_t=dqr_t+(1-dt*eps_gm)*dqr_tp 
       dqr_t=(dqr_t+dqr_tp)/(1.d0+eps_gm*dt)
! SC avoid developing a total charge
       dqr_t=dqr_t-sum(dqr_t)/nts_act
       dqr_tp=dqr_t
      ! Local Field
       if(Floc.eq."loc") then
         qx_t=qx_tp+dt*dqx_tp
         call DGEMV('N',nts_act,nts_act,dt,BEM_Qf,nts_act,potf_tp,one_i,&
                       zero,dqx_t,one_i)
         call DGEMV('N',nts_act,nts_act,-dt,BEM_Qw,nts_act,qx_tp,one_i,&
                       one,dqx_t,one_i)
! SC 17/8/2016: changed the following, 
!        dqx_t=dqx_t+(1.d0-dt*eps_gm)*dqx_tp 
         dqx_t=(dqx_t+dqx_tp)/(1.d0+eps_gm*dt)
! SC avoid developing a total charge
         dqx_t=dqx_t-sum(dqx_t)/nts_act
         dqx_tp=dqx_t
       endif

       return

      end subroutine


      subroutine prop_vv_ief_drl
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation Fprop=chr-ief velocity-verlet
! algorithm 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! Charge propagation with drude/lorentz and IEF equations 
       integer(i4b) :: its  

!      qr_t=qr_tp+dt*(1.d0-dt*0.5d0*eps_gm)*dqr_tp+dt*dt*0.5d0*fqr_tp
!      fqr_t=-matmul(BEM_Qw,qr_t)+matmul(BEM_Qf,pot_tp)
!      dqr_t=(1.d0-dt*0.5d0*eps_gm)*dqr_tp+0.5d0*dt*(fqr_t+fqr_tp)
!      dqr_t=dqr_t/(1.d0+dt*0.5d0*eps_gm)
! SC integrator from E. Vanden-Eijnden, G. Ciccotti CPL 429 (2006) 310–316
!      qr_t=qr_tp+dt*(1.d0-dt*0.5d0*eps_gm)*dqr_tp+dt*dt*0.5d0*fqr_tp
!      fqr_t=-matmul(BEM_Qw,qr_t)+matmul(BEM_Qf,pot_tp)
!      dqr_t=(1.d0-dt*eps_gm*(1.d0-dt*0.5*eps_gm))*dqr_tp+ &
!           0.5d0*dt*(fqr_t+fqr_tp)-eps_gm*dt*dt*0.5*fqr_tp
!      f1=dt*(1.d0-dt*0.5d0*eps_gm)
!      f2=dt*dt*0.5d0
!      f3=1.d0-dt*eps_gm*(1.d0-dt*0.5*eps_gm)
!      f4=0.5d0*dt
!      f5=eps_gm*f2
       qr_t=qr_tp+f1*dqr_tp+f2*fqr_tp

       fqr_t=-matmul(BEM_Qw,qr_t)+matmul(BEM_Qf,pot_tp)
       dqr_t=f3*dqr_tp+f4*(fqr_t+fqr_tp)-f5*fqr_tp
       fqr_tp=fqr_t
       dqr_tp=dqr_t
      ! Local Field
       if(Floc.eq."loc") then
        qx_t=qx_tp+f1*dqx_tp+f2*fqx_tp
        fqx_t=-matmul(BEM_Qw,qx_t)+matmul(BEM_Qf,potf_tp)
        dqx_t=f3*dqx_tp+f4*(fqx_t+fqx_tp)-f5*fqx_tp
        fqx_tp=fqx_t
        dqx_tp=dqx_t
       endif

       return

      end subroutine


      subroutine prop_ief_deb
!------------------------------------------------------------------------
! @brief Debye propagation Fprop=chr-ief 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! Charge propagation with debye and IEF equations 
       integer(i4b) :: its  

      ! Reaction Field
       qr_t=qr_tp-dt*matmul(BEM_R,qr_tp)+dt*matmul(BEM_Qt,pot_tp) &
                                    +matmul(BEM_Qd,pot_tp-pot_tp2)
      ! Local Field eq.47 JPCA 2015
       if(Floc.eq."loc") then
         qx_t=qx_tp-dt*matmul(BEM_R,qx_tp)+dt*matmul(BEM_Qt,potf_tp) &
                                     +matmul(BEM_Qd,potf_tp-potf_tp2)
       endif

       return

      end subroutine


      subroutine prop_ied_deb
!------------------------------------------------------------------------
! @brief  Debye propagation Fprop=chr-ied  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

      ! Charge propagation with dedye and IEF equations one taud 
       integer(i4b) :: its  

      ! Reaction Field
       qr_t=qr_tp-dt*taum1*qr_tp+dt*taum1*matmul(BEM_Q0,pot_tp) &
                                  +matmul(BEM_Qd,pot_tp-pot_tp2)
      ! Local Field eq.47 JPCA 2015
       if(Floc.eq."loc") then
         qx_t=qx_tp-dt*taum1*qx_tp+dt*taum1*matmul(BEM_Q0,potf_tp) &
                                   +matmul(BEM_Qd,potf_tp-potf_tp2)
       endif

       return

      end subroutine


      subroutine prop_dip_spho_deb(f_tp)
!------------------------------------------------------------------------
! @brief  Debye propagation of the spheroid reaction and local fields  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(IN):: f_tp(3)  
       real(dbl) :: f(3),fp(3),mp(3),mp2(3)  
       integer(i4b) :: i

       ! Equation (23) JCP 146, 064116 (2017)
       ! SP 05/07/17: this is easily genealizable to multipoles see JCP 146, 064116 (2017) 
       fr_t=zero 
       do i=1,3
         ! Determine previous field component along spheroid axis:
         fp=dot_product(fr_tp,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
         mp=dot_product(mu_tp,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
         mp2=dot_product(mu_tp2,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
         f=fp+matmul(MPL_Fd(:,:,i,1),mp-mp2)+dt*n_q*  &
             (matmul(MPL_Ft0(:,:,i),mp)-MPL_Taum1(i)*fp)
         if(Fdeb.eq."equ") f=matmul(MPL_F0(:,:,i,1),mu_tp) !Equilibrium RF for debug purposes
         fr_t=fr_t+f
       enddo
       ! Local Field
       if(Floc.eq."loc") then
         ! Equation (23) JCP 146, 064116 (2017)
         fx_t=zero
         do i=1,3
           fp=dot_product(fx_tp,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
           mp=dot_product(mu_tp,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
           mp2=dot_product(mu_tp2,sph_vrs(:,i,1))*sph_vrs(:,i,1) 
           f=fp+matmul(MPL_Fxd(:,:,i),mp-mp2)+dt*n_q*  &
               (matmul(MPL_Ftx0(:,:,i),mp)-MPL_Taum1(i)*fp)
           fx_t=fx_t+f
         enddo
       endif

       return

      end subroutine


      subroutine prop_dip_sphe_deb(f_tp)
!------------------------------------------------------------------------
! @brief  Debye propagation of the sphere reaction and local fields  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(IN):: f_tp(3)  
       real(dbl) :: f(3),fp(3),mt(3),mp(3)  
       integer(i4b) :: i

       fr_t=fr_tp+ONS_fd*(mu_tp-mu_tp2)+dt*n_q*ONS_taum1* &
           (ONS_f0*mu_tp-fr_tp)
       if(Fdeb.eq."equ") fr_t=ONS_f0*mu_tp
       ! Local Field
       if(Floc.eq."loc") then
        fx_t=fx_tp+ONS_fxd*(f_tp-f_tp2)+dt*n_q*ONS_tauxm1* &
           (ONS_fx0*f_tp-fx_tp)
       endif

       return

      end subroutine


      subroutine prop_dip_spho_drl(f_tp)
!------------------------------------------------------------------------
! @brief Drude-Lorentz propagation of the spheroid reaction and local fields  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(IN):: f_tp(3)  
       real(dbl):: f(3),fld(3),mr(3)
       integer(i4b):: i,j,k

       ! Equation ....  paper ....           
       fr_t=zero
       do i=1,nsph
         ! SP 07/07/17: Determine previous field component along spheroid axis:
         mr_t(:,i)=mr_tp(:,i)+f1*dmr_tp(:,i)+f2*fmr_tp(:,i)
         ! Determine total field acting on each dipole from:
         ! other particles
         fld=zero
         do k=1,nsph
           if (k.ne.i) then
             call do_field_from_dip(mr_tp(:,k),sph_centre(:,k),f, &
                                              sph_centre(:,i))
           fld=fld+f
           endif 
         enddo
         ! and the molecule
         call do_field_from_dip(mu_tp,mol_cc,f,sph_centre(:,i))
         fld=fld+f
         fmr_t=zero
         do j=1,3 
           mr(:)=dot_product(mr_t(:,i),sph_vrs(:,j,i))*sph_vrs(:,j,i) 
           fmr_t(:,i)=fmr_t(:,i)-MPL_Fw(j,i)*mr(:) &
                                +matmul(MPL_Ff(:,:,j,i),fld)
         enddo
         ! propagate dipole using vv
         dmr_t(:,i)=f3*dmr_tp(:,i)+f4*(fmr_t(:,i)+fmr_tp(:,i))- &
                                                f5*fmr_tp(:,i)
         call do_field_from_dip(mr_t(:,i),sph_centre(:,i),f,mol_cc)
         fr_t=fr_t+f
      enddo
       fmr_tp=fmr_t
       dmr_tp=dmr_t
       mr_tp=mr_t
       ! Local Field
       if(Floc.eq."loc") then
         ! Equation ....  paper ....           
         fx_t=zero
         do i=1,nsph
           ! SP 07/07/17: Determine previous field component along spheroid axis:
           mx_t(:,i)=mx_tp(:,i)+f1*dmx_tp(:,i)+f2*fmx_tp(:,i)
           fld=zero
           do k=1,nsph
             if (k.ne.i) then
               call do_field_from_dip(mx_tp(:,k),sph_centre(:,k),f, &
                                                sph_centre(:,i))
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
           call do_field_from_dip(mx_t(:,i),sph_centre(:,i),f,mol_cc)
           fx_t=fx_t+f
         enddo
         fmx_tp=fmx_t
         dmx_tp=dmx_t
         mx_tp=mx_t
       endif

       return

      end subroutine


      subroutine prop_dip_sphe_drl(f_tp)
!------------------------------------------------------------------------
! @brief  Drude-Lorentz propagation of the sphere reaction and local fields  
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(IN):: f_tp(3)  
       real(dbl):: f(3),fld(3)
       integer(i4b):: i,j,k

       ! Equation ....  paper ....           
       fr_t=zero
       do i=1,nsph
         mr_t(:,i)=mr_tp(:,i)+f1*dmr_tp(:,i)+f2*fmr_tp(:,i)
         fld=zero
         do k=1,nsph
           if (k.ne.i) then
             call do_field_from_dip(mr_tp(:,k),sph_centre(:,k),f, &
                                              sph_centre(:,i))
           fld=fld+f
           endif 
         enddo
         call do_field_from_dip(mu_tp,mol_cc,f,sph_centre(:,i))
         fld=fld+f
         fmr_t(:,i)=-ONS_fw*mr_t(:,i)+ONS_ff(i)*fld
         dmr_t(:,i)=f3*dmr_tp(:,i)+f4*(fmr_t(:,i)+fmr_tp(:,i))- &
                                                f5*fmr_tp(:,i)
         call do_field_from_dip(mr_t(:,i),sph_centre(:,i),f,mol_cc)
         fr_t=fr_t+f
       enddo
       fmr_tp=fmr_t
       dmr_tp=dmr_t
       mr_tp=mr_t
       if(Floc.eq."loc") then
         fx_t=zero
         do i=1,nsph
           mx_t(:,i)=mx_tp(:,i)+f1*dmx_tp(:,i)+f2*fmx_tp(:,i)
           fld=zero
           do k=1,nsph
             if (k.ne.i) then
               call do_field_from_dip(mx_tp(:,k),sph_centre(:,k),f, &
                                                sph_centre(:,i))
             fld=fld+f
             endif 
           enddo
           ! and the external field
           fld=fld+f_tp
           fmx_t=zero
           fmx_t(:,i)=-ONS_fw*mx_t(:,i)+ONS_ff(i)*fld
           dmx_t(:,i)=f3*dmx_tp(:,i)+f4*(fmx_t(:,i)+fmx_tp(:,i))-& 
                                                 f5*fmx_tp(:,i)
           call do_field_from_dip(mx_t(:,i),sph_centre(:,i),f,mol_cc)
           fx_t=fx_t+f
         enddo
         fmx_tp=fmx_t
         dmx_tp=dmx_t
         mx_tp=mx_t
       endif

       return

      end subroutine


      subroutine get_ons(f_med)
!------------------------------------------------------------------------
! @brief Fetch the current value of the onsager field from elsewhere   
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl):: f_med(3)

       f_med=fr_t

       return

      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Free-energy   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine do_gneq(c,mu_or_v,df_or_dq,f_or_q,f_or_q0,fact_d, &
                           n_coor_or_ts,sig)
!------------------------------------------------------------------------
! @brief Update non-equilibrium free energy   
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

        implicit none

        complex(cmp), intent(in) :: c(n_ci)
        integer(i4b),intent(in) :: n_coor_or_ts,sig
        real(dbl),intent(in) :: mu_or_v(n_coor_or_ts,n_ci,n_ci)
        real(dbl),intent(in) :: df_or_dq(n_coor_or_ts), &
                               f_or_q(n_coor_or_ts),&
                               f_or_q0(n_coor_or_ts),&
                               fact_d(n_coor_or_ts,n_coor_or_ts)
        real(dbl),allocatable :: v_avg(:)
!       real(dbl) :: de_a
        integer(i4b) :: its

        g_eq=0.d0
!       de_a=0.d0
        allocate(v_avg(n_coor_or_ts))
        do its=1,n_coor_or_ts
          v_avg(its)=dot_product(c,matmul(mu_or_v(its,:,:),c))
          g_neq1_part=g_neq1_part+sig*v_avg(its)*df_or_dq(its)
          g_eq=g_eq+sig*f_or_q(its)*v_avg(its)
!         de_a=de_a+sig*f_or_q0(its)*v_avg(its)
        enddo
        g_eq=0.5*g_eq
! SC 27/09/2016: corrected bug in expression of e_vac
!       e_vac=2.*g_eq_gs-de_a
        e_vac=-2.*g_eq+2.*g_eq_gs
!       g_neq1=g_neq_0+g_neq1_part-g_eq
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
        deallocate(v_avg)

        return

      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Test/Debug    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine do_ref(c)      
!------------------------------------------------------------------------
! @brief Test for Solvent/Nanoparticle reaction/local fields    
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       complex(cmp), intent(IN) :: c(n_ci)
       complex(cmp) :: refc, eps, E0
       integer(i4b) :: its  
       real(dbl):: dist,pos(3),dp,emol(3),rr,facd,fac0,tau

       select case (Ftest)
         ! Spherical Nanoparticle reaction field
         case ("n-r")
           call do_dip_from_coeff(c,dip,n_ci)
           if(Fprop.eq."dip") then
             pos(1)=sph_centre(1,1)-mol_cc(1)   
             pos(2)=sph_centre(2,1)-mol_cc(2) 
             pos(3)=sph_centre(3,1)-mol_cc(3)
             rr=sph_maj(1)
           else
             pos(1)=sfe_act(1)%x-mol_cc(1)   
             pos(2)=sfe_act(1)%y-mol_cc(2) 
             pos(3)=sfe_act(1)%z-mol_cc(3)
             rr=cts_act(1)%rsfe
           endif
           dist=sqrt(dot_product(pos,pos))
           dp=dot_product(dip,pos)
           emol(:)=(3*dp*pos(:)/dist**2-dip(:))/dist**3
           ref=emol(3)*rr**3
         ! Spherical Nanoparticle local field
         case ("n-l")
           if(Fprop.eq."dip") then
             rr=sph_maj(1)
           else
             rr=cts_act(1)%rsfe 
           endif
           call do_eps_drl
           !refc=eps_f*ui*exp(-ui*omega*t)
           E0=dcmplx(zero,0.5d0*sqrt(dot_product(fmax(:,1),fmax(:,1))))
           refc=eps_f*E0*exp(-ui*omega(1)*t)
           ref=real(refc+conjg(refc))*rr**3
           !write(*,*) eps_f,refc,ref, cts_act(1)%rsfe 
           !stop
         case ("s-r")
         !Corni JCPA 2015, fig.1
           if(Fprop.eq."dip") then
             rr=sph_maj(1)
           else
             rr=sqrt(cts_act(1)%x**2+cts_act(1)%y**2+cts_act(1)%z**2) 
           endif
           facd=(two*eps_d-two)/(two*eps_d+one)
           fac0=two*three*(eps_0-eps_d)/(two*eps_d-two)/(two*eps_0+one)
           tau=(two*eps_d+one)/(two*eps_0+one)*tau_deb
           ref=facd*(one-fac0*(exp(-t/tau)-one))/rr**3*mu_tp(3)
         case ("s-l")
           if(Fprop.eq."dip") then
             rr=sph_maj(1)
           else
             rr=sqrt(cts_act(1)%x**2+cts_act(1)%y**2+cts_act(1)%z**2) 
           endif
           call do_eps_deb
           !refc=eps_f*ui*exp(-ui*omega*t)
           E0=dcmplx(zero,0.5d0*sqrt(dot_product(fmax(:,1),fmax(:,1))))
           refc=eps_f*E0*exp(-ui*omega(1)*t)
           ref=real(refc+conjg(refc))
       end select

       return

      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Input/Output  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine out_mdm(i)
!------------------------------------------------------------------------
! @brief Output medium reaction/local field/dipoles    
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b),intent(IN) :: i  
       real(dbl):: fm(3)

       select case(Ftest)
       case ('n-r','n-l')
         write (file_med,'(i8,f12.2,4e22.10)') i,t,mu_mdm(:,1),ref
       case ('s-r')
         if(Fprop(1:3).eq."chr") call do_field_from_charges(qr_t,fr_t)
         write (file_med,'(i8,f12.2,4e22.10)') i,t,fr_t(:),ref
       case ('s-l')
         write (file_med,'(i8,f12.2,4e22.10)') i,t,fx_t(:),ref
       case default
         if(Fprop(1:3).eq."dip") then
           write (file_med,'(i8,f12.2,3e22.10)') i,t,fr_t(:)
         else
           call do_field_from_charges(qr_t,fr_t)
           write (file_med,'(i8,f12.2,9e22.10)')i,t,mu_mdm(:,1),&
                                                    fr_t(:),qtot,qtot0
         endif
       end select
      
       return

      end subroutine


      subroutine read_charges_gau
!------------------------------------------------------------------------
! @brief Read charges from gaussian "charges0.inp" output    
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

        integer(4) :: i,nts

        open(7,file="charges0.inp",status="old",err=10)
        write (6,*) "Initial charges read from charges0.inp"
          read(7,*) nts
          if(nts_act.eq.0.or.nts.eq.nts_act) then
            nts_act=nts
          else
            write(*,*) "Tesserae number conflict"
            stop
          endif
          qtot0=zero
          do i=1,nts_act 
            read(7,*) q0(i)
            qtot0=qtot0+q0(i)
          enddo
        close(7)

        return

10      write (6,*) "No file charges0.inp found,", &
                    "charges calculated for the initial state"
        return

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
      subroutine init_BEM_Q0
!------------------------------------------------------------------------
! @brief Create a new BEM_Q0=BEM_Qw^-1*BEM_Qf that should avoid
!               spurious charge dynamics for stationary states    
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       implicit none

       real(dbl),allocatable :: mat_t(:,:),mat_t2(:,:)
       integer(i4b) :: info,its
       integer(i4b),allocatable :: ipiv(:)

       allocate(mat_t(nts_act,nts_act))
       allocate(mat_t2(nts_act,nts_act))
       allocate(ipiv(nts_act))
       mat_t=BEM_Qw
       mat_t2=BEM_Qf
       call dgetrf(nts_act,nts_act,mat_t,nts_act,ipiv,info)
       call dgetrs('N',nts_act,nts_act,mat_t,nts_act,ipiv,mat_t2,nts_act,info)
       BEM_Q0=mat_t2
       deallocate (mat_t)
       deallocate (mat_t2)
       deallocate (ipiv)

       return

      end subroutine


      subroutine wrt_restart_mdm()
!------------------------------------------------------------------------
! @brief write restart 
!
! @date Created   : E. Coccia 28 Nov 2017
! Modified  :
!------------------------------------------------------------------------

       implicit none

       integer(i4b)     :: i

       open(778, file='restart_mdm')

      ! if (Fint.eq.'ons') then
          write(778,*) 'Dipoles for Onsager' 
          do i=1,3
             write(778,*) fr_t(1), fr_t(2), fr_t(3)
          enddo
          if (Floc.eq.'loc') then
             write(778,*) 'Dipoles for Onsager (local)' 
             do i=1,3
                write(778,*) fx_t(1), fx_t(2), fx_t(3)
             enddo
          endif
       !elseif (Fint.eq.'pcm') then
          write(778,*) 'Charges for PCM' 
          do i=1,nts_act
             write(778,*) qr_t(i)
          enddo
          if (Floc.eq.'loc') then
             write(778,*) 'Charges for PCM (local)' 
             do i=1,nts_act
                write(778,*) qx_t(i)
             enddo
          endif
       !endif

       close(778)

       return
 
      end subroutine wrt_restart_mdm

      end module
