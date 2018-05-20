      Module BEM_medium      
      use constants    
      use readio_medium
      use global_tdplas
      use pedra_friends
      use MathTools
#ifdef OMP
      use omp_lib
#endif
!      use, intrinsic :: iso_c_binding

      implicit none
! SP 25/06/17: variables starting with "BEM_" are public. 
      real(dbl), allocatable :: BEM_S(:,:),BEM_D(:,:)    !< Calderon S and D matrices 
      real(dbl), allocatable :: BEM_Sm1(:,:)             !< Onsager $S^{-1}$ matrix
      real(dbl), allocatable :: BEM_L(:),BEM_T(:,:)      !< $\Lambda$ and T eigenMatrices
      real(dbl), allocatable :: BEM_W2(:),BEM_Modes(:,:) !< BEM squared frequencies and Modes ($T*S^{1/2}$)
! SP 25/06/17: K0 and Kd are still common to 'deb' and 'drl' cases
      real(dbl), allocatable :: K0(:),Kd(:)              !< Diagonal $K_0$ and $K_d$ matrices
      real(dbl), allocatable :: K0x(:),Kdx(:)              !< Diagonal $K_0$ and $K_d$ matrices for local (x) field
      real(dbl), allocatable :: fact1(:),fact2(:)        !< Diagonal vectors for propagation matrices
      real(dbl), allocatable :: fact2x(:)                !< Diagonal vectors for propagation matrices for local (x) field
      real(dbl), allocatable :: Sm12T(:,:),TSm12(:,:)    !< $S^{-1/2}T$ and $T*S^{-1/2}$ matrices
      real(dbl), allocatable :: TSp12(:,:)               !< $T*S^{1/2}$ matrices
      real(dbl), allocatable :: BEM_Sm12(:,:),Sp12(:,:)  !< $S^{-1/2}$ and $S^{1/2}$ matrices
      real(dbl), allocatable :: BEM_Q0(:,:),BEM_Qd(:,:)  !< Static and Dyanamic BEM matrices $Q_0$ and $Q_d$
      real(dbl), allocatable :: BEM_Qt(:,:),BEM_R(:,:)   !< Debye propagation matrices $\tilde{Q}$ and $R$
      real(dbl), allocatable :: BEM_Qw(:,:),BEM_Qf(:,:)  !< Drude-Lorents propagation matrices $Q_\omega$ and $Q_f$
      real(dbl), allocatable :: BEM_Q0x(:,:),BEM_Qdx(:,:)  !< Static and Dyanamic BEM matrices $Q_0$ and $Q_d$ for local (x) field
      real(dbl), allocatable :: BEM_Qtx(:,:)						!< Debye propagation matrix $\tilde{Q}$ for local (x) field
      real(dbl), allocatable :: BEM_Qfx(:,:)						!< Drude-Lorents propagation matrix $Q_f$ for local (x) field
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
      complex(cmp) :: eps,eps_f                          !< drl complex eps(\omega) and (eps(\omega)-1)/(eps(\omega)+2) 
      real(dbl), allocatable :: lambda(:,:)              !< depolarizing factors (3,nsph)       
      complex(cmp), allocatable :: q_omega(:) !< Medium charges in frequency domain
      complex(cmp), allocatable :: Kdiag_omega(:)         !< Diagonal K matrix in frequency domain
      real(dbl), allocatable :: scrd3(:) ! Scratch vector dim=3
      save
      private
      public eps,eps_f,BEM_L,BEM_T,BEM_Sm1,ONS_ff,ONS_fw,              &
             BEM_Sm12,MPL_F0,MPL_Ft0,MPL_Fd,MPL_Fx0,MPL_Ftx0,MPL_Fxd,  &
             MPL_Tauxm1,MPL_Taum1,mat_f0,mat_fd,MPL_Ff,MPL_Fw,         &
             ONS_f0,ONS_fd,ONS_taum1,ONS_fx0,ONS_fxd,ONS_tauxm1,       &
             BEM_Qt,BEM_R,BEM_Qw,BEM_Qf,BEM_Qd,BEM_Q0,BEM_W2,BEM_Modes,&
             BEM_Qtx,BEM_Qfx,BEM_Qdx,BEM_Q0x,                          &
             do_BEM_prop,do_BEM_freq,do_BEM_quant,do_MPL_prop,         &
             do_eps_drl,do_eps_deb,do_charge_freq,                     &
             deallocate_BEM_public,deallocate_MPL_public

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  DRIVER  ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      subroutine do_BEM_prop     
!------------------------------------------------------------------------
! @brief BEM driver routine for propagation
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       !Cavity read/write and S D matrices 
       call init_BEM
       if(Fgamess.eq.'yes') then
         allocate(BEM_Qd(nts_act,nts_act))
         allocate(BEM_Q0(nts_act,nts_act))
         !Standard or Diagonal BEM           
         if(Fbem(1:4).eq.'stan') then
           write(6,*) "Standard BEM not implemented yet"
           stop
         elseif(Fbem(1:4).eq.'diag') then
           call init_BEM_diagonal
           call do_BEM_diagonal
         endif
         !Write out matrices for gamess                     
         call out_BEM_gamess
         call finalize_BEM
         stop
       endif
       if(Fprop(1:3).eq."chr") then
         allocate(BEM_Qd(nts_act,nts_act))
         allocate(BEM_Q0(nts_act,nts_act))
       endif
       if(Fprop(1:6).eq."chr-ie") then
         if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') allocate(BEM_Q0x(nts_act,nts_act),BEM_Qdx(nts_act,nts_act))
         !Standard or Diagonal BEM           
         if(Fbem(1:4).eq.'stan') then
           write(6,*) "Standard BEM not implemented yet"
           stop
         elseif(Fbem(1:4).eq.'diag') then
           call init_BEM_diagonal
           call do_BEM_diagonal
           !Save Modes for quantum BEM         
           if(Fmdm(1:1).eq."Q") then
             allocate(BEM_Modes(nts_act,nts_act))
             BEM_Modes=TSm12
           endif
         endif
         !Write out matrices for gamess                     
         if(Fwrite.eq."high") call out_BEM_gamess
         if(Fwrite.eq."high") call out_BEM_mat
         !Build propagation Matrices 
         if(Feps.eq."deb") then
           allocate(BEM_R(nts_act,nts_act))
           allocate(BEM_Qt(nts_act,nts_act))
           if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') allocate(BEM_Qtx(nts_act,nts_act))
           if(Fbem(1:4).eq.'stan') call do_propBEM_dia_deb ! 'dia' to be replace by 'std' 
           if(Fbem(1:4).eq.'diag') call do_propBEM_dia_deb
         elseif(Feps.eq."drl") then
           allocate(BEM_Qw(nts_act,nts_act))
           allocate(BEM_Qf(nts_act,nts_act))
           if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') allocate(BEM_Qfx(nts_act,nts_act))
           if(Fbem(1:4).eq.'stan') call do_propBEM_dia_drl ! 'dia' to be replace by 'std' 
           if(Fbem(1:4).eq.'diag') call do_propBEM_dia_drl
         endif
         !Write out propagation matrices         
         if(Fwrite.eq."high") call out_BEM_propmat  
       elseif(Fprop(1:7).eq."chr-ons") then
         allocate(BEM_Sm1(nts_act,nts_act))
         ! Form $S^{-1}$ matrix
         BEM_Sm1=inv(BEM_S)
         nsph=1
         if(Feps.eq."deb") then
           call do_propfact_ons_deb
         elseif(Feps.eq."drl") then
           allocate(ONS_ff(nsph))
           call do_propfact_ons_drl
         endif
         BEM_Q0=-ONS_f0*BEM_Sm1
         BEM_Qd=-ONS_fd*BEM_Sm1
       endif
       !Deallocate private arrays              
       call finalize_BEM

       return
 
      end subroutine


      subroutine do_BEM_freq(omega_list,n_omega)     
!------------------------------------------------------------------------
! @brief BEM driver routine for frequency calculation (old do_freq_mat) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(in):: omega_list(:)
       integer(i4b), intent(in):: n_omega
       real(dbl), allocatable :: pot(:)
       complex(cmp) :: mu_omega(3)
       integer(i4b):: i

       ! Cavity read/write and S D matrices 
       call init_BEM
       allocate(BEM_Qd(nts_act,nts_act))
       allocate(BEM_Q0(nts_act,nts_act))
       ! Using Diagonal BEM           
       call init_BEM_diagonal
       call do_BEM_diagonal
       ! Write out matrices                     
       call out_BEM_gamess
       ! Calculate potential on tesserae
       allocate(pot(nts_act))
       call do_pot_from_field(fmax(:,1),pot)
       allocate(Kdiag_omega(nts_act))
       allocate(q_omega(nts_act))
       open(7,file="dipole_freq.dat",status="unknown")
       write (7,*)"freq re(mux) re(muy) re(muz) im(mux) im(muy) im(muz)"
       do i=1,n_omega
        call do_charge_freq(omega_list(i),pot,mu_omega)
        write (7,'(7e15.6)') omega_list(i),real(mu_omega(:)),aimag(mu_omega(:))
       enddo
       close(7)
       deallocate(pot,q_omega,Kdiag_omega)
       !Deallocate private arrays              
       call finalize_BEM

       return
 
      end subroutine


      subroutine do_BEM_quant     
!------------------------------------------------------------------------
! @brief BEM driver routine for quantum BEM calculation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       ! Cavity read/write and S D matrices 
       call init_BEM
       allocate(BEM_Qd(nts_act,nts_act))
       allocate(BEM_Q0(nts_act,nts_act))
       ! Diagonal BEM           
       call init_BEM_diagonal
       call do_BEM_diagonal
       !Save Modes for quantum BEM         
       allocate(BEM_Modes(nts_act,nts_act))
       BEM_Modes=TSm12
       call finalize_BEM

       return
 
      end subroutine


      subroutine init_BEM     
!------------------------------------------------------------------------
! @brief BEM initialization routine: cavity read/write and S and D
! matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
     
       integer(i4b) :: its                    
      
       allocate(scrd3(3))
       sgn=one                 
       if(Fmdm(2:4).eq."nan") sgn=-one  
       if (FinitBEM.eq.'wri') then
       ! Write out geometric info and stop
         ! Build the cavity/nanoparticle surface
         if(Fsurf.eq.'fil') then 
           call read_cavity_full_file
         elseif(Fsurf.eq.'gms') then
           call read_gmsh_file
         else
           if(Fmdm(2:4).eq.'sol') call pedra_int('act')
           if(Fmdm(2:4).eq.'nan') call pedra_int('met')
         endif
         ! write out the cavity/nanoparticle surface
         call output_surf
         write(6,*) "Created output file with surface points"
         ! Build and write out Calderon SD matrices
         allocate(BEM_S(nts_act,nts_act))
         if (Fprop(1:7).ne.'chr-ons') allocate(BEM_D(nts_act,nts_act))
         call do_BEM_SD
         call write_BEM_SD
         write(6,*) "Matrixes S D have been written out"
         stop
       elseif (FinitBEM.eq.'rea') then
       !Read in geometric info and proceed
         call read_cavity_file
         allocate(BEM_S(nts_act,nts_act))
         if (Fprop(1:7).ne.'chr-ons') allocate(BEM_D(nts_act,nts_act))
         call read_BEM_SD     
         write(6,*) "BEM surface and Matrixes S D have been read in"
       endif 
       write(6,*) "BEM correctly initialized"

       return
 
      end subroutine


      subroutine finalize_BEM     
!------------------------------------------------------------------------
! @brief BEM finalized and deallocation routine 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       deallocate(scrd3)
       deallocate(BEM_S)
       if(allocated(BEM_D)) deallocate(BEM_D)
       if(allocated(fact1)) deallocate(fact1)
       if(allocated(fact2)) deallocate(fact2)
       if(allocated(K0)) deallocate(K0)
       if(allocated(Kd)) deallocate(Kd)
       if(allocated(fact2x)) deallocate(fact2x)
       if(allocated(K0x)) deallocate(K0x)
       if(allocated(Kdx)) deallocate(Kdx)
       if(allocated(Sp12)) deallocate(Sp12)
       if(allocated(Sm12T)) deallocate(Sm12T)
       if(allocated(TSm12)) deallocate(TSm12)
       if(allocated(TSp12)) deallocate(TSp12)

       return
 
      end subroutine


      subroutine deallocate_BEM_public     
!------------------------------------------------------------------------
! @brief BEM finalized and deallocation routine 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       if (Fprop(1:3).eq.'chr') then
         if(allocated(BEM_Qd)) deallocate(BEM_Qd)
         if(allocated(BEM_Q0)) deallocate(BEM_Q0)
         if(allocated(BEM_Qt)) deallocate(BEM_Qt)
         if(allocated(BEM_R)) deallocate(BEM_R)
         if(allocated(BEM_Qw)) deallocate(BEM_Qw)
         if(allocated(BEM_Qf)) deallocate(BEM_Qf)
         if(allocated(BEM_Qdx)) deallocate(BEM_Qdx)
         if(allocated(BEM_Q0x)) deallocate(BEM_Q0x)
         if(allocated(BEM_Qtx)) deallocate(BEM_Qtx)
         if(allocated(BEM_Qf)) deallocate(BEM_Qfx)
         if(allocated(BEM_L)) deallocate(BEM_L)
         if(allocated(BEM_W2)) deallocate(BEM_W2)
         if(allocated(BEM_T)) deallocate(BEM_T)
         if(allocated(BEM_Sm1)) deallocate(BEM_Sm1)
         if(allocated(BEM_Sm12)) deallocate(BEM_Sm12)
       endif

       return
 
      end subroutine


      subroutine do_MPL_prop         
!------------------------------------------------------------------------
! @brief Calculate propagation Onsager matrices from factors including
! depolarization 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl):: tmp(3),m
       integer(i4b):: i,j

       ! SP 05/07/17 Only one cavity!!! 
       if(Fmdm(2:4).eq."sol") nsph=1 
       call init_MPL
       if(MPL_ord.eq.1) then
       !SPHEROID
         if(Fshape.eq."spho") then   
           do i=1,nsph
           ! Determine minor axis versors for spheroids                          
             ! build a vector not parallel to major unit vector
             m=sph_maj(i)/sph_min(i)
             tmp=zero
             if(sph_vrs(1,1,i).ne.one) tmp(1)=one
             if(sph_vrs(1,1,i).eq.one) tmp(2)=one
             ! build the two minor unit vectors
             sph_vrs(:,2,i)=vprod(sph_vrs(:,1,i),tmp) 
             sph_vrs(:,2,i)=sph_vrs(:,2,i)/mdl(sph_vrs(:,2,i))
             sph_vrs(:,3,i)=vprod(sph_vrs(:,1,i),sph_vrs(:,2,i)) 
             sph_vrs(:,3,i)=sph_vrs(:,3,i)/mdl(sph_vrs(:,3,i))
            ! Determine depolarization factors: Osborn Phys. Rev. 67 (1945), 351.
             if(int(m*100).gt.100) then   ! Prolate spheroid
               lambda(1,i)=-4.d0*pi/(m*m-one)*(m/two/sqrt(m*m-one)* &
                     log((m+sqrt(m*m-one))/(m-sqrt(m*m-one)))-one)
             elseif(int(m*100).lt.100) then    ! Oblate spheroid
               m=1/m
               lambda(1,i)=-4.d0*pi*m*m/(m*m-one)*&
                       (one-one/sqrt(m*m-one)*asin(sqrt(m*m-one)/m))
             else 
               write(6,*) "This is a Sphere"
               lambda(1,i)=one/three
             endif
             lambda(2,i)=pt5*(one-lambda(1,i))
             lambda(3,i)=pt5*(one-lambda(1,i))
           enddo
           if(Feps.eq."deb") call do_propMPL_deb
           if(Feps.eq."drl") call do_propMPL_drl
           ! mat_f0 and mat_fd for scf, g_ and prop
           mat_f0=zero
           mat_fd=zero
           do j=1,nsph
             do i=1,3 
               mat_f0(:,:)=mat_f0(:,:)+MPL_F0(:,:,i,j)
               mat_fd(:,:)=mat_fd(:,:)+MPL_Fd(:,:,i,j)
             enddo
           enddo
         else                         
       !SPHERE  
         !Build propagation Matrices with factors  
           ! mat_f0 and mat_fd for scf, g_ and prop
           if(Feps.eq."drl")then 
             call do_propfact_ons_drl 
             do i=1,nsph
               ONS_ff=ONS_ff*sph_maj(i)**3
               ONS_f0=ONS_f0*sph_maj(i)**3
               ONS_fd=ONS_fd*sph_maj(i)**3
             enddo
           elseif(Feps.eq."deb") then 
             call do_propfact_ons_deb 
             do i=1,nsph
               ONS_f0=ONS_f0/sph_maj(i)**3
               ONS_fd=ONS_fd/sph_maj(i)**3
             enddo
           endif
           mat_fd=zero
           mat_f0=zero
           do i=1,3
             mat_fd(i,i)=ONS_fd
             mat_f0(i,i)=ONS_f0
           enddo
           write (6,*) "Onsager"
           write (6,*) "eps_0,eps_d",eps_0,eps_d
           write (6,*) "lambda",1/three     
           write (6,*) "f0",ONS_f0
           write (6,*) "fd",ONS_fd
           if(Feps.eq."deb")write (6,*) "tau",1./ONS_taum1
         endif
       else
         write(6,*) "Higher multipoles not implemented "
         stop
       endif
       call finalize_MPL
 
       return
  
      end subroutine


      subroutine init_MPL     
!------------------------------------------------------------------------
! @brief Allocate arrays for multipolar routines 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       allocate(mat_fd(3,3))
       allocate(mat_f0(3,3))
       if(Fshape.eq."spho") then
         allocate(lambda(3,nsph))
         allocate(MPL_F0(3,3,3,nsph))
         allocate(MPL_Fd(3,3,3,nsph))
         if(Feps.eq."drl") then
           allocate(MPL_Fw(3,nsph))
           allocate(MPL_Ff(3,3,3,nsph))
         elseif(Feps.eq."deb") then
           allocate(MPL_Ft0(3,3,3))
           allocate(MPL_Taum1(3))
           if(Floc.eq.'loc') then
             allocate(MPL_Fx0(3,3,3))
             allocate(MPL_Fxd(3,3,3))
             allocate(MPL_Ftx0(3,3,3))
             allocate(MPL_Tauxm1(3))
           endif
         endif
       else 
         allocate(ONS_ff(nsph))
       endif

       return
 
      end subroutine


      subroutine finalize_MPL     
!------------------------------------------------------------------------
! @brief Deallocate lambda 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       if(allocated(lambda)) deallocate(lambda)

       return
 
      end subroutine


      subroutine deallocate_MPL_public     
!------------------------------------------------------------------------
! @brief Deallocate MPL arrays 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

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
 
      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  CORE ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine do_BEM_SD
!------------------------------------------------------------------------
! @brief Compute Calderon's S and D matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl) :: temp
       integer(i4b) :: i,j

#ifdef OMP 
!$OMP PARALLEL 
!$OMP DO 
#endif
       do i=1,nts_act
        do j=1,nts_act
          call green_s(i,j,temp)
          BEM_S(i,j)=temp
          if (Fprop(1:7).ne.'chr-ons') then 
            call green_d(i,j,temp)
            BEM_D(i,j)=temp
          endif
        enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       return
 
      end subroutine


      subroutine green_d (i,j,value)
!------------------------------------------------------------------------
! @brief Calderon D matrix with Purisima Dii elements 
! SC: changed to a diagonal value of D_ii that should be more general
! than those for the sphere
!
! @date Created: S. Pipolo
! Modified: S. Corni 30/5/17
!------------------------------------------------------------------------

       integer(i4b), intent(in):: i,j
       real(dbl), intent(out) :: value
       real(dbl):: dist,sum_d
       integer(i4b) :: k

       if (i.ne.j) then
          scrd3(1)=(cts_act(i)%x-cts_act(j)%x)
          scrd3(2)=(cts_act(i)%y-cts_act(j)%y)
          scrd3(3)=(cts_act(i)%z-cts_act(j)%z)
          dist=sqrt(dot_product(scrd3,scrd3))
          value=dot_product(cts_act(j)%n,scrd3)/dist**3 
       else
          sum_d=0.d0
#ifdef OMP
!$OMP PARALLEL REDUCTION(+:sum_d)
!$OMP DO 
#endif
          do k=1,i-1
             scrd3(1)=(cts_act(i)%x-cts_act(k)%x)
             scrd3(2)=(cts_act(i)%y-cts_act(k)%y)
             scrd3(3)=(cts_act(i)%z-cts_act(k)%z)
             dist=sqrt(dot_product(scrd3,scrd3))
             sum_d=sum_d+dot_product(cts_act(k)%n,scrd3)/dist**3*cts_act(k)%area 
          enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:sum_d)
!$OMP DO 
#endif
          do k=i+1,nts_act
             scrd3(1)=(cts_act(i)%x-cts_act(k)%x)
             scrd3(2)=(cts_act(i)%y-cts_act(k)%y)
             scrd3(3)=(cts_act(i)%z-cts_act(k)%z)
             dist=sqrt(dot_product(scrd3,scrd3))
             sum_d=sum_d+dot_product(cts_act(k)%n,scrd3)/dist**3*cts_act(k)%area 
          enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
          sum_d=-(2.0*pi+sum_d)/cts_act(i)%area
          value=sum_d
          !value=-1.0694*sqrt(4.d0*pi*cts_act(i)%area)/(2.d0* &
          !       cts_act(i)%rsfe)/cts_act(i)%area
       endif

       return

      end subroutine


      subroutine green_s (i,j,value)
!------------------------------------------------------------------------
! @brief Calderon S matrix 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b), intent(in):: i,j
       real(dbl), intent(out) :: value
       real(dbl):: dist

       if (i.ne.j) then
         scrd3(1)=(cts_act(i)%x-cts_act(j)%x)
         scrd3(2)=(cts_act(i)%y-cts_act(j)%y)
         scrd3(3)=(cts_act(i)%z-cts_act(j)%z)
         dist=sqrt(dot_product(scrd3,scrd3))
         value=one/dist 
       else
         value=1.0694*sqrt(4.d0*pi/cts_act(i)%area)
       endif

       return

      end subroutine


      subroutine do_BEM_diagonal
!------------------------------------------------------------------------
! @brief Compute BEM matrices within diagonal approach 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       integer(i4b) :: i,j
       real(8), allocatable :: scr1(:,:),scr2(:,:),scr3(:,:)
       real(8), allocatable :: eigt(:,:),eigt_t(:,:)
       real(8), allocatable :: eigv(:)
       real(dbl) :: fac_eps0,fac_epsd

       allocate(scr1(nts_act,nts_act),scr2(nts_act,nts_act))
       allocate(scr3(nts_act,nts_act))
       allocate(eigv(nts_act))
       allocate(eigt(nts_act,nts_act),eigt_t(nts_act,nts_act))
       ! Form S^1/2 and S^-1/2
       ! Copy the matrix in the eigenvector matrix
       eigt = BEM_S
       call diag_mat(eigt,eigv,nts_act)
       if(Fwrite.eq."high")write(6,*) "S matrix diagonalized "
#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
          if(eigv(i).le.0.d0) then
            write(6,*) "WARNING:",i," eig of S is negative or zero!"
            write(6,*) "   I put it to 1e-8"
            eigv(i)=1.d-8
          endif
          scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
       enddo
#ifdef OMP
!$OMP ENDDO 
!$OMP END PARALLEL
#endif
       eigt_t=transpose(eigt)
       Sp12=matmul(scr1,eigt_t)                   

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
       enddo
#ifdef OMP
!$OMP ENDDO 
!$OMP END PARALLEL
#endif

       BEM_Sm12=matmul(scr1,eigt_t)                   
!      Form the S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 , and diagonalize it
       !S^-1/2 D A S^1/2

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=BEM_D(:,i)*cts_act(i)%area
       enddo
#ifdef OMP
!$OMP ENDDO 
!$OMP END PARALLEL
#endif

       scr3=matmul(BEM_Sm12,scr1)                   
       scr2=matmul(scr3,Sp12)                   
       !S^-1/2 D A S^1/2+S^1/2 A D* S^-1/2 and diagonalize

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do j=1,nts_act
        do i=1,nts_act
         BEM_T(i,j)=0.5*(scr2(i,j)+scr2(j,i))
        enddo
       enddo
#ifdef OMP
!$OMP ENDDO 
!$OMP END PARALLEL
#endif

       deallocate(scr2,scr3)
       call diag_mat(BEM_T,BEM_L,nts_act)
       if(Fwrite.eq."high")&
       write(6,*) "S^-1/2DAS^1/2+S^1/2AD*S^-1/2 matrix diagonalized"
       if (Feps.eq."deb") then
!       debye dielectric function  
         if(eps_0.ne.one) then
           fac_eps0=(eps_0+one)/(eps_0-one)
           K0(:)=(twp-sgn*BEM_L(:))/(twp*fac_eps0-sgn*BEM_L(:))
           ! GG: analogous to K_0 matrix in the case of local-field for solvent external medium
           if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') K0x(:)=-(twp+BEM_L(:))/(twp*fac_eps0-BEM_L(:))
         else
           K0(:)=zero
         endif
         if(eps_d.ne.one) then
           fac_epsd=(eps_d+one)/(eps_d-one)
           Kd(:)=(twp-sgn*BEM_L(:))/(twp*fac_epsd-sgn*BEM_L(:))
           ! GG: analogous to K_d matrix in the case of local-field for solvent external medium
           if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') Kdx(:)=-(twp+BEM_L(:))/(twp*fac_epsd-BEM_L(:))
         else
           Kd=zero
         endif
         ! SP: Need to check the signs of the second part for a debye medium localized in space  
         fact1(:)=((twp-sgn*BEM_L(:))*eps_0+twp+BEM_L(:))/ &
                 ((twp-sgn*BEM_L(:))*eps_d+twp+BEM_L(:))/tau_deb
         fact2(:)=K0(:)*fact1(:)
         ! GG: analogous to \tau K_0 matrix in the case of local-field for solvent external medium
         if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') fact2x(:)=K0x(:)*fact1(:)
       elseif (Feps.eq."drl") then       
!        Drude-Lorentz dielectric function
         Kd=zero 
         fact2(:)=(twp-sgn*BEM_L(:))*eps_A/(two*twp)  
! SC: the first eigenvector should be 0 for the NP
         if (Fmdm(2:4).eq.'nan') fact2(1)=0.d0
         ! SC: no spurious negative square frequencies


#ifdef OMP
!$OMP PARALLEL 
!$OMP DO 
#endif
         do i=1,nts_act
           if(fact2(i).lt.0.d0) then
             write(6,*) "WARNING: BEM_W2(",i,") is ", fact2(i)+eps_w0*eps_w0
             write(6,*) "   I put it to 1e-8"
             fact2(i)=1.d-8
             BEM_L(i)=-twp
           endif
         enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
         if (eps_w0.eq.zero) eps_w0=1.d-8
         BEM_W2(:)=fact2(:)+eps_w0*eps_w0  
         K0(:)=fact2(:)/BEM_W2(:)
         ! GG: analogous to K_f and K_0 matrices in the case of local-field for solvent external medium
         if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') then
          fact2x(:)=-(twp+BEM_L(:))*eps_A/(two*twp)
          K0x(:)=fact2x(:)/BEM_W2(:)
         endif
       endif
       if(Fwrite.eq."high") write(6,*) "Done BEM eigenmodes"
       Sm12T=matmul(BEM_Sm12,BEM_T)
       TSm12=transpose(Sm12T)
       TSp12=matmul(transpose(BEM_T),Sp12)
      ! SC 05/11/2016 write out the transition charges in pqr format
       if(Fwrite.eq."high") call output_charge_pqr
      ! Do BEM_Q0 and and BEM_Qd 

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=Sm12T(:,i)*K0(i) 
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       BEM_Q0=-matmul(scr1,TSm12) 

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=Sm12T(:,i)*Kd(i) 
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       BEM_Qd=-matmul(scr1,TSm12) 
       ! GG: analogous to Q_0 and Q_d matrices in the case of local-field for solvent external medium
       if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') then
        do i=1,nts_act
         scr1(:,i)=Sm12T(:,i)*K0x(i)
        enddo
        BEM_Q0x=-matmul(scr1,TSm12)
        do i=1,nts_act
          scr1(:,i)=Sm12T(:,i)*Kdx(i)
        enddo
        BEM_Qdx=-matmul(scr1,TSm12)
       endif
      ! Print matrices in output
       if(Fwrite.eq."high") call out_BEM_diagmat 
       deallocate(scr1,eigv,eigt,eigt_t)
       write(6,*) "Done BEM diagonal" 

       return
 
      end subroutine
 

      subroutine init_BEM_diagonal
!------------------------------------------------------------------------
! @brief Initialize diagonal BEM 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       allocate(fact1(nts_act),fact2(nts_act))
       allocate(Kd(nts_act),K0(nts_act))
       if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') then
        allocate(fact2x(nts_act))
        allocate(Kdx(nts_act),K0x(nts_act))
       endif
       allocate(BEM_L(nts_act))
       allocate(BEM_W2(nts_act))
       allocate(BEM_T(nts_act,nts_act))
       allocate(BEM_Sm12(nts_act,nts_act))
       allocate(Sp12(nts_act,nts_act))
       allocate(Sm12T(nts_act,nts_act))
       allocate(TSm12(nts_act,nts_act))
       allocate(TSp12(nts_act,nts_act))

       return

      end subroutine


      subroutine do_charge_freq(omega_a,pot,mu_omega)
!------------------------------------------------------------------------
! @brief Compute charges in the frequency domain 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl), intent(in):: omega_a
       real(dbl), intent(in) :: pot(:)
       complex(cmp), intent(out) :: mu_omega(3)
       real(dbl) :: a,b               
       integer(4) :: its
! SP 16/07/17: avoiding allocations loops
       !complex(cmp), allocatable :: q_omega(:),mu_omega(:)
       !complex(cmp), allocatable :: Kdiag_omega(:)

       sgn=-1
! SP 16/07/17: added eps_w0 to Kdiag_omega
       if(eps_w0.eq.zero) eps_w0=1.d-8 
! SP 16/07/17: changed the following to avoid divergence at low omega 
!      do its=1,nts_act
       Kdiag_omega(1)=zero                   

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do its=2,nts_act
        Kdiag_omega(its)=(twp-sgn*BEM_L(its))/ &
!        ((one-omega_a*(omega_a+ui*eps_gm)*2.d0/eps_A)*twp-sgn*BEM_L(its))
         ((one+(eps_w0**2-omega_a*(omega_a+ui*eps_gm))*two/eps_A)*twp-&
           sgn*BEM_L(its))
       !a=twp+BEM_L(its)+two*twp*(eps_w0**2-omega_a**2)/eps_A
       !b=two*twp*eps_gm*omega_a/eps_A
       !if(its.eq.1) write(6,*) a**2, b**2
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       q_omega=matmul(BEM_Sm12,pot)
       q_omega=matmul(transpose(BEM_T),q_omega)
       q_omega=Kdiag_omega*q_omega
       q_omega=matmul(BEM_T,q_omega)
       q_omega=-matmul(BEM_Sm12,q_omega)
       mu_omega=0.d0

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:mu_omega)
!$OMP DO 
#endif
       do its=1,nts_act
        mu_omega(1)=mu_omega(1)+q_omega(its)*(cts_act(its)%x)
        mu_omega(2)=mu_omega(2)+q_omega(its)*(cts_act(its)%y)
        mu_omega(3)=mu_omega(3)+q_omega(its)*(cts_act(its)%z)
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       return

      end subroutine


      subroutine do_propBEM_dia_deb
!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (debye) 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       integer(i4b) :: i
       real(8), allocatable :: scr1(:,:)

       allocate(scr1(nts_act,nts_act))
!      Form the \tilde{Q} and R for debye propagation

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO 
#endif
        do i=1,nts_act
          scr1(:,i)=Sm12T(:,i)*fact1(i) 
        enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

        BEM_R=matmul(scr1,TSp12)

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
        do i=1,nts_act
          scr1(:,i)=Sm12T(:,i)*fact2(i) 
        enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

        BEM_Qt=-matmul(scr1,TSm12)
        ! GG: analogous to \tilde{Q} matrix in the case of local-field for solvent external medium
        if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') then
         do i=1,nts_act
           scr1(:,i)=Sm12T(:,i)*fact2x(i)
         enddo
         BEM_Qtx=-matmul(scr1,TSm12)
        endif

        deallocate(scr1)

        return

      end subroutine


      subroutine do_propBEM_dia_drl
!------------------------------------------------------------------------
! @brief Propagation of matrices for diagonal BEM (drude-lorentz) 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       integer(i4b) :: i
       real(8), allocatable :: scr1(:,:)

       allocate(scr1(nts_act,nts_act))
!      Form the Q_w and Q_f for drude-lorentz propagation

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=Sm12T(:,i)*BEM_W2(i) 
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       BEM_Qw=matmul(scr1,TSp12)

#ifdef OMP
!$OMP PARALLEL 
!$OMP DO
#endif
       do i=1,nts_act
         scr1(:,i)=Sm12T(:,i)*fact2(i) 
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif

       BEM_Qf=-matmul(scr1,TSm12)
       if(Floc.eq.'loc'.and.Fmdm(2:4).eq.'sol') then
        do i=1,nts_act
          scr1(:,i)=Sm12T(:,i)*fact2x(i)
        enddo
        BEM_Qfx=-matmul(scr1,TSm12)
       endif
       deallocate(scr1)

       return

      end subroutine


      subroutine do_propMPL_deb   
!------------------------------------------------------------------------
! @brief Initialize factors for debye dipole propagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       real(dbl):: k1,f1,g(3,3),m
       integer(i4b):: i,j,k,l

       ! build geometric factors due to orientations      
       MPL_Fd=zero  
       MPL_F0=zero  
       MPL_Taum1=zero
       MPL_Ft0=zero
       if(Floc.eq."loc") then
         MPL_Fx0=zero
         MPL_Fxd=zero
         MPL_Tauxm1=zero
         MPL_Ftx0=zero
       endif
       do i=1,nsph
         m=sph_maj(i)/sph_min(i)
         do j=1,3 !Spheroids maj/min
           g=zero
           do l=1,3 !xyz
             do k=1,3 !xyz
               g(k,l)=sph_vrs(k,j,i)*sph_vrs(l,j,i) 
             enddo
           enddo
           ! Add depolarized reaction factors 3l/ab^2*(eps-1)/(eps+l/(1-l))
           ! Onsager JACS 58 (1936), Buckingham Trans.Far.Soc. 49 (1953), Abbott Trans.Far.Soc. 48 (1952), 
           k1=lambda(j,i)/(one-lambda(j,i))
           f1=three*lambda(j,i)/sph_maj(i)/sph_min(i)**2
           MPL_F0(:,:,j,i)=g(:,:)*f1*(eps_0-one)/(eps_0+k1)  
           MPL_Fd(:,:,j,i)=g(:,:)*f1*(eps_d-one)/(eps_d+k1)  
           ! f0/tau for propagation equations
           MPL_Ft0(:,:,j)=g(:,:)*f1*(eps_0-one)/(eps_d+k1)/tau_deb  
           ! SP 29/06/17: changed from tau_deb to 1/tau_deb
           MPL_Taum1(j)=(eps_0+k1)/(eps_d+k1)/tau_deb 
           if(Floc.eq."loc") then
             ! Add depolarized local factors eps/(eps+(1-eps)*l), 3eps/(2eps+1) for spheres
             ! Sihvola Journal of Nanomaterials 2007
             MPL_Fx0(:,:,j)=g(:,:)*eps_0/(eps_0+(one-eps_0)*lambda(j,i))
             MPL_Fxd(:,:,j)=g(:,:)*eps_d/(eps_d+(one-eps_d)*lambda(j,i))
             MPL_Ftx0(:,:,j)=g(:,:)*eps_0/(eps_d+(1-eps_d)*lambda(j,i))&
                                         /tau_deb
             MPL_Tauxm1(j)=(eps_0+(1-eps_0)*lambda(j,i))/ &
                           (eps_d+(1-eps_d)*lambda(j,i))/tau_deb 
           endif
         enddo
       enddo

       return

      end subroutine


      subroutine do_propfact_ons_deb   
!------------------------------------------------------------------------
! @brief Onsager propagation matrix for debye 
!
! @date Created: S. Pipolo
! Modified: G. Gil
!------------------------------------------------------------------------

       ONS_f0=(eps_0-one)/(eps_0+pt5) 
       ONS_fd=(eps_d-one)/(eps_d+pt5) 
       ONS_fx0=three*eps_0/(two*eps_0+one) 
       ONS_fxd=three*eps_d/(two*eps_d+one) 
       ONS_taum1=(two*eps_0+one)/(two*eps_d+one)/tau_deb

       return

      end subroutine


      subroutine do_propMPL_drl   
!------------------------------------------------------------------------
! @brief Initialize factors for drude-lorentz dipole propagation 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

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
       do i=1,nsph
         m=sph_maj(i)/sph_min(i)
         ! SP 4\pi\epsilon_0=1
         !f1=two*two*pi*sph_maj(i)*sph_min(i)**2
         f1=sph_maj(i)*sph_min(i)**2
         do j=1,3 !Spheroids maj/min
           MPL_Fw(j,i)=eps_A*lambda(j,i)+eps_w0*eps_w0
           g=zero
           do l=1,3 !xyz
             do k=1,3 !xyz
               g(k,l)=sph_vrs(k,j,i)*sph_vrs(l,j,i) 
             enddo
           enddo
           MPL_Ff(:,:,j,i)=g(:,:)*f1*eps_A/three  
           MPL_F0(:,:,j,i)=g(:,:)*f1/three*(eps_0-one)/(one+& 
                                     (eps_0-one)*lambda(j,i))  
           MPL_Fd(:,:,j,i)=g(:,:)*f1/three*(eps_d-one)/(one+& 
                                     (eps_d-one)*lambda(j,i))  
         enddo
       enddo

       return

      end subroutine


      subroutine do_propfact_ons_drl  
!------------------------------------------------------------------------
! @brief Onsager propagation matrix for drude-lorentz 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b)::i

       ! Form the ONS_fw=-1/tau+w_0^2 and ONS_ff=A/3                           
       ! SP 4\pi\epsilon_0=1
       ! May be extended to spheres with different \epsilon(\omega)
       ONS_fw=eps_A/three+eps_w0*eps_w0 
       ONS_f0=(eps_0-one)/(eps_0+two) 
       ONS_fd=(eps_d-one)/(eps_d+two) 
       ONS_ff(:)=eps_A/three 

       return

      end subroutine


      subroutine do_eps_drl      
!------------------------------------------------------------------------
! @brief Compute drl cmplx eps(\omega) and (eps(\omega)-1)/(eps(\omega)+2) 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       !eps_gm=eps_gm+f_vel/sfe_act(1)%r
       eps=dcmplx(eps_A,zero)/dcmplx(eps_w0**2-omega(1)**2,-omega(1)*eps_gm)
       eps=eps+onec
       eps_f=(eps-onec)/(eps+twoc)

       return

      end subroutine


      subroutine do_eps_deb      
!------------------------------------------------------------------------
! @brief Compute deb cmplx eps(\omega) and (3*eps(\omega))/(2*eps(\omega)+1)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       !eps_gm=eps_gm+f_vel/sfe_act(1)%r
       eps=dcmplx(eps_d,zero)+dcmplx(eps_0-eps_d,zero)/ &
                              dcmplx(one,-omega(1)*tau_deb)
       eps_f=(three*eps)/(two*eps+onec)

       return
 
      end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   INPUT/OUTPUT ROUTINES   !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine write_BEM_SD
!------------------------------------------------------------------------
! @brief Write out Calderon's D and S matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b) :: i,j

       open(7,file="mat_SD.inp",status="unknown")
       write(7,*) nts_act
       do j=1,nts_act
        do i=1,nts_act
          if (Fprop(1:7).eq.'chr-ons') then 
            write(7,'(2E26.16)')BEM_S(i,j)
          else
            write(7,'(2E26.16)')BEM_S(i,j),BEM_D(i,j)
          endif
        enddo
       enddo

       close(7)

       return

      end subroutine


      subroutine read_BEM_SD
!------------------------------------------------------------------------
! @brief Read Calderon's D and S matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b) :: i,j

       open(7,file="mat_SD.inp",status="old")
       read(7,*) nts_act
       do j=1,nts_act
        do i=1,nts_act
          if (Fprop(1:7).eq.'chr-ons') then 
            read(7,*) BEM_S(i,j)
          else
            read(7,*) BEM_S(i,j), BEM_D(i,j)
          endif
        enddo
       enddo

       close(7)

       return

      end subroutine


      subroutine out_BEM_mat 
!------------------------------------------------------------------------
! @brief Output BEM matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b):: i,j

       open(7,file="BEM_matrices.mat",status="unknown")
       write(7,*) "# Q_0 , Q_d "
       write(7,*) nts_act
       do j=1,nts_act
        do i=j,nts_act
         write(7,'(2E26.16)')BEM_Q0(i,j),BEM_Qd(i,j)
        enddo
       enddo
       close(7)
       write(6,*) "Written out the propagation BEM matrixes"

       return

      end subroutine


      subroutine out_BEM_propmat 
!------------------------------------------------------------------------
! @brief Output propagation matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b):: i,j

       open(7,file="BEM_propmat.mat",status="unknown")
       if(Feps.eq."deb")write(7,*) "# \tilde{Q} , R "
       if(Feps.eq."drl")write(7,*) "# Q_w , Q_t "
       write(7,*) nts_act
       do j=1,nts_act
        do i=j,nts_act
         if(Feps.eq."deb")write(7,'(2E26.16)')BEM_Qt(i,j),BEM_R(i,j)
         if(Feps.eq."drl")write(7,'(2E26.16)')BEM_Qw(i,j),BEM_Qt(i,j)
        enddo
       enddo
       close(7)
       write(6,*) "Written out the propagation BEM matrixes"

       return

      end subroutine


      subroutine out_BEM_gamess  
!------------------------------------------------------------------------
! @brief Output propagation matrices 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b):: i,j

       open(7,file="np_bem.mat",status="unknown")
       do j=1,nts_act
        do i=1,nts_act
         write(7,'(D20.12)') BEM_Q0(i,j)/cts_act(i)%area
        enddo
       enddo
       close(7)
       write(6,*) "Written out the static matrix for gamess"
       open(7,file="np_bem.mdy",status="unknown")
       do j=1,nts_act
        do i=1,nts_act
         write(7,'(D20.12)') BEM_Qd(i,j)/cts_act(i)%area
        enddo
       enddo
       close(7)
       write(6,*) "Written out the dynamic matrix for gamess"

       return

      end subroutine


      subroutine out_BEM_diagmat 
!------------------------------------------------------------------------
! @brief Output diagonal matrices and frequencies 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b):: i,j

       open(7,file="BEM_TSpm12.mat",status="unknown")
       write(7,*) "# T , S^1/2, S^-1/2 "
       write(7,*) nts_act
       do j=1,nts_act
        do i=j,nts_act
         write(7,'(3D26.16)')BEM_T(i,j),Sp12(i,j),BEM_Sm12(i,j)
        enddo
       enddo
       close(7)
       open(7,file="BEM_W2L.mat",status="unknown")
       write(7,*) "# \omega^2 (drude-lorentz) , \Lambda, K_0, K_d "
       write(7,*) nts_act
       do j=1,nts_act
         write(7,'(i10, 4D20.12)') j,BEM_W2(j),BEM_L(j),K0(j),Kd(j)
       enddo
       close(7)
       write(6,*) "Written out BEM diagonal matrices in BEM_TSpm12.mat"
       write(6,*) "Written out BEM squared frequencies in BEM_W2L.mat"

       return

      end subroutine


      subroutine output_surf
!------------------------------------------------------------------------
! @brief Output cavity files 
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer(i4b) :: i
! This routine creates the same files also created by Gaussian

       open(unit=7,file="cavity.inp",status="unknown",form="formatted")
       write (7,*) nts_act,nesf_act
       do i=1,nesf_act
         write (7,'(3F22.10)') sfe_act(i)%x,sfe_act(i)%y, &
                               sfe_act(i)%z
       enddo
       do i=1,nts_act
         write (7,'(4F22.10,D14.5)') cts_act(i)%x,cts_act(i)%y,cts_act(i)%z, &
                               cts_act(i)%area,cts_act(i)%rsfe
       enddo
       close(unit=7)
! SC 28/9/2016: write out cavity in GAMESS ineq format
! WHICH UNITS ARE ASSUMED IN GAMESS? CHECK!!
       open(unit=7,file="np_bem.cav",status="unknown",form="formatted")
       write (7,*) nts_act
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%x
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%y
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%z
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%area
       enddo
       close(unit=7)
     
       return 

      end subroutine


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
       area=sum(cts_act(:)%area)

       do i=1,10
        write (fname,'("charge_freq_",I0,".pqr")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) nts_act
        do its=1,nts_act
          write (7,'("ATOM ",I6," H    H  ",I6,3F11.3,F15.5,"  1.5")') &
               its,its,cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
               Sm12T(its,i+1)/cts_act(its)%area*area
        enddo
        close(unit=7)
        ! SP : eigenvectors
        write (fname,'("eigenvector_",I0,".pdb")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) "MODEL        1" 
        maxt=zero
        do its=1,nts_act
          if(sqrt(TSm12(i+1,its)*TSm12(i+1,its)).gt.maxt)&
             maxt=sqrt(TSm12(i+1,its)*TSm12(i+1,its))
        enddo
        do its=1,nts_act
          write (7,'("ATOM ",I6,"  H   HHH H ",I3,"    ",3F8.3,2F6.2)') &
               its,i,cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
               TSm12(i+1,its)/maxt*10.,TSm12(i+1,its)/maxt*10.
        enddo
        close(unit=7)
       enddo
       do i=nts_act-10,nts_act-1
        write (fname,'("eigenvector_",I0,".pdb")') i
        open(unit=7,file=fname,status="unknown", &
           form="formatted")
        write (7,*) "MODEL        1" 
        maxt=zero
        do its=1,nts_act
          if(sqrt(TSm12(i+1,its)*TSm12(i+1,its)).gt.maxt)&
             maxt=sqrt(TSm12(i+1,its)*TSm12(i+1,its))
        enddo
        do its=1,nts_act
          write (7,'("ATOM ",I6,"  H   HHH H ",I3,"    ",3F8.3,2F6.2)') &
               its,i,cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
               TSm12(i+1,its)/maxt*10.,TSm12(i+1,its)/maxt*10.
        enddo
        close(unit=7)
       enddo
       ! SC 31/10/2016: print out total charges associated with eigenvectors
       open(unit=7,file="charge_eigv.dat",status="unknown", &
            form="formatted")
       write (6,*) "Total charge associated to each eigenvector written"
       write (6,*) "   in file charge_eigv.dat"
       do i=1,nts_act
         write(7,'(i20, 2D20.12)') i,sum(Sm12T(:,i))
       enddo
       close(unit=7)

       return

      end subroutine

      end module
