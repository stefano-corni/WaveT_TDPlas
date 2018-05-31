!> Module that reads the input of the medium. 
!! It contains all public medium variables
      module readio_medium
      use constants           
      use global_tdplas       
      use pedra_friends
      implicit none
      save
!
! medium object/surface variables 
      integer(i4b) :: n_q                        !< stride in updating the medium-molecule interaction
      integer(i4b) :: nsph                       !< number of spheres/oids for both dipole propagation and building BEM surface 
      ! Dipole propagation
      real(dbl), allocatable :: sph_centre(:,:)  !< Centers of spheres/oids (:,nsph)
      real(dbl), allocatable :: sph_maj(:)       !< Principal axis modulus (nsph) of spheroids (spheres: sph_maj=sph_min)
      real(dbl), allocatable :: sph_min(:)       !< Secondary axis modulus (nsph) of spheroids (spheres: sph_maj=sph_min)
      real(dbl), allocatable :: sph_vrs(:,:,:)   !< versors of principal (:,1,:) and secondary axis of nsph spheroids (:,:,nsph)
      real(dbl) :: fr_0(3)                       !< Reaction field at time 0 defined with Finit_mdm, here because used in scf
      ! Charges propagation
      real(dbl), allocatable :: vts(:,:,:)       !<transition potentials on tesserae from cis
      real(dbl), allocatable :: vtsn(:)          !<nuclear potential on tesserae
      real(dbl), allocatable :: q0(:)            !< Charges at time 0 defined with Finit_mdm, here because used in scf
      ! Restart
      !real(dbl)              :: fr_i(3),fx_i(3)  !< restart Onsager 
      !real(dbl), allocatable :: qr_i(:),qx_i(:)  !< restart pcm  
! Dielectric function variables 
      real(dbl) :: eps_0,eps_d                   !< $\omega \rightarrow 0$ and $\omega \rightarrow \infty$ limits of $\epsilon(\omega)$
      real(dbl) :: tau_deb                       !< Debye's $\tau_D$
      real(dbl) :: eps_A,eps_gm,eps_w0,f_vel     !< Drude lorentz $\omega^2_p$, $\gamma$, $\omega_$, and fermi velocity $v_f$
! SCF variables
      integer(i4b) :: ncycmax !< maximum number of SCF cycles
      real(dbl) :: thrshld    !< SCF threshold on (i) eigenvalues 10^-thrshld (ii) eigenvectors 10^-(thrshld+2)
      real(dbl) :: mix_coef   !< SCF mixing ratio of old (1-mix_coef) and new (mix_coef) charges/field                
! Frequancy calculation specifica variables
      integer   :: n_omega                !< number of points of the discretized spectrum
      real(dbl) :: omega_ini,omega_end    !< initial and final value of $\omega$ for spectrum
! SP 220216: vtsn: 
      integer(i4b) :: MPL_ord             ! Order of the multipole expansion (not used for the moment)
      integer(i4b), parameter :: nsmax=10 !< Maximum number of speres/spheroids to read from input
!SP 270517:
      character(flg) :: Fint,      & !< Interaction type "ons" for -mu*F "pcm" for q*V
                        Finit_int, & !< SCF "sce" for self-consistent initialization of the System-medium interaction
                        Fprop,     & !< Propagation type "dip" for dipole/field "chr-ief", "chr-ied", "chr-ons" ... for charges
                        Floc,      & !< "loc" to incude the medium polarization "external field" in the local field.
                        !Fmdm,      & !< Medium type nan or sol should be defined here and not in readio.f90
                        Finit_mdm, & !< Medium at time 0: charges/field read from file "rea", zero "vac" or frozen "fro" () 
                        Fmdm_pol,  & !< Medium polarization given by apparant charges "chr" or dipole "dip"                 
                        Fbem,      & !< Type of BEM calculation "dia" for diagonal version "std" for standard version
                        FinitBEM,  & !< Read "rea" or write "wri" cavity and BEM matrices.
                        Fsurf,     & !< Surface to be built "bui", read from file "fil", or read from a gms file "gms" 
                        Fshape,    & !< shape of the medium 'sphe' for sphere and 'spho' for spheroid
                        Feps,      & !< Epsilon choice "deb" for Debye and "drl" for Drude-Lorentz
                        Fqbem,     & !< type of BEM quantization, only full diagonalization "diag-all" for the moment
                        Fwrite,    & !< Modulates the level of output "low" and "high"
                        Fmdm_relax,& !< Medium charges follow the quantum jump "rel" or not "non"
                        Fgamess,&    !< if 'yes' write out matrix for gamess calculations of states
                        Ftest,     & !< Test Flag: see below
                        Fdeb       !< Debug Flag: see below
                        !Fmdm_res     !< Medium restart
                                     !! 
      
! namelists user-friendly variables 
      real(dbl) :: interaction_stride
      real(dbl) :: spheres_number
      real(dbl) :: sphere_position_x(nsmax)            
      real(dbl) :: sphere_position_y(nsmax)            
      real(dbl) :: sphere_position_z(nsmax)            
      real(dbl) :: sphere_radius(nsmax)            
      real(dbl) :: spheroids_number
      real(dbl) :: spheroid_axis_x(nsmax)
      real(dbl) :: spheroid_axis_y(nsmax)
      real(dbl) :: spheroid_axis_z(nsmax)
      real(dbl) :: spheroid_position_x(nsmax)            
      real(dbl) :: spheroid_position_y(nsmax)            
      real(dbl) :: spheroid_position_z(nsmax)            
      real(dbl) :: spheroid_radius(nsmax)            
      real(dbl) :: scf_threshold,scf_mix_coeff,scf_max_cycles
      character(flg) :: interaction_type,propagation_type,medium_pol,&
                        bem_read_write,input_surface,medium_init,    &
                        debug_type,bem_type,local_field,medium_type, &
                        out_level,interaction_init,epsilon_omega,    &
                        test_type,medium_relax,gamess

      ! variables read from eps.inp in the case of the general dielectric function case (i.e., eps_omega = 'gen')
      integer(i4b) :: npts
      real(dbl), allocatable    :: omegas(:)          !< sampling frequencies for the complex dielectric function
      complex(cmp), allocatable :: eps_omegas(:)      !< complex dielectric function values for the sampling frequencies
      real(dbl), allocatable    :: re_deps_domegas(:) !< real part of the derivative of the dielectric function at the sampling frequencies

      private
      public read_medium,deallocate_medium,Fint,Feps,Fprop,          &
             nsph,sph_maj,sph_min,sph_centre,sph_vrs,                &
             eps_0,eps_d,tau_deb,eps_A,eps_gm,eps_w0,f_vel,          &
             npts,omegas,eps_omegas,re_deps_domegas,                 &
             vts,n_q,Fmdm_pol,                                       &
             MPL_ord,Fbem,Fshape,fr_0,q0,Floc,                       &
             Fdeb,vtsn,Finit_int,Fqbem,Ftest,                        &
             ncycmax,thrshld,mix_coef,                               &
             FinitBEM,Fsurf,Finit_mdm,read_medium_freq,              &
             read_medium_tdplas,n_omega,omega_ini,omega_end,         &
             Fwrite,Fmdm_relax,Fgamess!,Fmdm_res,fr_i,fx_i,qr_i,qx_i,         &
             !read_medium_restart 
!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  DRIVER  ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine read_medium
!------------------------------------------------------------------------
! @brief Driver routine for reading medium input 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       namelist /propagate/interaction_stride,interaction_init,        &
                         interaction_type,propagation_type,            &
                         scf_mix_coeff,scf_max_cycles,scf_threshold,   &
                         local_field,debug_type,out_level,test_type,   &
                         medium_relax
       namelist /medium/ medium_type,medium_init,medium_pol,bem_type,  &
                         bem_read_write                   
       namelist /surface/input_surface,spheres_number,spheroids_number,&
         sphere_position_x,sphere_position_y,sphere_position_z,        &
         spheroid_axis_x,spheroid_axis_y,spheroid_axis_z,              &
         spheroid_position_x,spheroid_position_y,spheroid_position_z,  &
         sphere_radius,spheroid_radius                                
       namelist /eps_function/epsilon_omega,eps_0,eps_d,eps_A, &
                         eps_gm,eps_w0,f_vel,tau_deb       
       call init_nml_all() 
       call init_nml_propagate() 
       read(*,nml=propagate) 
       call write_nml_propagate()
       call write_nml_interaction()
       read(*,nml=medium) 
       call write_nml_medium()
       read(*,nml=surface) 
       call write_nml_surface()
       if (Fmdm(2:4).eq.'nan') then
         call init_nml_nanoparticle() 
       elseif (Fmdm(2:4).eq.'sol') then
         call init_nml_solvent()
       endif
       read(*,nml=eps_function) 
       call write_nml_eps_function()
       call write_nml_all()

       return

      end subroutine


      subroutine read_medium_freq
!------------------------------------------------------------------------
! @brief Driver routine for reading medium input form main_freq 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       namelist /freq/ fmax,n_omega,omega_ini,omega_end,debug_type, &
                       out_level,test_type
! SP 13/07/18:: medium and surface not used but added for future implementations
       namelist /medium/ medium_type,medium_init,medium_pol,bem_type,  &
                         bem_read_write                   
       namelist /surface/input_surface,spheres_number,spheroids_number,&
         sphere_position_x,sphere_position_y,sphere_position_z,        &
         spheroid_axis_x,spheroid_axis_y,spheroid_axis_z,              &
         spheroid_position_x,spheroid_position_y,spheroid_position_z,  &
         sphere_radius,spheroid_radius                                
       namelist /eps_function/epsilon_omega,eps_0,eps_d,eps_A, &
                         eps_gm,eps_w0,f_vel,tau_deb       
       call init_nml_all() 
       call init_nml_freq()
       read(*,nml=freq) 
       write(*,nml=freq)
       read(*,nml=medium) 
       write(*,nml=medium)
       call write_nml_medium() 
       read(*,nml=surface) 
       call write_nml_surface() 
       read(*,nml=eps_function) 
       write(*,nml=eps_function)
       call write_nml_eps_function()
       call write_nml_all() 

       return

      end subroutine


      subroutine read_medium_tdplas
!------------------------------------------------------------------------
! @brief Driver routine for main_tdplas 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------


       !namelist /tdplas/ debug
       namelist /medium/ medium_type,medium_init,medium_pol,bem_type,  &
                         bem_read_write,out_level,debug_type,test_type 
       namelist /surface/input_surface,spheres_number,spheroids_number,&
         sphere_position_x,sphere_position_y,sphere_position_z,        &
         spheroid_axis_x,spheroid_axis_y,spheroid_axis_z,              &
         spheroid_position_x,spheroid_position_y,spheroid_position_z,  &
         sphere_radius,spheroid_radius                                
       namelist /eps_function/ epsilon_omega,eps_0,eps_d,eps_A,&
                               eps_gm,eps_w0,f_vel,tau_deb
       namelist /out_matrix/ gamess
       call init_nml_all() 
       call init_nml_tdplas() 
       !read(*,nml=tdplas)
       read(*,nml=medium) 
       call write_nml_medium() 
       read(*,nml=surface) 
       call write_nml_surface() 
       read(*,nml=eps_function) 
       call write_nml_eps_function()
       read(*,nml=out_matrix,end=10) 
10     call write_nml_all()

       return
      end subroutine read_medium_tdplas
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!  INITIALIZATION  ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine init_nml_all 
!------------------------------------------------------------------------
! @brief Initialize variables for all mains to safe values 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------


       ! Output and debug
       out_level="low"
       debug_type="non"
       test_type="non"
       ! Sphere and Spheroid
       spheres_number=0 
       sphere_radius=zero
       sphere_position_x=zero
       sphere_position_y=zero
       sphere_position_z=zero
       spheroids_number=0 
       spheroid_axis_x=zero
       spheroid_axis_y=zero
       spheroid_axis_z=zero
       spheroid_position_x=zero
       spheroid_position_y=zero
       spheroid_position_z=zero
       spheroid_radius=zero
       ! SCF
       scf_threshold=10
       scf_mix_coeff=0.2
       scf_max_cycles=600
       ! Various Medium
       medium_init='fro'
       MPL_ord=1 ! SP 29/06/17: multipole order (not used)
       nts_act=0 

       return

      end subroutine


      subroutine init_nml_propagate()
!------------------------------------------------------------------------
! @brief Initialize variables for propagation main (will be tdplas) 
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       medium_init='fro'
       medium_type='nan'
       medium_pol='chr'
       medium_relax="non"
       interaction_stride=1
       interaction_type='pcm'
       interaction_init='non-scf'
       bem_type='diag'
       bem_read_write='rea'
       input_surface='fil'
       propagation_type='ief'
       local_field='loc'

       return

      end subroutine 


      subroutine init_nml_nanoparticle()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist nanoparticle 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  : SP 10/07/17
! @param epsilon_omega,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel
!------------------------------------------------------------------------

       epsilon_omega='drl'
       tau_deb=1000.
       eps_0=1000. 
       eps_d=1.      
       eps_A=0.110224
       eps_gm=0.000757576
       eps_w0=0. 
       f_vel=0. 

       return

      end subroutine init_nml_nanoparticle


      subroutine init_nml_solvent()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist solvent 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  : SP 10/07/17
! @param epsilon_omega,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel
!------------------------------------------------------------------------

       epsilon_omega='deb'
       tau_deb=1000.
       eps_0=35.688
       eps_d=1.806874
       eps_A=0.110224
       eps_gm=0.000757576
       eps_w0=0.0
       f_vel=0.0

       return

      end subroutine init_nml_solvent


      subroutine init_nml_freq()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist freq 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  : SP 10/07/17
!------------------------------------------------------------------------

       ! SP: No propagation: Fprop set to other than "dip" or "chr" 
       Fprop="non"
       medium_type='nan'
       medium_pol='chr'
       bem_type='diag'
       bem_read_write='rea'
       input_surface='fil'
       epsilon_omega='drl'
       eps_0=1000.   
       eps_d=1.      
       eps_A=0.110224
       eps_gm=0.000757576
       eps_w0=zero
       f_vel=zero
       fmax=0.0001
       n_omega=5000 
       omega_ini=0.01
       omega_end= 0.35
       ! The following are not 
       interaction_type='pcm'

       return

      end subroutine init_nml_freq


      subroutine init_nml_tdplas()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist tdplas 
!
! @date Created   : E. Coccia 16 May 2017
! Modified  : SP 14/07/17
!------------------------------------------------------------------------

       ! SP: No propagation: Fprop set to other than "dip" or "chr" 
       Fprop="non"
       medium_type='nan'
       medium_pol='chr'
       bem_read_write='rea'
       bem_type='diag'
       input_surface='bui'
       epsilon_omega='deb'
       eps_0=35.688
       eps_d=1.806874
       tau_deb=1000.
       eps_A=0.110224
       eps_w0=zero
       eps_gm=0.000757576
       f_vel=zero
       gamess='no '

       return

      end subroutine init_nml_tdplas
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  VARIABLE DEFINITION ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine write_nml_all 
!------------------------------------------------------------------------
! @brief Write variables for all mains 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       ! Output level
       select case (out_level)
        case ('high','High','HIGH')
         Fwrite='high'
         write(*,*) 'Complete output written.'
        case ('low','Low','LOW')
         Fwrite='low'
         write(*,*) 'Essential output written.'
        case default
         write(*,*) "Please specify one of the following output levels"
         write(*,*) "- high"
         write(*,*) "- low"
         stop
       end select
       ! Test  
       select case(test_type)
       case ('n-r','N-r','n-R','N-R')
        Ftest='n-r'
        Floc='non'
        write(6,*) "TEST: Nanoparticle reaction field"
       case ('n-l','N-l','n-L','N-L')
        Ftest='n-l'
        Ffld='snd'
        Floc='loc'
        write(6,*) "TEST: Nanoparticle local field"
       case ('s-r','S-r','s-R','S-R')
        Ftest='s-r'
        Floc='non'
        write(6,*) "TEST: Solvent reaction field"
       case ('s-l','S-l','s-L','S-L')
        Ftest='s-l'
        Ffld='snd'
        Floc='loc'
        write(6,*) "TEST: Solvent local field"
       case ('QMT','Qmt','qmt')
        Ftest='qmt'
        write(6,*) "TEST: dipolar quantum coupling calculation"
       case default
        Ftest='non'
       end select
       ! Debug
       select case(debug_type)
       case ('equ','Equ','EQU')
        Fdeb='equ'
        write(6,*) "DEBUG: Equilibrium reaction field calculation"
       case ('vmu','Vmu','VMU','VMu')
        Fdeb='vmu'
        write(6,*) "DEBUG: Potentials calculated from Dipoles "
       case ('off','Off','OFF')
        Fdeb='off'
        write(6,*) "DEBUG: Molecule - Medium interaction turned off"
       case default
        Fdeb='non'
       end select
       select case(gamess)
       case ('yes','Yes','YES')
        Fgamess='yes'
       case default
        Fgamess='no '
       end select
       return

      end subroutine


      subroutine write_nml_propagate()
!------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       ! propagation_type refers to which quantity is propagated by equations of motions
       !   dip: only the dipolar (i.e., Onsager) reaction/local field/dipole is propagated
       !   ief,ied,ons: the apparent charges are propagated, within different schemes
       select case (propagation_type)
        case ('ief','Ief','IEF')
         write(6,*) 'Apparent charges are propagated'
         Fprop='chr-ief'
        case ('ied','Ied','IED')
         write(6,*) 'Apparent charges are propagated'
         Fprop='chr-ied'
        case ('ons','Ons','ONS')
       ! COSMO propagation PCM, Onsager Nanoparticle
         write(6,*) 'Apparent charges are propagated'
         Fprop='chr-ons'
        case ('dip','Dip','DIP')
       ! Dipole propagation no charges involved 
         write(6,*) 'Only the dipolar reaction field is propagated'
         Fprop='dip'
        case default
         write(*,*) "Error, specify the propagation type "
         stop
       end select
       ! Local Field: include or not external field into local field.
       select case(local_field)
       case ('loc','Loc','LOC')
        Floc='loc'
        write(6,*) "Local field effects are included"
       case default
        write(6,*) "Local field effects are NOT included"
       end select
! SP 270917: added when merging to newer master 
       select case(medium_relax)
          case ('rel', 'Rel', 'REL', 'REl')
           write(*,*) 'Medium charges follow the quantum jump'
           !np_relax=.true.
           Fmdm_relax="rel"
          case ('non', 'NoN', 'NON', 'NOn')
           write(*,*) 'Medium charges do not follow the quantum jump'
           !np_relax=.false.
           Fmdm_relax="non"
          case default
           Fmdm_relax="non"
           !np_relax=.true.
       end select
! EC 281117: added restart for medium
       !select case(medium_res)
       !   case ('y','Y')
       !    write(*,*) 'Restart for medium'
       !    Fmdm_res='Yesr'
       !    call read_medium_restart()
       !   case ('n','N')
       !    Fmdm_res='Nonr' 
       !end select

       return

      end subroutine


      subroutine write_nml_tdplas()
!------------------------------------------------------------------------
! @brief Write variables in the namelist tdplas and put conditions 
!
! @date Created   : E. Coccia 16 May 2017
! Modified  :
! @param epsilon_omega,eps_0,eps_d,tau_deb,eps_A,eps_gm,
!        eps_w0,f_vel,input_surface,xr,yr,zr,rr,nsph  
!------------------------------------------------------------------------
       return
      end subroutine  write_nml_tdplas


      subroutine write_nml_eps_function()
!------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia, G. Gil
!------------------------------------------------------------------------

       integer(i4b) :: i

       select case (epsilon_omega)
         case ('deb','Deb','DEB')
           Feps='deb'
         case ('drl','Drl','DRL')
           Feps='drl'
         case ('gen','Gen','GEN','gral','Gral','GRAL')
           Feps='gen'
           open(1,file='eps.inp')
           read(1,*) npts
           allocate(omegas(npts),eps_omegas(npts),re_deps_domegas(npts))
           do i=1, npts
            read(1,*) omegas(i), eps_omegas(i)
            if(i.eq.1) cycle
            re_deps_domegas(i-1) = real(eps_omegas(i)-eps_omegas(i-1))/(omegas(i)-omegas(i-1))
           enddo
           re_deps_domegas(npts) = zero
           close(1)
         case default
           write(*,*) "Error, specify eps(omega) type DEB, DRL or GEN"
           stop
       end select        

       return

      end subroutine


      subroutine write_nml_interaction()
!------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       n_q=interaction_stride
       write(*,*) 'Frequency of updating the interaction potential', n_q
       select case(interaction_init)
! SC & SP: determine the starting state & charges of the simulation:
!      SCF_ES: !      Default: no self-consistent optimization, initial charges
!               are calculated for the ground state
        case('SCF_ES')
         ! self consistent optimization of the states in the RF of the 
         ! state defined by ci_ini.inp, charges for such state
         Finit_int='sce'
         thrshld=scf_threshold
         ncycmax=scf_max_cycles
         mix_coef=scf_mix_coeff
         write(6,*) "Do SCF calculation for the state in ci_ini.inp"
         write(*,*) 'SCF threshold', thrshld
         write(*,*) 'Mixing coefficient', mix_coef
         write(*,*) 'Max number of SCF cycles', ncycmax
        case default
         Finit_int='nsc'
         write(6,*) "Use GS RF as coded in the ci_energy.inp values"
       end select
       select case (interaction_type)
       ! read interaction type: PCM or Onsager
       ! interaction_type refers to the coupling between the molecule and the reaction field
       !   ons: the reaction field is considered constant in space and the coupling is
       !   -mu*F_RF
       !   pcm: the reaction potential is used, the coupling is \int dr rho(r) V_RF(r)=sum_i q_i*V_rho(r_i)
        case ('ons','Ons','ONS')
         Fint='ons'
         write(*,*) 'The coupling is of Onsager type: dip \cdot Fld'
        case ('pcm','Pcm','PCM')
         Fint='pcm'
         write(*,*) 'The coupling is of PCM type: pot \cdot chr'
         if (Fprop.eq.'dip') then
           write(*,*) &
           "Error: PCM incompatible with dipole/field propagation"
           stop
         endif
        case default
         write(*,*) "Error, specify interaction type ONS or PCM"
         stop
       end select

       return

      end subroutine


      subroutine write_nml_medium()
!------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

! PER STEFANO:
! SP 14/07/17: medium_type is the same of medium in maedium.f90 except for the vacuum case
!              the medium type could be specified at this level and the activation of the 
!              call to the read_medium should be due to the "-medium medium.inp" as an 
!              input file passed at the command line level, see:
!              https://genomeek.wordpress.com/2012/02/09/fortran-command-line-arguments/  
       select case (medium_type)
        case ('sol','Sol','SOL')
          write(*,*) "Solvent as external medium"
          Fmdm='Csol'
        case ('qso','Qso','QSO')
          write(*,*) "Quantum Solvent as external medium"
          Fmdm='Csol'
        case ('nan','Nan','NAN')
          write(*,*) "Nanoparticle as external medium"
          Fmdm='Cnan'
        case ('Qna','qna','QNA')
          write(*,*) "Quantum Nanoparticle as external medium"
          Fmdm='Qnan'
        case default
          write(*,*) "Error: Choose a medium type"
          stop 
       end select
       select case (medium_pol)
        case ('chr','Chr','CHR')
          write(*,*) "Medium polarization described by apparent charges"
          Fmdm_pol='chr'
        case ('dip','Dip','DIP')
          write(*,*) "Medium polarization described by -mu*F term"
          Fmdm_pol='dip'
        case default
          write(*,*) "Error: Choose a medium type"
          stop 
       end select
       if (Fmdm_pol.eq.'chr') then
       ! if this is a run with an explicit boundary (fprop.neq.dip),
       ! then
       ! we need to know if the boundary and the matrix are read from outside
       ! or are produced here and just written out, with no real propagation
         select case (bem_type)
          case ('diag','Diag','DIAG')
           write(6,*) 'Diagonal BEM furmulation'
           Fbem="diag"
          case ('stan','Stan','STAN')
           write(6,*) 'Standard BEM furmulation not implemented yet'
           stop
          case default
           write(*,*) "Error, specify a BEM type "
           stop
         end select
         if (Fmdm(1:1).eq.'C') then
           write (6,*) "This is a Classical BEM run"
         elseif (Fmdm(1:1).eq.'Q') then
! SP 220617: Adding quantum coupling                            
           write (6,*) "This is a Quantum BEM run"
           ! This is the only option by now
           Fqbem='diag-all' ! couple with all modes but only one occupied
         endif
         select case(bem_read_write)
         case ('rea','Rea','REA')
          FinitBEM='rea'
          if(Fprop(1:3).eq."chr") call read_gau_out_medium
          write(6,*) "This is full run reading matrix and boundary"
         case ('wri','Wri', 'WRI')
          FinitBEM='wri'
          write(6,*) "This run just writes matrices and boundary"
         case default
          write(*,*) "Error, specify if boundary data & matrices", &
             "are read (read) or made and written out (write) "
          stop
         end select
       endif
       if (Fprop(1:3).eq.'chr'.or.Fprop(1:3).eq.'dip') then
         ! Medium initialization for propagation: how to set q0 and fr_0
         select case(medium_init)
         case ('vac','VAC','Vac') ! q0=zero, fr_0=zero
          Finit_mdm='vac'
          write(6,*) "Medium polarization set initially to zero"
         case ('fro','Fro','FRO') ! q0=matmul(Q_0,pot_0), fr_0=ONS_f0*mu_0
          write(6,*) "Medium polarization frozen to molecule in its GS"
          Finit_mdm='fro'
         case ('rea','Rea','REA') ! read from file    
          write(6,*) "Medium reaction field/charges read from",& 
                                           " charges0.inp file"
          Finit_mdm='rea'
         case default
          write(6,*) "Please choose how to initialize charges/RField"
          stop
         end select
       endif     

       return

      end subroutine 


!------------------------------------------------------------------------
! SP 14/07/17 calculations should probably go in a different module. which one?
!             Probably pedra_firends....
      subroutine write_nml_surface()
!------------------------------------------------------------------------
! @brief Write variables for surface/medium object 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       integer(i4b)::i,j

       if (spheres_number.gt.0) nsph=spheres_number
       if (spheroids_number.gt.0) nsph=spheroids_number
       if (Fprop(1:3).eq.'dip') then
       ! fill sph_centre,sph_min,sph_maj,sph_vrs that define the medium object for calculations/propagation
         allocate(sph_maj(nsph))
         allocate(sph_min(nsph))
         allocate(sph_vrs(3,3,nsph))
         allocate(sph_centre(3,nsph))
         if(spheres_number.gt.0) then
           write(6,*) 'This is a Spherical Onsager run in solution'
           ! sphere major axis = minor axis = radius           
           Fshape='sphe' ! Sphere
           if(Fmdm(2:4).eq.'sol')write(*,*)'Spherical cavity'
           if(Fmdm(2:4).eq.'nan')write(*,*)'Spherical nanoparticle' 
           do i=1,nsph
             sph_centre(1,i)=sphere_position_x(i)
             sph_centre(2,i)=sphere_position_y(i)
             sph_centre(3,i)=sphere_position_z(i)
             sph_min(i)=sphere_radius(i)
             sph_maj(i)=sphere_radius(i)
             sph_vrs(:,:,i)=zero
             write(*,*) 'Radius (a.u.) ', sph_maj(i)
           enddo
         else
           write(6,*) 'This is a Spheroidal Onsager run in solution'
           Fshape='spho' ! Spheroid
           ! calculate the major axis modulus and unit vector
           do i=1,nsph
             sph_centre(1,i)=spheroid_position_x(i)
             sph_centre(2,i)=spheroid_position_y(i)
             sph_centre(3,i)=spheroid_position_z(i)
             sph_min(i)=spheroid_radius(i)
             sph_vrs(1,1,i)=spheroid_axis_x(i)
             sph_vrs(2,1,i)=spheroid_axis_y(i)
             sph_vrs(3,1,i)=spheroid_axis_z(i)
             sph_maj(i)=sqrt(sph_vrs(1,1,i)**2+ &
                             sph_vrs(2,1,i)**2+ &    
                             sph_vrs(3,1,i)**2)     
             if(sph_maj(i).eq.zero) then
               write(6,*) 'ERROR: please provide spheroid axis'
               stop
             endif
             sph_vrs(:,1,i)=sph_vrs(:,1,i)/sph_maj(i)
             if(Fmdm(2:4).eq.'sol')write(*,*)'Spheroidal cavity'
             if(Fmdm(2:4).eq.'nan')write(*,*)'Spheroidal nanoparticle'
             write(*,*) 'Principal axis (a.u)', sph_maj(i)
             write(*,'(a,3F10.5)') 'Principal direction (a.u.) ', &
                                        (sph_vrs(j,1,i),j=1,3)
             write(*,*) 'Secondary axis (a.u)', sph_min(i)
             write(*,'(A,3F10.5)') 'Centre ', (sph_centre(j,i),j=1,3)
           enddo
         endif
       else
       ! fill sph_centre,sph_min,sph_maj,sph_vrs that define the medium object for calculations/propagation
         select case(input_surface)
         case ('fil','FIL','Fil')
          Fsurf='fil'
          write(6,*) "Surface read from file cavity.inp"
         case ('gms','GMS','Gms')
          Fsurf='gms'
          write(6,*) "Surface read from file surface_msh.inp"
         case ('bui','Bui','BUI')
          Fsurf='bui'
          write(6,*) "Building surface from spheres."
          if (spheres_number.le.0) then
             write(*,*) 'ERROR: number of spheres is zero or negative '
             stop
          endif
          if (spheres_number.gt.nsmax) then
             write(*,*) 'ERROR: nsph is larger than', nsmax
             stop
          endif
          ! Build surface from spheres
          call read_act(sphere_position_x,& 
                        sphere_position_y,&
                        sphere_position_z,&
                        sphere_radius,nsph,nsmax)
         case default
          write(6,*) "Please choose: build or read surface?"
          stop
         end select
       endif

       return

      end subroutine 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  READ/WRITE ROUTINES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      subroutine read_sph_fromfile
!------------------------------------------------------------------------
! @brief Read spheres/oids parameters from file 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       integer(i4b) :: i,j,its
       real(dbl)  :: scr       

       open(7,file="sph.inp",status="old")
       read(7,*) scr 
       if(scr.ne.nsph) then
         write(6,*) "Wrong number of spheres/oids"
         stop
       endif
       nsph=scr
       allocate(sph_centre(3,nsph),sph_min(nsph),sph_vrs(3,3,nsph))
       do i=1,nsph
         read(7,*) (sph_centre(j,i),j=1,3),sph_min, &
                   (sph_vrs(j,1,i),j=1,3) 
       enddo

       close(7)

       return

      end subroutine


      subroutine read_gau_out_medium
!------------------------------------------------------------------------
! @brief Read transition potentials on tesserae 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       integer(i4b) :: i,j,its,nts
       real(dbl)  :: scr       

       open(7,file="ci_pot.inp",status="old")
       read(7,*) nts
       if(nts_act.eq.0.or.nts.eq.nts_act) then
         nts_act=nts
       else
         write(*,*) "Tesserae number conflict"
         stop
       endif
       allocate (vts(nts_act,n_ci,n_ci))
       allocate (vtsn(nts_act))
       ! V00
       read(7,*) 
       do its=1,nts_act
        read(7,*) vts(its,1,1),scr,vtsn(its)
        vts(its,1,1)=vts(its,1,1)+vtsn(its)
       enddo
       !V0j
       do j=2,n_ci_read
         read(7,*) 
         if (j.le.n_ci) then
          do its=1,nts_act
          read(7,*) vts(its,1,j)
          enddo
          vts(:,j,1)=vts(:,1,j)
         else
          do its=1,nts_act
           read(7,*)
          enddo
         endif
       enddo
       !Vij
       do i=2,n_ci_read
        do j=2,i   
         read(7,*) 
         if (i.le.n_ci.and.j.le.n_ci) then
          do its=1,nts_act
           read(7,*) vts(its,i,j)             
          enddo
          vts(:,j,i)=vts(:,i,j)
         else
          do its=1,nts_act
           read(7,*) 
          enddo
         endif
        enddo
        ! add nuclear potential
        if (i.le.n_ci) then
         do its=1,nts_act
          vts(its,i,i)=vts(its,i,i)+vtsn(its)
         enddo
        endif
       enddo
       close(7)

       return

      end subroutine


      subroutine output_surf
!------------------------------------------------------------------------
! @brief Output surface.xyz file 
!      
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------

       integer :: i

       open(unit=7,file="surface.xyz",status="unknown",form="formatted")
        write (7,*) nts_act
        do i=1,nts_act
          write (7,'(3F22.10)') cts_act(i)%x,cts_act(i)%y,cts_act(i)%z
        enddo
       close(unit=7)
      end subroutine


      subroutine deallocate_medium
!------------------------------------------------------------------------
! @brief Deallocate medium arrays 
!      
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------

       if(allocated(q0)) deallocate(q0)
       if(allocated(vts)) deallocate(vts)
       if(allocated(vtsn)) deallocate(vtsn)
       if(allocated(sph_maj)) deallocate(sph_maj)
       if(allocated(sph_min)) deallocate(sph_min)
       if(allocated(sph_vrs)) deallocate(sph_vrs)
       if(allocated(sph_centre)) deallocate(sph_centre)

       return

      end subroutine

     
!      subroutine read_medium_restart() 
!------------------------------------------------------------------------
! @brief Read restart 
!
! @date Created   : E. Coccia 28 Nov 2017
! Modified  :
!------------------------------------------------------------------------
       
!       implicit none

!       integer(i4b)     :: i
!       character(3)     :: cdum 
!       logical          :: exist

!       inquire(file='restart_mdm', exist=exist)
!       if (exist) then
!          open(779, file='restart_mdm', status="old")
!       else
!          write(*,*) 'ERROR:  file restart_mdm is missing'
!          stop
!       endif

       !if (Fint.eq.'ons') then
!          read(779,*) cdum 
!          do i=1,3
!             read(779,*) fr_i(1), fr_i(2), fr_i(3)
!          enddo
!          if (Floc.eq.'loc') then
!             read(779,*) cdum 
!             do i=1,3
!                read(779,*) fx_i(1), fx_i(2), fx_i(3)
!             enddo
!          endif
       !elseif (Fint.eq.'pcm') then
!          read(779,*) cdum 
!          do i=1,nts_act
!             read(779,*) qr_i(i)
!          enddo
!          if (Floc.eq.'loc') then
!             read(779,*) cdum 
!             do i=1,nts_act
!                read(779,*) qx_i(i)
!             enddo
!          endif
       !endif

!       close(779)

!       return

!      end subroutine read_medium_restart

 end module
