      module readio
      use constants       
      use random
#ifdef MPI
      use mpi
#endif
      implicit none
      save
!
      integer(i4b) :: n_f,n_ci,n_ci_read,n_step,n_out,ncit 
      !integer(i4b) :: imar !imar=0 Markvian, imar=1, nonMarkovian
      integer(i4b) :: i_sp=0,i_nr=0,i_de=0 !counters for quantum jump occurrences
      integer(i4b) :: nrnd !the time step for Euler-Maruyama is dt/nrnd
! SP 17/07/17: Changed to char flags
      integer(i4b) :: nr_typ !input integer for type of decay for the internal conversion
      integer(i4b) :: idep   !input integer for the dephasing operator
      integer(i4b) :: tdis   !input integer for Euler tdis=0, Matthews tdis=1 
! SP270917: added for merging to newer master
      integer(i4b)              :: npulse    !number of pulses
      integer(i4b), allocatable :: irel(:,:) !mapping for intermediate relaxations  
      integer(i4b), allocatable :: ik(:,:)   !mapping for intermediate relaxations 
      integer(i4b)              :: restart_i ! step for restart
      integer(i4b)              :: n_restart ! frequency of restart writing
      integer(i4b)              :: n_jump    ! number of quantum jumps along a simulation
      integer(i4b)              :: diff_step ! effective number of steps for restart
      integer(i4b)              :: restart_seed  ! seed for restart
      real(dbl), allocatable    :: mut_np2(:,:) !squared dipole from NP
      real(dbl)                 :: tdelay(npulsemax), pshift(npulsemax)  ! time delay and phase shift with two pulses
!
      !real(dbl), allocatable    :: c_i(:),e_ci(:)  ! energy from cis
      real(dbl), allocatable    :: e_ci(:)  ! energy from cis
      complex(cmp), allocatable :: c_i(:),c_i_prev(:),c_i_prev2(:) ! coefficients from cis
      real(dbl)                 :: mu_i_prev(3),mu_i_prev2(3),mu_i_prev3(3),mu_i_prev4(3),mu_i_prev5(3)
      real(dbl), allocatable    :: mut(:,:,:) !transition dipoles from cis
      real(dbl), allocatable    :: nr_gam(:), de_gam(:) !decay rates for nonradiative and dephasing events
      real(dbl), allocatable    :: sp_gam(:) !decay rate for spontaneous emission  
      real(dbl), allocatable    :: sp_fact(:) !multiplicative factor for the decay rate for spontaneous emission
      real(dbl), allocatable    :: tmom2(:) !square transition moments i->0
      real(dbl), allocatable    :: delta(:) !phases randomly added during the propagation 
      real(dbl), allocatable    :: de_gam1(:) !combined decay rates for idep=1
      real(dbl), allocatable    :: tomega(:)  !generalized transition frequencies 
      real(dbl), allocatable    :: ion_rate(:) !ionization rates
      real(dbl)                 :: restart_t  ! time for restart
      real(dbl)                 :: dt,tau(2),start,krnd
! SP 17/07/17: Changed to char flags
      !logical :: dis !turns on the dissipation
      !logical :: qjump ! =.true. quantum jump, =.false. stochastic propagation
      !logical :: ernd=.false. !add normal number to E: E -> E + krnd*rnd()
! Global flags        
      character(flg) :: Fdis_rel  !< Flag for decay for internal conversion, relaxation via dipole "dip" or matrix "mat"
      character(flg) :: Fdis_deph !< Flag for dephasing operator: exp(i delta_i)|i><i| "exp" or |i><i|-|0><0| "i-0" 
      character(flg) :: Fdis !< Flag for dissipation type: 
                             !! Markovian     : quantum jumps "mar-qjump", Euler-Maruyama "mar-EuMar", Leimkuhler-Matthews "mar-LeiMa"
                             !! Non-Markovian : quantum jumps "nma-qjump", Continuous stochastic propagator "nma-cstoc"
                             !! Random        : random energy term "ernd" 
      ! qjump works for the Markovian case
      ! Stochastic propagation from:
      ! Appl. Math. Res. Express vol. 2013 34-56 (2013)
      ! IMA J. Numer. Anal. vol. 36 13-79 (2016)
      real(dbl) :: t_mid,sigma(npulsemax),dir_ft(3),fmax(3,npulsemax),omega(npulsemax),mol_cc(3)
      character(flg) :: Ffld !< Field type 
      character(flg) :: Fmdm !< Flag for medium type, this will be defined in readio_medium after separation
      character(flg) :: Frad !< Flag for radiative damping 
      character(flg) :: Fres !< Flag for restart
      character(flg) :: Fful !< Flag for relaxation matrix
      character(flg) :: Fexp !< Flag for propagation of the energy term
      character(flg) :: Fsim !< Flag for the dynamics length in case of restart
      character(flg) :: Fabs !< Flag for the dynamics in presence of an absorber
      character(flg) :: Fbin !< Flag for writing output files with binary format
      character(flg) :: Fopt !< Flag for using OMP-optimized matrix/vector multiplication 
      character(flg) :: Fwrt !< Flag for SSE output
! Flags read from input file
      character(flg) :: medium,radiative,dissipative,lsim,absorber,binary,out_sse
      character(flg) :: dis_prop
      character(flg) :: restart 
      character(flg) :: propa
      character(flg) :: full ! full or only |e> -> |0> relaxation
      integer(i4b) :: iseed  ! seed for random number generator
      integer(i4b) :: nexc   ! number of excited states
      integer(i4b) :: nrel   ! number of relaxation channels
      integer(i4b) :: nf     ! number of relaxation channels (including |e> -> |0> terms) 
      integer(i4b) :: i,nspectra
! kind of surrounding medium and shape of the impulse
!     Fmdm=sol: solvent
!     Fmdm=nan: nanoparticle
!     Fmdm=vac: no medium
!
!     Ffld=gau: gaussian impulse
!SC
!     Frad: wheter or not to apply radiative damping
!     TO BE COMPLETED, SEE PROPAGATE.F90      \
      real(dbl) :: eps_A,eps_gm,eps_w0,f_vel
      real(dbl) :: r
      
      private
      public read_input,n_ci,n_ci_read,n_step,dt,           &
             Ffld,t_mid,sigma,omega,fmax,restart,           & 
             Fmdm,mol_cc,tau,start,c_i,e_ci,mut,            &
             Frad,n_out,iseed,n_f,dir_ft,full,              &
! SP 17/07/17: Changed to char flags
             tdis,nr_gam,de_gam,sp_gam,tmom2,nexc,delta,    &
             deallocate_dis,i_sp,i_nr,i_de,nrnd,sp_fact,    &
!             nr_typ,idep,imar,de_gam1,krnd,ernd       
             de_gam1,krnd,Fdis,Fdis_deph,Fdis_rel,nf,irel,  &
             npulse,tdelay,pshift,nrel,Fful,  &
             Fexp,Fres,restart_t,restart_i,n_restart,       &
             c_i_prev,c_i_prev2,mu_i_prev,mu_i_prev2,       &
             mu_i_prev3,mu_i_prev4,mu_i_prev5,restart_seed, &
             n_jump,Fsim,diff_step,mpibcast_readio,         &
             mpibcast_e_dip,mpibcast_sse,mpibcast_restart,  &
             nspectra,Fabs,ion_rate,mpibcast_ion_rate,Fbin, &
             ncit,Fopt,ik,Fwrt 
             
!
      contains
!
      subroutine read_input
!------------------------------------------------------------------------
! @brief Read input namelists 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------


       !integer(i4b):: i,nspectra
       !character(3) :: medium,radiative,dissipative
       !character(5) :: dis_prop 
     
       !Molecular parameters 
       namelist /general/n_ci_read,n_ci,mol_cc,n_f,medium,restart,full,& 
                         dt,n_step,n_out,propa,n_restart,lsim,absorber,&
                         binary,ncit
       !External field paramaters
       namelist /field/ Ffld,t_mid,sigma,omega,radiative,iseed,fmax, &
                        npulse,tdelay,pshift
       !Stochastic Schroedinger equation
       namelist /sse/ dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd,out_sse
       !Namelist spectra
       namelist /spectra/ start,tau,dir_ft

       write(*,*) ''
       write(*,*) '*****************************************************'
       write(*,*) '*****************************************************'
       write(*,*) '**                                                 **' 
       write(*,*) '**                    WaveT                        **'
       write(*,*) '**     Evolution of molecular wavefunctions        **'
       write(*,*) '**   under external electromagnetic perturbations  **'
       write(*,*) '**                                                 **'
       write(*,*) '**                     by                          **'
       write(*,*) '**                Stefano Corni                    **'
       write(*,*) '**                Silvio Pipolo                    **'
       write(*,*) '**               Emanuele Coccia                   **'
       write(*,*) '**                 Gabriel Gil                     **'
       write(*,*) '**                Jacopo Fregoni                   **'
       write(*,*) '**                                                 **'
       write(*,*) '*****************************************************'
       write(*,*) '*****************************************************'    
       write(*,*) ''

       !Namelist mol 
       call init_nml_general()
       read(*,nml=general)
       call write_nml_general()

       !Namelist field
       call init_nml_field()
       read(*,nml=field)
       if (npulse.gt.npulsemax) then
          write(*,*) 'Error: npulse', npulse,'larger than', npulsemax
#ifdef MPI
       call mpi_finalize(ierr_mpi)
#endif
          stop 
       endif
       call write_nml_field() 

       !Namelist sse
       call init_nml_sse()
       read(*,nml=sse)
       call write_nml_sse()

       !read gaussian output for CIS propagation
       call read_gau_out

       !Namelist spectra
       call init_nml_spectra()
       read(*,nml=spectra) 
       call write_nml_spectra() 

       if (Fdis(1:5).ne."nodis") call read_dis_params

       if (Fres.eq.'Yesr') call read_restart()

       if (Fdis(1:5).ne.'nodis'.or.Fexp.ne.'exp') then 
           Fabs(1:3)='non'
           write(*,*) 'Absorber switched off with SSE or full Euler'
       endif

       if (Fabs(1:3).eq.'abs') call read_ion_rate

       return

      end subroutine read_input
!
!
      subroutine read_gau_out
!------------------------------------------------------------------------
! @brief Read input files 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------

       implicit none

       integer(i4b) :: i,j
       real(dbl)    :: rtmp 
       character(4) :: junk

!    read ci energy
       open(7,file="ci_energy.inp",status="old")
       ! add ground state
       n_ci=n_ci+1
       n_ci_read=n_ci_read+1
       allocate (e_ci(n_ci))
       e_ci(1)=0.d0
       write (6,*) "Excitation energies read from input, in Hartree"
       do i=2,n_ci
         read(7,*) junk,junk,junk,e_ci(i)
         e_ci(i)=e_ci(i)*ev_to_au
         write(6,*) i-1,e_ci(i)
       enddo
       write (6,*) 
       close(7)
!    read transition dipoles (also for pcm: useful for analysis)
       open(7,file="ci_mut.inp",status="old")
       allocate (mut(3,n_ci,n_ci))
! SC 31/05/2016: now the GS data from gaussian is in debye and same format as others
!       read(7,*) junk,mut(1,1,1),mut(2,1,1),mut(3,1,1)
!       mut(:,1,1)=mut(:,1,1)*debye_to_au
       do i=1,n_ci_read
!         read(7,*) junk,mut(1,1,i),mut(2,1,i),mut(3,1,i)
        if (i.le.n_ci) then
         read(7,*)junk,junk,junk,junk,mut(1,1,i),mut(2,1,i),mut(3,1,i)
         mut(:,i,1)=mut(:,1,i)
        else
         read(7,*)
        endif
       enddo
       do i=2,n_ci_read
         do j=2,i   
          if (i.le.n_ci.and.j.le.n_ci) then
           read(7,*)junk,junk,junk,junk,mut(1,i,j),mut(2,i,j),mut(3,i,j)
           mut(:,j,i)=mut(:,i,j)
          else
           read(7,*)
          endif
         enddo
       enddo
!       write(6,*) "mut"
!       do i=1,n_ci
!        do j=1,n_ci
!         write(6,'(2i4,3f8.4)') i,j,mut(:,i,j)
!        enddo
!       enddo
       close(7)


!    read initial coefficients for the dynamics using the Slater determinants instead of the CIS_0 states.
       allocate (c_i(n_ci))
       if (Fres.eq.'Nonr') then
          open(7,file="ci_ini.inp",status="old")
          do i=1,n_ci
             read(7,*) rtmp 
             c_i(i) = dcmplx(rtmp,0.d0)   
          enddo
          c_i=c_i/sqrt(dot_product(c_i,c_i))

       !elseif (Fres.eq.'Yesr') then
       !   call checkfile('restart',7)
       !   read(7,*) junk 
       !   read(7,*) restart_t,restart_i,diff_step
       !   allocate(c_i_prev(n_ci),c_i_prev2(n_ci))
       !   write(*,*) ''
       !   write(*,*) 'Restart from time', restart_t
       !   write(*,*) ''
       !   read(7,*) junk 
       !   do i=1,n_ci
       !      read(7,*) c_i(i)
       !   enddo
       !   read(7,*) junk
       !   do i=1,n_ci
       !      read(7,*) c_i_prev(i)
       !   enddo
       !   read(7,*) junk
       !   do i=1,n_ci
       !      read(7,*) c_i_prev2(i)
       !   enddo
       !   if (Fdis(1:5).ne.'nodis') then 
       !      read(7,*) junk 
       !      read(7,*) restart_seed
       !      iseed=restart_seed
       !      write(*,*) ''
       !      write(*,*) 'Used iseed:', iseed
       !      write(*,*) ''
       !      read(7,*) junk
       !      read(7,*) n_jump
       !   endif
       !   read(7,*) junk 
       !   read(7,*) mu_i_prev(1),mu_i_prev(2),mu_i_prev(3)
       !   read(7,*) mu_i_prev2(1),mu_i_prev2(2),mu_i_prev2(3)
       !   read(7,*) mu_i_prev3(1),mu_i_prev3(2),mu_i_prev3(3)
       !   read(7,*) mu_i_prev4(1),mu_i_prev4(2),mu_i_prev4(3)
       !   read(7,*) mu_i_prev5(1),mu_i_prev5(2),mu_i_prev5(3)


         close(7) 

       endif

       return

      end subroutine

      subroutine read_restart()
!------------------------------------------------------------------------
! @brief Read restart file
!
! @date Created   : E. Coccia 18 Apr 2018
! Modified        :
!------------------------------------------------------------------------

          implicit none

          integer(i4b)  :: i,ii 
          character(4)  :: junk
          character(32) :: filename

#ifndef MPI
          myrank=0
#endif

          if (myrank+1.lt.10) then
             write(filename,'("restart",I1)') myrank+1 
          elseif (myrank+1.lt.100) then
             write(filename,'("restart",I2)') myrank+1
          elseif (myrank+1.lt.1000) then   
             write(filename,'("restart",I3)') myrank+1
          elseif (myrank+1.lt.10000) then
             write(filename,'("restart",I4)') myrank+1  
          elseif (myrank+1.lt.100000) then
             write(filename,'("restart",I5)') myrank+1
          endif

          ii=7+myrank

          call checkfile(filename,ii)
          read(ii,*) junk
          read(ii,*) restart_t,restart_i,diff_step
          allocate(c_i_prev(n_ci),c_i_prev2(n_ci))
          write(*,*) ''
          write(*,*) 'Restart from time', restart_t
          write(*,*) ''
          read(ii,*) junk
          do i=1,n_ci
             read(ii,*) c_i(i)
          enddo
          read(ii,*) junk
          do i=1,n_ci
             read(ii,*) c_i_prev(i)
          enddo
          read(ii,*) junk
          do i=1,n_ci
             read(ii,*) c_i_prev2(i)
          enddo
          if (Fdis(1:5).ne.'nodis') then
             read(ii,*) junk
             read(ii,*) restart_seed
             iseed=restart_seed
             write(*,*) ''
             write(*,*) 'Used iseed:', iseed
             write(*,*) ''
             read(ii,*) junk
             read(ii,*) n_jump
          endif
          read(ii,*) junk
          read(ii,*) mu_i_prev(1),mu_i_prev(2),mu_i_prev(3)
          read(ii,*) mu_i_prev2(1),mu_i_prev2(2),mu_i_prev2(3)
          read(ii,*) mu_i_prev3(1),mu_i_prev3(2),mu_i_prev3(3)
          read(ii,*) mu_i_prev4(1),mu_i_prev4(2),mu_i_prev4(3)
          read(ii,*) mu_i_prev5(1),mu_i_prev5(2),mu_i_prev5(3)

         close(ii)

         return

      end subroutine read_restart

!
      subroutine read_dis_params()
!------------------------------------------------------------------------
! @brief Read nonradiative and dephasing rates 
! Define spontaneous emission coefs from Einstein coefficients 
! Read phases randomly added during the propagation
!
! @date Created   : E. Coccia 21 Dec 2016
! Modified  :
!------------------------------------------------------------------------
     
       implicit none
        integer     :: i,j,k,idum,ierr0,ierr1,ierr2,ierr3,ierr4,err,kk   
        real(dbl)   :: term   
        real(dbl)   :: rdum
        real(dbl)   :: rx,ry,rz
        real(dbl)   :: ix,iy,iz

       open(8,file='nr_rate.inp',status="old",iostat=ierr0,err=100)
       open(9,file='de_rate.inp',status="old",iostat=ierr1,err=101)
       open(10,file='de_phase.inp',status="old",iostat=ierr2,err=102)
       open(11,file='sp_rate.inp',status="old",iostat=ierr3,err=103)

       nexc = n_ci - 1
       if (Fful.eq.'Yesf') then
          nrel = nexc*(nexc-1)/2
          allocate(irel(nrel,2))
          allocate(ik(nexc,nexc))
          call map_nrel(nexc,nrel,irel)
          ik=0
          kk=nexc
          do j=nexc,1,-1
             do k=j-1,1,-1
                kk=kk+1
                ik(k,j)=kk
             enddo
          enddo
       elseif (Fful.eq.'Nonf') then
          nrel = 0 
       endif
       nf=nexc+nrel 

       if (idep.ne.0.and.idep.ne.1) then
          write(*,*) 'Invalid value for idep, must be 0 or 1'
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif

       if (nr_typ.ne.0.and.nr_typ.ne.1) then
          write(*,*) 'Invalid value for nr_typ, must be 0 or 1'
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif 

       allocate(nr_gam(nf))
       if (idep.eq.0) then
          allocate(de_gam(nexc+1))
          allocate(delta(nexc+1))
       elseif (idep.eq.1) then
          allocate(de_gam(nexc))
          allocate(de_gam1(nexc+1))
       endif
       allocate(sp_gam(nf))
       allocate(tmom2(nf))
       allocate(sp_fact(nf))
       allocate(tomega(nf))

! Transition frequencies
       do i=1,nexc
          tomega(i) = e_ci(i+1)
       enddo
       if (Fful.eq.'Yesf') then 
          k=nexc
          do i=nexc,1,-1
             do j=i-1,1,-1
                k=k+1
                tomega(k) = abs(e_ci(i+1) - e_ci(j+1))
             enddo
          enddo
       endif

! Read dephasing terms
       do i=1,nexc
          read(9,*) idum, de_gam(i)
       enddo
       if (idep.eq.0) read(9,*) idum, de_gam(nexc+1)
! Read relaxation terms
       do i=1,nf
          read(8,*) idum, nr_gam(i)
          read(11,*) idum, sp_fact(i)
       enddo

! The sigma_z operator for dephasing has an extra factor 2
! If idep.eq.1 S_alpha = \sum_beta M(alpha, beta) |beta><beta|
! if alpha.eq. beta then M(alpha,beta) = -1
! otherwise M(alpha,beta) = 1
       if (idep.eq.1) then 
          de_gam=0.5d0*de_gam
          call define_gamma(de_gam,de_gam1,n_ci)
       endif
! Define spontaneous emission terms 
       term=4.d0/(3.d0*clight**3)
       do i=1,nf
          sp_gam(i) = sp_fact(i)*term*tomega(i)**3
       enddo

! Define tmom2 =  <i|x|j>**2 + <i|y|j>**2 + <i|z|j>**2  
       do i=1,nexc
          tmom2(i) = mut(1,1,i+1)**2 + mut(2,1,i+1)**2 + mut(3,1,i+1)**2
       enddo
       if (Fful.eq.'Yesf') then 
          k=nexc
          do i=nexc,1,-1
             do j=i-1,1,-1 
                k=k+1
                tmom2(k) = mut(1,i+1,j+1)**2 + mut(2,i+1,j+1)**2 + mut(3,i+1,j+1)**2
             enddo
          enddo
       endif

       if (Fmdm.eq.'nan') then
          open(7,file="ci_mut_np.inp",status="old",iostat=ierr4,err=104)
          allocate(mut_np2(nf,3))
          do i=1,nf
             read(7,*) rdum, rx, ry, rz, ix, iy, iz, rdum
             mut_np2(i,1) = rx**2 + ix**2
             mut_np2(i,2) = ry**2 + iy**2
             mut_np2(i,3) = rz**2 + iz**2
          enddo
         close(7)
         do i=1,nf
            tmom2(i) = tmom2(i) + mut_np2(i,1) + mut_np2(i,2) + mut_np2(i,3)
         enddo

         write(*,*)"Contribution to the transition dipole <i|mu|j> ", &
                   "from NP"

104      if (ierr4.ne.0) then
           write(*,*)
           write(*,*) 'Dissipation and NP:'
           write(*,*) 'An error occurred during reading ci_mut_np.inp'
           write(*,*)
#ifdef MPI
           call mpi_finalize(ierr_mpi)
#endif
           stop
         endif
       endif
! Read phase randomly added during the propagation
       if (idep.eq.0) then
          do i=1,nexc+1
             read(10,*) idum, delta(i)
          enddo 
       endif

100    if (ierr0.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading nr_rate.inp'
          write(*,*)
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif  

101    if (ierr1.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading de_rate.inp'
          write(*,*)
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif  
          stop
       endif  
     
102    if (ierr2.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading de_phase.inp'
          write(*,*)
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif

103    if (ierr3.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading sp_rate.inp'
          write(*,*)
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif

       close(8)
       close(9)
       close(10) 
       close(11)

       return

      end subroutine read_dis_params

      subroutine map_nrel(nexc,nrel,irel)
!------------------------------------------------------------------------
! @brief Map state pairs for intermediate relaxations 
!
! @date Created   : E. Coccia 10 Oct 2017
! Modified  :
!------------------------------------------------------------------------
     
       implicit none

       integer(i4b)                  :: i,j,k,kk
       integer(i4b),  intent(in)     :: nrel,nexc
       integer(i4b),  intent(inout)  :: irel(nrel,2)


       irel=0
       do j=nexc,1,-1
          do k=j-1,1,-1
             kk=kk+1
             irel(kk,1)=j
             irel(kk,2)=k
          enddo
       enddo 

       return

      end subroutine map_nrel

      subroutine deallocate_dis()
!------------------------------------------------------------------------
! @brief Deallocate arrays for the dissipation 
!
! @date Created   : E. Coccia 22 Dec 2016
! Modified  :
!------------------------------------------------------------------------

#ifndef MPI
       myrank=0
#endif

       deallocate(nr_gam)
       deallocate(de_gam)
       if (idep.eq.1) deallocate(de_gam1)
       deallocate(sp_gam)
       deallocate(sp_fact)
       deallocate(tmom2)
       deallocate(tomega)
       if (idep.eq.0) deallocate(delta)
       if (myrank.eq.0) then
          if (Fdis.ne."nodis".and.Fmdm.eq.'nan') deallocate(mut_np2)
          if (Fful.eq.'Yesf') then
             deallocate(irel)
          endif
       endif
       if (Fful.eq.'Yesf') deallocate(ik)
       if (Fabs(1:3).eq.'abs') deallocate(ion_rate)

       return

      end subroutine deallocate_dis

      subroutine define_gamma(de_gam,de_gam1,nci)
!------------------------------------------------------------------------
! @brief Move from "physical" gammas to gammas 
! for keeping the population constant
! for each trajectory  
!
! @date Created   : E. Coccia 24 Mar 2017
! Modified  :
! @param de_gam(:), de_gam1(:)  
!------------------------------------------------------------------------

       implicit none
       integer                :: i
       integer, intent(in)    :: nci
       real(dbl), intent(in)    :: de_gam(nci-1)
       real(dbl), intent(inout) :: de_gam1(nci)

       de_gam1(1) = 0.d0
       do i=1,nexc
          de_gam1(i+1) = de_gam(i) - de_gam1(1) 
       enddo

       return

      end subroutine define_gamma

      subroutine init_nml_general()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist general 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_ci_read,n_ci,mol_cc,n_f,medium,restart,full,propa,n_restart 
!        lsim,absorber,binary,ncit 
!------------------------------------------------------------------------

       ! Time step in the propagation (a.u.)
       dt=0.05
       ! Number of steps
       n_step=10000
       ! Frequency in writing output files
       n_out=1
       ! Molecular center of charge
       mol_cc=0.d0
       ! Integer to be appended in all output files
       n_f=1
       ! Vacuum calculation
       medium='vac'
       ! Restart
       restart='n'
       ! Full or only |e> -> |0> relaxation
       full='n'
       ! Propagation of the energy term
       propa='e'
       ! Frequency in writing restart file
       n_restart=10
       ! Dynamics length in case of restart ('n' n_step from input, 'y'
       ! n_step = n_step - restart_step)
       lsim='n'
       ! Ionization rates
       absorber='n'
       ! Binary output
       binary='n'
       ! Threshold value for doing matmul or explicit loop in prop()
       ncit=150

       return

      end subroutine init_nml_general

      subroutine init_nml_field()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist field 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dt,n_step,n_out,Ffld,t_mid,sigma,omega,radiative,iseed,fmax 
!------------------------------------------------------------------------

       ! Type of field
       Ffld='gau'
       ! Center of the pulse
       t_mid=200
       ! Sigma for the pulse
       sigma=0.d0
       ! Frequency (no oscillation)
       omega=0.d0
       ! Stochastic field
       radiative='non'
       ! ieed for radiative and/or dissipation
       iseed=12345678
       ! Amplitude of the external field
       fmax=0.d0
       ! Number of pulses
       npulse=1
       ! Time delay
       tdelay=0.d0
       ! Phase shift
       pshift=0.d0

       return

      end subroutine init_nml_field

      subroutine init_nml_spectra()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist spectra 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param start,tau,dir_ft 
!------------------------------------------------------------------------

       ! Parameter for computing spectra
       nspectra=1
       ! Parameter for computing spectra
       tau(:)=zero
       ! Direction fo the field (no field)
       dir_ft=0.d0

       return

      end subroutine init_nml_spectra

      subroutine init_nml_sse()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist sse 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd 
!------------------------------------------------------------------------

       ! Do dissipation 
       dissipative='non'
       ! Type of the dephasing operator
       idep=1
       ! Type of propagator
       dis_prop='qjump'
       ! If dis_prop='euler', steps for accumulating the Wiener process
       nrnd=1
       ! If dis_prop='euler', use the Euler-Maruyama algorithm
       tdis=0
       ! Type of nonradiative relaxation
       nr_typ=1
       ! Factor is dissipation='ernd'
       krnd=1.d0
       ! Output level
       out_sse='n'

       return

      end subroutine init_nml_sse


      subroutine write_nml_general()
!------------------------------------------------------------------------
! @brief Write variables in the namelist general and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_ci_read,n_ci,mol_cc,n_f,medium 
!------------------------------------------------------------------------

       if (n_ci.le.0) then
          write(*,*) 'ERROR: number of excited states is wrong', n_ci
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif 
       if (n_ci_read.le.0) then
          write(*,*) 'ERROR: number of excited states is wrong', n_ci
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif
       write(6,*) "Number to append to dat file", n_f
       write (*,*) "Number of CIS states to be read",n_ci_read
       write (*,*) "Number of CIS states to be used",n_ci
       write(*,*) 'Molecular center of charge'
       write(*,*) mol_cc 
       write(*,*) ''
       write (*,*) "Time step (in au), number of steps, stride",dt, &
                 n_step,n_out
       write(*,*) ''
       select case (medium)
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
          write(*,*) "No external medium, vacuum calculation"
          Fmdm='vac'
       end select
       select case (restart)
        case ('n', 'N')
          write(*,*) 'No restart'
          Fres='Nonr'
        case ('y', 'Y')
          write(*,*) 'Restart from a previous calculation'
          Fres='Yesr'
          write(*,*) 'Frequency in writing restart file:', n_restart
!EC 15/12/17: Restart only for vacuum calculations
          if (Fmdm.ne.'vac') Fres='Nonr'
       end select
       select case (full)
        case ('n', 'N') 
          write(*,*) 'Only |k> -> |0> relaxation'
          Fful='Nonf'
        case ('y', 'Y')
          write(*,*) 'Full relaxation'
          Fful='Yesf'
       end select
       select case (propa)
        case ('e', 'E')
          write(*,*) 'Energy term is propagated analytically'
          Fexp='exp'
        case ('n','N')
          write(*,*) 'Energy term is propagated via second-order Euler'
          Fexp='non' 
       end select
       select case (lsim)
        case ('y', 'Y')
         write(*,*) 'Effective number of steps is n_step - restart_step'
         Fsim='y'
        case ('n', 'N')
         Fsim='n'
       end select
       select case (absorber)
        case ('y', 'Y')
         write(*,*) 'Absorber in dynamics'
         Fabs(1:3)='abs'
        case ('n', 'N')
         Fabs(1:3)='non'
       end select
       select case (binary)
        case ('y','Y')
         write(*,*) 'Output files in binary format'
         Fbin(1:3)='bin'
        case ('n','N')
         Fbin(1:3)='non'
       end select
       if (n_ci.gt.ncit.and.nthreads.gt.1) then
           write(*,*) 'Explicit loops are used for matrix/vector'
           write(*,*) 'for propagation, if OMP is switched on.'
           Fopt(1:3)='omp'
       else
          write(*,*) 'Matmul is used in the propagation.'
          Fopt(1:3)='non'
       endif 
       write(*,*) ''

       return

      end subroutine write_nml_general

      subroutine write_nml_field()
!------------------------------------------------------------------------
! @brief Write variables in the namelist field and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dt,n_step,n_out,Ffld,t_mid,sigma,omega,radiative,iseed,fmax 
!------------------------------------------------------------------------

       integer :: i

       if (npulse.lt.1) then
           write(*,*) 'ERROR: number of pulses in input '
           write(*,*) 'is less than zero'
#ifdef MPI
           call mpi_finalize(ierr_mpi)
#endif
           stop
       endif

       write (*,*) "Time shape of the perturbing field",Ffld
       write (*,*) "time at the center of the pulse (au):",t_mid
       write (*,*) "Width of the pulse (time au):",sigma(1)
       write (*,*) "Frequency (au):",omega(1)
       write (*,*) "Maximum E field",fmax(:,1)
       !SC
       select case (radiative)
        case ('rad','Rad','RAD')
         Frad='arl'
        case default
         Frad='non'
       end select
       do i=2,npulse
          write(*,*) 'Pulse no.', i
          write(*,*) 'Frequency (au) =', omega(i)
          write(*,*) 'Sigma pulse =', sigma(i)
          write(*,*) 'Time delay (au) between pulse',i-1,'and pulse', i, '=', tdelay(i-1)
          write(*,*) 'Phase shift between pulse',i-1,'and pulse', i, '=', pshift(i-1)
          write(*,*) ''
       enddo
       write(*,*) ''

       return

      end subroutine write_nml_field

      subroutine write_nml_spectra()
!------------------------------------------------------------------------
! @brief Write variables in the namelist spectra and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param start,tau,dir_ft 
!------------------------------------------------------------------------
 
       if (medium.ne.'vac') then 
          nspectra=2
       endif

       write(*,*) 'Starting point for FT calculation', start
       write(*,*) 'Artificial damping', (tau(i),i=1,nspectra)
       write(*,*) 'Direction along which the field is oriented', dir_ft
       write(*,*) ''

       return

      end subroutine write_nml_spectra

      subroutine write_nml_sse()
!------------------------------------------------------------------------
! @brief Write variables in the namelist sse and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd 
!------------------------------------------------------------------------

       select case (dissipative)
        case ('mar', 'Mar', 'MAR')
          select case (idep)
           case (0)
            Fdis_deph="exp"
            write(*,*) 'Use sqrt(gamma_i) exp(i delta_i)|i><i|(i=0,nci)'
           case (1)
            Fdis_deph="i-0"
            write(*,*) 'Use sqrt(gamma_i) (|i><i|-|0><0|)'
          end select
          write(*,*) 'Markovian dissipation'
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             Fdis="mar-qjump"
             write(*,*) 'Quantum jump algorithm'
           case ('Euler', 'EUler', 'EULER', 'euler')
             write(*,*) 'Continuous stochastic propagator'
             write(*,*) 'Time step for the Brownian motion is:', dt/nrnd
             select case (tdis)
              case (0)
                write(*,*) 'Euler-Maruyama algorithm'
                Fdis="mar-EuMar"
              case (1)
                write(*,*) 'Leimkuhler-Matthews algorithm'
                Fdis="mar-LeiMa"
             end select
          end select
          select case (nr_typ)
           case (0)
            Fdis_rel="dip"
            write(*,*) 'Internal conversion relaxation via dipole'
           case (1)
            Fdis_rel="mat"
            write(*,*) 'Internal conversion relaxation via given matrix'
          end select
        case ('nma', 'NMa', 'NMA', 'Nma')
          write(*,*) 'NonMarkovian dissipation'
          read(*,*) dis_prop
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             Fdis="nma-qjump"
             write(*,*) 'Quantum jump algorithm'
           case default
             Fdis="nma-cstoc"
             write(*,*) 'Continuous stochastic propagator'
          end select
        case ('rnd', 'RND', 'Rnd', 'RNd')
          Fdis="ernd"
          write(*,*) "Normal random term added to CI energies"
        case default
          Fdis="nodis"
          write(*,*) 'No dissipation'
       end select
       select case (out_sse)
        case ('y','Y')
          write(*,*) 'SSE quantum jumps written in output'
          Fwrt(1:3)='yes'
        case ('n','N')
          write(*,*) 'SSE quantum jumps not written'
          Fwrt(1:3)='non' 
       end select
       write(*,*) ''

       return

      end subroutine write_nml_sse

      subroutine checkfile(filename,channel)
!------------------------------------------------------------------------
! @brief Check file existence 
!
! @date Created   : E. Coccia 20 Nov 2017
! Modified  :
!------------------------------------------------------------------------
        implicit none

        character*12,   intent(in)  :: filename
        integer(i4b),  intent(in)  :: channel   
        logical                    :: exist

        inquire(file=filename, exist=exist)
        if (exist) then
           open(channel, file=filename, status="old")
        else
           write(*,*) 'ERROR:  file', filename,'is missing'       
#ifdef MPI
           call mpi_finalize(ierr_mpi)
#endif
           stop 
        endif    

        return 

      end subroutine checkfile 

      subroutine mpibcast_readio()
!------------------------------------------------------------------------
! @brief MPI broadcast of input variables 
!
! @date Created   : E. Coccia 20 Apr 2018
! Modified  :
!------------------------------------------------------------------------

#ifdef MPI
       !n_f and iseed generated ad hoc for each process
       n_f=myrank+1
       call seed_random_number_sc(iseed+myrank)
       call random_number(r)
       iseed=floor(r*1000000000)

       !Broadcast to other processes
       call mpi_bcast(n_ci,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(n_ci_read, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(n_step,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(n_out,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(n_restart, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(npulse,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(idep,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(nrnd,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(tdis,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)    
       call mpi_bcast(nr_typ,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

       call mpi_bcast(dt,        1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(t_mid,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(krnd,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(start,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(tau,       2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mol_cc,    3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(dir_ft,    3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(sigma,     npulsemax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(omega,     npulsemax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(pshift,    npulsemax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(tdelay,    npulsemax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(fmax,      3*npulsemax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)

       call mpi_bcast(propa,       flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(lsim,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(medium,      flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(restart,     flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(full,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Ffld,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(radiative,   flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(dissipative, flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(dis_prop,    flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fmdm,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fres,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fful,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fexp,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fsim,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Frad,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fdis,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fdis_rel,    flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fdis_deph,   flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fabs,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fbin,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(Fopt,        flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
#endif

       return 

      end subroutine mpibcast_readio 

      subroutine mpibcast_e_dip()
!------------------------------------------------------------------------
! @brief MPI broadcast of energies, dipoles and initial coefs 
!
! @date Created   : E. Coccia 23 Apr 2018
! Modified  :
!------------------------------------------------------------------------

#ifdef MPI
       if (myrank.ne.0) then
          allocate(e_ci(n_ci))
          allocate(mut(3,n_ci,n_ci))
          allocate(c_i(n_ci))
       endif

       call mpi_bcast(e_ci,      n_ci,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mut,       3*n_ci*n_ci,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)

       call mpi_bcast(c_i,       2*n_ci,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr_mpi)
#endif

       return

      end subroutine mpibcast_e_dip

      subroutine mpibcast_sse()
!------------------------------------------------------------------------
! @brief MPI broadcast of SSE relaxation and dephasing rates 
!
! @date Created   : E. Coccia 23 Apr 2018
! Modified  :
!------------------------------------------------------------------------

#ifdef MPI
       call mpi_bcast(nexc,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(nf,         1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(nrel,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

       if (myrank.ne.0) then
          allocate(nr_gam(nf))
          if (idep.eq.0) then
             allocate(de_gam(nexc+1))
             allocate(delta(nexc+1))
          elseif (idep.eq.1) then
             allocate(de_gam(nexc))
             allocate(de_gam1(nexc+1))
          endif
          allocate(sp_gam(nf))
          allocate(tmom2(nf))
          allocate(sp_fact(nf))
          allocate(tomega(nf))
          if (Fful.eq.'Yesf') allocate(ik(nexc,nexc))
       endif

       call mpi_bcast(nr_gam,     nf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       if (idep.eq.0) then
          call mpi_bcast(de_gam,  nexc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
          call mpi_bcast(delta,   nexc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 
       elseif (idep.eq.1) then
          call mpi_bcast(de_gam,  nexc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
          call mpi_bcast(de_gam1, nexc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 
       endif
       call mpi_bcast(sp_gam,     nf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(sp_fact,    nf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(tmom2,      nf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 
       call mpi_bcast(tomega,     nf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 
       if (Fful.eq.'Yesf') call mpi_bcast(ik,nexc*nexc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi) 

       call mpi_bcast(Fwrt,       flg,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)

#endif

       return

      end subroutine mpibcast_sse


      subroutine mpibcast_restart()
!------------------------------------------------------------------------
! @brief MPI broadcast of restart variables 
!
! @date Created   : E. Coccia 23 Apr 2018
! Modified  :
!------------------------------------------------------------------------

#ifdef MPI

       if (myrank.ne.0) then
          allocate(c_i(n_ci))
          allocate(c_i_prev(n_ci))
          allocate(c_i_prev2(n_ci)) 
       endif

       call mpi_bcast(restart_i,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi) 
       call mpi_bcast(diff_step,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)  
       call mpi_bcast(n_jump,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi) 
       call mpi_bcast(restart_seed, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)  
       call mpi_bcast(iseed,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi) 

       call mpi_bcast(restart_t,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 
       call mpi_bcast(mu_i_prev,    3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mu_i_prev2,   3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mu_i_prev3,   3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mu_i_prev4,   3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(mu_i_prev5,   3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)

       call mpi_bcast(c_i,          2*n_ci,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(c_i_prev,     2*n_ci,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr_mpi)
       call mpi_bcast(c_i_prev2,    2*n_ci,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr_mpi)

#endif

       return

      end subroutine mpibcast_restart 
   
      subroutine mpibcast_ion_rate
!------------------------------------------------------------------------
! @brief Broadcast ion_rate 
!
! @date Created   : E. Coccia 31 May 2018
! Modified  :
!------------------------------------------------------------------------

#ifdef MPI

       if (myrank.ne.0) allocate(ion_rate(n_ci))

       call mpi_bcast(ion_rate,n_ci,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 

#endif

       return

      end subroutine mpibcast_ion_rate

 
      subroutine read_ion_rate()
!------------------------------------------------------------------------
! @brief Read ionization rates 
!
! @date Created   : E. Coccia 31 May 2018
! Modified  :
!------------------------------------------------------------------------

       implicit none

       integer(i4b) :: ierr5,err,idum,i
       real(dbl)    :: rdum

       open(12,file='ion_rate.dat',status="old",iostat=ierr5,err=105)

! Read ionization rates
! IJQC, 116, 1120 (2016), Eq. 14
! Obtained from calculation using Light code
! J. Chem. Phys., 139, 164121 (2013)
       allocate(ion_rate(n_ci))
       do i=1,n_ci
          read(12,*) idum,rdum, ion_rate(i)
       enddo
       close(12)

105    if (ierr5.ne.0) then
          write(*,*)
          write(*,*) 'Ionization:'
          write(*,*) 'An error occurred during reading ion_rate.dat'
          write(*,*)
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif

       return

      end subroutine read_ion_rate


    end module readio
