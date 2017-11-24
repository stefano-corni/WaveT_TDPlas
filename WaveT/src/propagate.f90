      module propagate
      use constants   
      use global_wavet
      use readio
      !use pedra_friends  
      use spectra
      use random
      use sse 

      implicit none

      integer(i4b)                :: ijump=0
      real(dbl),allocatable       :: f(:,:)
      complex(cmp), allocatable   :: c(:),c_prev(:),c_prev2(:),h_rnd(:,:), h_rnd2(:,:)
      real(dbl),     allocatable  :: h_int(:,:), h_dis(:)
      real(dbl),     allocatable  :: pjump(:)
      real(dbl)                   :: f_prev(3),f_prev2(3)
      real(dbl)                   :: mu_prev(3),mu_prev2(3),mu_prev3(3),&
                                    mu_prev4(3), mu_prev5(3)
      real(dbl),     allocatable  :: w(:), w_prev(:)
      real(dbl)                   :: eps
      logical                     :: first=.true.

! SC mu_a is the dipole moment at current step,
!    int_rad is the classical radiated power at current step
!    int_rad_int is the integral of the classical radiated power at current step
      real(dbl) :: mu_a(3),int_rad,int_rad_int
      integer(i4b) :: file_c=10,file_e=8,file_mu=9 
      integer(i4b) :: file_d=14
      save
      private
      public create_field, prop
!
      contains
!
      subroutine prop
!------------------------------------------------------------------------
! @brief Propogate C(t) using a second
! order Euler algorithm 
! 
! 
! @date Created   : 
! Modified  : E. Coccia Dec-Apr 2017
!------------------------------------------------------------------------

       implicit none
       integer(i4b)                :: i,j,k
       complex(cmp), allocatable   :: ccexp(:) !SC 31/10/17: added to store exp(-ui*e(:)*dt), used in propagation
       character(20)               :: name_f


! OPEN FILES
       if (Fres.eq.'Yesr') then
          write(name_f,'(a4,i0,a4)') "c_t_",n_f,".dat"
          open (file_c,file=name_f,status="unknown",access="append")
          write(name_f,'(a4,i0,a4)') "e_t_",n_f,".dat"
          open (file_e,file=name_f,status="unknown",access="append")
          write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
          open (file_mu,file=name_f,status="unknown",access="append")
          write(name_f,'(a4,i0,a4)') "d_t_",n_f,".dat"
          open (file_d,file=name_f,status="unknown",access="append")
       elseif (Fres.eq.'Nonr') then
          write(name_f,'(a4,i0,a4)') "c_t_",n_f,".dat"
          open (file_c,file=name_f,status="unknown")
          write(name_f,'(a4,i0,a4)') "e_t_",n_f,".dat"
          open (file_e,file=name_f,status="unknown")
          write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
          open (file_mu,file=name_f,status="unknown")
          write(name_f,'(a4,i0,a4)') "d_t_",n_f,".dat"
          open (file_d,file=name_f,status="unknown")
       endif
! ALLOCATING
       allocate (c(n_ci))
       allocate (c_prev(n_ci))
       allocate (c_prev2(n_ci))
       allocate (h_int(n_ci,n_ci))
       if (Fexp.eq."exp") then 
          allocate (ccexp(n_ci))
          ccexp=exp(-ui*dt*e_ci)
       endif
! SP 17/07/17: new flags
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          allocate(h_dis(n_ci))
          if (Fdis(5:9).eq."qjump") then
             allocate (pjump(2*nf+nexc+1))
          else
             allocate (h_rnd(n_ci,n_ci))
             allocate (h_rnd2(n_ci,n_ci))
             allocate (w(3*n_ci), w_prev(3*n_ci))
          endif
       endif

! STEP ZERO: build interaction matrices to do a first evolution
! Different initialization in case of restart
       if (Fres.eq.'Nonr') then
          c_prev2=c_i
          c_prev=c_i
          !f_prev2=f(:,2)
          f_prev2=f(:,1)
          f_prev=f(:,1)
          mu_prev=0.d0
          mu_prev2=0.d0
          mu_prev3=0.d0
          mu_prev4=0.d0
          mu_prev5=0.d0
          n_jump=0 
       elseif (Fres.eq.'Yesr') then
          c_prev2=c_i_prev2
          c_prev=c_i_prev
          !f_prev2=f(:,restart_i-1)
          !f_prev=f(:,restart_i)
          mu_prev=mu_i_prev
          mu_prev2=mu_i_prev2
          mu_prev3=mu_i_prev3
          mu_prev4=mu_i_prev4
          mu_prev5=mu_i_prev5 
       endif
       c=c_i
! SP: 17/07/17: changed the following f_prev2 <= f_prev
       h_int=zero  
       int_rad_int=0.d0
       if(Frad.eq."arl".or.Fdis.ne."nodis") &
                               call seed_random_number_sc(iseed)
       if (Fmdm(1:3).ne."vac") then
           call init_medium(c_prev,f_prev,h_int)
           if (Fres.eq.'Nonr') then
              i=1
              call prop_medium(i,c_prev,f_prev,h_int)
           endif
       endif
       if (Fres.eq.'Nonr') then
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
          call add_int_vac(f_prev,h_int)
! SP 16/07/17: added call to output at step 0 to have full output in outfiles
          call out_header
          call output(1,c,f_prev,h_int)
       elseif (Fres.eq.'Yesr') then
          if (Fdis.ne."nodis") call random_seq(restart_i)
       endif

! EC 20/12/16
! Dissipation according to the Markovian SSE (eq 25 J. Phys: Condens.
! Matter vol. 24 (2012) 273201)
! OR
! Dissipation according to the non-Markovian SSE (eq 24 J. Phys:
! Condens. Matter vol. 24 (2012) 273201) 
! Add a random fluctuation for the stochastic propagation (.not.qjump)
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          call define_h_dis(h_dis,n_ci)
          if (Fdis(5:9).ne."qjump") then 
             call rnd_noise(w,w_prev,n_ci,first)
             first=.false.
             call add_h_rnd(h_rnd,n_ci,w,w_prev) 
             call add_h_rnd2(h_rnd2,n_ci)
          endif
       endif

       if (Fexp.eq.'exp') then
! Energy term is propagated analytically
! Interaction term via second-order Euler
          call exp_euler_prop(ccexp,n_ci)
       elseif (Fexp.eq.'non') then
! Both energy and intercation terms are
! porpagated via second-order Euler
          call full_euler_prop(n_ci)
       endif


! DEALLOCATION AND CLOSING
       deallocate(c,c_prev,c_prev2,h_int)
       deallocate(c_i)
       if (Fres.eq."Yesr") deallocate(c_i_prev,c_i_prev2)
       if (Fexp.eq."exp") deallocate(ccexp)
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          call deallocate_dis()
          deallocate(h_dis)
          if (Fdis(5:9).ne."qjump") then
             deallocate(h_rnd)
             deallocate(h_rnd2)
             deallocate(w)
             deallocate(w_prev)
          endif 
          if (Fdis(5:9).eq."qjump") then
             write(*,*)
             write(*,*) 'Total number of quantum jumps',i_sp+i_nr+i_de,&
                        '(',  real(i_sp+i_nr+i_de)/real(n_step),')'
             write(*,*) 'Spontaneous emission quantum jumps', i_sp,    &
                        '(',  real(i_sp)/real(n_step),')'
             write(*,*) 'Nonradiarive relaxation quantum jumps', i_nr, &
                        '(', real(i_nr)/real(n_step),')'
             write(*,*) 'Dephasing relaxation quantum jumps', i_de,    &
                        '(', real(i_de)/real(n_step),')'
             write(*,*)
          endif
       endif
       close (file_c)
       close (file_e)
       close (file_mu)
       close (file_d)

       if(Fmdm(1:3).ne.'vac') call finalize_medium

       return

      end subroutine prop
!
      subroutine create_field
!------------------------------------------------------------------------
! @brief Create electric field 
! 
! 
! @date Created   : 
! Modified  : 
!------------------------------------------------------------------------

       implicit none

       integer(i4b) :: i,i_max,n_tot
       real(dbl) :: t_a,ti,tf,arg
       character(15) :: name_f

       if (Fres.eq.'Nonr') then
          n_tot=n_step
       elseif (Fres.eq.'Yesr') then
          n_tot=n_step+restart_i
       endif

       allocate (f(3,n_tot))

       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (7,file=name_f,status="unknown")
        f(:,:)=0.d0
        select case (Ffld)
        case ("mdg")
        ! Gaussian modulated sinusoid: exp(-(t-t0)^2/s^2) * sin(wt) 
         select case (npulse)
          case (1)
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
              if (mod(i,n_out).eq.0) &
              write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          case (2)
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
              f(:,i)=f(:,i) + fmax(:)*exp(-(t_a-(t_mid+tdelay))**2 &
                     /(sigma1**2))*sin(omega1*t_a+pshift) 
              if (mod(i,n_out).eq.0) &
              write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
         end select 
        case ("mds")
        ! Cosine^2 modulated sinusoid: 1/2* cos^2(pi(t-t0)/(2t0)) * sin(wt) 
        !          f=0 for t>t0
         i_max=int(t_mid/dt)
         do i=1,2*i_max
          t_a=dt*(dble(i)-1)
          f(:,i)=fmax(:)*cos(pi*(t_a-t_mid)/(2*t_mid))**2/2.d0* &
                 sin(omega*t_a)
         enddo
         do i=2*i_max+1,n_tot
          t_a=dt*(i-1)
          f(:,i)=0.
         enddo
        case ("pip")
        ! Pi pulse: cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)) 
         i_max=int(t_mid/dt)
         do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=0.
          if (abs(t_a-t_mid).lt.sigma) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(2*sigma)))**2* &
               cos(omega*(t_a-t_mid))
          endif
         enddo
        case ("sin")
        ! Sinusoid:  sin(wt) 
         do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=fmax(:)*sin(omega*t_a)
         enddo
        case ("snd")
        ! Linearly modulated (up to t0) Sinusoid:
        !         0 < t < t0 : t/to* sin(wt) 
        !             t > t0 :       sin(wt) 
         do i=1,n_tot
          t_a=dt*(dble(i)-1)
          if (t_a.gt.t_mid) then
            f(:,i)=fmax(:)*sin(omega*t_a)
          else
            if (t_a.gt.zero) f(:,i)=fmax(:)*t_a/t_mid*sin(omega*t_a)
          endif            
         enddo
        case ("gau")
        ! Gaussian pulse: exp(-(t-t0)^2/s^2) 
          select case (npulse)
          case (1)
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-pt5*(t_a-t_mid)**2/(sigma**2))
              if (mod(i,n_out).eq.0) &
                write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          case(2)
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-pt5*(t_a-t_mid)**2/(sigma**2))
              f(:,i)=f(:,i)+fmax(:)*exp(-pt5*(t_a-(t_mid+tdelay))**2/(sigma**2))
              if (mod(i,n_out).eq.0) &
                write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          end select
! SP 270817: the following (commented) is probably needed for spectra  
         !do i=1,n_tot
         ! t_a=dt*(dble(i)-1)
         ! arg=-pt5*(t_a-t_mid)**2/(sigma**2)
         ! if(arg.lt.-50.d0) then 
         !   f(:,i)=zero
         ! else
         !   f(:,i)=fmax(:)*exp(arg)
         ! endif
         !enddo
        case ("css")
        ! Cos^2 pulse (only half a period): cos^2(pi*(t-t0)/(s)) 
         ti=t_mid-sigma/two
         tf=t_mid+sigma/two
         do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=zero
          if (t_a.gt.ti.and.t_a.le.tf) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(sigma)))**2
          endif
         enddo
        case default
         write(*,*)  "Error: wrong field type !"
         stop
        end select
        ! write out field 
        do i=1,n_tot
         t_a=dt*(i-1)
         if (mod(i,n_out).eq.0) &
           write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
        enddo

        close(7)

        return

      end subroutine create_field

!
      subroutine do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
!------------------------------------------------------------------------
! @brief Compute C^T mu C and save previous dipoles 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------

       implicit none

       complex(cmp), intent(IN) :: c(n_ci)
       real(dbl)::mu_prev(3),mu_prev2(3),mu_prev3(3),mu_prev4(3), &
                mu_prev5(3)

       mu_a(1)=dot_product(c,matmul(mut(1,:,:),c))
       mu_a(2)=dot_product(c,matmul(mut(2,:,:),c))
       mu_a(3)=dot_product(c,matmul(mut(3,:,:),c))
! SC save previous mu for radiative damping
       mu_prev5=mu_prev4
       mu_prev4=mu_prev3
       mu_prev3=mu_prev2
       mu_prev2=mu_prev
       mu_prev=mu_a

       return
 
      end subroutine do_mu

      subroutine output(i,c,f_prev,h_int)     
!------------------------------------------------------------------------
! @brief Write output files 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------

       implicit none

       integer(i4b), intent(IN) :: i
       complex(cmp), intent(IN) :: c(n_ci)
       real(dbl), intent(IN) :: h_int(n_ci,n_ci)
       real(dbl), intent(IN) :: f_prev(3)
       real(dbl) :: e_a,e_vac,t,g_neq_t,g_neq2_t,g_eq_t,f_med(3)
       character(4000) :: fmt_ci,fmt_ci2
       integer(i4b)    :: itmp,j

       t=(i-1)*dt 
       e_a=dot_product(c,e_ci*c+matmul(h_int,c))
! SC 07/02/16: added printing of g_neq, g_eq 
       if(Fmdm(1:3).ne.'vac') then 
          g_eq_t=e_a
          g_neq_t=e_a
          g_neq2_t=e_a
          e_vac=e_a
          call get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)
          write (file_e,'(i8,f14.4,7e20.8)') i,t,e_a,e_vac, &
                   g_eq_t,g_neq2_t,g_neq_t,int_rad,int_rad_int
       else
          write (file_e,'(i8,f14.4,3e22.10)') i,t,e_a,int_rad,int_rad_int
       endif
       itmp = n_ci*(n_ci-1)/2
       write (fmt_ci,'("(i8,f14.4,",I0,"e13.5)")') n_ci
       write (fmt_ci2,'("(i8,f14.4,",I0,"2e13.5)")') 2*itmp
! SP 10/07/17: commented the following, strange error 
       call wrt_decoherence(i,t,file_d,fmt_ci2,c,n_ci)
       write (file_c,fmt_ci) i,t,real(c(:)*conjg(c(:)))
       write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
       j=int(dble(i)/dble(n_out))
       if(j.lt.1) j=1
       Sdip(:,1,j)=mu_a(:)
! SP 270817: using get_* functions to communicate with TDPlas
       if(Fmdm(1:3).ne."vac") call get_medium_dip(Sdip(:,2,j))
       Sfld(:,j)=f(:,i)

       return

      end subroutine output


      subroutine add_int_vac(f_prev,h_int)
!------------------------------------------------------------------------
! @brief Create the field term of the hamiltonian 
!
! @date Created   : 
! Modified  : E. Coccia 22/11/2017
!------------------------------------------------------------------------

       implicit none

       real(dbl), intent(IN) :: f_prev(3)
       real(dbl), intent(INOUT) :: h_int(n_ci,n_ci)
! SC 16/02/2016: changed to - sign, 

       h_int(:,:)=h_int(:,:)-mut(1,:,:)*f_prev(1)-             &
                   mut(2,:,:)*f_prev(2)-mut(3,:,:)*f_prev(3)
       return
 
      end subroutine add_int_vac

      subroutine add_int_rad(mu_prev,mu_prev2,mu_prev3,mu_prev4, &
                                                   mu_prev5,h_int)
!------------------------------------------------------------------------
! @brief Calculate the Aharonov Lorentz radiative damping 
!
! @date Created   : S. Corni 
! Modified  : E. Coccia 22/11/2017
!------------------------------------------------------------------------

       implicit none

       real(dbl),intent(in) :: mu_prev(3),mu_prev2(3),mu_prev3(3), &
                 mu_prev4(3), mu_prev5(3)
       real(dbl), intent(INOUT) :: h_int(n_ci,n_ci)
       real(dbl) :: d3_mu(3),d2_mu(3),d_mu(3),d2_mod_mu,coeff,scoeff, &
                  force(3),de


!       d3_mu=(mu_prev-3.*mu_prev2+3.*mu_prev3-mu_prev4)/(dt*dt*dt)
       d3_mu=(2.5*mu_prev-9.*mu_prev2+12.*mu_prev3-7.*mu_prev4+ &
              1.5*mu_prev5)/(dt*dt*dt)
! SC 08/06/2016: i'm confused on the right formula for istantaneous radiated intensity:
!  Novotny-Hecht is pro. to (d^2/dt^2 |mu|)^2 (eq. 8.70), however in Jackson
! for Larmor one has pro. to (d^2/dt^2 mu \cdot d^2/dt^2 mu). Since Novotny-Hech
! in eq. 8.82 seems to use jacson definition, I also use it here
!       d2_mod_mu=(2.*sqrt(dot_product(mu_prev,mu_prev)) &
!                  -5.*sqrt(dot_product(mu_prev2,mu_prev2))+ &
!                  4.*sqrt(dot_product(mu_prev3,mu_prev3)) &
!               -sqrt(dot_product(mu_prev4,mu_prev4)))/(dt*dt)
       d2_mu=(2.*mu_prev-5.*mu_prev2+4.*mu_prev3-mu_prev4)/(dt*dt)
       d_mu=(1.5*mu_prev-2.*mu_prev2+0.5*mu_prev3)/dt
!SC coefficient 1/(6 pi eps0 c^3) in atomic units
       coeff=2.d0/3.d0/137.036**3.
       scoeff=coeff
!SC: added random force             
!       d3_mu=d3_mu*coeff
       force(1)=random_normal()*scoeff
       force(2)=random_normal()*scoeff
       force(3)=random_normal()*scoeff
       d3_mu=d3_mu*coeff
       de=dot_product(d_mu,force)
       if(de.lt.0) d3_mu=d3_mu+force
!       write (6,*) d3_mu
!SC Instantenous emitted intensity (from Novotny Hech eq. 8.70)
!       int_rad=coeff*d2_mod_mu*d2_mod_mu
       int_rad=coeff*dot_product(d2_mu,d2_mu)+de/dt
       int_rad_int=int_rad_int+int_rad*dt
       h_int(:,:)=h_int(:,:)-mut(1,:,:)*d3_mu(1)-             &
                   mut(2,:,:)*d3_mu(2)-mut(3,:,:)*d3_mu(3)

       return

      end subroutine add_int_rad

      subroutine out_header
!------------------------------------------------------------------------
! @brief Write headers to output files 
!
! @date Created   : S. Corni 
! Modified  : E. Coccia 22/11/2017
!------------------------------------------------------------------------

       implicit none

       write(file_c,'(2a)')'# istep   time (au)', &
                           '    Population 0,1,...,n_ci'

       write(file_e,'(8a)') '#   istep time (au)',' <H(t)>-E_gs(0)', &
              ' DE_vac(t)',' DG_eq(t)',' DG_neq(t)',  '  Const', &
              '  Rad. Int', '  Rad. Ene'   
      
       write(file_mu,'(5a)') '#   istep time (au)',' dipole-x ', &
              ' dipole-y ',' dipole-z '
       
       write(file_d,'(3a)')'# istep   time (au)', & 
                           '    Coherence (0,i=1,n_ci)', &
                           '    Coherence (j=1,n_ci,k=j+1,n_ci)'

       return
    
      end subroutine out_header

      subroutine wrt_decoherence(i,t,int1,char2,c,nci)
!------------------------------------------------------------------------
! @brief Print the tridiagional C*_iC_j (i.ne.j) matrix 
! corresponding to the decoherence 
! 
! @date Created   : E. Coccia 3 Feb 2017
! Modified  :
!------------------------------------------------------------------------
  
        implicit none

        integer(i4b),    intent(in)    :: i, nci
        integer(i4b),    intent(in)    :: int1 
        character(4000), intent(in)    :: char2 
        real(dbl),       intent(in)    :: t
        complex(cmp),    intent(in)    :: c(nci)
        complex(cmp)                   :: cc(nci,nci)
        integer(i4b)                   :: j, k
        complex(cmp)                   :: tmp

        tmp=cmplx(0.d0,0.d0)

        do k=2,nci
           tmp = conjg(c(k))*c(1)
           cc(1,k) = tmp
        enddo  

        if (nci.gt.2) then
           do j=2,nci
              do k=j+1,nci
                 tmp = conjg(c(k))*c(j)
                 cc(j,k) = tmp
              enddo
           enddo
        endif

        if (nci.gt.2) then
           write(int1,char2) i, t, ((cc(j,k), k=j+1,nci),j=1,nci-1)
        else 
           write(int1,char2) i, t, cc(1,2)
        endif
    
        return

      end subroutine wrt_decoherence


      subroutine exp_euler_prop(ccexp,nci)
!------------------------------------------------------------------------
! @brief Energy term is propagated analytically
! Interaction term via second-order Euler 
! 
! @date Created   : E. Coccia 15 Nov 2017
! Modified  :
!------------------------------------------------------------------------

        implicit none

        integer(i4b),  intent(in)  :: nci
        complex(cmp),  intent(in)  :: ccexp(nci) 

        integer(i4b)               :: i,j,istart,iend
        real(dbl)                  :: t 
        complex(cmp)               :: dis(nci)

! Initialization only without restart
       if (Fres.eq.'Nonr') then
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
! SC: 31/10/17 modified the propagation with ccexp 
          c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev))
          if (Fdis.eq."ernd") then
             do j=1,n_ci
                c(j) = c(j) - ccexp(j)*ui*dt*krnd*random_normal()*c_prev(j)
             enddo
          endif
          if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
             dis=disp(h_dis,c_prev,nci)
             c=c-ccexp*dt*dis
             if (Fdis(5:9).eq."EuMar") then
          ! Euler-Maruyama
                c=c-ccexp*(ui*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev))               
             elseif (Fdis(5:9).eq."LeiMa") then
          ! Leimkuhler-Matthews
                c=c-ccexp*(ui*0.5d0*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev))
             endif
          endif
          c=c/sqrt(dot_product(c,c))
          c_prev=c

! SP 16/07/17: added call to medium propagation at step 2 to have full output
          if (Fmdm(1:3).ne."vac") then
             i=2
             call prop_medium(i,c_prev,f_prev,h_int)
          endif
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
! SP 16/07/17: heder called at step 1                                        
          !call out_header
          if (mod(2,n_out).eq.0) call output(2,c,f_prev,h_int)
       endif

       if (Fres.eq.'Nonr') then
          istart=3
          iend=n_step
       elseif (Fres.eq.'Yesr') then
          istart=restart_i+1
          iend=n_step+istart-1
       endif

! PROPAGATION CYCLE: starts the propagation at timestep 3
! without restar, at timestep restart_i+1 otherwise
! Markovian dissipation (quantum jump) 
       if (Fdis(5:9).eq."qjump") then
          !do i=3,n_step
          do i=istart,iend 
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            dis=disp(h_dis,c_prev,nci)
            if (i.eq.ijump+1) then
               c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev)-dt*dis)
            else 
               c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev)-2.d0*dt*dis)
            endif 
! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c_prev,n_ci,pjump)
            call random_number(eps)  
            if (dtot.gt.eps)  then
               call quan_jump(c,c_prev,n_ci,pjump)
               ijump=i
               n_jump=n_jump+1
               write(*,*) 'Quantum jump at step:', i, (i-1)*dt 
               c_prev=c
            else
                c=c/sqrt(dot_product(c,c))
                c_prev2=c_prev
                c_prev=c
            endif
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            endif
          enddo
! Markovian dissipation (Euler-Maruyama) 
       elseif (Fdis(1:3).eq."mar") then
          !do i=3,n_step
          do i=istart,iend
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! Dissipation by a continuous stochastic propagation
            call rnd_noise(w,w_prev,n_ci,first)
            call add_h_rnd(h_rnd,n_ci,w,w_prev)
            dis=disp(h_dis,c_prev,nci)
            if (Fdis(5:9).eq."EuMar") then
            ! Euler-Maruyama 
                c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev)-dt* &
                dis-ui*sqrt(dt)*matmul(h_rnd,c_prev)- &
                dt*matmul(h_rnd2,c_prev))
            elseif (Fdis(5:9).eq."LeiMa") then
                c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev)-dt* &
                dis-ui*0.5d0*sqrt(dt)* &
                matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev))
            endif
            !c=c/sqrt(dot_product(c,c))
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5) 
            endif
          enddo
       elseif (Fdis.eq."nodis".or.Fdis.eq."ernd") then
! No dissipation in the propagation 
          !do i=3,n_step
          do i=istart,iend
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! SC 31/10/17: modified propagation by adding the exp term
            c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev))
            if (Fdis.eq."ernd") then
               do j=1,n_ci
                  c(j) = c(j) - 2.d0*ccexp(j)*ui*dt*krnd*random_normal()*c_prev(j)
               enddo
            endif
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5) 
            endif
          enddo
       endif 

       return

      end subroutine exp_euler_prop


      subroutine full_euler_prop(nci)
!------------------------------------------------------------------------
! @brief Energy and interaction terms are propagated
! via second-order Euler 
! 
! @date Created   : E. Coccia 15 Nov 2017
! Modified  :
!------------------------------------------------------------------------      

        implicit none

        integer(i4b), intent(in)  :: nci

        integer(i4b)              :: i,j,istart,iend
        real(dbl)                 :: t
        complex(cmp)              :: dis(nci)

!Initialization only without restart
       if (Fres.eq.'Nonr') then
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
          c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
          if (Fdis.eq."ernd") then
             do j=1,n_ci
                c(j) = c(j) - ui*dt*krnd*random_normal()*c_prev(j)
             enddo
          endif
          if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
             dis=disp(h_dis,c_prev,nci)
             c=c-dt*dis
             if (Fdis(5:9).eq."EuMar") then
          ! Euler-Maruyama
                c=c-ui*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)               
             elseif (Fdis(5:9).eq."LeiMa") then
          ! Leimkuhler-Matthews
                c=c-ui*0.5d0*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)
             endif
          endif
          c=c/sqrt(dot_product(c,c))
          c_prev=c

! SP 16/07/17: added call to medium propagation at step 2 to have full
! output
          if (Fmdm(1:3).ne."vac") then
             i=2
             call prop_medium(i,c_prev,f_prev,h_int)
          endif
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
          if (mod(2,n_out).eq.0) call output(2,c,f_prev,h_int)
       endif

       if (Fres.eq.'Nonr') then
          istart=3
          iend=n_step
       elseif (Fres.eq.'Yesr') then
          istart=restart_i+1
          iend=n_step+istart-1
       endif 

! PROPAGATION CYCLE: starts the propagation at timestep 3
! Markovian dissipation (quantum jump) -> qjump
       if (Fdis(5:9).eq."qjump") then
          !do i=3,n_step
          do i=istart,iend
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                               mu_prev4,mu_prev5,h_int)
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            dis=disp(h_dis,c_prev,nci)
            if (i.eq.ijump+1) then
               c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt*dis
            else
               c=c_prev2-2.d0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-2.d0*dt*dis
            endif
! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c_prev,n_ci,pjump)
            call random_number(eps)
            if (dtot.gt.eps)  then
               call quan_jump(c,c_prev,n_ci,pjump)
               ijump=i
               n_jump=n_jump+1
               write(*,*) 'Quantum jump at step:', i, (i-1)*dt
               c_prev=c
            else
                c=c/sqrt(dot_product(c,c))
                c_prev2=c_prev
                c_prev=c
            endif
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5) 
            endif
          enddo
! Markovian dissipation (Euler-Maruyama) 
       elseif (Fdis(1:3).eq."mar") then
          !do i=3,n_step
          do i=istart,iend
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! Dissipation by a continuous stochastic propagation
            call rnd_noise(w,w_prev,n_ci,first)
            call add_h_rnd(h_rnd,n_ci,w,w_prev)
            dis=disp(h_dis,c_prev,nci)
            if (Fdis(5:9).eq."EuMar") then
            ! Euler-Maruyama 
              c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt* &
                dis-ui*sqrt(dt)*matmul(h_rnd,c_prev)- &
                dt*matmul(h_rnd2,c_prev)
            elseif (Fdis(5:9).eq."LeiMa") then
              c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt* &
                dis-ui*0.5d0*sqrt(dt)* & 
                matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)
            endif
            !c=c/sqrt(dot_product(c,c))
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5) 
            endif
          enddo
       elseif (Fdis.eq."nodis".or.Fdis.eq."ernd") then
! No dissipation in the propagation 
          !do i=3,n_step
          do i=istart,iend
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
            c=c_prev2-2.d0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
            if (Fdis.eq."ernd") then
               do j=1,n_ci
                  c(j) = c(j) - 2.d0*ui*dt*krnd*random_normal()*c_prev(j)
               enddo
            endif
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,n_ci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5) 
            endif 
          enddo
       endif

       return
      
      end subroutine full_euler_prop

      subroutine wrt_restart(i,t,c,c_prev,c_prev2,nci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
!------------------------------------------------------------------------
! @brief Write restart file 
! 
! @date Created   : E. Coccia 21 Nov 2017
! Modified  :
!------------------------------------------------------------------------      
  
       implicit none

       integer(i4b),  intent(in)  :: i,nci,iseed
       real(dbl),     intent(in)  :: t
       real(dbl),     intent(in)  :: mu_prev(3),mu_prev2(3),mu_prev3(3)
       real(dbl),     intent(in)  :: mu_prev4(3),mu_prev5(3)
       complex(cmp),  intent(in)  :: c(nci), c_prev(nci), c_prev2(nci)

       integer(i4b)               :: j

       open(777,file='restart')
       rewind(777)
       
       write(777,*) 'Restart time in au, Restart step'
       write(777,*) t,i
       write(777,*) 'Coefficients'
       do j=1,nci
          write(777,*) c(j)
       enddo
       write(777,*) 'Coefficients -1'
       do j=1,nci
          write(777,*) c_prev(j)
       enddo
       write(777,*) 'Coefficients -2'
       do j=1,nci
          write(777,*) c_prev2(j)
       enddo
       if (Fdis(1:5).ne.'nodis') then
          write(777,*) 'Seed'
          write(777,*) iseed
          write(777,*) 'Number of quantum jumps'
          write(777,*) n_jump
       endif
       write(777,*) 'Dipoles'
       write(777,*) mu_prev(1), mu_prev(2), mu_prev(3)
       write(777,*) mu_prev2(1), mu_prev2(2), mu_prev2(3)
       write(777,*) mu_prev3(1), mu_prev3(2), mu_prev3(3)
       write(777,*) mu_prev4(1), mu_prev4(2), mu_prev4(3)
       write(777,*) mu_prev5(1), mu_prev5(2), mu_prev5(3)


       close(777)
 
       !flush(6)

       return

      end subroutine wrt_restart 

      end module


