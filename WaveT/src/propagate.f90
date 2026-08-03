module propagate
      use constants
      use readio
      use spectra
      use random
      use dissipation
      use scf
      use interface_classic
      use initialise
      use WTMathTools
#ifdef OMP
      use omp_lib
#endif
#ifdef MPI
      use mpi
#endif

      implicit none

      integer(i4b)                :: ijump=0
      real(dbl),     allocatable  :: f(:,:)
      real(dbl),     allocatable  :: avec(:,:)
      complex(cmp),  allocatable  :: trans_p(:,:,:) ! Added by Manuel Sanchez 2026-04-22
      complex(cmp),  allocatable  :: h_int_vg(:,:) ! Added by Manuel Sanchez 2026-04-22
      complex(cmp),  allocatable  :: c(:),c_prev(:),c_prev2(:),h_rnd(:,:), h_rnd2(:,:)
      real(dbl),     allocatable  :: h_int(:,:), h_dis(:), gamma_sum(:), Rp(:,:), Rn(:,:), gamma_nr(:,:)
      real(dbl),     allocatable  :: pjump(:)
      real(dbl)                   :: f_prev(3),f_prev2(3)
      real(dbl)                   :: mu_prev(3),mu_prev2(3),mu_prev3(3),&
                                    mu_prev4(3), mu_prev5(3)
      complex(cmp)                :: m_prev(3),m_prev2(3),m_prev3(3),&
                                     m_prev4(3), m_prev5(3) ! MM - test
      real(dbl),     allocatable  :: w(:), w_prev(:)
      real(dbl)                   :: eps
      logical                     :: first=.true.

! SC mu_a is the dipole moment at current step,
!    int_rad is the classical radiated power at current step
!    int_rad_int is the integral of the classical radiated power at current step
      real(dbl) :: int_rad,int_rad_int,mu_a(3),sm,mu_a_esa(3)
      complex(cmp) :: m_a(3)
      integer(i4b) :: file_c=10,file_e=8,file_mu=9,file_m=611 !MM 
      save
      private
      public create_field, create_vector_potential, build_trans_p_from_dipoles, build_h_int_vg, prop, create_2d_map, print_time
!
      contains
!
!------------------------------------------------------------------------
! @brief Propogate C(t) using a second
! order Euler algorithm 
! 
! 
! @date Created   : 
! Modified  : E. Coccia Dec-Apr 2017
! Modified  : Manuel Sanchez 02/04/2026
!------------------------------------------------------------------------
      subroutine prop

       implicit none
       integer(i4b)                :: i,j,k
       complex(cmp), allocatable   :: ccexp(:) !SC 31/10/17: added to store exp(-ui*e(:)*dt), used in propagation
       character(20)               :: name_e,name_c,name_d,name_mu, &
                                      name_m ! MM
       !MR
       real :: start, finish

! OPEN FILES
       write(name_c,'(a4,i0,a4)') "c_t_",n_f,".dat"
       write(name_e,'(a4,i0,a4)') "e_t_",n_f,".dat"
       write(name_mu,'(a5,i0,a4)') "mu_t_",n_f,".dat"
       if (Fmag.eq.'mag') then
          write(name_m,'(a4,i0,a4)') "m_t_",n_f,".dat"
       endif
       if (Fres.eq.'Yesr') then
          if (Fbin.ne.'bin') then
             open (file_c,file=name_c,status="unknown",access="append")
             open (file_e,file=name_e,status="unknown",access="append")
             open (file_mu,file=name_mu,status="unknown",access="append")
             if (Fmag.eq.'mag') then !MM
              open(file_m,file=name_m,status="unknown",access="append")
             endif
          else
             open (file_c,file=name_c,status="unknown",access="append",form="unformatted")   
             open (file_e,file=name_e,status="unknown",access="append",form="unformatted")
             open (file_mu,file=name_mu,status="unknown",access="append",form="unformatted")  
             if (Fmag.eq.'mag') then !MM
              open(file_m,file=name_m,status="unknown",access="append",form="unformatted")
             endif
          endif
       elseif (Fres.eq.'Nonr') then
             if (Fbin.ne.'bin') then
                open (file_c,file=name_c,status="unknown")
                open (file_e,file=name_e,status="unknown")
                open (file_mu,file=name_mu,status="unknown")
                if (Fmag.eq.'mag') then !MM
                 open(file_m,file=name_m,status="unknown")
                endif
             else
                open(file_c,file=name_c,status="unknown",form="unformatted")
                open(file_e,file=name_e,status="unknown",form="unformatted")
              open(file_mu,file=name_mu,status="unknown",form="unformatted")
                if (Fmag.eq.'mag') then !MM
                   open(file_m,file=name_m,status="unknown",form="unformatted")
                endif
             endif
       else
          !write(name_c,'(a9,i0,a1,i0,a4)') "mu_t_esa_",nmap,"_",n_f,".dat"
          write(name_mu,'(a5,i0,a1,i0,a4)') "mu_t_",nmap,"_",n_f,".dat"
          if (Fbin.ne.'bin') then
                !open (file_c,file=name_c,status="unknown")
                open (file_mu,file=name_mu,status="unknown")
          else
             !open(file_c,file=name_c,status="unknown",form="unformatted")
             open(file_mu,file=name_mu,status="unknown",form="unformatted")
          endif
       endif
! ALLOCATING
       allocate (c(nstates))
       allocate (c_prev(nstates))
       allocate (c_prev2(nstates))
       allocate (h_int(nstates,nstates))
       allocate (h_int_vg(nstates,nstates)) ! Added by Manuel Sanchez 2026-04-22
       h_int_vg = zeroc ! Added by Manuel Sanchez 2026-04-22
       if (Fexp.eq."exp") then 
          allocate (ccexp(nstates))
          ccexp=exp(-ui*dt*energies)
          if (Fabs.eq.'abs') ccexp=ccexp*exp(-ion_rate*dt/2.d0) 
       endif
! SP 17/07/17: new flags
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          allocate(h_dis(nstates))
          allocate(gamma_sum(nstates))
          allocate(Rp(nstates,nstates))
          allocate(Rn(nstates,nstates))
          allocate(gamma_nr(nstates,nstates))
          if (Fdis(5:9).eq."qjump") then
             allocate (pjump(2*nf+nexc+1))
          else
             allocate (h_rnd(nstates,nstates))
             allocate (h_rnd2(nstates,nstates))
             allocate (w(3*nstates), w_prev(3*nstates))
          endif
       endif

! STEP ZERO: build interaction matrices to do a first evolution
! Different initialization in case of restart
       if (Fres.eq.'Nonr') then
          c=coeff0
          c_prev2=coeff0
          c_prev=coeff0
          f_prev2=f(:,1)
          f_prev=f(:,1)
          ! SP 18/05/20: shouldn't we initialize mu_prev* to mu at time zero?
          mu_prev=0.d0
          mu_prev2=0.d0
          mu_prev3=0.d0
          mu_prev4=0.d0
          mu_prev5=0.d0
          n_jump=0 
       elseif (Fres.eq.'Yesr') then
          c=c_i_t
          c_prev2=c_i_prev2
          c_prev=c_i_prev
          f_prev2=f(:,restart_i-1)
          f_prev=f(:,restart_i)
          mu_prev=mu_i_prev
          mu_prev2=mu_i_prev2
          mu_prev3=mu_i_prev3
          mu_prev4=mu_i_prev4
          mu_prev5=mu_i_prev5 
          if (Fmag.eq.'mag') then
             m_prev=m_i_prev
             m_prev2=m_i_prev2
             m_prev3=m_i_prev3
             m_prev4=m_i_prev4
             m_prev5=m_i_prev5
          endif
       endif
       h_int=zero  
       int_rad_int=0.d0
       if(Frad.eq."arl".or.Fdis.ne."nodis") &
                               call seed_random_number_sc(iseed)
       !> Initialize medium for propagation
       if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") then
           call init_env_prop(c_prev,mu_prev,f_prev,h_int)
       endif
       if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
          if (allocated(trans_p)) deallocate(trans_p) 
          allocate(trans_p(3,nstates,nstates)) 
          call build_trans_p_from_dipoles(trans_p)
       endif 
       if (Fres.eq.'Nonr') then
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)          
          if (Fmag.eq.'mag') then ! MM
             call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
          endif
          if (gauge.ne.'vg') call add_int_vac(f_prev,h_int) ! Added by Manuel Sanchez 2026-05-03
! SP 16/07/17: added call to output at step 0 to have full output in outfiles
          if (Fbin.ne.'bin'.and.twod.eq.'no') call out_header
          call output(1,c,f_prev,h_int)
       elseif (Fres.eq.'Yesr') then
          if (Fdis.ne."nodis") call random_seq(restart_i)
       endif
       ! SP 19/06/26: added call for ropagation at step 1
       if (Fres.ne.'Yesr') then
         if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") then
           i=1
           !write(*,*) "mu_prev ", mu_prev
           !stop
           call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
         endif
       endif
! EC 20/12/16
! Dissipation according to the Markovian SSE (eq 25 J. Phys: Condens.
! Matter vol. 24 (2012) 273201)
! OR
! Dissipation according to the non-Markovian SSE (eq 24 J. Phys:
! Condens. Matter vol. 24 (2012) 273201) 
! Add a random fluctuation for the stochastic propagation (.not.qjump)
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          call define_h_dis(h_dis,nstates)
          if (Fdis(5:9).eq."EuMar".or.Fdis(5:9).eq."RuKu4".or.Fdis(5:9).eq."HeuSt") then
             call build_gamma_nr_matrix(gamma_nr,nstates)
             call sqrt_gamma_nr_matrix(gamma_nr,Rp,nstates)
             call build_gamma_sum_from_gamma_nr(gamma_nr,gamma_sum,nstates)
          endif
          !if (Fdis(5:9).ne."qjump") then 
          !   call rnd_noise(w,w_prev,nstates,first)
          !   first=.false.
          !   call add_h_rnd(h_rnd,nstates,w,w_prev) 
          !   call add_h_rnd2(h_rnd2,nstates)
          !endif
       endif


       if (n_step.gt.1) then
         if (Fexp.eq.'exp') then
! Energy term is propagated analytically
! Interaction term via second-order Euler
           call exp_euler_prop(ccexp,nstates)
         elseif (Fexp.eq.'non') then
! Both energy and intercation terms are
! porpagated via second-order Euler
           call full_euler_prop(nstates)
         endif
       endif


! DEALLOCATION AND CLOSING
       deallocate(c,c_prev,c_prev2,h_int,f)
       deallocate(h_int_vg) ! Added by Manuel Sanchez 2026-04-22
       if (allocated(avec)) deallocate(avec)
       if (allocated(trans_p)) deallocate(trans_p) ! Added by Manuel Sanchez 2026-04-22
       if (Fres.eq."Yesr") deallocate(c_i_prev,c_i_t,c_i_prev2)
       if (Fexp.eq."exp") deallocate(ccexp)
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then

          !call deallocate_dis()
          deallocate(h_dis)
          deallocate(gamma_sum)
          deallocate(Rp)
          deallocate(Rn)
          deallocate(gamma_nr)

          if (Fdis(5:9).ne."qjump") then
             deallocate(h_rnd)
             deallocate(h_rnd2)
             deallocate(w)
             deallocate(w_prev)
          else
             deallocate(pjump)
          endif

#ifndef MPI
          if (Fdis(5:9).eq."qjump") then
             write(*,*)
             write(*,*) 'Total number of quantum jumps',i_sp+i_nr+i_de,&
                        '(',  real(i_sp+i_nr+i_de)/real(n_step),')'
             write(*,*) 'Spontaneous emission quantum jumps', i_sp,    &
                        '(',  real(i_sp)/real(n_step),')'
             write(*,*) 'Nonradiative relaxation quantum jumps', i_nr, &
                        '(', real(i_nr)/real(n_step),')'
             write(*,*) 'Dephasing relaxation quantum jumps', i_de,    &
                        '(', real(i_de)/real(n_step),')'
             write(*,*)
          endif
#endif
       endif
       if (twod.eq.'no') then
         close (file_c)
         close (file_e)
         close (file_mu)
         close (file_m)
       else
         !close(file_c)
         close(file_mu)
       endif

       if(Fmdm.ne.'vac') call finalize_medium

       return

      end subroutine prop
!
!------------------------------------------------------------------------
! @brief Create electric field 
! 
! 
! @date Created   : 
! Modified  : E. Coccia 16 Jan 2018
!------------------------------------------------------------------------
      subroutine create_field(f00)

       implicit none
       real(dbl), intent(out) :: f00(3)

       integer(i4b) :: i,j,i_max,n_tot
       real(dbl) :: t_a,ti,tf,arg
       character(15) :: name_f

       if (Fres.eq.'Nonr') then
          n_tot=n_step
       elseif (Fres.eq.'Yesr') then
          if (Fsim.eq.'y') then
             n_tot=diff_step+restart_i
          elseif (Fsim.eq.'n') then 
             n_tot=n_step+restart_i
          endif
       endif

       allocate (f(3,n_tot))
       if (twod.eq.'no') then
#ifndef MPI
          myrank=0
          write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
          if (Fbin.ne.'bin') then
             open (7,file=name_f,status="unknown")
          else
          open (7,file=name_f,status="unknown",form="unformatted")
       endif  
#endif
#ifdef MPI
       if (myrank.eq.0) then
          write(name_f,'(a9)') "field.dat"
             if (Fbin.ne.'bin') then
                open (7,file=name_f,status="unknown")
             else
                open (7,file=name_f,status="unknown",form="unformatted")
             endif
       endif
#endif
        endif

        f(:,:)=0.d0
        if (Flig.eq.'lin') then
           select case (Ffld)
           case("tdg")
        ! Gaussian modulated cosine with phase for twod
           do i=1,n_tot
               t_a=dt*(i-1)
               f(:,i) = fmax(:,1)*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))* &
                       cos(omega(1)*(t_a-t_mid)-pshift(1))
               do j=2,npulse
                  f(:,i) = f(:,i) + fmax(:,j)*                                   &
                  exp(-pt5*(t_a-(t_mid+sum(tdelay(1:j-1))))**2/(sigma(j)**2))*   &
                  cos(omega(j)*(t_a-(t_mid+sum(tdelay(1:j-1))))-pshift(j))
               enddo
           enddo
           case ("mdg")
        ! Gaussian modulated sinusoid: exp(-(t-t0)^2/s^2) * sin(wt) 
             do i=1,n_tot
                t_a=dt*(i-1)
                f(:,i) = fmax(:,1)*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))* & 
                         sin(omega(1)*t_a)
                do j=2,npulse
                   f(:,i) = f(:,i) + fmax(:,j)*                   &
                          exp(-pt5*(t_a-(t_mid+sum(tdelay(1:j-1))))**2/(sigma(j)**2))*   &
                          sin(omega(j)*t_a+sum(pshift(1:j-1)))
                enddo 
             enddo
            case ("snc")             
        ! Sinc apodized pulse: sin(Tt/2)/(Tt/2) * sin(wt) *
        ! exp(-Gamma*abs(t))         
             do i=1,n_tot
                t_a=dt*(i-1)      
               f(:,i) = fmax(:,1) * &
                        (sin(0.5*sigma(1)*(t_a - t_mid + 0.5*dt)) / &
                        (0.5*sigma(1)*(t_a - t_mid + 0.5*dt))) * &
                        sin(omega(1)*t_a) * &
                        exp(-(2.d0/t_ap)*abs(t_a - t_mid + 0.5*dt))
                do j=2,npulse
                  f(:,i) = f(:,i) + fmax(:,j) * &
                           (sin(0.5*sigma(j)*(t_a -(t_mid + sum(tdelay(1:j-1))) + 0.5*dt)) / &
                           (0.5*sigma(j)*(t_a - (t_mid + sum(tdelay(1:j-1))) + 0.5*dt))) * &
                           sin(omega(j)*t_a + sum(pshift(1:j-1))) * &
                           exp(-(2.d0/t_ap)*abs(t_a - (t_mid + sum(tdelay(1:j-1))) + 0.5*dt))
                enddo 
             enddo
           case ("mds")
        ! Cosine^2 modulated sinusoid: 1/2* cos^2(pi(t-t0)/(2t0)) * sin(wt) 
        !          f=0 for t>t0
            i_max=int(t_mid/dt)
            if (2*i_max.gt.n_tot) then
              write(*,*) 'ERROR: 2*t_mid/dt must be smaller than', n_tot
#ifdef MPI
              call mpi_finalize(ierr_mpi)  
#endif
              stop
            endif
            do i=1,2*i_max
               t_a=dt*(dble(i)-1)
               f(:,i)=fmax(:,1)*cos(pi*(t_a-t_mid)/(2*t_mid))**2* &
                     sin(omega(1)*t_a)
            enddo
            do i=2*i_max+1,n_tot
               t_a=dt*(i-1)
               f(:,i)=0.
            enddo
         !do j=2,npulse
         !   i_max=int((t_mid+sum(tdelay(1:j-1)))/dt)
         !   do i=1,2*i_max
         !      t_a=dt*(dble(i)-1)
         !      f(:,i)=f(:,i)+fmax(:,j)*cos(pi*(t_a-(t_mid+sum(tdelay(1:j-1))))/ &
         !             (2.d0*(t_mid+sum(tdelay(1:j-1)))))**2/2.d0* &
         !             sin(omega(j)*t_a+sum(pshift(1:j-1)))
         !   enddo
         !   do i=2*i_max+1,n_tot
         !      t_a=dt*(i-1)
         !      f(:,i)=0.
         !   enddo
         !enddo
         case ("pip")
        ! Pi pulse: cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)) 
           do i=1,n_tot
              t_a=dt*(dble(i)-1)
              f(:,i)=0.d0
               if (abs(t_a-t_mid).lt.sigma(1)) then
                  f(:,i)=fmax(:,1)*(cos(pi*(t_a-t_mid)/(2*sigma(1))))**2* &
                  cos(omega(1)*(t_a-t_mid))
               endif
           enddo
           do j=2,npulse
              do i=1,n_tot
                 t_a=dt*(dble(i)-1)
                 if (abs(t_a-(t_mid+sum(tdelay(1:j-1)))).lt.sigma(j)) then
                    f(:,i)=f(:,i)+fmax(:,j)*(cos(pi*(t_a-(t_mid+sum(tdelay(1:j-1))))/ &
                    (2*sigma(j))))**2* &
                    cos(omega(j)*(t_a-(t_mid+sum(tdelay(1:j-1))))+sum(pshift(1:j-1)))
                 endif
              enddo
           enddo
           case ("sin")
        ! Sinusoid:  sin(wt) 
           do i=1,n_tot
              t_a=dt*(dble(i)-1)
              f(:,i)=fmax(:,1)*sin(omega(1)*t_a)
           enddo
         
         case ("snd")
        ! Linearly modulated (up to t0) Sinusoid:
        !         0 < t < t0 : t/to* sin(wt) 
        !             t > t0 :       sin(wt) 
           do i=1,n_tot
            t_a=dt*(dble(i)-1)
            if (t_a.gt.t_mid) then
               f(:,i)=fmax(:,1)*sin(omega(1)*t_a)
            else
               if (t_a.gt.zero) f(:,i)=fmax(:,1)*t_a/t_mid*sin(omega(1)*t_a)
            endif            
         enddo
        case ("gau")
        ! Gaussian pulse: exp(-(t-t0)^2/s^2) 
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:,1)*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))
              do j=2,npulse
                 f(:,i)=f(:,i) + fmax(:,j)*exp(-pt5*(t_a- &
                 (t_mid+sum(tdelay(1:j-1))))**2/(sigma(j)**2)) 
              enddo
           enddo
        !EC 150424
        case ("tra")
        ! Trapezoidal pulse
           do i=1,n_tot
              if (i.gt.pini.and.i.lt.pfin) f(:,i)=fmax(:,1)
           enddo  
           ! Linear Gaussian pulse: exp(-(t-t0)^2/s^2)
        case ("lga")
           do i=1,n_tot
              t_a=dt*(i-1)
              f(:,i)=fmax(:,1)*t_a*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))
              do j=2,npulse
                 f(:,i)=f(:,i) + fmax(:,j)*t_a*exp(-pt5*(t_a- &
                 (t_mid+sum(tdelay(1:j-1))))**2/(sigma(j)**2))
              enddo
           enddo
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
         ti=t_mid-sigma(1)/two
         tf=t_mid+sigma(1)/two
         do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=zero
          if (t_a.gt.ti.and.t_a.le.tf) then
            f(:,i)=fmax(:,1)*(cos(pi*(t_a-t_mid)/(sigma(1))))**2
          endif
         enddo
        case ("sta")
         ! static field
          do i=1,n_tot
           f(:,i)=fmax(:,1)
          enddo
        case default
         write(*,*)  "Error: wrong field type !"
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
        end select

        elseif (Flig.eq.'cir') then
           if (e_dir(1).ne.0) then
           ! Light in the yz plane
              do i=1,n_tot
                 t_a=dt*(i-1)
                 f(2,i)=f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*cos(omega(1)*t_a)
                 f(3,i)=-f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*sin(omega(1)*t_a)
              enddo
           elseif (e_dir(2).ne.0) then
           ! Light in the xz plane
              do i=1,n_tot
                 t_a=dt*(i-1)
                 f(1,i)=f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*cos(omega(1)*t_a)
                 f(3,i)=-f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*sin(omega(1)*t_a)
              enddo
           elseif (e_dir(3).ne.0) then
           ! Light in the xy plane
              do i=1,n_tot
                 t_a=dt*(i-1)
                 f(1,i)=f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*cos(omega(1)*t_a)
                 f(2,i)=-f0*exp(-pt5*(t_a-t_mid)**2/(sigma(1)**2))*sin(omega(1)*t_a)
              enddo
           endif
        endif


        if (myrank.eq.0) then
        ! write out field 
          if (twod.eq.'no') then
             if (Fbin.ne.'bin') then
                do i=1,n_tot
                   t_a=dt*(i-1)
                   if (mod(i,n_out).eq.0) &
                       write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
                enddo
             else
                do i=1,n_tot
                   t_a=dt*(i-1)
                   if (mod(i,n_out).eq.0) write (7) t_a,f(:,i)
                enddo
             endif 
          endif
        endif
       

        close(7)
        f00(:)=f(:,1) 
        return

      end subroutine create_field

!------------------------------------------------------------------------
! @brief Create vector potential from electric field
! 
! 
! @date Created   : 
! Modified  : M. Sanchez 21 Apr 2026
!------------------------------------------------------------------------
      subroutine create_vector_potential

       implicit none

       integer(i4b) :: i, n_tot
       real(dbl)    :: t_a
       character(25) :: name_a

       if (.not.allocated(f)) then
          write(*,*) 'ERROR: electric field not available. Call create_field first.'
#ifdef MPI
          call mpi_finalize(ierr_mpi)
#endif
          stop
       endif

       n_tot = size(f,2)

       if (allocated(avec)) deallocate(avec)
       allocate(avec(3,n_tot))
       avec(:,1)= -pt5*dt*f(:,1) !zero
       do i=2,n_tot
          ! Length/velocity-gauge convention: E(t) = -dA(t)/dt
          avec(:,i)=avec(:,i-1)-pt5*dt*(f(:,i-1)+f(:,i))
       enddo

       if (twod.eq.'no') then
#ifndef MPI
          myrank=0
          write(name_a,'(a17,i0,a4)') "vector_potential",n_f,".dat"
          if (Fbin.ne.'bin') then
             open (77,file=name_a,status="unknown")
          else
             open (77,file=name_a,status="unknown",form="unformatted")
          endif
#endif
#ifdef MPI
          if (myrank.eq.0) then
             write(name_a,'(a20)') "vector_potential.dat"
             if (Fbin.ne.'bin') then
                open (77,file=name_a,status="unknown")
             else
                open (77,file=name_a,status="unknown",form="unformatted")
             endif
          endif
#endif
          if (myrank.eq.0) then
             if (Fbin.ne.'bin') then
                do i=1,n_tot
                   t_a=dt*(i-1)
                   if (mod(i,n_out).eq.0) &
                       write (77,'(f12.2,3e22.10e3)') t_a,avec(:,i)
                enddo
             else
                do i=1,n_tot
                   t_a=dt*(i-1)
                   if (mod(i,n_out).eq.0) write (77) t_a,avec(:,i)
                enddo
             endif
          endif
          close(77)
       endif

       return
 
      end subroutine create_vector_potential

!
!------------------------------------------------------------------------
! @brief Build momentum-like transition matrix trans_p from trans_dipoles
!        using i*p_ab = E_ab*mu_ab, with E_ab = energies(a)-energies(b).
!        trans_p is antisymmetric in (a,b) when trans_dipoles is symmetric.
!
! @date Created   : Manuel Sanchez 22/04/2026
! Modified  : Manuel Sanchez 03/05/2026 (upper-triangle build + antisymmetry)
!------------------------------------------------------------------------
      subroutine build_trans_p_from_dipoles(trans_p)

       implicit none

       complex(cmp), intent(out) :: trans_p(3,nstates,nstates)
       integer(i4b)              :: icart, ia, ib
       real(dbl)                 :: eab

       trans_p = zeroc
       do ia=1,nstates-1
          do ib=ia+1,nstates
             eab = energies(ib)-energies(ia)
             do icart=1,3
                trans_p(icart,ia,ib)=(-ui)*eab*trans_dipoles(icart,ia,ib) ! Added by Manuel Sanchez 2026-05-03
                trans_p(icart,ib,ia)=-trans_p(icart,ia,ib) ! Added by Manuel Sanchez 2026-05-03
             enddo
          enddo
       enddo

       return
 
      end subroutine build_trans_p_from_dipoles

!------------------------------------------------------------------------
! @brief Compute C^T mu C and save previous dipoles 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------
      subroutine do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)

       implicit none

       complex(cmp), intent(IN) :: c(nstates)
       real(dbl)                :: mu_prev(3),mu_prev2(3),mu_prev3(3),mu_prev4(3), &
                                   mu_prev5(3)
       complex(cmp)             :: ctmp(nstates),c_esa(nstates-1)
       integer(i4b)             :: j,k

#ifdef OMP
       if (Fopt.eq.'omp') then
          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_dipoles(1,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          mu_a(1)=dot_product(c,ctmp)

          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_dipoles(2,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          mu_a(2)=dot_product(c,ctmp)

          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_dipoles(3,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          mu_a(3)=dot_product(c,ctmp)
       else
          mu_a(1)=dot_product(c,matmul(trans_dipoles(1,:,:),c))
          mu_a(2)=dot_product(c,matmul(trans_dipoles(2,:,:),c))
          mu_a(3)=dot_product(c,matmul(trans_dipoles(3,:,:),c))
       endif
#endif
#ifndef OMP
       mu_a(1)=dot_product(c,matmul(trans_dipoles(1,:,:),c))
       mu_a(2)=dot_product(c,matmul(trans_dipoles(2,:,:),c))
       mu_a(3)=dot_product(c,matmul(trans_dipoles(3,:,:),c))
#endif
       if (twod.eq.'yes') then
           c_esa = c(2:)
           mu_a_esa(1)=dot_product(c_esa,matmul(trans_dipoles(1,2:,2:),c_esa))
           mu_a_esa(2)=dot_product(c_esa,matmul(trans_dipoles(2,2:,2:),c_esa))
           mu_a_esa(3)=dot_product(c_esa,matmul(trans_dipoles(3,2:,2:),c_esa))
       endif
! SC save previous mu for radiative damping
       mu_prev5=mu_prev4
       mu_prev4=mu_prev3
       mu_prev3=mu_prev2
       mu_prev2=mu_prev
       mu_prev=mu_a

       return
 
      end subroutine do_mu

!------------------------------------------------------------------------
! @brief print time in a separate file for 2D calculation
!
! @date Created   : G.Dall'Osto 30/04/2025
! Modified  :
!------------------------------------------------------------------------
      subroutine print_time
           integer(i4b)                :: i

           if (Fbin.ne.'bin') then
               open(22,file="time.dat",status="unknown")
               do i=1,n_step
                  write(22,*) (i-1)*dt
               enddo
           else
            open(22,file="time.dat",status="unknown",form="unformatted")
               do i=1,n_step
                  write(22) (i-1)*dt
               enddo
           endif

           close(22)

      end subroutine print_time

!------------------------------------------------------------------------
! @brief create dipoles along the ks=-k1+k2+k3 and ks=+k1-k2+k3 direction 
!        for 2D calculation
!
! @date Created   : G.Dall'Osto 30/04/2025
! Modified  :
!------------------------------------------------------------------------
      subroutine create_2d_map
           integer(i4b)                :: j, i, idum
           real(dbl)                   :: rdum, mu_read(3),mu_esa_read(3)
           character(20)               :: name_mu,name_c,name_esa,name_all
           complex(cmp),allocatable    :: mu_phase(:,:), mu_esa_phase(:,:)

      allocate(mu_phase(n_step,2), mu_esa_phase(n_step,2))
      mu_phase = 0
      mu_esa_phase = 0
      write(name_esa,'(a7,i0,a4)') "mu_esa_",n_f,".dat"
      write(name_all,'(a7,i0,a4)') "mu_all_",n_f,".dat"
      do j=1,12
        !write(name_c,'(a9,i0,a1,i0,a4)') "mu_t_esa_",j,"_",n_f,".dat"
        write(name_mu,'(a5,i0,a1,i0,a4)') "mu_t_",j,"_",n_f,".dat"
        if (Fbin.ne.'bin') then
          !open (file_c,file=name_c,status="unknown")
          open (file_mu,file=name_mu,status="unknown")
          !read(file_mu,*)
           !read(file_c,*)
          do i=1,n_step
             read(file_mu,*) idum, rdum, mu_read
             !read(file_c,*) idum, rdum, mu_esa_read
             mu_phase(i,1) = mu_phase(i,1) + (mu_read(1) + mu_read(2) + mu_read(3))*mat_c_inv(j,1)
             mu_phase(i,2) = mu_phase(i,2) + (mu_read(1) + mu_read(2) + mu_read(3))*mat_c_inv(j,2)
             !mu_esa_phase(i,1) = mu_esa_phase(i,1) + (mu_esa_read(1) + mu_esa_read(2) + mu_esa_read(3))*mat_c_inv(j,1)
             !mu_esa_phase(i,2) = mu_esa_phase(i,2) + (mu_esa_read(1) + mu_esa_read(2) + mu_esa_read(3))*mat_c_inv(j,2)
          enddo
          close(file_c, status='delete')
          close(file_mu, status='delete')
        else
          open(file_c,file=name_c,status="unknown",form="unformatted")
          open(file_mu,file=name_mu,status="unknown",form="unformatted")
          do i=1,n_step
             read(file_mu) idum, rdum, mu_read
             read(file_c) idum, rdum, mu_esa_read
             mu_phase(i,1) = mu_phase(i,1) + (mu_read(1) + mu_read(2) + mu_read(3))*mat_c_inv(j,1)
             mu_phase(i,2) = mu_phase(i,2) + (mu_read(1) + mu_read(2) + mu_read(3))*mat_c_inv(j,2)
             !mu_esa_phase(i,1) = mu_esa_phase(i,1) + (mu_esa_read(1) + mu_esa_read(2) + mu_esa_read(3))*mat_c_inv(j,1)
             !mu_esa_phase(i,2) = mu_esa_phase(i,2) + (mu_esa_read(1) + mu_esa_read(2) + mu_esa_read(3))*mat_c_inv(j,2)
          enddo
          !close(file_c, status='delete')
          close(file_mu, status='delete')
        endif
      enddo

      if (Fbin.ne.'bin') then
          open(20,file=name_all,status="unknown")
          !open(21,file=name_esa,status="unknown")
          do i=1,n_step
             write(20,*) real(mu_phase(i,1)), aimag(mu_phase(i,1)), &
                         real(mu_phase(i,2)), aimag(mu_phase(i,2))
             !write(21,*) real(mu_esa_phase(i,1)), aimag(mu_esa_phase(i,1)),&
             !            real(mu_esa_phase(i,2)), aimag(mu_esa_phase(i,2))
          enddo
      else
          open(20,file=name_all,status="unknown",form="unformatted")
          !open(21,file=name_esa,status="unknown",form="unformatted")
          do i=1,n_step
             write(20) mu_phase(i,:)
          !   write(21) mu_esa_phase(i,:)
         enddo
      endif
      deallocate(mu_phase, mu_esa_phase)
      close(20)
      !close(21)

      end subroutine

!------------------------------------------------------------------------
! @brief Compute C^T m C and save previous dipoles
!
! @date Created   : M. Monti 21/07/2022 MM
! Modified  : 
!------------------------------------------------------------------------
      subroutine do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)

       implicit none

       complex(cmp), intent(IN) :: c(nstates)
       complex(cmp)             :: m_prev(3),m_prev2(3),m_prev3(3),m_prev4(3), &
                                   m_prev5(3)
       complex(cmp)             :: ctmp(nstates)
       complex(cmp)             :: trans_mag_cmp(3,nstates,nstates)
       integer(i4b)             :: j,k

       trans_mag_cmp=cmplx(0.d0,trans_mag)
#ifdef OMP
       if (Fopt.eq.'omp') then
          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp)
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_mag_cmp(1,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          m_a(1)=dot_product(c,ctmp)

          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp)
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_mag_cmp(2,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          m_a(2)=dot_product(c,ctmp)

          ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp)
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ trans_mag_cmp(3,k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          m_a(3)=dot_product(c,ctmp)
       else
          m_a(1)=dot_product(c,matmul(trans_mag_cmp(1,:,:),c))
          m_a(2)=dot_product(c,matmul(trans_mag_cmp(2,:,:),c))
          m_a(3)=dot_product(c,matmul(trans_mag_cmp(3,:,:),c))
       endif
#endif
#ifndef OMP
       m_a(1)=dot_product(c,matmul(trans_mag_cmp(1,:,:),c))
       m_a(2)=dot_product(c,matmul(trans_mag_cmp(2,:,:),c))
       m_a(3)=dot_product(c,matmul(trans_mag_cmp(3,:,:),c))

#endif
! SC save previous mu for radiative damping
       m_prev5=m_prev4
       m_prev4=m_prev3
       m_prev3=m_prev2
       m_prev2=m_prev
       m_prev=m_a

!EC: scalar product between electric and magnetic dipole
       sm = dot_product(mu_a,m_a)

       return

      end subroutine do_m
!------------------------------------------------------------------------
! @brief Write output files 
!
! @date Created   : 
! Modified  : E. Coccia 20/11/2017
!------------------------------------------------------------------------
      subroutine output(i,c,f_prev,h_int)

       implicit none

       integer(i4b),    intent(IN) :: i
       complex(cmp),    intent(IN) :: c(nstates)
       real(dbl),       intent(IN) :: h_int(nstates,nstates)
       real(dbl),       intent(IN) :: f_prev(3)
       real(dbl)                   :: e_a,e_vac,t,g_neq_t,g_neq2_t,g_eq_t,f_med(3)
       character(4000)             :: fmt_ci
       integer(i4b)                :: itmp,j,k
       complex(cmp)                :: ctmp(nstates)

       t=(i-1)*dt
       if (twod.eq.'no') then 
#ifdef OMP
         if (Fopt.eq.'omp') then
            ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
          do k=1,nstates
             do j=1,nstates
                ctmp(k)=ctmp(k)+ h_int(k,j)*c(j)
             enddo
          enddo
!$OMP END PARALLEL
          e_a=dot_product(c,energies*c+ctmp)
         else      
          e_a=dot_product(c,energies*c+matmul(h_int,c))
         endif
#endif

#ifndef OMP
         e_a=dot_product(c,energies*c+matmul(h_int,c))
#endif

! SC 07/02/16: added printing of g_neq, g_eq 
         if(Fmdm.ne.'vac') then 
           g_eq_t=e_a
           g_neq_t=e_a
           g_neq2_t=e_a
           e_vac=e_a
           call get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)
           if (Fbin.ne.'bin') then
              write (file_e,'(i8,f14.4,7e20.8)') i,t,e_a,e_vac, &
                    g_eq_t,g_neq2_t,g_neq_t,int_rad,int_rad_int
           else
              write (file_e) i,t,e_a,e_vac, &
                    g_eq_t,g_neq2_t,g_neq_t,int_rad,int_rad_int
           endif
         else
           if (Fbin.ne.'bin') then
              write (file_e,'(i8,f14.4,3e22.10)') i,t,e_a,int_rad,int_rad_int
           else
              write (file_e) i,t,e_a,int_rad,int_rad_int
           endif 
         endif

         if (Fbin.ne.'bin') then
            write (fmt_ci,'("(i8,f14.4,",I0,"e17.8E3)")') 2*nstates
            write (file_c,fmt_ci) i,t,c(:)
            write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
            if (Fmag.eq.'mag') then
               write (file_m,'(i8,f14.4,7e22.10)') i,t,dble(m_a(1)),aimag(m_a(1)),&
                    dble(m_a(2)),aimag(m_a(2)),dble(m_a(3)),aimag(m_a(3)), sm
            endif
         else
            write (file_c) i,t,c(:)
            write (file_mu) i,t,mu_a(:)
            if (Fmag.eq.'mag') then
                write (file_m) i,t,dble(m_a(1)),aimag(m_a(1)),&
                    dble(m_a(2)),aimag(m_a(2)),dble(m_a(3)),aimag(m_a(3)),sm
            endif
         endif
       else
          if (Fbin.ne.'bin') then
             write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
             ! mut without ground state contribution
             !write (file_c,'(i8,f14.4,3e22.10)') i,t,mu_a_esa(:)
          else
             write(file_mu) i,t,mu_a(:)
             !write(file_c) i,t,mu_a_esa(:)
          endif
       endif
       j=int(dble(i)/dble(n_out))
       if(j.lt.1) j=1
       Sdip(:,1,j)=mu_a(:)
       if (Fmag.eq.'mag') Smag(:,1,j)=m_a(:)
! SP 270817: using get_* functions to communicate with TDPlas
       if(Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call get_medium_dip(Sdip(:,2,j))
       Sfld(:,j)=f(:,i)

       return

      end subroutine output


!------------------------------------------------------------------------
! @brief Create the field term of the hamiltonian 
!
! @date Created   : 
! Modified  : E. Coccia 22/11/2017
! Modified  : M. Monti 19/07/2022 MM
!------------------------------------------------------------------------
      subroutine add_int_vac(f_prev,h_int)

       implicit none

       real(dbl), intent(IN)    :: f_prev(3)
       real(dbl), intent(INOUT) :: h_int(nstates,nstates)
       real(dbl)                :: vec_prod(3),half_alpha ! e_dir vector f_prev -> vec_prod
       integer(i4b)             :: i,j

       half_alpha = 1.d0/(2.d0*clight)

! SC 16/02/2016: changed to - sign, 

       h_int(:,:)=h_int(:,:)-trans_dipoles(1,:,:)*f_prev(1)-             &
                 trans_dipoles(2,:,:)*f_prev(2)-trans_dipoles(3,:,:)*f_prev(3)
        !if (Fmag.eq.'mag') then
        !   call cross(e_dir,f_prev,vec_prod)
        !   h_int(:,:)=h_int(:,:)-half_alpha*trans_mag(1,:,:)*vec_prod(1)-  &
        !              half_alpha*trans_mag(2,:,:)*vec_prod(2)- &
        !              half_alpha*trans_mag(3,:,:)*vec_prod(3)
        !endif

       return
 
      end subroutine add_int_vac

!------------------------------------------------------------------------
! @brief Create the velocity-gauge interaction term from vector potential.
!        Uses trans_p matrix and updates h_int with -p·A.
!
! @date Created   : Manuel Sanchez 22/04/2026
! Modified  :
!------------------------------------------------------------------------
      subroutine add_int_vac_vg(a_prev,h_int_vg)

       implicit none

       real(dbl), intent(IN)    :: a_prev(3)
       complex(cmp), intent(INOUT) :: h_int_vg(nstates,nstates) ! Added by Manuel Sanchez 2026-04-22

       h_int_vg(:,:)=h_int_vg(:,:)-trans_p(1,:,:)*a_prev(1)-             & ! Added by Manuel Sanchez 2026-04-22
                 trans_p(2,:,:)*a_prev(2)-trans_p(3,:,:)*a_prev(3) ! Added by Manuel Sanchez 2026-04-22

       return
 
      end subroutine add_int_vac_vg
!------------------------------------------------------------------------
! @brief Build vg interaction matrix independently from lg interaction.
!        Uses only vacuum vg term (no medium contribution).
!
! @date Created   : Manuel Sanchez 22/04/2026
! Modified  :
!------------------------------------------------------------------------
      subroutine build_h_int_vg(step_idx,h_int_vg) ! Added by Manuel Sanchez 2026-04-22

       implicit none

       integer(i4b), intent(IN)      :: step_idx ! Added by Manuel Sanchez 2026-04-22
       complex(cmp), intent(INOUT)   :: h_int_vg(nstates,nstates) ! Added by Manuel Sanchez 2026-04-22

       h_int_vg=zeroc ! Added by Manuel Sanchez 2026-04-22
       call add_int_vac_vg(avec(:,step_idx),h_int_vg) ! Added by Manuel Sanchez 2026-04-22

       return

      end subroutine build_h_int_vg ! Added by Manuel Sanchez 2026-04-22
!------------------------------------------------------------------------
      subroutine cross(e_dir,f_prev,vec_prod) !MM
        implicit none

        real(dbl), intent(IN)  :: f_prev(3), e_dir(3)
        real(dbl), intent(OUT) :: vec_prod(3)

        vec_prod(1)=e_dir(2)*f_prev(3)-e_dir(3)*f_prev(2)
        vec_prod(2)=e_dir(3)*f_prev(1)-e_dir(1)*f_prev(3)
        vec_prod(3)=e_dir(1)*f_prev(2)-e_dir(2)*f_prev(1)

        return
      end subroutine cross
!------------------------------------------------------------------------
! @brief Calculate the Aharonov Lorentz radiative damping 
!
! @date Created   : S. Corni 
! Modified  : E. Coccia 22/11/2017
!------------------------------------------------------------------------
      subroutine add_int_rad(mu_prev,mu_prev2,mu_prev3,mu_prev4, &
                                                   mu_prev5,h_int)

       implicit none

       integer(i4b)                :: i,j
       real(dbl),    intent(in)    :: mu_prev(3),mu_prev2(3),mu_prev3(3), &
                 mu_prev4(3), mu_prev5(3)
       real(dbl),    intent(INOUT) :: h_int(nstates,nstates)
       real(dbl)                   :: d3_mu(3),d2_mu(3),d_mu(3),d2_mod_mu,coeff,scoeff, &
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

       h_int(:,:)=h_int(:,:)-trans_dipoles(1,:,:)*d3_mu(1)-             &
                 trans_dipoles(2,:,:)*d3_mu(2)-trans_dipoles(3,:,:)*d3_mu(3)

       return

      end subroutine add_int_rad

!------------------------------------------------------------------------
! @brief Write headers to output files 
!
! @date Created   : S. Corni 
! Modified  : E. Coccia 22/11/2017
!------------------------------------------------------------------------
      subroutine out_header

       implicit none

       write(file_c,'(2a)')'# istep   time (au)', &
                           '    Re(C_0) Im(C_0), Re(C_1) Im(C_1),...'

       write(file_e,'(8a)') '#   istep time (au)',' <H(t)>-E_gs(0)', &
              ' DE_vac(t)',' DG_eq(t)',' DG_neq(t)',  '  Const', &
              '  Rad. Int', '  Rad. Ene'   
      
       write(file_mu,'(4a)') '#   istep time (au)',' dipole-x ', &
              ' dipole-y ',' dipole-z '
 
       if (Fmag.eq.'mag') then   
          write(file_m,'(5a)') '#   istep time (au)',' mag_dipole-x ', &
                ' mag_dipole-y ',' mag_dipole-z ', ' dot_product(m,mu) ' ! MM
       endif

       return
    
      end subroutine out_header

!------------------------------------------------------------------------
! @brief Normalize coefficients for EuMar, RuKu4, and HeuSt (same rule).
!
! Fixed vs stochastic classification:
! - A state i is "fixed" if the i-th row of Rn is entirely zeros.
! - Otherwise it is "stochastic".
!
! Populations:
! - p_fixed = sum_{i fixed} |c(i)|^2
! - p_stch  = sum_{i stochastic} |c(i)|^2
!
! Normalization rule (EuMar, RuKu4, HeuSt):
! - fixed coefficients unchanged
! - stochastic coefficients multiplied by sqrt((1 - p_fixed)/p_stch)
! @date Created   : Manuel Sanchez 02/04/2026
! Modified  : Manuel Sanchez 02/04/2026
!------------------------------------------------------------------------
      subroutine normalize_c_eumar(c,Rn,nci)

        implicit none

        integer(i4b), intent(in)     :: nci
        complex(cmp), intent(inout)  :: c(nci)
        real(dbl),    intent(in)      :: Rn(nci,nci)

        integer(i4b) :: i
        logical       :: is_fixed(nci)
        real(dbl)     :: p_fixed, p_stch, scale, one_minus_pf

        do i=1,nci
          ! If Rp(i,*) is all zeros then Rp(i,*)*random_normal() gives Rn(:,i)=0 exactly.
          is_fixed(i) = all(Rn(i,:) == zero)
        enddo

        p_fixed = zero
        p_stch  = zero
        do i=1,nci
          if (is_fixed(i)) then
            p_fixed = p_fixed + abs(c(i))**2
          else
            p_stch  = p_stch  + abs(c(i))**2
          endif
        enddo

        if (p_stch > zero) then
          one_minus_pf = one - p_fixed!max(zero, one - p_fixed)
          scale = sqrt(one_minus_pf/p_stch)
          if (one_minus_pf < zero) then
            c = c / sqrt(dot_product(c,c))
          else
            do i=1,nci
              if (.not.is_fixed(i)) c(i) = c(i) * scale
            enddo
          endif
        endif

        return

      end subroutine normalize_c_eumar

!------------------------------------------------------------------------
! @brief One Markovian step: RK4 deterministic (-i Hn - 0.5 Gamma) c
!        then Euler-Maruyama noise -i*sqrt(dt)*Rn on the RK4 state.
!        Hn is diag(energies)+h_int; Gamma is diagonal (gamma_sum).
! @date Created   : Manuel Sanchez 03/04/2026
!------------------------------------------------------------------------
      subroutine mar_ruku4_apply(c_new, c_old, nci)

        implicit none
        integer(i4b), intent(in)    :: nci
        complex(cmp), intent(in)    :: c_old(nci)
        complex(cmp), intent(out)   :: c_new(nci)
        complex(cmp)                :: k1(nci), k2(nci), k3(nci), k4(nci)
        complex(cmp)                :: cst(nci), c_det(nci)
        real(dbl), parameter       :: one_sixth = 1.d0/6.d0

        k1 = (-ui)*(energies*c_old + matmul(h_int,c_old)) &
             - 0.5d0*gamma_sum(1:nci)*c_old
        cst = c_old + 0.5d0*dt*k1
        k2 = (-ui)*(energies*cst + matmul(h_int,cst)) &
             - 0.5d0*gamma_sum(1:nci)*cst
        cst = c_old + 0.5d0*dt*k2
        k3 = (-ui)*(energies*cst + matmul(h_int,cst)) &
             - 0.5d0*gamma_sum(1:nci)*cst
        cst = c_old + dt*k3
        k4 = (-ui)*(energies*cst + matmul(h_int,cst)) &
             - 0.5d0*gamma_sum(1:nci)*cst
        c_det = c_old + dt*one_sixth*(k1 + 2*k2 + 2*k3 + k4)
        c_new = c_det - ui*sqrt(dt)*matmul(Rn,c_det)

        return
      end subroutine mar_ruku4_apply

!------------------------------------------------------------------------
! @brief HeuSt step: predictor (-i Hn-0.5*Gamma)*dt and -i*sqrt(dt)*Rn on c,
!        same noise Rn on c_tilde; corrector averages drift and diffusion.
! @date Created   : Manuel Sanchez 03/04/2026
!------------------------------------------------------------------------
      subroutine mar_heust_apply(c_new, c_old, nci)

        implicit none
        integer(i4b), intent(in)    :: nci
        complex(cmp), intent(in)    :: c_old(nci)
        complex(cmp), intent(out)   :: c_new(nci)
        complex(cmp)                :: k(nci), f1d(nci), f1s(nci)
        complex(cmp)                :: f2d(nci), f2s(nci), ctil(nci)

        k = (-ui)*(energies*c_old + matmul(h_int,c_old)) &
            - 0.5d0*gamma_sum(1:nci)*c_old
        f1d = dt * k
        f1s = (-ui)*sqrt(dt)*matmul(Rn,c_old)
        ctil = c_old + f1d + f1s
        k = (-ui)*(energies*ctil + matmul(h_int,ctil)) &
            - 0.5d0*gamma_sum(1:nci)*ctil
        f2d = dt * k
        f2s = (-ui)*sqrt(dt)*matmul(Rn,ctil)
        c_new = c_old + 0.5d0*(f1d + f2d) + 0.5d0*(f1s + f2s)

        return
      end subroutine mar_heust_apply

!------------------------------------------------------------------------
! @brief Energy term is propagated analytically
! Interaction term via second-order Euler 
! 
! @date Created   : E. Coccia 15 Nov 2017
! Modified  : Manuel Sanchez 02/04/2026
!------------------------------------------------------------------------
      subroutine exp_euler_prop(ccexp,nci)

        implicit none

        integer(i4b),  intent(in)  :: nci
        complex(cmp),  intent(in)  :: ccexp(nci) 

        integer(i4b)               :: i,j,istart,iend,k
        real(dbl)                  :: t 
        complex(cmp)               :: dis(nci),ctmp(nci)

! Initialization only without restart
       if (Fres.eq.'Nonr') then
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
! SC: 31/10/17 modified the propagation with ccexp 
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             c=ccexp*(c_prev-ui*dt*matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
          else
             c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev))
          endif
          if (Fdis.eq."ernd") then
             do j=1,nstates
                c(j) = c(j) - ccexp(j)*ui*dt*krnd*random_normal()*c_prev(j)
             enddo
          endif
          if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
             dis=disp(h_dis,c_prev,nci)
             c=c-ccexp*dt*dis
             if (Fdis(5:9).eq."EuMar") then
               ! Euler-Maruyama stochastic step
                call build_rp_random_matrix(Rp,Rn,nstates)
                if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                   c=c-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev) ! Added by Manuel Sanchez 2026-04-22
                else
                   c=c-ui*dt*(energies*c_prev+matmul(h_int,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev)               
                endif
             elseif (Fdis(5:9).eq."RuKu4") then
                call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_ruku4_apply(c,c_prev,nci)
             elseif (Fdis(5:9).eq."HeuSt") then
                call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_heust_apply(c,c_prev,nci)
             endif
          endif
          if (Fdis(5:9).eq."EuMar".or.Fdis(5:9).eq."RuKu4".or.Fdis(5:9).eq."HeuSt") then
             call normalize_c_eumar(c,Rn,nci)
          else
             c=c/sqrt(dot_product(c,c))
          endif
          c_prev=c

! SP 16/07/17: added call to medium propagation at step 2 to have full output
          f_prev=f(:,2)
          h_int=zero
          if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") then
            i=2
            call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
          endif
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             call build_h_int_vg(2,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
          else
             call add_int_vac(f_prev,h_int)
          endif ! Added by Manuel Sanchez 2026-05-03
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
          if (Fmag.eq.'mag') then ! MM
             call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
          endif
          if (n_out.eq.1) call output(2,c,f_prev,h_int)
       endif

       if (Fres.eq.'Nonr') then
          istart=3
          iend=n_step
       elseif (Fres.eq.'Yesr') then
          istart=restart_i+1
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             call build_h_int_vg(restart_i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
          else
             call add_int_vac(f_prev,h_int)
          endif ! Added by Manuel Sanchez 2026-05-03
          if (Fsim.eq.'y') then 
             iend=diff_step+restart_i
          elseif (Fsim.eq.'n') then
             iend=n_step+istart-1
          endif
       endif

! PROPAGATION CYCLE: starts the propagation at timestep 3
! without restart, at timestep restart_i+1 otherwise
! Markovian dissipation (quantum jump) 
       if (Fdis(5:9).eq."qjump") then
          !do i=3,n_step
          do i=istart,iend 
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            if (prop_type.eq."mix") then
              dis=disp(h_dis,c_prev2,nci)
            elseif (prop_type.eq."full") then
              dis=disp(h_dis,c_prev,nci)
            endif
#ifndef OMP
            if (i.eq.ijump+1) then
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=ccexp*(c_prev-ui*dt*matmul(h_int_vg,c_prev)-dt*dis) ! Added by Manuel Sanchez 2026-04-22
               else
                  c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev)-dt*dis)
               endif
            else 
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int_vg,c_prev)-2.d0*dt*dis) ! Added by Manuel Sanchez 2026-04-22
               else
                  c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev)-2.d0*dt*dis)
               endif
            endif 
#endif
#ifdef OMP
            if (Fopt.eq.'omp') then
               ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
               do k=1,nci
                  do j=1,nci
                     ctmp(k)=ctmp(k)+ h_int(k,j)*c_prev(j)
                  enddo
               enddo
!$OMP END PARALLEL
               if (i.eq.ijump+1) then
                  c=ccexp*(c_prev-ui*dt*ctmp-dt*dis)
               else
                  c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*ctmp-2.d0*dt*dis)
               endif
            else
               if (i.eq.ijump+1) then
                  if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                     c=ccexp*(c_prev-ui*dt*matmul(h_int_vg,c_prev)-dt*dis) ! Added by Manuel Sanchez 2026-04-22
                  else
                     c=ccexp*(c_prev-ui*dt*matmul(h_int,c_prev)-dt*dis)
                  endif
               else
                  if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                     c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int_vg,c_prev)-2.d0*dt*dis) ! Added by Manuel Sanchez 2026-04-22
                  else
                     c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev)-2.d0*dt*dis)
                  endif
               endif
            endif
#endif
! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c_prev,nstates,pjump)
            call random_number(eps) 
            if (dtot.gt.eps)  then
               call quan_jump(c,c_prev,nstates,pjump)
               ijump=i
               n_jump=n_jump+1
#ifndef MPI
              if (Fwrt.eq.'yes') write(*,*) 'Quantum jump at step:', i, (i-1)*dt 
#endif 
               c_prev=c
            else
                c=c/sqrt(dot_product(c,c))
                c_prev2=c_prev
                c_prev=c
            endif

            f_prev=f(:,i)
            h_int=zero

            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, & 
                                               mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            if (mod(i-1,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend)
            endif
          enddo
! Markovian dissipation: EuMar; RuKu4 (RK4 drift + EM); HeuSt (Heun-type + averaged noise)
       elseif (Fdis(1:3).eq."mar") then
          !do i=3,n_step
          do i=istart,iend
! Dissipation by a continuous stochastic propagation
            ! call rnd_noise(w,w_prev,nstates,first)
            ! call add_h_rnd(h_rnd,nstates,w,w_prev)
            dis=disp(h_dis,c_prev,nci)
            if (Fdis(5:9).eq."EuMar") then
            call build_rp_random_matrix(Rp,Rn,nstates)
            ! Euler-Maruyama stochastic step
                if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                   c=c-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev) ! Added by Manuel Sanchez 2026-04-22
                else
                   c=c-ui*dt*(energies*c_prev+matmul(h_int,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev)
                endif
            elseif (Fdis(5:9).eq."RuKu4") then
            call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_ruku4_apply(c,c_prev,nci)
            elseif (Fdis(5:9).eq."HeuSt") then
            call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_heust_apply(c,c_prev,nci)
            endif
            if (Fdis(5:9).eq."EuMar".or.Fdis(5:9).eq."RuKu4".or.Fdis(5:9).eq."HeuSt") then
               call normalize_c_eumar(c,Rn,nci)
            else
               c=c/sqrt(dot_product(c,c))
            endif
            c_prev=c

            f_prev=f(:,i)
            h_int=zero

            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                  mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            ! Correct if n_out is different from 1 it print also t=0, t=n etc
            if (mod(i-1,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend) 
            endif
          enddo
       elseif (Fdis.eq."nodis".or.Fdis.eq."ernd") then
! No dissipation in the propagation 
          !do i=3,n_step
          do i=istart,iend
! SC 31/10/17: modified propagation by adding the exp term
#ifndef OMP
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
            else
               c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev))
            endif
#endif
#ifdef OMP
            if (Fopt.eq.'omp') then
               ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
               do k=1,nci 
                  do j=1,nci 
                     ctmp(k)=ctmp(k)+ h_int(k,j)*c_prev(j)
                  enddo
               enddo
!$OMP END PARALLEL
               c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*ctmp)  
            else
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
               else
                  c=ccexp*(ccexp*c_prev2-2.d0*ui*dt*matmul(h_int,c_prev))
               endif
            endif
#endif
            if (Fdis.eq."ernd") then
               do j=1,nstates
                  c(j) = c(j) - 2.d0*ccexp(j)*ui*dt*krnd*random_normal()*c_prev(j)
               enddo
            endif
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
            c_prev=c

            f_prev=f(:,i)
            h_int=zero
            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                               mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            if (mod(i-1,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend) 
            endif
          enddo
       endif 

       return

      end subroutine exp_euler_prop

!------------------------------------------------------------------------
! @brief Energy and interaction terms are propagated
! via second-order Euler 
! 
! @date Created   : E. Coccia 15 Nov 2017
! Modified  : Manuel Sanchez 02/04/2026
!------------------------------------------------------------------------      
      subroutine full_euler_prop(nci)

        implicit none

        integer(i4b), intent(in)  :: nci

        integer(i4b)              :: i,j,istart,iend,k
        real(dbl)                 :: t
        complex(cmp)              :: dis(nci),ctmp(nci)

!Initialization only without restart
       if (Fres.eq.'Nonr') then
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             c=c_prev-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
          else
             c=c_prev-ui*dt*(energies*c_prev+matmul(h_int,c_prev))
          endif
          if (Fdis.eq."ernd") then
             do j=1,nstates
                c(j) = c(j) - ui*dt*krnd*random_normal()*c_prev(j)
             enddo
          endif
          if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
             dis=disp(h_dis,c_prev,nci)
             c=c-dt*dis
             if (Fdis(5:9).eq."EuMar") then
         call build_rp_random_matrix(Rp,Rn,nstates)
          ! Euler-Maruyama stochastic step
                if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                   c=c-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev) ! Added by Manuel Sanchez 2026-04-22
                else
                   c=c-ui*dt*(energies*c_prev+matmul(h_int,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev)               
                endif
             elseif (Fdis(5:9).eq."RuKu4") then
         call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_ruku4_apply(c,c_prev,nci)
             elseif (Fdis(5:9).eq."HeuSt") then
         call build_rp_random_matrix(Rp,Rn,nstates)
                call mar_heust_apply(c,c_prev,nci)
             endif
          endif
          if (Fdis(5:9).eq."EuMar".or.Fdis(5:9).eq."RuKu4".or.Fdis(5:9).eq."HeuSt") then
             call normalize_c_eumar(c,Rn,nci)
          else
             c=c/sqrt(dot_product(c,c))
          endif
          c_prev=c

! SP 16/07/17: added call to medium propagation at step 2 to have full
! output
          f_prev=f(:,2)
          h_int=zero

          if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") then
             i=2
             call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
          endif
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             call build_h_int_vg(2,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
          else
             call add_int_vac(f_prev,h_int)
          endif ! Added by Manuel Sanchez 2026-05-03
          call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
          if (Fmag.eq.'mag') then ! MM
             call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
          endif
          if (mod(2,n_out).eq.0) call output(2,c,f_prev,h_int)
       endif

       if (Fres.eq.'Nonr') then
          istart=3
          iend=n_step
       elseif (Fres.eq.'Yesr') then
          istart=restart_i+1
          if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
             call build_h_int_vg(restart_i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
          else
             call add_int_vac(f_prev,h_int)
          endif ! Added by Manuel Sanchez 2026-05-03
          if (Fsim.eq.'y') then 
             iend=diff_step+restart_i
          elseif (Fsim.eq.'n') then
             iend=n_step+istart-1
          endif 
       endif 


! PROPAGATION CYCLE: starts the propagation at timestep 3
! Markovian dissipation (quantum jump) -> qjump
       if (Fdis(5:9).eq."qjump") then
          !do i=3,n_step
          do i=istart,iend
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            if (prop_type.eq."mix") then
              dis=disp(h_dis,c_prev2,nci)
            elseif (prop_type.eq."full") then
              dis=disp(h_dis,c_prev,nci)
            endif
#ifndef OMP
            if (i.eq.ijump+1) then
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=c_prev-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev))-dt*dis ! Added by Manuel Sanchez 2026-04-22
               else
                  c=c_prev-ui*dt*(energies*c_prev+matmul(h_int,c_prev))-dt*dis
               endif
            else
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev))-2.d0*dt*dis ! Added by Manuel Sanchez 2026-04-22
               else
                  c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int,c_prev))-2.d0*dt*dis
               endif
            endif
#endif
#ifdef OMP
            if (Fopt.eq.'omp') then
               ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
               do k=1,nci
                  do j=1,nci
                     ctmp(k)=ctmp(k)+ h_int(k,j)*c_prev(j)
                  enddo
               enddo
!$OMP END PARALLEL
               if (i.eq.ijump+1) then
                   c=c_prev-ui*dt*(energies*c_prev+ctmp)-dt*dis
               else
                   c=c_prev2-2.d0*ui*dt*(energies*c_prev+ctmp)-2.d0*dt*dis
               endif
            else
               if (i.eq.ijump+1) then
                   if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                      c=c_prev-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev))-dt*dis ! Added by Manuel Sanchez 2026-04-22
                   else
                      c=c_prev-ui*dt*(energies*c_prev+matmul(h_int,c_prev))-dt*dis
                   endif
               else
                   if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                      c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev))-2.d0*dt*dis ! Added by Manuel Sanchez 2026-04-22
                   else
                      c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int,c_prev))-2.d0*dt*dis
                   endif
               endif
            endif 
#endif


! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c_prev,nstates,pjump)
            call random_number(eps)
            if (dtot.gt.eps)  then
               call quan_jump(c,c_prev,nstates,pjump)
               ijump=i
               n_jump=n_jump+1
#ifndef MPI
              if (Fwrt.eq.'yes') write(*,*) 'Quantum jump at step:', i, (i-1)*dt
#endif
               c_prev=c
            else
                c=c/sqrt(dot_product(c,c))
                c_prev2=c_prev
                c_prev=c
            endif

            f_prev=f(:,i)
            h_int=zero
            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend) 
            endif
          enddo
! Markovian dissipation: EuMar; RuKu4 (RK4 drift + EM); HeuSt (Heun-type + averaged noise)
       elseif (Fdis(1:3).eq."mar") then
          !do i=3,n_step
          do i=istart,iend
! Dissipation by a continuous stochastic propagation
            ! call rnd_noise(w,w_prev,nstates,first)
            ! call add_h_rnd(h_rnd,nstates,w,w_prev)
            dis=disp(h_dis,c_prev,nci)
            if (Fdis(5:9).eq."EuMar") then
            call build_rp_random_matrix(Rp,Rn,nstates)
            ! Euler-Maruyama stochastic step
              if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                 c=c-ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev) ! Added by Manuel Sanchez 2026-04-22
              else
                 c=c-ui*dt*(energies*c_prev+matmul(h_int,c_prev)) - 0.5*dt*gamma_sum*c_prev - ui*sqrt(dt)*matmul(Rn,c_prev)
              endif
            elseif (Fdis(5:9).eq."RuKu4") then
            call build_rp_random_matrix(Rp,Rn,nstates)
              call mar_ruku4_apply(c,c_prev,nci)
            elseif (Fdis(5:9).eq."HeuSt") then
            call build_rp_random_matrix(Rp,Rn,nstates)
              call mar_heust_apply(c,c_prev,nci)
            endif
            if (Fdis(5:9).eq."EuMar".or.Fdis(5:9).eq."RuKu4".or.Fdis(5:9).eq."HeuSt") then
               call normalize_c_eumar(c,Rn,nci)
            else
               c=c/sqrt(dot_product(c,c))
            endif
            c_prev=c

            f_prev=f(:,i)
            h_int=zero

            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend) 
            endif
          enddo
       elseif (Fdis.eq."nodis".or.Fdis.eq."ernd") then
! No dissipation in the propagation 
          !do i=3,n_step
          do i=istart,iend
#ifndef OMP
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
            else
               c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int,c_prev))
            endif
#endif
#ifdef OMP
            if (Fopt.eq.'omp') then
               ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
               do k=1,nci
                  do j=1,nci
                     ctmp(k)=ctmp(k)+ h_int(k,j)*c_prev(j)
                  enddo
               enddo
!$OMP END PARALLEL
               c=c_prev2-2.d0*ui*dt*(energies*c_prev+ctmp)
            else
               if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
                  c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int_vg,c_prev)) ! Added by Manuel Sanchez 2026-04-22
               else
                  c=c_prev2-2.d0*ui*dt*(energies*c_prev+matmul(h_int,c_prev))
               endif
            endif
#endif

            if (Fdis.eq."ernd") then
               do j=1,nstates
                  c(j) = c(j) - 2.d0*ui*dt*krnd*random_normal()*c_prev(j)
               enddo
            endif
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
            c_prev=c

            f_prev=f(:,i)
            h_int=zero
            if (Fmdm.ne."vac".and.this_Finit_int.ne."qmt") call prop_medium(i,c_prev,mu_prev,f_prev,h_int)
            if (gauge.eq.'vg') then ! Added by Manuel Sanchez 2026-04-22
               call build_h_int_vg(i,h_int_vg) ! Added by Manuel Sanchez 2026-04-22
            else
               call add_int_vac(f_prev,h_int)
            endif ! Added by Manuel Sanchez 2026-05-03
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)

            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (Fmag.eq.'mag') then ! MM
               call do_m(c,m_prev,m_prev2,m_prev3,m_prev4,m_prev5)
            endif
            if (mod(i,n_out).eq.0) call output(i,c,f_prev,h_int)
            ! Restart
            if (mod(i,n_restart).eq.0) then
               t=(i-1)*dt
               call wrt_restart(i,t,c,c_prev,c_prev2,nstates,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend) 
            endif 
          enddo
       endif

       return
      
      end subroutine full_euler_prop

!------------------------------------------------------------------------
! @brief Write restart file 
! 
! @date Created   : E. Coccia 21 Nov 2017
! Modified  :
!------------------------------------------------------------------------      
      subroutine wrt_restart(i,t,c,c_prev,c_prev2,nci,iseed,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5,iend)
  
       implicit none

       integer(i4b),  intent(in)  :: i,nci,iseed,iend
       real(dbl),     intent(in)  :: t
       real(dbl),     intent(in)  :: mu_prev(3),mu_prev2(3),mu_prev3(3)
       real(dbl),     intent(in)  :: mu_prev4(3),mu_prev5(3)
       complex(cmp),  intent(in)  :: c(nci), c_prev(nci), c_prev2(nci)
     
       integer(i4b)               :: j,ii
       character(20)              :: filename

#ifndef MPI
       myrank=0
#endif

       ii=777+myrank
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

       !open(ii,file='restart')
       open(ii,file=filename)
       rewind(ii)
 
       write(ii,*) 'Restart time in au, Restart step'
       write(ii,*) t,i,iend-i,size(c)
       write(ii,*) 'Coefficients'
       do j=1,nci
          write(ii,*) c(j)
       enddo
       write(ii,*) 'Coefficients -1'
       do j=1,nci
          write(ii,*) c_prev(j)
       enddo
       write(ii,*) 'Coefficients -2'
       do j=1,nci
          write(ii,*) c_prev2(j)
       enddo
       if (Fdis.ne.'nodis') then
          write(ii,*) 'Seed'
          write(ii,*) iseed
          write(ii,*) 'Number of quantum jumps'
          write(ii,*) n_jump
       endif
       write(ii,*) 'Dipoles'
       write(ii,*) mu_prev(1), mu_prev(2), mu_prev(3)
       write(ii,*) mu_prev2(1), mu_prev2(2), mu_prev2(3)
       write(ii,*) mu_prev3(1), mu_prev3(2), mu_prev3(3)
       write(ii,*) mu_prev4(1), mu_prev4(2), mu_prev4(3)
       write(ii,*) mu_prev5(1), mu_prev5(2), mu_prev5(3)
       if (Fmag.eq.'mag') then
          write(ii,*) 'Magnetic Dipoles'
          write(ii,*) m_prev(1), m_prev(2), m_prev(3)
          write(ii,*) m_prev2(1), m_prev2(2), m_prev2(3)
          write(ii,*) m_prev3(1), m_prev3(2), m_prev3(3)
          write(ii,*) m_prev4(1), m_prev4(2), m_prev4(3)
          write(ii,*) m_prev5(1), m_prev5(2), m_prev5(3)
       endif

       close(ii)
 
       !flush(6)
      
       return

      end subroutine wrt_restart 

      end module
