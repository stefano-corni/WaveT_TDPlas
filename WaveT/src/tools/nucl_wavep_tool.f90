!------------------------------------------------------------------------
! @brief Time evolution of the nuclear wave packet
! 
! 
! @date Created   : E. Coccia 26 Feb 2018
! Modified  :       E. Coccia 18 Dec 2018 
!------------------------------------------------------------------------
program nuclear_wp

  use constants


  implicit none

  integer(i4b)               ::  i,j,k,l,m
  integer(i4b)               ::  ne,nv,nstep,ntot,nx,nout
  integer(i4b), allocatable  ::  istep(:) 
  real(dbl)                  ::  x,hk,hl,xmin,xmax,dx
  real(dbl)                  ::  expp,cc,sq,tmp,fact
  real(dbl)                  ::  w(20),deq(20),mu(20)
  real(dbl),    allocatable  ::  tstep(:),hv(:,:),rep(:,:),imp(:,:),npe(:)
  real(dbl)                  ::  ctmp,ccc
  real(dbl),    allocatable  ::  rc(:),ic(:)          
  complex(cmp), allocatable  ::  c(:,:),cpe(:,:)
  character*32               ::  filename,str
  character*200              ::  cdum
  character*4000             ::  fmt_ci



  namelist /nuclwp/ nstep,ne,nv,filename,dx,xmin,xmax,w,deq,mu,nout 

  w=0.d0
  deq=0.d0
  mu=0.d0

  read(*,nml=nuclwp)

  mu=mu*amu_to_au

  write(*,*) '**********************************************'
  write(*,*) '*                                            *'
  write(*,*) '*               Tool for the                 *'
  write(*,*) '*              time evolution                *'
  write(*,*) '*                  of the                    *'
  write(*,*) '*            nuclear wavepacket              *'
  write(*,*) '*                                            *'
  write(*,*) '*                   by                       *'
  write(*,*) '*              Stefano Corni                 *'
  write(*,*) '*                  and                       *'
  write(*,*) '*             Emanuele Coccia                *'
  write(*,*) '*                                            *'
  write(*,*) '**********************************************'

  write(*,*) ''
  write(*,*) 'File of the coefficients', filename
  write(*,*) 'Number of steps', nstep
  write(*,*) 'Number of electronic states', ne
  write(*,*) 'Total number of vibrational states per electronic state', nv
  write(*,*) 'Minimum x value (au)', xmin
  write(*,*) 'Maximum x value (au)', xmax
  write(*,*) 'Spatial step (au)', dx

  if (nstep.le.0) then
     write(*,*) 'ERROR: nstep must be positive',nstep
     stop
  endif

  if (ne.le.0) then
     write(*,*) 'ERROR: ne must be positive',ne
     stop
  endif

  do i=1,ne
     write(*,*) 'Frequency (cm-1) of harmonic oscillator for state', i, w(i)
  enddo
  do i=1,ne-1
     write(*,*) 'Displacement (bohr) for harmonic oscillators for state', i+1, deq(i+1)
  enddo
  write(*,*) ''

  if (nv.le.0) then
     write(*,*) 'ERROR: nv must be positive',nv
     stop
  endif

  if (xmax.le.xmin) then
     write(*,*) 'ERROR: xmax must be larger than xmin',xmax,xmin 
     stop
  endif

  if (dx.le.0.d0) then
     write(*,*) 'ERROR: dx must be positive',dx
     stop
  endif

  ntot=ne*nv
  nx=int((xmax-xmin)/dx)

  allocate(istep(nstep))
  allocate(tstep(nstep))
  allocate(hv(ntot,nx))
  allocate(cpe(nx,nstep))
  allocate(rep(nx,nstep))
  allocate(imp(nx,nstep))
  allocate(npe(nstep))
  allocate(c(ntot,nstep))
  allocate(rc(ntot),ic(ntot)) 

  w=w*cm_to_au

  open(10,file=filename,status="unknown")

  read(10,*) cdum
  write (fmt_ci,'("(i8,f14.4,",I0,"e17.8E3)")') 2*ntot
  do i=1,nstep
     read(10,fmt_ci) istep(i), tstep(i), (rc(j), ic(j), j=1,ntot)
     do j=1,ntot
        c(j,i) = dcmplx(rc(j),ic(j))
     enddo
  enddo

  deallocate(rc,ic)

  close(10)

  do i=1,nx
     l=0
     do k=1,ne
        x = xmin + (i-1)*dx + deq(k)
        x = dsqrt(mu(k)*w(k))*x
        expp=dexp(-0.5d0*x**2)
        sq=dsqrt(mu(k)*w(k)/pi)
        cc=dsqrt(sq)
        cc=cc*expp  
        do j=1,nv
           l=l+1
           call hermite(j-1,x,hv(l,i))
           tmp=1.d0/dsqrt(2.d0**(j-1)*fact(j-1))
           hv(l,i)=tmp*cc*hv(l,i)
        enddo
     enddo
  enddo

  rep=0.d0
  imp=0.d0
  do i=1,nstep
     !Do not include vibrational ground state of the electronic ground state
     do j=2,ntot
        do m=1,nx
           rep(m,i) = rep(m,i) + dble(c(j,i))*hv(j,m)
           imp(m,i) = imp(m,i) + aimag(c(j,i))*hv(j,m)
        enddo
     enddo
  enddo

  cpe=dcmplx(0.d0,0.d0)
  do i=1,nstep
     do j=2,ntot
        do k=2,ntot
           ccc=conjg(c(j,i))*c(k,i)
           do m=1,nx
              cpe(m,i) = cpe(m,i) + ccc*hv(j,m)*hv(k,m)
           enddo
        enddo
     enddo 
  enddo

  do i=1,nstep
     npe(i) = sum(cpe(:,i))
  enddo

  !Square modulus of the wave function
  !First step always printed out
  open(11,file='1m.dat')
  write(11,*) '#xstep   |Psi(x)|^2'
  do m=1,nx
     x = xmin + (m-1)*dx
     write(11,*)  x, (rep(m,1)**2+imp(m,1)**2)/npe(1)
  enddo
  close(11)
  do i=1,nstep
     if (mod(i,nout).eq.0) then
         write(str,*) i
         open(11+i,file=trim(str)//'m.dat')
         write(11+i,*) '#xstep   |Psi(x)|^2'
         do m=1,nx
            x = xmin + (m-1)*dx
            !write(11+i,*)  x, (rep(m,i)**2+imp(m,i)**2)/npe(i)
            write(11+i,*)  x, (imp(m,i)**2)/npe(i)
         enddo
         close(11+i)
     endif
  enddo


  deallocate(istep)
  deallocate(tstep)
  deallocate(hv)
  deallocate(cpe)
  deallocate(rep)
  deallocate(imp)
  deallocate(npe)
  deallocate(c)

  write(*,*) ''
  write(*,*) 'End of simulation'
  write(*,*) ''

  stop

end program nuclear_wp


!------------------------------------------------------------------------
!   @brief Computes the value of the Hermite polynomial of degree nn             
!   at a given point               
!   nn = degree of the polynomial                                   
!   x  = point at which the calculation is done 
!   y  = value of the polynomial in x                                
!
!   @date Created  :  E. Coccia 8 Sep 2017 
!   Modified   :
!------------------------------------------------------------------------
subroutine hermite(nn,x,y)
  
        use constants

        implicit none


        integer(i4b), intent(in)  :: nn
        real(dbl),    intent(in)  :: x
        real(dbl),    intent(out) :: y
        integer(i4b)              :: k
        real(dbl)                 :: yp,dk,ym

        y = 1.d0
        if (nn.eq.0) return

        y = 2.d0*x
        if (nn.eq.1) return

        yp=1.d0
        do 1 k=2,nn
           dk = dble(k-1)
           ym = y
           y  = 2.d0*x*y-2.d0*dk*yp
           yp = ym
1       continue

        return

end subroutine hermite

!------------------------------------------------------------------------
! @brief Function computing the factorial n!
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
real(dbl) function fact(n)

       use constants

       integer(i4b), intent(in) :: n
       integer(i4b)             :: k
       real(dbl)                :: s

       s=1.d0

       do k=1,n
          s=s*dble(k)
       enddo

       fact=s

       if (n.eq.0) fact=1.d0

       if (fact.gt.plus_inf) fact=plus_inf

       return

end function fact


