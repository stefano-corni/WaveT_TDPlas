program nuclear_wp

!------------------------------------------------------------------------
! @brief Time evolution of the nuclear wave packet
! Pe(x,t) = \sum_ve,we Cve*(t)Cwe(t) Xve(x)Xwe(x)
 
! 
! @date Created   : E. Coccia 26 Feb 2018
! Modified  :
!------------------------------------------------------------------------

  use constants


  implicit none

  integer(i4b)               ::  i,j,k,l,m
  integer(i4b)               ::  ne,nv,nstep,ntot,nx,nout
  integer(i4b), allocatable  ::  istep(:) 
  real(dbl)                  ::  x,hk,hl,xmin,xmax,dx,mu
  real(dbl)                  ::  expp,cc,sq,tmp,fact
  real(dbl)                  ::  w(20),deq(20)
  real(dbl),    allocatable  ::  tstep(:),hv(:,:),rep(:,:),imp(:,:),npe(:)
  real(dbl)                  ::  ctmp,ccc
  complex(cmp), allocatable  ::  c(:,:),cpe(:,:)
  character*32               ::  filename,str
  character*200              ::  cdum



  namelist /nuclwp/ nstep,ne,nv,filename,dx,xmin,xmax,w,deq,mu,nout 

  w=0.d0
  deq=0.d0

  read(*,nml=nuclwp)


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
  write(*,*) 'Reduced mass (au)', mu
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

  w=w*cm_to_au

  open(10,file=filename)

  read(10,*) cdum
  do i=1,nstep
     read(10,*) istep(i), tstep(i), (c(j,i), j=1,ntot) 
  enddo

  close(10)

  do i=1,nx
     l=0
     do k=1,ne
        x = xmin + (i-1)*dx + deq(k)
        x = dsqrt(mu*w(k))*x
        expp=dexp(-0.5d0*x**2)
        sq=dsqrt(mu*w(k)/pi)
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

  !pe=0.d0
  rep=0.d0
  imp=0.d0
  do i=1,nstep
     do j=1,ntot
        do m=1,nx
           rep(m,i) = rep(m,i) + dble(c(j,i))*hv(j,m)
           imp(m,i) = imp(m,i) + aimag(c(j,i))*hv(j,m)
        enddo
     enddo
  enddo

  cpe=dcmplx(0.d0,0.d0)
  do i=1,nstep
     do j=1,ntot
        do k=1,ntot
           ccc=conjg(c(j,i))*c(k,i)
           do m=1,nx
              cpe(m,i) = cpe(m,i) + cc*hv(j,m)*hv(k,m)
           enddo
        enddo
     enddo 
  enddo

  do i=1,nstep
     npe(i) = sum(cpe(:,i))
  enddo

  !Population
  !do i=1,nstep
  !   if (mod(i,nout).eq.0) then
  !       write(str,*) i
  !       open(11+i,file=trim(str)//'.dat')
  !       write(11+i,*) '#xstep   P(x)'
  !       do m=1,nx
  !          x = xmin + (m-1)*dx
  !          write(11+i,*)  x, dsqrt(rep)/npe(i)
  !       enddo
  !       close(11+i)
  !   endif
  !enddo

  !do i=1,nstep
  !   write(999,*) i, sum(rep(:,i)**2)
  !enddo

  !Real part of the wave function
  !First step always printed out
  open(11,file='1r.dat')
  write(11,*) '#xstep   Re(x)'
  do m=1,nx
     x = xmin + (m-1)*dx
     write(11,*)  x, rep(m,1)**2/sum(rep(:,1)**2)
  enddo
  close(11)
  do i=1,nstep
     if (mod(i,nout).eq.0) then
         write(str,*) i
         open(11+i,file=trim(str)//'r.dat')
         write(11+i,*) '#xstep   Re(x)'
         do m=1,nx
            x = xmin + (m-1)*dx
            write(11+i,*)  x, rep(m,i)**2/sum(rep(:,i)**2)
         enddo
         close(11+i)
     endif
  enddo

  !Imaginary part of the wave function
  !First step always printed out
  open(11,file='1i.dat')
  write(11,*) '#xstep   Im(x)'
  do m=1,nx
     x = xmin + (m-1)*dx
     write(11,*)  x, imp(m,1)/npe(1)
  enddo
  close(11)
  do i=1,nstep
     if (mod(i,nout).eq.0) then
         write(str,*) i
         open(11+i,file=trim(str)//'i.dat')
         write(11+i,*) '#xstep   Im(x)'
         do m=1,nx
            x = xmin + (m-1)*dx
            write(11+i,*)  x, imp(m,i)/npe(i)
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

  stop

end program nuclear_wp



subroutine hermite(nn,x,y)
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

real(dbl) function fact(n)
!------------------------------------------------------------------------
! @brief Function computing the factorial n!
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

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

end function


