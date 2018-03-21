program twodfft

!------------------------------------------------------------------------
! @brief Energy term is propagated analytically
! Interaction term via second-order Euler 
! 
! @date Created   : E. Coccia 22 Jan 2018
! Modified  :
!------------------------------------------------------------------------

 use constants 
 use, intrinsic :: iso_c_binding

 implicit none

 integer(i4b)               :: dim2,i,j,plan,maxd3,dscan1s,nfout
 integer(i4b)               :: jj,n3,idum,m,dscan1,dmid,tdim3
 integer(i4b)               :: nhann1,nhann3
 integer(i4b), allocatable  :: dim1(:),dim3(:),dim3e(:),ntot(:)
 real(dbl)                  :: rdum,dw1,dw3,dt1,dt3,modr,modi,w1,w3
 real(dbl)                  :: w1min,w1max,w3min,w3max,tmid,modd
 real(dbl)                  :: dir(3),mu(3),sigma3
 real(dbl)                  :: lo_fmax(3),lo_sigma,lo_w,lo_tmid
 real(dbl),    allocatable  :: dinp(:,:,:),inp(:,:),f(:,:),wf1(:),wf3(:) 
 complex(cmp), allocatable  :: doutp(:,:),inp1(:,:,:)
 character*3                :: field
 character*20               :: filename,detector
 character*100              :: cdum

 include 'fftw3.f03'

 namelist /fft/ dscan1,dim2,dt1,dt3,dir,sigma3,w1min,w1max,w3min,w3max,&
                tmid,detector,field,lo_fmax,lo_sigma,lo_w,lo_tmid,nfout

 read(*,nml=fft)

 write(*,*) '**********************************************'
 write(*,*) '*                                            *'
 write(*,*) '*               Tool for a                   *'
 write(*,*) '*                 2D FFT                     *'
 write(*,*) '*                                            *'
 write(*,*) '*                   by                       *'
 write(*,*) '*              Stefano Corni                 *'
 write(*,*) '*                  and                       *'
 write(*,*) '*             Emanuele Coccia                *'
 write(*,*) '*                                            *'
 write(*,*) '**********************************************'

 allocate(dim1(dscan1))
 allocate(dim3(dscan1))
 allocate(dim3e(dscan1))
 allocate(ntot(dscan1))

 open(9,file='t1_t3.dat')
 do m=1,dscan1
    read(9,*) dim1(m),dim3(m)
 enddo
 close(9)

 if (nfout.le.0) then
    write(*,*) 'ERROR: nfout must be positive'
    stop
 endif

 if (dscan1.lt.0) then
    write(*,*) 'ERROR: dscan1 must be positive'
    stop
 endif 

 if (dim2.lt.0) then
    write(*,*) 'ERROR: dim2 must be positive'
    stop
 endif

 if (dt1.le.0) then
    write(*,*) 'ERROR: dt1 must be larger than zero'
    stop
 endif

 if (dt3.le.0) then
    write(*,*) 'ERROR: dt3 must be larger than zero'
    stop
 endif

 if (sigma3.lt.0) then
    write(*,*) 'ERROR: sigma3 must be positive'
    stop
 endif 

 if (dir(1).le.0.and.dir(2).le.0.and.dir(3).le.0) then
    write(*,*) 'ERROR: One component of dir_ft must be larger than zero'
    stop
 endif


 write(*,*) ''
 write(*,*) 'Number of points for t1', dscan1
 write(*,*) 'Center of the first pulse', tmid
 write(*,*) 'Time step along the first interval (au)', dt1
 write(*,*) 'Time step along the third interval (au)', dt3
 write(*,*) 'Direction for FFT', dir(1),dir(2),dir(3)
 write(*,*) 'Sigma of the enevelope of the third pulse (au)', sigma3
 write(*,*) 'Frequency range for w1 (au)', w1min,w1max
 write(*,*) 'Frequency range for w3 (au)', w3min,w3max 
 write(*,*) 'Detector', detector 
 write(*,*) 'Output frequency', nfout
 if (detector(1:10).eq.'heterodyne') then
    write(*,*) 'Local oscillator field (lof) shape', field
    write(*,*) 'Maximum lof', lo_fmax(:)
    write(*,*) 'Center of lof', lo_tmid
    write(*,*) 'Lof frequency', lo_w
    write(*,*) 'Lof sigma', lo_sigma
 endif
 write(*,*) ''


 n3=int((3.d0*sigma3)/dt3)
 dmid=int(tmid/dt3)

 do m=1,dscan1
   dim3e(m)=dim3(m)+n3
 enddo
 if(mod(dscan1,2).gt.0) dscan1=dscan1-1

 do m=1,dscan1
    ntot(m)=dmid+dim1(m)+dim2+dim3(m)
 enddo

 maxd3=maxval(ntot(:))

 if(mod(maxd3,2).gt.0) maxd3=maxd3-1

 nhann1=int(two*dscan1/dt1)
 nhann3=int(two*maxd3/dt3)

 allocate (wf1(dscan1))
 allocate (wf3(maxd3))
 allocate (dinp(maxd3,dscan1,3))
 allocate (inp(maxd3,dscan1))
 allocate (inp1(maxd3,dscan1,1))
 allocate (doutp(int(dble(maxd3)/two)+1,dscan1))

 dw1=2*pi/(dble(dscan1)*dt1)
 dw3=2*pi/(dble(maxd3)*dt3)

 open(8,file='2d_spectrum.dat',status="unknown")

 dinp(:,:,:)=0.d0
 wf1=0.d0
 wf3=0.d0
 inp=0.d0

 do m=1,dscan1
    wf1(m)=0.5d0*(1.d0 - cos(2.d0*pi*m/nhann1))
 enddo

 do i=1,maxd3
    wf3(i)=0.5d0*(1.d0 - cos(2.d0*pi*i/nhann3))
 enddo

 do m=1,dscan1

    if (m.lt.10) then
        WRITE(filename,'(a,i1.1,a)') "mu_t_",m,".dat"
    elseif (m.ge.10.and.m.lt.100) then
        WRITE(filename,'(a,i2.2,a)') "mu_t_",m,".dat"
    elseif (m.ge.100.and.m.lt.1000) then
        WRITE(filename,'(a,i3.3,a)') "mu_t_",m,".dat"
    elseif (m.ge.1000.and.m.lt.10000) then
        WRITE(filename,'(a,i4.4,a)') "mu_t_",m,".dat"
    elseif (m.ge.10000.and.m.lt.100000) then
        WRITE(filename,'(a,i5.5,a)') "mu_t_",m,".dat"
    endif

    open(7,file=filename)

    jj=0
    tdim3 = dmid+dim1(m)+dim2
    read(7,*) cdum
    if (n3.ge.tdim3) then
       do i=1,maxd3
          read(7,*) idum, rdum, mu(:)
          dinp(i,m,1) = mu(1)
          dinp(i,m,2) = mu(2)
          dinp(i,m,3) = mu(3) 
       enddo
    elseif (n3.lt.tdim3) then
       do i=1,ntot(m)
          read(7,*) idum, rdum, mu(:) 
          if (i.gt.(tdim3-n3)) then
             jj=jj+1 
             dinp(jj,m,1) = mu(1)
             dinp(jj,m,2) = mu(2)
             dinp(jj,m,3) = mu(3) 
          endif
       enddo
    endif

    close(7)

 enddo

 doutp=cmplx(0.d0,0.d0)

 dir(:)=dir(:)/sqrt(dot_product(dir,dir))

 do j=1,dscan1
    do i=1,maxd3
       inp(i,j)=dot_product(dinp(i,j,:),dir(:))
    enddo
 enddo

 do m=1,dscan1
    do i=1,maxd3
       inp(i,m)=inp(i,m)*wf1(m)*wf3(i)
     enddo
 enddo

 !Homodyne detection
 if (detector(1:8).eq.'homodyne') then
    !inp(:,:) = inp(:,:)**2
 elseif (detector(1:10).eq.'heterodyne') then
    allocate(f(3,maxd3))
    f(:,:)=0.d0
    call lo_field(dt3,field,maxd3,lo_fmax,lo_sigma,lo_w,lo_tmid,f)
    !do i=1,maxd3
    !  inp(:,i) = inp(:,i)*f()
    !enddo
    write(*,*) 'Not implemented yet'
    stop
 endif

 dscan1s=dscan1

 do j=1,dscan1
    do i=1,maxd3
       inp1(i,j,1) = inp(i,j)
    enddo
 enddo

 !call rlft3(inp1,doutp,maxd3,dscan1,1,1)

 call dfftw_plan_dft_r2c_2d(plan,maxd3,dscan1,inp,doutp,FFTW_ESTIMATE)
 call dfftw_execute_dft_r2c(plan,inp,doutp)
 call dfftw_destroy_plan(plan)

 !The signal is in [-w1,w3] quadrant
 !A(-w) = A*(w)
 !doutp=conjg(doutp)

 do i=1,dscan1s
    w1=(i-1)*dw1
    do j=1,int(maxd3/two)
       w3=(j-1)*dw3
       if (w1.ge.w1min.and.w1.le.w1max) then
          if (w3.ge.w3min.and.w3.le.w3max) then
             modi=-aimag(doutp(j,i))
             modr=real(doutp(j,i))
             modd=dsqrt(modr**2+modi**2)
             if (mod(i,nfout).eq.0.and.mod(j,nfout).eq.0) then
                 write(8,'(6e20.10)') w1/ev_to_au, w3/ev_to_au, modr, modi, modd 
             endif
          endif
       endif
    enddo 
    if (w1.ge.w1min.and.w1.le.w1max) write(8,'(2x)') 
 enddo

 close(8)

 deallocate(dim1)
 deallocate(dim3)
 deallocate(dim3e)
 deallocate(ntot)
 deallocate(dinp)
 deallocate(doutp)
 deallocate(inp)
 deallocate(inp1)
 deallocate(wf1)
 deallocate(wf3)
 if (detector(1:10).eq.'heterodyne') deallocate(f)

 write(*,*) 'End of the calculation'
 
 stop


end program twodfft

subroutine lo_field(dt,field,n_tot,lo_fmax,lo_sigma,lo_w,lo_tmid,f)

!------------------------------------------------------------------------
! @brief Local oscillator field for heterodyne detection 
! 
! @date Created   : E. Coccia 13 Feb 2018
! Modified  :
!------------------------------------------------------------------------

 use constants

 implicit none

 integer(i4b), intent(in)     :: n_tot
 character*3,  intent(in)     :: field
 real(dbl),    intent(in)     :: dt,lo_fmax(3),lo_sigma,lo_w,lo_tmid 
 real(dbl),    intent(inout)  :: f(3,n_tot) 

 integer(i4b) :: i,j,i_max
 real(dbl)    :: t_a,ti,tf,arg

 select case (field)
    case ("mdg")
    ! Gaussian modulated sinusoid: exp(-(t-t0)^2/s^2) * sin(wt) 
       do i=1,n_tot
          t_a=dt*(i-1)
          f(:,i) = lo_fmax(:)*exp(-pt5*(t_a-lo_tmid)**2/(lo_sigma**2))*& 
                   sin(lo_w*t_a)
       enddo
    case ("mds")
    ! Cosine^2 modulated sinusoid: 1/2* cos^2(pi(t-t0)/(2t0)) * sin(wt) 
    !          f=0 for t>t0
       i_max=int(lo_tmid/dt)
       if (2*i_max.gt.n_tot) then
          write(*,*) 'ERROR: 2*t_mid/dt must be smaller than', n_tot
          stop
       endif
       do i=1,2*i_max
          t_a=dt*(dble(i)-1)
          f(:,i)=lo_fmax(:)*cos(pi*(t_a-lo_tmid)/(2*lo_tmid))**2/2.d0*&
                 sin(lo_w*t_a)
       enddo
       do i=2*i_max+1,n_tot
          t_a=dt*(i-1)
          f(:,i)=0.
       enddo
    case ("pip")
    ! Pi pulse: cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)) 
       do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=0.d0
          if (abs(t_a-lo_tmid).lt.lo_sigma) then
             f(:,i)=lo_fmax(:)*(cos(pi*(t_a-lo_tmid)/(2*lo_sigma)))**2*&
             cos(lo_w*(t_a-lo_tmid))
          endif
       enddo
    case ("gau")
    ! Gaussian pulse: exp(-(t-t0)^2/s^2) 
       do i=1,n_tot
          t_a=dt*(i-1)
          f(:,i)=lo_fmax(:)*exp(-pt5*(t_a-lo_tmid)**2/(lo_sigma**2))
       enddo
    case ("css")
    ! Cos^2 pulse (only half a period): cos^2(pi*(t-t0)/(s)) 
       ti=lo_tmid-lo_sigma/two
       tf=lo_tmid+lo_sigma/two
       do i=1,n_tot
          t_a=dt*(dble(i)-1)
          f(:,i)=zero
          if (t_a.gt.ti.and.t_a.le.tf) then
             f(:,i)=lo_fmax(:)*(cos(pi*(t_a-lo_tmid)/(lo_sigma)))**2
          endif
       enddo
    case default
       write(*,*)  "Error: wrong field type"
       stop
    end select

    return

end subroutine lo_field

SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)

use constants

INTEGER(i4b)   isign,nn1,nn2,nn3
COMPLEX(cmp)   data(nn1/2,nn2,nn3),speq(nn2,nn3)
INTEGER(i4b)   i1,i2,i3,j1,j2,j3,nn(3)
real(dbl)      theta,wi,wpi,wpr,wr,wtemp
COMPLEX(cmp)   c1,c2,h1,h2,w

c1=cmplx(0.5,0.0)
c2=cmplx(0.0,-0.5*isign)
theta=6.28318530717959d0/dble(isign*nn1)
wpr=-2.0d0*sin(0.5d0*theta)**2
wpi=sin(theta)
nn(1)=nn1/2
nn(2)=nn2
nn(3)=nn3

if(isign.eq.1)then
  call fourn(data,nn,3,isign)
  do i3=1,nn3
     do i2=1,nn2
        speq(i2,i3)=data(1,i2,i3)
     enddo
   enddo
endif

do i3=1,nn3
   j3=1
   if (i3.ne.1) j3=nn3-i3+2
   wr=1.0d0
   wi=0.0d0
   do i1=1,nn1/4+1
     j1=nn1/2-i1+2
     do i2=1,nn2
        j2=1
        if (i2.ne.1) j2=nn2-i2+2
           if (i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
           else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
        endif
     enddo
     wtemp=wr
     wr=wr*wpr-wi*wpi+wr
     wi=wi*wpr+wtemp*wpi+wi
     w=cmplx(sngl(wr),sngl(wi))
   enddo
enddo
if(isign.eq.-1)then
  call fourn(data,nn,3,isign)
endif

return

END


SUBROUTINE fourn(data,nn,ndim,isign)

use constants

INTEGER(i4b)  isign,ndim,nn(ndim)
REAL(cmp)     data(*)
INTEGER(i4b)  i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2
INTEGER(i4b)  ip3,k1,k2,n,nprev,nrem,ntot
REAL(dbl)     tempi,tempr
real(dbl)     theta,wi,wpi,wpr,wr,wtemp

ntot=1
do idim=1,ndim
   ntot=ntot*nn(idim)
enddo

nprev=1
do idim=1,ndim
   n=nn(idim)
   nrem=ntot/(n*nprev)
   ip1=2*nprev
   ip2=ip1*n
   ip3=ip2*nrem
   i2rev=1
   do i2=1,ip2,ip1
      if (i2.lt.i2rev)then
        do i1=i2,i2+ip1-2,2
           do i3=i1,ip3,ip2
              i3rev=i2rev+i3-i2
              tempr=data(i3)
              tempi=data(i3+1)
              data(i3)=data(i3rev)
              data(i3+1)=data(i3rev+1)
              data(i3rev)=tempr
              data(i3rev+1)=tempi
           enddo
        enddo
      endif
      ibit=ip2/2
1     if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
         i2rev=i2rev-ibit
         ibit=ibit/2
      goto 1
      endif
      i2rev=i2rev+ibit
   enddo
   ifp1=ip1
2  if (ifp1.lt.ip2)then
      ifp2=2*ifp1
      theta=isign*6.28318530717959d0/(ifp2/ip1)
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do i3=1,ifp1,ip1
         do i1=i3,i3+ip1-2,2
            do i2=i1,ip3,ifp2
               k1=i2
               k2=k1+ifp1
               tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
               tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
               data(k2)=data(k1)-tempr
               data(k2+1)=data(k1+1)-tempi
               data(k1)=data(k1)+tempr
               data(k1+1)=data(k1+1)+tempi
            enddo
         enddo
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
      enddo
      ifp1=ifp2
   goto 2
   endif
   nprev=n*nprev
enddo

return

END

