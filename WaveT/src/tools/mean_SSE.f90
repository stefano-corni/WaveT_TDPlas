program mean_for_PostProc
        implicit none
 integer                 :: nrep,nstates,nsteps,neff,nfreq
 complex*16              :: tmp
 complex*16, allocatable :: coeff(:,:,:), cmean(:,:),corr(:,:), coeff_read(:)
 integer, allocatable    :: step(:)
 real(8), allocatable    :: t(:),rdum(:),idum(:),mix(:,:)
 integer                 :: m,j,k,i,kk
 character(4000)         :: fmt_c, fmt_r, filename, mixing
 logical                 :: binary

  read(*,*) nsteps  !number of time steps
  read(*,*) nrep    !number of trajs
  read(*,*) nstates !total number of states
  read(*,*) binary  !write files as formatted or unformatted
  read(*,*) nfreq   !how many times write the output files
  read(*,*) mixing 

  neff=nsteps/nfreq

  allocate(step(neff),coeff(neff,nrep,nstates))
  allocate(t(neff))
  allocate(corr(neff,nstates*(nstates+1)/2),cmean(neff,nstates))
  allocate(rdum(nstates),idum(nstates))

  do m=1,nrep
     if (m.lt.10) then
        WRITE(filename,'(a,i1.1,a)') "c_t_",m,".dat"
     elseif (m.ge.10.and.m.lt.100) then
        WRITE(filename,'(a,i2.2,a)') "c_t_",m,".dat"
     elseif (m.ge.100.and.m.lt.1000) then
        WRITE(filename,'(a,i3.3,a)') "c_t_",m,".dat"
     elseif (m.ge.1000.and.m.lt.10000) then
        WRITE(filename,'(a,i4.4,a)') "c_t_",m,".dat"
     elseif (m.ge.10000.and.m.lt.100000) then
        WRITE(filename,'(a,i5.5,a)') "c_t_",m,".dat"
     endif
     open(20+m,file=filename)
     !write(*,*) filename
  enddo
  if (mixing.eq."yes") then
       allocate(mix(nstates,nstates),coeff_read(nstates))
       open(19, file="mixing_wavet.dat", status="old")
       do i=1,nstates
           read(19,*) (mix(i,k), k=1,nstates)
       enddo
       close(19)
  endif

  !write(*,*) 'uno'
  !stop
  coeff(:,:,:)=0
  do m=1,nrep
     read(20+m,*) 
     kk=0
     do j=1,nsteps
        if (mod(j,nfreq).eq.0) then
           kk=kk+1
           read(20+m,*) step(kk), t(kk), (rdum(k),idum(k), k=1,nstates)
           if (mixing.eq."yes") then
               do k=1,nstates
                 coeff_read(k)=dcmplx(rdum(k),idum(k))
               enddo
               do k=1,nstates
                      do i=1,nstates
                           coeff(kk,m,k)=coeff(kk,m,k)+coeff_read(i)*mix(i,k)
                      enddo
               enddo
           else
               do k=1,nstates
                 coeff(kk,m,k)=dcmplx(rdum(k),idum(k))
               enddo
           endif
        endif
     enddo
  enddo

  deallocate(rdum)
  deallocate(idum) 
  !write(*,*) 'due'

  do j=1,neff
     do i=1,nstates
        tmp=0.d0
        do m=1,nrep
           tmp = tmp + coeff(j,m,i)
        enddo
        tmp=tmp/real(nrep)
        cmean(j,i)=tmp
     enddo
  enddo

  write(*,*) 'tre'

  do j=1,neff
     kk=0
     do i=1,nstates
        do k=i,nstates
           kk=kk+1
           tmp=0.d0
           do m=1,nrep
              tmp=tmp+conjg(coeff(j,m,i))*coeff(j,m,k)
           enddo
           tmp=tmp/real(nrep)
           corr(j,kk)=tmp-conjg(cmean(j,i))*cmean(j,k)
        enddo
     enddo
  enddo

  deallocate(coeff)
  write(*,*) 'quattro'

  open(10,file='c_t_avg.dat')
  open(11,file='corr_t_avg.dat')

  write(10,*) 'Mean coefficients, same structure of coefficient file'
  if (.not.binary) then
     write (fmt_c,'("(i8,f14.4,",I0,"e17.8E3)")') 2*nstates
     write (fmt_r,'("(i8,f14.4,",I0,"e17.8E3)")') 2*(nstates*(nstates+1)/2)
     do j=1,neff
        write(10, fmt_c) step(j), t(j),(real(cmean(j,k)),aimag(cmean(j,k)), k=1,nstates)
        write(11, fmt_r) step(j), t(j),(real(corr(j,i)),aimag(corr(j,i)), i=1,(nstates*(nstates+1)/2)) 
     enddo
  else
     do j=1,neff
        write(10) step(j), t(j),(real(cmean(j,k)),aimag(cmean(j,k)),k=1,nstates)
        write(11) step(j), t(j),(real(corr(j,i)),aimag(corr(j,i)),i=1,(nstates*(nstates+1)/2))
     enddo
  endif
  close(10)
  close(11)

  write(*,*) 'cinque'

  deallocate(step)
  deallocate(t)
  deallocate(corr)
  deallocate(cmean)

  if (mixing.eq."yes") deallocate(mix,coeff_read)

  write(*,*) 'End of simulation'
 
  stop 

end program
