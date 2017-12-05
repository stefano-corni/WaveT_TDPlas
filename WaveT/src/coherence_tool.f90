program decoherence 

 implicit none

 integer                 :: nrep,nstates,nsteps,npair
 real(8), allocatable    :: rdum(:,:,:),pop(:,:),perr(:,:)
 real(8), allocatable    :: cor(:,:),coi(:,:),co(:,:)
 real(8), allocatable    :: rerr(:,:), ierr(:,:)
 integer, allocatable    :: i(:)
 real(8), allocatable    :: t(:),l1(:),l1_err(:)
 real(8), allocatable    :: idum(:,:,:),ferr(:,:)

 integer       :: j,k,m 
 character(30) :: filename, quant


 namelist /coherence/ nstates,nrep,nsteps,quant

 read(*,nml=coherence)

 write(*,*) '**********************************************'
 write(*,*) '*                                            *'
 write(*,*) '*               Tool for a                   *'
 write(*,*) '*         global and quantitative            *'
 write(*,*) '*                analysis                    *'
 write(*,*) '*                   of                       *'
 write(*,*) '*            quantum coherence               *'
 write(*,*) '*                                            *'
 write(*,*) '*                   by                       *'
 write(*,*) '*              Stefano Corni                 *'
 write(*,*) '*                  and                       *'
 write(*,*) '*             Emanuele Coccia                *'
 write(*,*) '*                                            *'
 write(*,*) '**********************************************'

 write(*,*) ''
 write(*,*) 'Number of steps:', nsteps
 write(*,*) 'Number of trajectories:', nrep
 write(*,*) 'Number of states', nstates
 write(*,*) 'Quantifier:'
 if (quant.eq.'l1') then
    write(*,*) 'l1-norm' 
 else
    write(*,*) 'Other quantifiers not implemented yet'
    stop
 endif
 write(*,*) '' 

 if (nsteps.lt.1) then
    write(*,*) 'ERROR: nsteps must be larger than zero' 
    stop
 endif
 if (nrep.lt.1) then
    write(*,*) 'ERROR: nrep must be larger than zero'
    stop
 endif
 if (nstates.lt.2) then
    write(*,*) 'ERROR: nstates must be larger than zero'
    stop
 endif

 npair=nstates*(nstates-1)/2 

 allocate(i(nsteps))
 allocate(t(nsteps))
 if (quant.ne.'l1') then
    allocate(pop(nsteps,nstates),perr(nsteps,nstates))
    allocate(rdum(nsteps,nstates,nrep))
 elseif (quant.eq.'l1') then
    allocate(cor(nsteps,npair),coi(nsteps,npair))
    allocate(l1(nsteps),co(nsteps,npair))
    allocate(rdum(nsteps,npair,nrep))
    allocate(idum(nsteps,npair,nrep))
    allocate(rerr(nsteps,npair),ierr(nsteps,npair))
    allocate(l1_err(nsteps))
    allocate(ferr(nsteps,npair))
 endif

 if (quant(1:2).ne.'l1') then

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
    enddo

    !Read populations
    do m=1,nrep
       do j=1,nsteps
          read(20+m,*) i(j), t(j), (rdum(j,k,m), k=1,nstates)
       enddo
    enddo

    ! Population 
    pop=0.d0
    do m=1,nrep
       do k=1,nstates
          do j=1,nsteps
             pop(j,k) = pop(j,k) + rdum(j,k,m)
          enddo
       enddo
    enddo
    pop=pop/real(nrep)

 ! Error 
    perr=0.d0
    do m=1,nrep
       do k=1,nstates
          do j=1,nsteps
             perr(j,k) = perr(j,k) + (rdum(j,k,m) - pop(j,k))**2
          enddo
       enddo 
    enddo
    if (nrep.gt.1) then
       perr=perr/real(nrep-1.d0)
    else
       perr=perr/real(nrep-1.d0)
    endif
    perr=sqrt(perr)

  elseif (quant(1:2).eq.'l1') then

     do m=1,nrep
       if (m.lt.10) then
          WRITE(filename,'(a,i1.1,a)') "d_t_",m,".dat"
       elseif (m.ge.10.and.m.lt.100) then
          WRITE(filename,'(a,i2.2,a)') "d_t_",m,".dat"
       elseif (m.ge.100.and.m.lt.1000) then
          WRITE(filename,'(a,i3.3,a)') "d_t_",m,".dat"
       elseif (m.ge.1000.and.m.lt.10000) then
          WRITE(filename,'(a,i4.4,a)') "d_t_",m,".dat"
       elseif (m.ge.10000.and.m.lt.100000) then
          WRITE(filename,'(a,i5.5,a)') "d_t_",m,".dat"
       endif
       open(20+m,file=filename)
     enddo


     !Read coherences
     do m=1,nrep
        do j=1,nsteps
          read(20+m,*) i(j),t(j), (rdum(j,k,m), idum(j,k,m), k=1,npair)
        enddo
     enddo


     ! Mean
     cor=0.d0
     coi=0.d0
     do m=1,nrep
        do k=1,npair
           do j=1,nsteps
              cor(j,k) = cor(j,k) + rdum(j,k,m)
              coi(j,k) = coi(j,k) + idum(j,k,m)
           enddo
        enddo
     enddo
     cor=cor/real(nrep)
     coi=coi/real(nrep)

     ! Error 
     rerr=0.d0
     ierr=0.d0 
     do m=1,nrep
        do k=1,npair
           do j=1,nsteps
              rerr(j,k) = rerr(j,k) + (rdum(j,k,m) - cor(j,k))**2
              ierr(j,k) = ierr(j,k) + (idum(j,k,m) - coi(j,k))**2 
           enddo
        enddo
     enddo
     if (nrep.gt.1) then
        rerr=rerr/real(nrep-1.d0)
        ierr=ierr/real(nrep-1.d0)
     else
        rerr=rerr/real(nrep)
        ierr=ierr/real(nrep)
     endif
     rerr=sqrt(rerr)
     ierr=sqrt(ierr)

     ferr=0.d0
     do k=1,npair
        do j=1,nsteps
           co(j,k) = sqrt(cor(j,k)**2 + coi(j,k)**2)
           if (co(j,k).ne.0.d0) then
              ferr(j,k) = sqrt( (cor(j,k)/co(j,k)*rerr(j,k))**2 + (coi(j,k)/co(j,k)*ierr(j,k))**2 )
              ferr(j,k) = ferr(j,k)/sqrt(real(nrep)) 
           endif
        enddo
     enddo


     write(*,*) ''
     write(*,*) 'Coherence quantifier: l1-norm'
     write(*,*) 'C(rho) = \sum_i.ne.j |rho_ij|'
     write(*,*) ''

     l1=0.d0
     do k=1,npair
        do j=1,nsteps
           l1(j) = l1(j) + co(j,k)
        enddo
     enddo
     l1=2.d0*l1
 

     l1_err=0.d0
     do k=1,npair 
        do j=1,nsteps
           l1_err(j) = l1_err(j) + ferr(j,k)**2
        enddo 
     enddo
     l1_err=2.d0*sqrt(l1_err)
 
     open(11,file='l1_norm')
     write(11,*) '# step     time(au)    C_l1(rho(t))     error(t)'
     do j=1,nsteps
        write(11,'(I8,3(E12.5))') i(j),t(j),l1(j),l1_err(j) 
     enddo

  endif

  deallocate(i) 
  deallocate(t)
  deallocate(rdum)
  if (quant.ne.'l1') then
     deallocate(perr)
     deallocate(pop)
  elseif (quant.eq.'l1') then
     deallocate(idum)
     deallocate(cor)
     deallocate(coi)
     deallocate(co)
     deallocate(l1)
     deallocate(rerr)
     deallocate(ierr)
     deallocate(ferr)
     deallocate(l1_err)
  endif

  close(11)

  do m=1,nrep
     close(20+m)
  enddo

  write(*,*) 'Calculation ended'

  stop

end program decoherence 
