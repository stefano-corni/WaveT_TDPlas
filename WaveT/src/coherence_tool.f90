program decoherence 

 use constants 

 implicit none

 integer(i4b)                 :: nrep,nstates,nsteps,npair,ngs
 real(dbl),    allocatable    :: rdum(:,:,:),pop(:,:),perr(:,:)
 real(dbl),    allocatable    :: cor(:,:),coi(:,:),co(:,:)
 real(dbl),    allocatable    :: rerr(:,:),ierr(:,:),rpop(:,:,:)
 integer(i4b), allocatable    :: i(:)
 real(dbl),    allocatable    :: t(:),l1(:),l1_err(:)
 real(dbl),    allocatable    :: idum(:,:,:),ferr(:,:)
 real(dbl),    allocatable    :: trp2(:),trp2_err(:)
 real(dbl),    allocatable    :: rho(:,:,:),rho2(:,:,:)
 real(dbl),    allocatable    :: leps(:), l1_norm_eps(:) 
 real(dbl),    allocatable    :: leps_err(:),l1_norm_eps_err(:)
 real(dbl),    allocatable    :: rho_err(:,:,:),rho2_err(:,:,:)


 integer(i4b)       :: j,k,m,ijunk,l
 real(dbl)          :: rjunk,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmperr  
 character(30)      :: filename, quant


 namelist /coherence/ nstates,nrep,nsteps,quant,ngs

 ngs=1
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
 write(*,*) '*                  and                       *'
 write(*,*) '*             Filippo Troiani                *'
 write(*,*) '*                                            *'
 write(*,*) '**********************************************'

 write(*,*) ''
 write(*,*) 'Number of steps:', nsteps
 write(*,*) 'Number of trajectories:', nrep
 write(*,*) 'Number of states', nstates
 write(*,*) 'Quantifier:'
 if (quant.eq.'l1') then
    write(*,*) 'l1-norm' 
 elseif (quant.eq.'le') then
    write(*,*) 'linear entropy'
 else
    write(*,*) 'Other quantifiers not implemented yet'
    stop
 endif
 if (ngs.gt.1) then
    write(*,*) 'Degenerate ground state:', ngs
 elseif (ngs.le.0) then
    write(*,*) 'ERROR: ngs must be large than zero'
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
 allocate(rho(nsteps,nstates,nstates))
 allocate(rho2(nsteps,nstates,nstates))
 allocate(rpop(nsteps,nstates,nrep))
 allocate(rdum(nsteps,npair,nrep))
 allocate(idum(nsteps,npair,nrep))
 allocate(pop(nsteps,nstates),perr(nsteps,nstates))
 allocate(rerr(nsteps,npair),ierr(nsteps,npair))
 allocate(ferr(nsteps,npair)) 
 allocate(cor(nsteps,npair),coi(nsteps,npair))
 allocate(co(nsteps,npair))
 allocate(leps(nsteps))
 allocate(leps_err(nsteps))
 allocate(trp2(nsteps),trp2_err(nsteps))
 allocate(rho_err(nsteps,nstates,nstates))
 allocate(rho2_err(nsteps,nstates,nstates))
 allocate(l1(nsteps))
 allocate(l1_err(nsteps))
 allocate(l1_norm_eps(nsteps))
 allocate(l1_norm_eps_err(nsteps))

 ! Population files
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
       read(20+m,*) i(j), t(j), (rpop(j,k,m), k=1,nstates)
    enddo
 enddo

 !Average population 
 pop=0.d0
 do m=1,nrep
    do k=1,nstates
       do j=1,nsteps
          pop(j,k) = pop(j,k) + rpop(j,k,m)
       enddo
    enddo
 enddo
 pop=pop/dble(nrep)

 !Building density matrix
 rho=0.d0 
 do k=1,nstates
    do j=1,nsteps
       rho(j,k,k) = pop(j,k)
    enddo
 enddo

 !Error on population 
 perr=0.d0
 do m=1,nrep
    do k=1,nstates
       do j=1,nsteps
          perr(j,k) = perr(j,k) + (rpop(j,k,m) - pop(j,k))**2
       enddo
    enddo
 enddo
 if (nrep.gt.1) then
    perr=perr/dble(nrep-1.d0)
 else
    perr=perr/dble(nrep)
 endif
 perr=dsqrt(perr)
 perr=perr/dsqrt(dble(nrep))


 !Coherences
 do m=1,nrep
    close(20+m)
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
       read(20+m,*) ijunk,rjunk, (rdum(j,k,m), idum(j,k,m), k=1,npair)
    enddo
 enddo

 do m=1,nrep
     close(20+m)
 enddo

 !Average coherence 
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
 cor=cor/dble(nrep)
 coi=coi/dble(nrep)

 !Error on coherence 
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
    rerr=rerr/dble(nrep-1.d0)
    ierr=ierr/dble(nrep-1.d0)
  else
    rerr=rerr/dble(nrep)
    ierr=ierr/dble(nrep)
  endif
  rerr=dsqrt(rerr)
  ierr=dsqrt(ierr)
  rerr=rerr/dsqrt(dble(nrep))
  ierr=ierr/dsqrt(dble(nrep))

  ferr=0.d0
  do k=1,npair
     do j=1,nsteps
        co(j,k) = dsqrt(cor(j,k)**2 + coi(j,k)**2)
        if (co(j,k).ne.0.d0) then
           ferr(j,k) = dsqrt( (cor(j,k)/co(j,k)*rerr(j,k))**2 + (coi(j,k)/co(j,k)*ierr(j,k))**2 )
        endif
     enddo
  enddo

  !Coherence from WaveT: (1,i=2,n_ci) (j=2,n_ci,k=j+1,n_ci)
  !Buiding density matrix
  k=0
  do m=1,nstates
     do l=m+1,nstates
        k=k+1
        do j=1,nsteps
           rho(j,l,m) = co(j,k) 
           rho(j,m,l) = co(j,k)
        enddo
     enddo
  enddo

  !Epsilon = Tr(Pe*rho*Pe)
  !Error on epsilon 
  leps_err=0.d0
  do k=ngs+1,nstates
     do j=1,nsteps
        leps(j) = leps(j) + rho(j,k,k)
        leps_err(j) = leps_err(j) + perr(j,k)**2
     enddo
  enddo
  leps_err=dsqrt(leps_err)

  open(11,file='eps')
  write(11,*) '# step     time(au)    eps(rho(t))     error(t)'
  do j=1,nsteps
     write(11,'(I8,3(E14.7))') i(j),t(j),leps(j),leps_err(j)
  enddo
  close(11)


  if (quant(1:2).eq.'le') then

    !Building rho^2
    do j=1,nsteps
       rho2(j,:,:)=matmul(rho(j,:,:),rho(j,:,:))
    enddo

    !Error on density matrix
    do k=1,nstates
       do j=1,nsteps
          rho_err(j,k,k) = perr(j,k)
       enddo
    enddo

    k=0
    do m=1,nstates
       do l=m+1,nstates
          k=k+1
          do j=1,nsteps
             rho_err(j,l,m) = ferr(j,k)
             rho_err(j,m,l) = ferr(j,k)
          enddo
       enddo
    enddo

    do j=1,nsteps
       rho2_err(j,:,:)=matmul(rho_err(j,:,:),rho_err(j,:,:))
    enddo

    !Building Tr(rho^2)
    trp2=0.d0
    do k=1,nstates
       do j=1,nsteps
          trp2(j) = trp2(j) + rho2(j,k,k)
       enddo
    enddo
 
    do j=1,nsteps
       if (trp2(j).gt.1.d0) trp2(j)=1.d0 
    enddo

    write(*,*) ''
    write(*,*) 'Coherence quantifier: linear entropy'
    write(*,*) 'S(rho) = 1 - Tr(rho^2)'
    write(*,*) ''


    trp2_err=0.d0
    do k=1,nstates
       do m=1,nstates
          do j=1,nsteps
             trp2_err(j) = trp2_err(j) + (2.d0*rho(j,k,k)*rho_err(j,k,k))**2 
             if (m.ne.k) trp2_err(j) = trp2_err(j) + 2.d0*(rho(j,m,k)*rho_err(j,m,k))**2
          enddo
       enddo
    enddo
    trp2_err=dsqrt(trp2_err)


    open(11,file='lin_entropy')
    write(11,*) '# step     time(au)    S(rho(t))     error(t)'
    do j=1,nsteps
       write(11,'(I8,3(E14.7))') i(j),t(j),1.d0-trp2(j),trp2_err(j)
    enddo
    close(11)

    open(11,file='rho2_eps')
    write(11,*) '# step     time(au)    Tr(rho2(eps))     error(t)'
    tmp4=dble(nstates)/dble(nstates-1)
    do j=1,nsteps
       tmp5 = 1.d0 - 2.d0*leps(j) + leps(j)**2*tmp4
       tmperr = 2.d0*(-1.d0 + leps(j)*tmp4)*leps_err(j)
       tmperr=tmperr**2
       tmperr=dsqrt(tmperr)
       write(11,'(I8,3(E14.7))') i(j),t(j),tmp5,tmperr
    enddo
    close(11)


  elseif (quant(1:2).eq.'l1') then

     tmp=dble(ngs-1)
     tmp1=dble(nstates-ngs-1)
     tmp3=tmp1+1.d0

     !Maximum l1-norm value for a given epsilon (Tr(Pe*rho*Pe)) 
     do j=1,nsteps
        tmp2=1.d0-leps(j)
        l1_norm_eps(j) = tmp*tmp2 + tmp1*leps(j) + 2.d0*dsqrt(ngs*tmp3*leps(j)*tmp2)
     enddo

     do j=1,nsteps
        tmp2=1.d0-leps(j)
        if (leps_err(j).ne.0.d0) then
           l1_norm_eps_err(j) = leps_err(j)*(-tmp + tmp1 +(ngs*tmp3)*(1.d0-2.d0*leps(j))/dsqrt(ngs*tmp3*leps(j)*tmp2))
           l1_norm_eps_err(j) = l1_norm_eps_err(j)**2 
        endif
     enddo
     l1_norm_eps_err=dsqrt(l1_norm_eps_err)


     open(11,file='l1_norm_eps')
     write(11,*) '# step     time(au)    l1_norm_eps(rho(t))     error(t)'
     do j=1,nsteps
        write(11,'(I8,3(E14.7))') i(j),t(j),l1_norm_eps(j),l1_norm_eps_err(j)
     enddo
     close(11)


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
     l1_err=2.d0*dsqrt(l1_err)
 
     open(11,file='l1_norm')
     write(11,*) '# step     time(au)    C_l1(eps)     error(t)    C_l1(rho(t)  error(t) '
     do j=1,nsteps
        if (l1_norm_eps(j).ne.0.d0) then
           tmperr=(l1_err(j)/l1_norm_eps(j))**2 + (l1(j)/(l1_norm_eps(j))**2*l1_norm_eps_err(j))**2
           tmperr=dsqrt(tmperr)
        else
           tmperr=0.d0
        endif
        if (l1_norm_eps(j).ne.0.d0) then
           write(11,'(I8,5(E14.7))') i(j),t(j),l1(j)/l1_norm_eps(j),tmperr,l1(j), l1_err(j) 
        else
           write(11,'(I8,5(E14.7))') i(j),t(j),l1(j),tmperr,l1(j), l1_err(j)
        endif
     enddo

  endif

  deallocate(i) 
  deallocate(t)
  deallocate(rdum)
  deallocate(idum)
  deallocate(rpop)
  deallocate(rho)
  deallocate(rho2)
  deallocate(perr)
  deallocate(pop)
  deallocate(rerr)
  deallocate(ierr)
  deallocate(ferr)
  deallocate(cor)
  deallocate(coi)
  deallocate(co)
  deallocate(leps)
  deallocate(leps_err)
  deallocate(trp2)
  deallocate(trp2_err) 
  deallocate(rho_err)
  deallocate(rho2_err)
  deallocate(l1)
  deallocate(l1_err)
  deallocate(l1_norm_eps)
  deallocate(l1_norm_eps_err)

  close(11)

  write(*,*) 'Calculation ended'

  stop

end program decoherence 
