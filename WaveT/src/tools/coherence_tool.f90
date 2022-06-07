program decoherence 

 use constants 

 implicit none

 integer(i4b)                 :: nrep,nstates,nsteps,npair,ngs
 real(dbl),    allocatable    :: rdum(:,:,:),pop(:,:),perr(:,:)
 real(dbl),    allocatable    :: cor(:,:),coi(:,:),co(:,:)
 real(dbl),    allocatable    :: rerr(:,:),ierr(:,:)
 real(dbl),    allocatable    :: rc(:), ic(:)
 integer(i4b), allocatable    :: i(:)
 real(dbl),    allocatable    :: t(:),l1(:),l1_err(:)
 real(dbl),    allocatable    :: ferr(:,:)
 real(dbl),    allocatable    :: trp2(:),trp2_err(:)
 real(dbl),    allocatable    :: rho(:,:,:)
 real(dbl),    allocatable    :: leps(:), l1_norm_eps(:) 
 real(dbl),    allocatable    :: leps_err(:),l1_norm_eps_err(:)
 real(dbl),    allocatable    :: rho_err(:,:,:)
 real(dbl),    allocatable    :: opop(:)
 complex(cmp), allocatable    :: cc(:,:,:)
 complex(cmp), allocatable    :: c(:,:,:)


 integer(i4b)       :: j,k,m,ijunk,l
 integer            :: st,current,rate
 real(dbl)          :: rjunk,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmperr  
 character(30)      :: filename, quant
 character(4000)    :: dum,fmt_ci,fmt_ci2
 character(1)       :: read_bin


 namelist /coherence/ nstates,nrep,nsteps,quant,ngs,read_bin

 call system_clock(st,rate)

 ngs=1
 read_bin='n'

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
 if (read_bin.eq.'n') then
    write(*,*) 'Use formatted coefficient files'
 elseif (read_bin.eq.'y') then 
    write(*,*) 'use unformatted coefficient files' 
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
 allocate(opop(nsteps))
 allocate(rho(nsteps,nstates,nstates))
 allocate(c(nsteps,nstates,nrep))
 allocate(rc(nstates),ic(nstates))
 allocate(rdum(nsteps,npair,nrep))
 allocate(cc(nsteps,npair,nrep))
 allocate(pop(nsteps,nstates),perr(nsteps,nstates))
 allocate(rerr(nsteps,npair),ierr(nsteps,npair))
 allocate(ferr(nsteps,npair)) 
 allocate(cor(nsteps,npair),coi(nsteps,npair))
 allocate(co(nsteps,npair))
 allocate(leps(nsteps))
 allocate(leps_err(nsteps))
 allocate(trp2(nsteps),trp2_err(nsteps))
 allocate(rho_err(nsteps,nstates,nstates))
 allocate(l1(nsteps))
 allocate(l1_err(nsteps))
 allocate(l1_norm_eps(nsteps))
 allocate(l1_norm_eps_err(nsteps))

 ! Coefficient files
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

 if (read_bin(1:1).ne.'y') then
    do m=1,nrep
       open (20+m,file=filename,status="unknown")
       read(20+m,*) dum
       write (fmt_ci,'("(i8,f14.4,",I0,"e17.8E3)")') 2*nstates
       do j=1,nsteps
          read(20+m,fmt_ci) i(j), t(j), (rc(k), ic(k), k=1,nstates)
          do k=1,nstates
             c(j,k,m) = dcmplx(rc(k),ic(k))
          enddo
       enddo
    enddo
 else
    do m=1,nrep
       open (20+m,file=filename,status="unknown",form="unformatted")
       do j=1,nsteps
          read(20+m) i(j),t(j),(rc(k),ic(k),k=1,nstates)
          do k=1,nstates
             c(j,k,m) = dcmplx(rc(k),ic(k))
          enddo   
       enddo
    enddo                                                                                       
 endif   

 !Average population 
 pop=0.d0
 do m=1,nrep
    do k=1,nstates
       do j=1,nsteps
          pop(j,k) = pop(j,k) + conjg(c(j,k,m))*real(c(j,k,m))
       enddo
    enddo
 enddo
 pop=pop/dble(nrep)

 do j=1,nsteps
    opop(j) = pop(j,1)
 enddo

 pop(:,1)=0.d0
 do j=1,nsteps
    pop(j,1) = 1.d0 - sum(pop(j,2:nstates))
 enddo

 !Building density matrix
 !Diagonal elements
 rho=0.d0
 do k=1,nstates
    do j=1,nsteps
       rho(j,k,k) = pop(j,k)
    enddo
 enddo

 !Error on population 
 perr=0.d0
 do m=1,nrep
    do k=2,nstates
       do j=1,nsteps
          perr(j,k) = perr(j,k) + (conjg(c(j,k,m))*real(c(j,k,m))  - pop(j,k))**2
       enddo
    enddo
 enddo
 do m=1,nrep
    do j=1,nsteps
       perr(j,1) = perr(j,1) + (conjg(c(j,1,m))*real(c(j,1,m)) - opop(j))**2
    enddo
 enddo
 if (nrep.gt.1) then
    perr=perr/dble(nrep-1.d0)
 else
    perr=perr/dble(nrep)
 endif
 perr=dsqrt(perr)
 perr=perr/dsqrt(dble(nrep))


 !Define coherences 
 do m=1,nrep
    do j=1,nsteps
       call compute_coherence(c(j,:,m),cc(j,:,m),nstates,npair)
    enddo
 enddo

 !Average coherence 
 cor=0.d0
 coi=0.d0
 do m=1,nrep
    do k=1,npair
       do j=1,nsteps
          cor(j,k) = cor(j,k) + real(cc(j,k,m))
          coi(j,k) = coi(j,k) + aimag(cc(j,k,m))
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
          rerr(j,k) = rerr(j,k) + (real(cc(j,k,m)) - cor(j,k))**2
          ierr(j,k) = ierr(j,k) + (aimag(cc(j,k,m)) - coi(j,k))**2
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
  !Off-diagonal elements
  k=0
  do m=1,nstates
     do l=m+1,nstates
        k=k+1
        do j=1,nsteps
           rho(j,l,m) = co(j,k) 
           rho(j,m,l) = co(j,k)
           if (j.eq.1.or.m.eq.1) then
              rho(j,l,m) = rho(j,l,m)*dsqrt(pop(j,1)/opop(j))
              rho(j,m,l) = rho(j,l,m)
           endif
           !rho(j,l,m) = cmplx(cor(j,k),coi(j,k))
           !rho(j,m,l) = conjg(rho(j,l,m))!cmplx(cor(j,k),-coi(j,k))
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

    !Building Tr(rho^2)
    trp2=0.d0
    do k=1,nstates
       do j=1,nsteps
          trp2(j) = trp2(j) + pop(j,k)**2
       enddo
    enddo    
    do l=1,nstates
       do k=1,nstates
          do j=1,nsteps
             if (l.ne.k) then
                trp2(j) = trp2(j) + rho(j,k,l)**2
             endif
          enddo
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
             !if (m.ne.k) trp2_err(j) = trp2_err(j) + 2.d0*(rho(j,m,k)*rho_err(j,m,k))**2
             if (m.ne.k) trp2_err(j) = trp2_err(j) + (2.d0*rho(j,m,k)*rho_err(j,m,k))**2
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
  deallocate(opop)
  deallocate(rdum)
  deallocate(rho)
  deallocate(perr)
  deallocate(pop)
  deallocate(rerr)
  deallocate(ierr)
  deallocate(ferr)
  deallocate(cor)
  deallocate(coi)
  deallocate(cc)
  deallocate(co)
  deallocate(leps)
  deallocate(leps_err)
  deallocate(trp2)
  deallocate(trp2_err) 
  deallocate(rho_err)
  deallocate(l1)
  deallocate(l1_err)
  deallocate(l1_norm_eps)
  deallocate(l1_norm_eps_err)

  close(11)

  call system_clock(current)
  write(*,*) ''
  write(6,'("Done , total elapsed time", &
           F10.3,"s")') real(current-st)/real(rate)

  stop

end program decoherence 


!------------------------------------------------------------------------
! @brief Compute C*_iC_j (i.ne.j)  
! 
! @date Created   : E. Coccia 22 Mar 2019
! Modified  :
!------------------------------------------------------------------------
subroutine compute_coherence(c,cc,nci,npair)

        use constants

        implicit none

        integer(i4b),    intent(in)    :: nci,npair
        complex(cmp),    intent(in)    :: c(nci)
        complex(cmp),    intent(out)   :: cc(npair)
        integer(i4b)                   :: j,k,kk
        complex(cmp)                   :: tmp

        !tmp=dcmplx(0.d0,0.d0)

        kk=0
        do j=1,nci
           do k=j+1,nci
              kk=kk+1
              cc(kk) = dconjg(c(k))*c(j)
           enddo
        enddo

        !do k=2,nci
        !   tmp = dconjg(c(k))*c(1)
        !   cc(1,k) = tmp
        !enddo

        !if (nci.gt.2) then
        !   kk=0
        !   do j=2,nci
        !      do k=j+1,nci
        !         kk=kk+1
        !         tmp = dconjg(c(k))*c(j)
        !         cc(j,kk) = tmp
        !      enddo
        !   enddo
        !endif

        return

end subroutine compute_coherence

