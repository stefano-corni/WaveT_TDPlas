PROGRAM write_wavet 
  
#ifdef OMP
  use omp_lib
#endif
  implicit none


  real*8,  allocatable :: dipmatx(:,:) , dipmaty(:,:) , dipmatz(:,:)
  real*8,  allocatable :: lmatx(:,:) , lmaty(:,:) , lmatz(:,:)
  real*8,  allocatable :: potmat(:,:,:), potmat_nuc(:,:,:), potmut(:,:), pot0(:), pot_nuc0(:)
  real*8,  allocatable :: exciten(:), tddfteig(:,:,:), tddfteigl(:,:,:)
  real*8,  allocatable :: dmut(:,:),norm(:), norml(:), lmut(:,:)
  real*8               :: dip0(3),dip_nuc(3),l0(3)

  real*8,  PARAMETER   ::  EVAU = 27.21139628D0 


  integer,  allocatable   :: cont(:,:)
  integer   :: ia,ib,i,j,a,itoten,iener,isym,is,js,jj,ints,kk,dd
  integer   :: kvirt,kocc,ntoten,ntotmo,nener,nsym,nthreads,nts,ntotentr,threads

!MPI_VAR
!  integer   :: ierr, nprocs, myid

  real      :: start, finish
  real      :: start_omp, finish_omp
  integer   :: rate, st, current

  character*20  :: cdum
  
  logical :: cdspectrum
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIER_e_LEO!!!!!!!!!!!!!!
  logical :: NP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  call MPI_INIT(ierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
!  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

#ifdef OMP
   nthreads=omp_get_max_threads( )
#endif
#ifndef OMP
   nthreads=1
#endif

#ifdef OMP
   write(*,*) '*************OMP*****************'
   write(*,'(a)')
   write(*,*) 'OMP parallelization on'
   write(*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs ( )
   write(*,'(a,i8)') 'The number of threads available    = ', nthreads 
   write(*,'(a)')
   write(*,*) '*************OMP*****************'
#endif

  call system_clock(st,rate)
  call CPU_TIME(start)
  open(80,file='eig.dat', form='unformatted')
  read(80) nsym,kocc,kvirt,ntoten,ntotmo,cdspectrum,NP
  write(*,*) 'sono arrivato dopo eig.dat'
  if (NP) then
  read(80) nts
  write(*,*) 'num tessere', nts
  endif
  ntotentr = ((ntoten+1)*(ntoten+2))/2
  if (cdspectrum) open(81,file='eig_l.dat')
  allocate(dipmatx(ntotmo,ntotmo))
  allocate(dipmaty(ntotmo,ntotmo))
  allocate(dipmatz(ntotmo,ntotmo))
  allocate(cont(ntoten+1,ntoten+1))
  write(*,*) 'sono arrivato dopo l allocazione di dipmat, questo ntotmo',ntotmo
  if (cdspectrum) then
     allocate(lmatx(ntotmo,ntotmo))
     allocate(lmaty(ntotmo,ntotmo))
     allocate(lmatz(ntotmo,ntotmo))
     lmatx=0.d0
     lmaty=0.d0
     lmatz=0.d0
  allocate(tddfteigl(kvirt,kocc,ntoten))
  allocate(lmut(3,ntotentr))
  allocate(norml(ntoten))
  endif
  write(*,*) 'numero orbitali virtuali', kvirt, 'num orb occ', kocc, 'num tot stati',ntoten
  
  dipmatx=0.d0
  dipmaty=0.d0
  dipmatz=0.d0
  
  cont=0.d0
  kk=0.d0
  do is=1,ntoten+1
     do js=is,ntoten+1
        kk=kk+1
        cont(is,js)=kk
        cont(js,is)=cont(is,js)
     enddo
 enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!matrice negli stati molecolari

  if (NP) then 
     allocate(potmat(ntotmo,ntotmo,nts))
     allocate(potmut(ntotentr,nts))
     allocate(pot0(nts))
     allocate(pot_nuc0(nts))
  endif 
  potmat=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(tddfteig(kvirt,kocc,ntoten))
  allocate(exciten(ntoten))
  allocate(norm(ntoten))
  allocate(dmut(3,ntotentr))
  
  itoten=0
  do isym=1,nsym
     read(80) nener
     do iener=1,nener
        itoten=itoten+1
        do i=1,kocc
           do j=1,kvirt
              read(80) tddfteig(j,i,itoten)
           enddo
        enddo
     enddo
  enddo
  close(80)

  if (cdspectrum) then
     itoten=0
     do isym=1,nsym
        do iener=1,nener
           itoten=itoten+1
           do i=1,kocc
              do j=1,kvirt
                 read(81,*) tddfteigl(j,i,itoten)
              enddo
           enddo
        enddo
     enddo
  close(81)
  endif

  open(70,file='ene.dat')
  do i=1,ntoten
     read(70,*) exciten(i)
  enddo

  itoten=0
  do isym=1,nsym
     do iener=1,nener
        itoten=itoten+1
        norm(itoten)=0.d0
        start_omp = omp_get_Wtime()
#ifdef OMP 
!$OMP PARALLEL DO REDUCTION(+:norm)
#endif 
        do i=1,kocc
           do j=1,kvirt
               norm(itoten)=norm(itoten)+tddfteig(j,i,itoten)**2
           enddo
        enddo
#ifdef OMP 
!$OMP END PARALLEL DO
#endif
     enddo
  enddo
   
  if (cdspectrum) then
     itoten=0
     do isym=1,nsym
        do iener=1,nener
           itoten=itoten+1
           norml(itoten)=0.d0
#ifdef OMP
!$OMP PARALLEL DO REDUCTION(+:norml)
#endif
           do i=1,kocc
              do j=1,kvirt
                 norml(itoten)=norml(itoten)+tddfteigl(j,i,itoten)**2
              enddo
           enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
        enddo
     enddo
  endif


  open(72,file='dip_nuc.dat')
  read(72,*) cdum
  read(72,*) dip_nuc(1),dip_nuc(2),dip_nuc(3)
  close(72)

  open(90,file='dipmat.dat')
  do i=1,ntotmo
     do j=1,ntotmo
        read(90,*) dipmatx(j,i),dipmaty(j,i),dipmatz(j,i)
     enddo
  enddo
  close(90)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (NP) then 
     open(92, file='potmat.dat')
   do ints=1,nts
     do i=1,ntotmo
        do j=1,ntotmo
              read(92,*) potmat(j,i,ints)
           enddo
        enddo
     enddo
     close(92)
  endif

  if (NP) then
     open(93, file='potmat_nuc.dat')
     do ints=1,nts
!        do i=1,ntotmo
!           do j=1,ntotmo
!                 read(93,*) potmat_nuc(j,i,ints)
         read(93,*) pot_nuc0(ints)
!              enddo
!           enddo
        enddo
     close(93)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
  if (cdspectrum) then
     open(91,file='lmat.dat')
     do i=1,ntotmo
        do j=1,ntotmo
           read(91,*) lmatx(j,i),lmaty(j,i),lmatz(j,i)
        enddo
     enddo
     close(91)
  endif

  !GS-GS
  dip0=0.d0
#ifdef OMP
!$OMP PARALLEL DO REDUCTION(+:dip0)
#endif
  do j=1,kocc
     dip0(1) = dip0(1) + dipmatx(j,j) 
     dip0(2) = dip0(2) + dipmaty(j,j)
     dip0(3) = dip0(3) + dipmatz(j,j)
  enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif

  dmut=0.d0
  dmut(:,1)=2.d0*dip0(:)


  if (NP) then 
  pot0=0.d0
#ifdef OMP
!$OMP  PARALLEL DO REDUCTION(+:pot0)
#endif
  do j=1,kocc
     do ints=1,nts
        pot0(ints) = pot0(ints) + potmat(j,j,ints)
     enddo
  enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
  
!   pot_nuc0=0.d0
  
!#ifdef OMP
!!$OMP PARALLEL DO REDUCTION(+:pot_nuc0)
!#endif
!  do j=1,kocc
!     do ints=1,nts
!        pot_nuc0(ints) = pot_nuc0(ints) + potmat_nuc(1,1,ints)
!     enddo
!  enddo
!#ifdef OMP
!!$OMP END PARALLEL DO
!#endif

  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (NP) then
     potmut = 0.d0
     !pot_nuc0 = 0.d0
     potmut(1,:)=2.d0*pot0(:)
     !potmut_nuc(1,:)=2.d0*pot_nuc0(:)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (cdspectrum) then
     l0=0.d0
#ifdef OMP
!$OMP PARALLEL DO REDUCTION(+:l0)
#endif
     do j=1,kocc
        l0(1) = l0(1) + lmatx(j,j)
        l0(2) = l0(2) + lmaty(j,j)
        l0(3) = l0(3) + lmatz(j,j)
     enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
     lmut(:,1)=2.d0*l0(:)
  endif
  
  
  !GS-EXC
#ifdef OMP
!$OMP PARALLEL DO  
#endif
  do is=1,ntoten
      kk = cont(1,is+1)
!     kk = dd
      do i=1,kocc
        do ia=1,kvirt
!#ifdef OMP
!            !$OMP ATOMIC
!#endif
           dmut(1,kk)  = dmut(1,kk) + tddfteig(ia,i,is)*dipmatx(i,kocc+ia)
           dmut(2,kk)  = dmut(2,kk) + tddfteig(ia,i,is)*dipmaty(i,kocc+ia) 
           dmut(3,kk)  = dmut(3,kk) + tddfteig(ia,i,is)*dipmatz(i,kocc+ia)
        enddo
     enddo
  enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
write(*,*) "Sono dopo dmut"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (NP) then
#ifdef OMP
!$OMP PARALLEL DO
#endif 
  do ints=1,nts
     do is=1,ntoten
     
      kk = cont(1,is+1)
!     dd = kk
        do i=1,kocc
           do ia=1,kvirt
!#ifdef OMP
!              !$OMP ATOMIC
!#endif
              potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*potmat(i,kocc+ia,ints)
           enddo
        enddo
      enddo
   enddo
#ifdef OMP 
!$OMP END PARALLEL DO
#endif
 
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (cdspectrum) then
#ifdef OMP
!$OMP PARALLEL DO 
#endif
    do is=1,ntoten 
      kk = cont(1,is+1)
!     dd = kk
        do i=1,kocc
           do ia=1,kvirt
!#ifdef OMP
!              !$OMP ATOMIC 
!#endif
              lmut(1,kk) = lmut(1,kk) + tddfteigl(ia,i,is)*lmatx(i,kocc+ia)
              lmut(2,kk) = lmut(2,kk) + tddfteigl(ia,i,is)*lmaty(i,kocc+ia)
              lmut(3,kk) = lmut(3,kk) + tddfteigl(ia,i,is)*lmatz(i,kocc+ia)
           enddo
        enddo
     enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
  endif
  

 

  itoten=0
  do isym=1,nsym
     do iener=1,nener
        itoten=itoten+1
       do i=1,kocc
           do j=1,kvirt
               tddfteig(j,i,itoten)=tddfteig(j,i,itoten)/sqrt(norm(itoten))
           enddo
       enddo
     enddo
  enddo

  if (cdspectrum) then
     itoten=0
     do isym=1,nsym
        do iener=1,nener
           itoten=itoten+1
           do i=1,kocc
              do j=1,kvirt
                 tddfteigl(j,i,itoten)=tddfteigl(j,i,itoten)/sqrt(norml(itoten))
              enddo
           enddo
        enddo
     enddo
  endif

!EXC-EXC
#ifdef OMP
!$OMP PARALLEL DO 
#endif
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
!       kk = dd
        do i=1,kocc
           do ia=1,kvirt
!#ifdef OMP            
!              !$OMP ATOMIC
!#endif
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(1) + dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(2) + dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(3) + dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              !do ib=1,ia-1
              !   dmut(1,is+1,js+1) = dmut(1,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ia,kocc+ib)
              !   dmut(2,is+1,js+1) = dmut(2,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ia,kocc+ib)
              !   dmut(3,is+1,js+1) = dmut(3,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ia,kocc+ib)
              !enddo
              !do ib=ia+1,kvirt
              !   dmut(1,is+1,js+1) = dmut(1,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ia,kocc+ib)
              !   dmut(2,is+1,js+1) = dmut(2,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ia,kocc+ib)
              !   dmut(3,is+1,js+1) = dmut(3,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ia,kocc+ib)
              !enddo
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
!#ifdef OMP
!                !$OMP ATOMIC
!#endif
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ia,kocc+ib)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ia,kocc+ib)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ia,kocc+ib)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
!#ifdef OMP
!                !$OMP ATOMIC 
!#endif
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

              !do j=1,i-1
              !   dmut(1,is+1,js+1) = dmut(1,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
              !   dmut(2,is+1,js+1) = dmut(2,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
              !   dmut(3,is+1,js+1) = dmut(3,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              !enddo
              !do j=i+1,kocc
              !   dmut(1,is+1,js+1) = dmut(1,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
              !   dmut(2,is+1,js+1) = dmut(2,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
              !   dmut(3,is+1,js+1) = dmut(3,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              !enddo
           enddo
        enddo
     enddo
  enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif 
!EXC-EXC
write(*,*) "Sono dopo gli stati eccitati"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (NP) then
#ifdef OMP
!$OMP PARALLEL DO
#endif
     do ints=1,nts
        do is=1,ntoten
           do js=is,ntoten
              kk = cont(is+1,js+1)
!             dd = kk
              do i=1,kocc
                 do ia=1,kvirt
!#ifdef OMP
!                       !$OMP ATOMIC
!#endif                      
                       potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*pot0(ints) + potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
                       do ib=1,kvirt
                          if (ib.eq.ia) cycle
!#ifdef OMP
!                          !$OMP ATOMIC
!#endif
                          potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*potmat(kocc+ia,kocc+ib,ints)

                       enddo
                       do j=1,kocc
                          if (j.eq.i) cycle
!#ifdef OMP
!                          !$OMP ATOMIC
!#endif 
                          potmut(kk,ints) = potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*potmat(j,i,ints)
                       enddo
                    enddo
                 enddo
           enddo
        enddo
     enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (cdspectrum) then
#ifdef OMP
!$OMP PARALLEL DO
#endif
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
!         dd = kk
          do i=1,kocc
             do ia=1,kvirt
!#ifdef OMP
!               !$OMP ATOMIC  
!#endif
                lmut(1,kk) = lmut(1,kk) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(1) + lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                lmut(2,kk) = lmut(2,kk) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(2) + lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                lmut(3,kk) = lmut(3,kk) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(3) + lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                do ib=1,kvirt
                   if (ib.eq.ia) cycle
!#ifdef OMP
!                  !$OMP ATOMIC
!#endif
                   lmut(1,kk) = lmut(1,kk) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmatx(kocc+ia,kocc+ib)
                   lmut(2,kk) = lmut(2,kk) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmaty(kocc+ia,kocc+ib)
                   lmut(3,kk) = lmut(3,kk) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmatz(kocc+ia,kocc+ib)
                enddo
                
                do j=1,kocc
                   if (j.eq.i) cycle
!#ifdef OMP
!                   !$OMP ATOMIC 
!#endif
                   lmut(1,kk) = lmut(1,kk) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmatx(j,i)
                   lmut(2,kk) = lmut(2,kk) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmaty(j,i)
                   lmut(3,kk) = lmut(3,kk) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmatz(j,i)
                enddo
             enddo
          enddo
       enddo
    enddo
#ifdef OMP
!$OMP END PARALLEL DO
#endif
  finish_omp = omp_get_Wtime()
  endif

  deallocate(tddfteig)
  deallocate(dipmatx)
  deallocate(dipmaty)
  deallocate(dipmatz)

  if (NP) then
    deallocate(potmat)
  endif 

  if (cdspectrum) then
     deallocate(lmatx)
     deallocate(lmaty)
     deallocate(lmatz)
     deallocate(tddfteigl)
  endif


  exciten=exciten*EVAU
  write(*,*) "Sto per scrivere"

  open(24,file='ci_energy.inp')
  do i=1,ntoten
     write(24,'("Root", I5, " : ", F12.5)') i, exciten(i)
  enddo
  close(24)

  deallocate(exciten)

  open(15,file='ci_mut.inp')


  dmut(:,cont(1,1))=dmut(:,cont(1,1))-dip_nuc(:)
  do i=2,ntoten+1
     dd=cont(i,i)
     dmut(1,dd)=dmut(1,dd)-dip_nuc(1)!*norm(i-1)
     dmut(2,dd)=dmut(2,dd)-dip_nuc(2)!*norm(i-1)
     dmut(3,dd)=dmut(3,dd)-dip_nuc(3)!*norm(i-1)
  enddo

  deallocate(norm)

  !do ints=1,nts
  !   potmut(1,1,ints)=potmut(1,1,ints)-pot_nuc0(ints)
  !   do i=2,ntoten+1
  !      potmut(i,i,ints)=potmut(i,i,ints)-pot_nuc0(ints)!*norm(i-1)
  !   enddo 
  !enddo

  do i=1,ntoten+1
     dd=cont(1,i)
     write(15,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  0,  i-1, dmut(1,dd),dmut(2,dd),dmut(3,dd)
  enddo

  do i=2,ntoten+1
     do j=2,i
        dd=cont(i,j)
        write(15,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  i-1, j-1, dmut(1,dd),dmut(2,dd),dmut(3,dd)
     enddo
  enddo
  close(15)
  deallocate(dmut)

  if (cdspectrum) then
     open(16,file='ci_lt.inp')

     do i=1,ntoten+1
        dd=cont(1,i)
        write(16,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  0,  i-1,lmut(1,dd),lmut(2,dd),lmut(3,dd)
     enddo

     do i=2,ntoten+1
        do j=i,ntoten+1
           dd=cont(i,j)
           write(16,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  i-1, j-1,lmut(1,dd),lmut(2,dd),lmut(3,dd)
        enddo
     enddo
     close(16)

     deallocate(lmut)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NANOPARTICLE by Pier e Leo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (NP) then 
     open(17, file = 'ci_pot.inp')
     write(17,*) nts
     i=0
     j=0
     write(17,*) i,j
     do ints=1,nts
        write(17,*) potmut(1,ints), 0, pot_nuc0(ints)  
     enddo

     do i=2,ntoten+1
        dd=cont(1,i)
        write(17,*) 0, i-1
        do ints=1,nts
           write(17,*) potmut(dd,ints)
        enddo
     enddo

      do i=2,ntoten+1
           do j=2,i
              dd=cont(i,j)
              write(17,*)  i-1, j-1
              do ints=1,nts
                 write(17,*) potmut(dd,ints)
              enddo
           enddo
       enddo
       close(17)
       deallocate(potmut)
       deallocate(pot_nuc0)
  endif

  deallocate(cont)


 write(*,*) "Ho finito di scrivere"


  WRITE(*,98) ' ********************************* '
  WRITE(*,96) ' * REGULAR TERMIMATION OF TRANSM * '
  WRITE(*,95) ' ********************************* '


  call CPU_TIME(finish)
  call system_clock(current)

  write(*,*) 'Elapsed time: ', real(current-st)/real(rate)
  write(*,*) 'CPU Elapsed time: ', finish-start
  write(*,*) 'OpenMP time: ', finish_omp-start_omp
98 FORMAT(//// , 30X , A     )
96 FORMAT(     30X , A     )
95 FORMAT(     30X , A ,//// )

  write(*,*) 'TERMINATION'

  !call MPI_FINALIZE(ierr)
  stop

END PROGRAM write_wavet 
