PROGRAM write_wavet 

  use omp_lib

  implicit none


  real*8,  allocatable :: dipmatx(:,:) , dipmaty(:,:) , dipmatz(:,:)
  real*8,  allocatable :: lmatx(:,:) , lmaty(:,:) , lmatz(:,:)
  real*8,  allocatable :: exciten(:), tddfteig(:,:,:), tddfteigl(:,:,:)
  real*8,  allocatable :: dmut(:,:,:),norm(:), norml(:), lmut(:,:,:)
  real*8               :: dip0(3),dip_nuc(3),l0(3)

  real*8,  PARAMETER   ::  EVAU = 27.21139628D0 

  integer   :: ia,ib,i,j,a,itoten,iener,isym,is,js,jj
  integer   :: kvirt,kocc,ntoten,ntotmo,nener,nsym,nthreads

!MPI_VAR
!  integer   :: ierr, nprocs, myid

  real      :: start, finish
  real      :: start_omp, finish_omp
  integer   :: rate, st, current

  character*20  :: cdum

  logical :: cdspectrum

!  call MPI_INIT(ierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
!  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

!#ifdef OMP
!   nthreads=omp_get_max_threads( )
!#endif
!#ifndef OMP
!   nthreads=1
!#endif

!#ifdef OMP
!   write(*,*) '*************OMP*****************'
!   write(*,'(a)')
!   write(*,*) 'OMP parallelization on'
!   write(*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs ( )
!   write(*,'(a,i8)') 'The number of threads available    = ', nthreads 
!   write(*,'(a)')
!   write(*,*) '*************OMP*****************'
!#endif

  call system_clock(st,rate)
  call CPU_TIME(start)
  open(80,file='eig.dat')
  read(80,*) nsym,kocc,kvirt,ntoten,ntotmo,cdspectrum 

  if (cdspectrum) open(81,file='eig_l.dat')
  allocate(dipmatx(ntotmo,ntotmo))
  allocate(dipmaty(ntotmo,ntotmo))
  allocate(dipmatz(ntotmo,ntotmo))
  if (cdspectrum) then
     allocate(lmatx(ntotmo,ntotmo))
     allocate(lmaty(ntotmo,ntotmo))
     allocate(lmatz(ntotmo,ntotmo))
     lmatx=0.d0
     lmaty=0.d0
     lmatz=0.d0
  allocate(tddfteigl(kvirt,kocc,ntoten))
  allocate(norml(ntoten))
  allocate(lmut(3,ntoten+1,ntoten+1))
  endif
  dipmatx=0.d0
  dipmaty=0.d0
  dipmatz=0.d0

  allocate(tddfteig(kvirt,kocc,ntoten))
  allocate(exciten(ntoten))
  allocate(norm(ntoten))
  allocate(dmut(3,ntoten+1,ntoten+1))

  itoten=0
  do isym=1,nsym
     read(80,*) nener
     do iener=1,nener
        itoten=itoten+1
        do i=1,kocc
           do j=1,kvirt
              read(80,*) tddfteig(j,i,itoten)
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
        !$OMP PARALLEL DO REDUCTION(+:norm) 
        do i=1,kocc
           do j=1,kvirt
               norm(itoten)=norm(itoten)+tddfteig(j,i,itoten)**2
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  enddo 
  if (cdspectrum) then
     itoten=0
     do isym=1,nsym
        do iener=1,nener
           itoten=itoten+1
           norml(itoten)=0.d0
           !$OMP PARALLEL DO REDUCTION(+:norml)
           do i=1,kocc
              do j=1,kvirt
                 norml(itoten)=norml(itoten)+tddfteigl(j,i,itoten)**2
              enddo
           enddo
           !$OMP END PARALLEL DO
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
  !$OMP  PARALLEL DO REDUCTION(+:dip0)
  do j=1,kocc
     dip0(1) = dip0(1) + dipmatx(j,j) 
     dip0(2) = dip0(2) + dipmaty(j,j)
     dip0(3) = dip0(3) + dipmatz(j,j)
  enddo
  !$OMP  END PARALLEL DO
  dmut=0.d0
  dmut(:,1,1)=2.d0*dip0(:)

  if (cdspectrum) then
     l0=0.d0
     !$OMP  PARALLEL DO REDUCTION(+:l0)
     do j=1,kocc
        l0(1) = l0(1) + lmatx(j,j)
        l0(2) = l0(2) + lmaty(j,j)
        l0(3) = l0(3) + lmatz(j,j)
     enddo
     !$OMP  END PARALLEL DO
     lmut(:,1,1)=2.d0*l0(:)
  endif

  !GS-EXC
  !$OMP  PARALLEL DO REDUCTION(+:dmut)
  do is=1,ntoten 
     do i=1,kocc
        do ia=1,kvirt
           dmut(1,1,is+1) = dmut(1,1,is+1) + tddfteig(ia,i,is)*dipmatx(i,kocc+ia)
           dmut(2,1,is+1) = dmut(2,1,is+1) + tddfteig(ia,i,is)*dipmaty(i,kocc+ia) 
           dmut(3,1,is+1) = dmut(3,1,is+1) + tddfteig(ia,i,is)*dipmatz(i,kocc+ia)
        enddo
     enddo
     dmut(:,is+1,1)=dmut(:,1,is+1)
  enddo
  !$OMP  END PARALLEL DO
  if (cdspectrum) then
     !$OMP  PARALLEL DO REDUCTION(+:lmut)
     do is=1,ntoten 
        do i=1,kocc
           do ia=1,kvirt
              lmut(1,1,is+1) = lmut(1,1,is+1) + tddfteigl(ia,i,is)*lmatx(i,kocc+ia)
              lmut(2,1,is+1) = lmut(2,1,is+1) + tddfteigl(ia,i,is)*lmaty(i,kocc+ia)
              lmut(3,1,is+1) = lmut(3,1,is+1) + tddfteigl(ia,i,is)*lmatz(i,kocc+ia)
           enddo
        enddo
        lmut(:,is+1,1)=-lmut(:,1,is+1) !QUESTIONE MENO?
     enddo
     !$OMP  END PARALLEL DO
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
 !$OMP PARALLEL DO REDUCTION(+:dmut)
  do is=1,ntoten
     do js=is,ntoten
        do i=1,kocc
           do ia=1,kvirt   
              dmut(1,is+1,js+1) = dmut(1,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(1) + dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,is+1,js+1) = dmut(2,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(2) + dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,is+1,js+1) = dmut(3,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(3) + dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
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
                 dmut(1,is+1,js+1) = dmut(1,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ia,kocc+ib)
                 dmut(2,is+1,js+1) = dmut(2,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ia,kocc+ib)
                 dmut(3,is+1,js+1) = dmut(3,is+1,js+1) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ia,kocc+ib)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,is+1,js+1) = dmut(1,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,is+1,js+1) = dmut(2,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,is+1,js+1) = dmut(3,is+1,js+1) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
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
        dmut(:,js+1,is+1) = dmut(:,is+1,js+1)
     enddo
  enddo
 !$OMP END PARALLEL DO 
!EXC-EXC

 if (cdspectrum) then
  !$OMP PARALLEL DO REDUCTION(+:lmut)
    do is=1,ntoten
       do js=is,ntoten
          do i=1,kocc
             do ia=1,kvirt   
                lmut(1,is+1,js+1) = lmut(1,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(1) + lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                lmut(2,is+1,js+1) = lmut(2,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(2) + lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                lmut(3,is+1,js+1) = lmut(3,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ia,i,js)*(2.d0*l0(3) + lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
            
                do ib=1,kvirt
                   if (ib.eq.ia) cycle
                   lmut(1,is+1,js+1) = lmut(1,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmatx(kocc+ia,kocc+ib)
                   lmut(2,is+1,js+1) = lmut(2,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmaty(kocc+ia,kocc+ib)
                   lmut(3,is+1,js+1) = lmut(3,is+1,js+1) + tddfteigl(ia,i,is)*tddfteigl(ib,i,js)*lmatz(kocc+ia,kocc+ib)
                enddo
                do j=1,kocc
                   if (j.eq.i) cycle 
                   lmut(1,is+1,js+1) = lmut(1,is+1,js+1) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmatx(j,i)
                   lmut(2,is+1,js+1) = lmut(2,is+1,js+1) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmaty(j,i)
                   lmut(3,is+1,js+1) = lmut(3,is+1,js+1) - tddfteigl(ia,i,is)*tddfteigl(ia,j,js)*lmatz(j,i)
                enddo
             enddo
          enddo
          lmut(:,js+1,is+1) = -lmut(:,is+1,js+1)
       enddo
    enddo
  !$OMP END PARALLEL DO
  finish_omp = omp_get_Wtime()
  endif

  deallocate(tddfteig)
  deallocate(dipmatx)
  deallocate(dipmaty)
  deallocate(dipmatz)
  if (cdspectrum) then
     deallocate(lmatx)
     deallocate(lmaty)
     deallocate(lmatz)
     deallocate(tddfteigl)
  endif


  exciten=exciten*EVAU

  open(24,file='ci_energy.inp')
  do i=1,ntoten
     write(24,'("Root", I5, " : ", F12.5)') i, exciten(i)
  enddo
  close(24)

  deallocate(exciten)

  open(15,file='ci_mut.inp')

  dmut(:,1,1)=dmut(:,1,1)-dip_nuc(:)
  do i=2,ntoten+1
     dmut(1,i,i)=dmut(1,i,i)-dip_nuc(1)!*norm(i-1)
     dmut(2,i,i)=dmut(2,i,i)-dip_nuc(2)!*norm(i-1)
     dmut(3,i,i)=dmut(3,i,i)-dip_nuc(3)!*norm(i-1)
  enddo

  deallocate(norm)

  do i=1,ntoten+1
     write(15,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  0,  i-1, dmut(1,1,i),dmut(2,1,i),dmut(3,1,i)
  enddo

  do i=2,ntoten+1
     do j=2,i
        write(15,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  i-1, j-1, dmut(1,i,j),dmut(2,i,j),dmut(3,i,j)
     enddo
  enddo
  close(15)
  deallocate(dmut)

  if (cdspectrum) then
     open(16,file='ci_lt.inp')

     do i=1,ntoten+1
        write(16,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  0,  i-1,lmut(1,1,i),lmut(2,1,i),lmut(3,1,i)
     enddo

     do i=2,ntoten+1
        do j=2,i
           write(16,'("States", I5, " and", I5,F17.8,F17.8,F17.8)')  i-1, j-1,lmut(1,i,j),lmut(2,i,j),lmut(3,i,j)
        enddo
     enddo
     close(16)

     deallocate(lmut)
  endif

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

!  call MPI_FINALIZE(ierr)
  stop

END PROGRAM write_wavet 
