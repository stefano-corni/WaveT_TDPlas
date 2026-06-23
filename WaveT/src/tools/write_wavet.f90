PROGRAM write_wavet 


  use omp_lib
!#ifdef MPI
!  use mpi
!#endif

  implicit none
  
  include 'mpif.h'


  integer              :: ierr, rank, size, nts_loc, frac2, nts_loc_0
#ifdef MPI
  integer, dimension(MPI_STATUS_SIZE) :: status 
#endif   
  real*8, allocatable  :: local_potmut(:,:), local_pot0(:),local_potmut1(:,:)


  real*8,  allocatable :: dipmatx(:,:) , dipmaty(:,:) , dipmatz(:,:), sqrtepsdiff(:,:)
  real*8,  allocatable :: lmatx(:,:) , lmaty(:,:) , lmatz(:,:)
  real*8,  allocatable :: potmat(:,:,:), potmat_nuc(:,:,:), potmut1(:,:),potmut(:,:), pot0(:), pot_nuc0(:), potmut_nuc(:)
  real*8,  allocatable :: exciten(:), tddfteig(:,:,:),tddfteigl(:,:,:),tddfteig1(:,:,:)
  real*8,  allocatable :: dmut(:,:),norm(:), norml(:), lmut(:,:),dmut1(:,:),lmut1(:,:)
  real*8               :: dip0(3),dip_nuc(3),l0(3)

  real*8,  PARAMETER   ::  EVAU = 27.21139628D0 


  integer,  allocatable   :: cont(:,:)
  integer   :: ia,ib,i,j,a,itoten,iener,isym,is,js,jj,ints,kk,dd,ll,request,dummy
  integer   :: kvirt,kocc,ntoten,ntotmo,nener,nsym,nthreads,nts,ntotentr,threads


  real      :: start, finish
  real      :: start_omp, finish_omp
  integer   :: rate, st, current

  character*20  :: cdum
  
  logical :: cdspectrum,tda,hybrid,dipo
  
  logical :: NP

  call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)


#ifndef MPI
   rank=0
#endif 


   nthreads=omp_get_max_threads( )
if (rank==0) then
   write(*,*) '*************OMP*****************'
   write(*,'(a)')
   write(*,*) 'OMP parallelization on'
   write(*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs ( )
   write(*,'(a,i8)') 'The number of threads available    = ', nthreads 
   write(*,'(a)')
   write(*,*) '*************OMP*****************'
endif
  call system_clock(st,rate)
  call CPU_TIME(start)
  open(80,file='eig.dat', form='unformatted')
  read(80) nsym,kocc,kvirt,ntoten,ntotmo,cdspectrum,NP
  write(*,*) 'sono arrivato dopo eig.dat'
  allocate(cont(ntoten+1,ntoten+1))
  if (NP) then
  read(80) nts
  endif
  
  ntotentr = ((ntoten+1)*(ntoten+2))/2
  
  if (NP) then
  nts_loc=nts/size
  frac2=mod(nts,size)
  nts_loc_0 = nts_loc + frac2
  write(*,*) 'nts, nts_loc e nts_loc_0', nts, nts_loc, nts_loc_0
  endif

  open(unit=14,file='info_qm.dat')
  read(14,*) tda,hybrid,dipo
  close(14)

  if (cdspectrum) open(81,file='eig_l.dat')
  if (dipo) then
     allocate(dipmatx(ntotmo,ntotmo))
     allocate(dipmaty(ntotmo,ntotmo))
     allocate(dipmatz(ntotmo,ntotmo))
     dipmatx=0.d0
     dipmaty=0.d0
     dipmatz=0.d0
  endif   
  if (cdspectrum) then
     allocate(lmatx(ntotmo,ntotmo))
     allocate(lmaty(ntotmo,ntotmo))
     allocate(lmatz(ntotmo,ntotmo))
     lmatx=0.d0
     lmaty=0.d0
     lmatz=0.d0
  allocate(lmut(3,ntotentr))
  allocate(lmut1(3,ntotentr))
  allocate(norml(ntoten))
  endif
  write(*,*) 'numero orbitali virtuali', kvirt, 'num orb occ', kocc, 'num tot stati',ntoten
  

  
  cont=0.d0
  kk=0.d0
  do is=1,ntoten+1
     do js=is,ntoten+1
        kk=kk+1
        cont(is,js)=kk
        cont(js,is)=cont(is,js)
     enddo
 enddo


  if (NP) then
   if (rank==0) then 
  write(*,*) "ciaoo" 
     allocate(potmut(ntotentr,nts))
  write(*,*) "ciaoo1" 
     allocate(potmut1(ntotentr,nts))
   endif
  write(*,*) "ciaoo2" 
     allocate(pot0(nts))
  write(*,*) "ciaoo3" 
     allocate(potmat(ntotmo,ntotmo,nts))
  write(*,*) "ciao04" 
     allocate(potmut_nuc(nts))
  write(*,*) "ciaoo5" 
     allocate(pot_nuc0(nts))
     potmat=0.d0
     pot0=0.d0
!     potmat_nuc=0.d0
  write(*,*) "ciaoo6" 
     allocate (local_potmut(ntotentr,nts_loc))
  write(*,*) "ciaoo7" 
     allocate (local_potmut1(ntotentr,nts_loc))
  write(*,*) "ciaoo8" 
  endif 

  allocate(tddfteig(kvirt,kocc,ntoten))
  allocate(tddfteig1(kvirt,kocc,ntoten))
  allocate(exciten(ntoten))
  allocate(norm(ntoten))
  if (dipo) then
    allocate(dmut(3,ntotentr))
    allocate(dmut1(3,ntotentr))
  endif
  allocate(sqrtepsdiff(kvirt,kocc))
 
  open(13, file='std_epsilons.dat')
  read(13,*) sqrtepsdiff
  close(13)
 
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
  
  write(*,*) 'sono arrivato dopo tddfteig'


  open(70,file='ene.dat')
  do i=1,ntoten
     read(70,*) exciten(i)
  enddo


  if (hybrid) then
  itoten=0
  open(82, file='eig_hybrid.dat',form='unformatted')
  do isym=1,nsym
     read(82) nener
     do iener=1,nener
        itoten=itoten+1
        do i=1,kocc
           do j=1,kvirt
              read(82) tddfteig1(j,i,itoten)
           enddo
        enddo
     enddo
  enddo
  close(82)
  elseif(.not.tda.and..not.hybrid) then 

   itoten=0
   do isym=1,nsym
      do iener=1,nener
         itoten=itoten+1
         do i=1,kocc
           do j=1,kvirt
             tddfteig1(j,i,itoten)=tddfteig(j,i,itoten)*sqrt(exciten(itoten))/sqrtepsdiff(j,i)
             tddfteig1(j,i,itoten)=tddfteig1(j,i,itoten)*sqrt(exciten(itoten))/sqrtepsdiff(j,i)
           enddo
         enddo
     enddo
   enddo
   endif


  if (dipo) then 
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
  endif

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

     open(93, file='potmat_nuc.dat')
     do ints=1,nts
        read(93,*) potmut_nuc(ints)
     enddo
     close(93)
  endif

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

  if (dipo) then
  dip0=0.d0
!$OMP PARALLEL DO REDUCTION(+:dip0)
  do j=1,kocc
     dip0(1) = dip0(1) + dipmatx(j,j) 
     dip0(2) = dip0(2) + dipmaty(j,j)
     dip0(3) = dip0(3) + dipmatz(j,j)
  enddo
!$OMP END PARALLEL DO

   dmut=0.d0
   dmut(:,1)=2.d0*dip0(:)
   dmut(:,1)=2.d0*dip0(:)-dip_nuc(:)
  endif
write(*,*) 'rank1', rank


if (NP) then 
if (rank.gt.0) then 
 allocate(local_pot0(nts_loc))
  local_pot0=0.d0
!$OMP PARALLEL DO  
  do ints=1,nts_loc
     do j=1,kocc
         local_pot0(ints) = local_pot0(ints) + potmat(j,j,ints + nts_loc_0+nts_loc*(rank-1))
     enddo
  enddo
!$OMP END PARALLEL DO
  call MPI_Send(local_pot0, nts_loc, MPI_DOUBLE_PRECISION, 0, 120+rank, MPI_COMM_WORLD,ierr) 
endif
if (rank==0) then
 pot0=0.d0
 allocate(local_pot0(nts_loc))
 local_pot0=0.d0
!$OMP PARALLEL DO  
 do ints=1,nts_loc_0
     do j=1,kocc
         pot0(ints) = pot0(ints) + potmat(j,j,ints)
     enddo
  enddo
!$OMP END PARALLEL DO
 do i=1,size-1     
  call MPI_Recv(local_pot0, nts_loc, MPI_DOUBLE_PRECISION, i, 120+i, MPI_COMM_WORLD,status,ierr)
     write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)
     do ints=1,nts_loc
         ll = ints+nts_loc_0+nts_loc*(i-1)
         pot0(ll) = local_pot0(ints) 
     enddo
 enddo 
endif
if (tda) then
  call MPI_Bcast(pot0,nts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
endif
!pot_nuc0=0.d0
deallocate (local_pot0)
endif

if (NP) then
 if (rank==0) then
    potmut = 0.d0
!    potmut_nuc = 0.d0
    potmut(1,:)=2.d0*pot0(:)
  endif
!    potmut_nuc(:)=potmat_nuc(1,1,:)
endif

 
  if (cdspectrum) then
     l0=0.d0
!$OMP PARALLEL DO REDUCTION(+:l0)
     do j=1,kocc
        l0(1) = l0(1) + lmatx(j,j)
        l0(2) = l0(2) + lmaty(j,j)
        l0(3) = l0(3) + lmatz(j,j)
     enddo
!$OMP END PARALLEL DO
      lmut=0.d0
      lmut(:,1)=2.d0*l0(:)
  endif



  !GS-EXC
 if (dipo) then
   if (rank==0.or.rank==1) then
!$OMP PARALLEL DO  
  do is=1,ntoten
      kk = cont(1,is+1)
      do i=1,kocc
        do ia=1,kvirt
           dmut(1,kk)  = dmut(1,kk) + tddfteig(ia,i,is)*dipmatx(i,kocc+ia)
           dmut(2,kk)  = dmut(2,kk) + tddfteig(ia,i,is)*dipmaty(i,kocc+ia) 
           dmut(3,kk)  = dmut(3,kk) + tddfteig(ia,i,is)*dipmatz(i,kocc+ia)
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
 write(*,*) "Sono dopo dmut"

  dmut(:,:) = sqrt(2.d0)*dmut(:,:)


  dmut(:,1) = dmut(:,1)/sqrt(2.d0)

  endif
endif


if (cdspectrum) then

if (rank==0.or.rank==1) then
  if (tda) then
!$OMP PARALLEL DO 
    do is=1,ntoten
      kk = cont(1,is+1)
        do i=1,kocc
           do ia=1,kvirt
              lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*lmatx(i,kocc+ia)
              lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*lmaty(i,kocc+ia)
              lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*lmatz(i,kocc+ia)
           enddo
        enddo
     enddo
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO 
    do is=1,ntoten
      kk = cont(1,is+1)
        do i=1,kocc
           do ia=1,kvirt
              lmut(1,kk) = lmut(1,kk) + tddfteig1(ia,i,is)*lmatx(i,kocc+ia)
              lmut(2,kk) = lmut(2,kk) + tddfteig1(ia,i,is)*lmaty(i,kocc+ia)
              lmut(3,kk) = lmut(3,kk) + tddfteig1(ia,i,is)*lmatz(i,kocc+ia)
           enddo
        enddo
     enddo
!$OMP END PARALLEL DO
endif

lmut(:,:) = sqrt(2.d0)*lmut(:,:)

lmut(:,1) = lmut(:,1)/sqrt(2.d0)
endif
endif

write(*,*) 'rank2', rank

if (NP) then
if (rank.gt.0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
!$OMP PARALLEL DO
  do ints=1,nts_loc
    do is=1,ntoten
        kk = cont(1,is+1)
        do i=1,kocc
          do ia=1,kvirt
           local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*&
                       potmat(i,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))
          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO
  write(*,*) 'nei rank',rank, 'nts_loc', nts_loc
  call MPI_Send(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 120+rank, MPI_COMM_WORLD,ierr) 
else
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
!$OMP PARALLEL DO
 do ints=1,nts_loc_0
     do is=1,ntoten
       kk = cont(1,is+1)
       do i=1,kocc
         do ia=1,kvirt
            potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*potmat(i,kocc+ia,ints)
         enddo
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO
 do i=1,size-1     
   write(*,*) 'ciao tocca a ', i
   call MPI_Recv(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 120+i, MPI_COMM_WORLD,status,ierr)
     write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)
     do ints=1,nts_loc
         ll = ints + nts_loc_0+nts_loc*(i-1) 
         do is = 1,ntoten
            kk=cont(1,is+1)
            potmut(kk,ll) = local_potmut(kk,ints) 
         enddo
     enddo
 enddo 


potmut(:,:) = sqrt(2.d0)*potmut(:,:)
potmut(1,:) = potmut(1,:)/sqrt(2.d0)

endif

endif

write(*,*) 'ho fatto potmut1', rank


if (dipo) then
 if (hybrid) then 
  if (size.gt.1) then 
   if (rank==1) then

!EXC-EXC
!$OMP PARALLEL DO 
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
!EXC-EXC
call MPI_Send(dmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 0, 191, MPI_COMM_WORLD,ierr) 
 
elseif (rank==0) then

 dmut1 = dmut

!$OMP PARALLEL DO
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 dmut1(1,kk) = dmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatx(j,i)
                 dmut1(2,kk) = dmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmaty(j,i)
                 dmut1(3,kk) = dmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
call MPI_Recv(dmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 1, 191, MPI_COMM_WORLD,status,ierr)
write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)

dmut=0.5*(dmut1 + dmut)

endif

else !entro qui se la size è 1 

 dmut1 = dmut

!$OMP PARALLEL DO 
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 dmut1(1,kk) = dmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatx(j,i)
                 dmut1(2,kk) = dmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmaty(j,i)
                 dmut1(3,kk) = dmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO

dmut=0.5*(dmut1 + dmut)

endif !finisce l'if per la size 

elseif(.not.tda.and..not.hybrid) then 

!   itoten=0
!  do isym=1,nsym
!     do iener=1,nener
!        itoten=itoten+1
!        do i=1,kocc
!           do j=1,kvirt
!             tddfteig1(j,i,itoten)=tddfteig(j,i,itoten)*sqrt(exciten(itoten))/sqrtepsdiff(j,i)
!             tddfteig1(j,i,itoten)=tddfteig1(j,i,itoten)*sqrt(exciten(itoten))/sqrtepsdiff(j,i)
!           enddo
!         enddo
!     enddo
!  enddo

if (size.gt.1) then 
 if (rank==1) then

!EXC-EXC
!$OMP PARALLEL DO 
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
!EXC-EXC


call MPI_Send(dmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 0, 191, MPI_COMM_WORLD,ierr) 

 elseif (rank==0) then

 dmut1 = dmut
!$OMP PARALLEL DO
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 dmut1(1,kk) = dmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatx(j,i)
                 dmut1(2,kk) = dmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmaty(j,i)
                 dmut1(3,kk) = dmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO

call MPI_Recv(dmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 1, 191, MPI_COMM_WORLD,status,ierr)
write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)

dmut=0.5*(dmut1 + dmut)

endif

else !entro qui se la size è 1

 dmut1 = dmut

 
!$OMP PARALLEL DO 
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
 
 
!$OMP PARALLEL DO
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*( dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut1(1,kk) = dmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut1(2,kk) = dmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut1(3,kk) = dmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 dmut1(1,kk) = dmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatx(j,i)
                 dmut1(2,kk) = dmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmaty(j,i)
                 dmut1(3,kk) = dmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
dmut=0.5*(dmut1 + dmut)

endif !finisce l'if della size

else 

if (rank==0) then
!EXC-EXC
!$OMP PARALLEL DO 
  do is=1,ntoten
     do js=is,ntoten
        kk = cont(is+1,js+1)
        do i=1,kocc
           do ia=1,kvirt
              dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(1) + dipmatx(kocc+ia,kocc+ia)-dipmatx(i,i))
              dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(2) + dipmaty(kocc+ia,kocc+ia)-dipmaty(i,i))
              dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*dip0(3) + dipmatz(kocc+ia,kocc+ia)-dipmatz(i,i))
              do ib=1,kvirt
                 if (ib.eq.ia) cycle
                 dmut(1,kk) = dmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatx(kocc+ib,kocc+ia)
                 dmut(2,kk) = dmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmaty(kocc+ib,kocc+ia)
                 dmut(3,kk) = dmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*dipmatz(kocc+ib,kocc+ia)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle 
                 dmut(1,kk) = dmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatx(j,i)
                 dmut(2,kk) = dmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmaty(j,i)
                 dmut(3,kk) = dmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*dipmatz(j,i)
              enddo

           enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
!EXC-EXC
  endif
 endif
endif 
 call MPI_Barrier(MPI_COMM_WORLD, ierr)
 
 if (NP) then

if(tda) then

 if (rank.gt.0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
!$OMP PARALLEL DO
  do ints=1,nts_loc
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                  (2.d0*pot0(ints + nts_loc_0+nts_loc*(rank-1)) + potmat(kocc+ia,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))-&
                    potmat(i,i,ints + nts_loc_0+nts_loc*(rank-1)))
                   do ib=1,kvirt
                       if (ib.eq.ia) cycle
                         local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                                 potmat(kocc+ia,kocc+ib,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                    do j=1,kocc
                       if (j.eq.i) cycle
                       local_potmut(kk,ints) = local_potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*&
                               potmat(j,i,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                enddo
             enddo
          enddo
      enddo
   enddo
!$OMP END PARALLEL DO
 call MPI_Send(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 120+rank, MPI_COMM_WORLD,ierr) 
 endif

if (rank==0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
!$OMP PARALLEL DO
  do ints=1,nts_loc_0
     do is=1,ntoten
        do js=is,ntoten
          kk = cont(is+1,js+1)
          do i=1,kocc
            do ia=1,kvirt
                 potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                         (2.d0*pot0(ints) + potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
             do ib=1,kvirt
                 if (ib.eq.ia) cycle
                  potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                          potmat(kocc+ia,kocc+ib,ints)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 potmut(kk,ints) = potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*potmat(j,i,ints)
              enddo
            enddo
          enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
 
 do i=1,size-1     
     call MPI_Recv(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 120+i, MPI_COMM_WORLD,status,ierr)
     write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)
   do ints=1,nts_loc
     ll = ints + nts_loc_0+nts_loc*(i-1) 
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            potmut(kk,ll) = local_potmut(kk,ints) 
        enddo
     enddo
   enddo
  enddo 
 endif

elseif(hybrid) then

 if (rank.gt.0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
  local_potmut1=local_potmut
!$OMP PARALLEL DO
  do ints=1,nts_loc
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                  (potmat(kocc+ia,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))-&
                  potmat(i,i,ints + nts_loc_0+nts_loc*(rank-1)))
                  local_potmut1(kk,ints) = local_potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*&
                  (potmat(kocc+ia,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))-&
                  potmat(i,i,ints + nts_loc_0+nts_loc*(rank-1)))
                   do ib=1,kvirt
                       if (ib.eq.ia) cycle
                         local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                                 potmat(kocc+ia,kocc+ib,ints + nts_loc_0+nts_loc*(rank-1))
                         local_potmut1(kk,ints) = local_potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*&
                                 potmat(kocc+ia,kocc+ib,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                    do j=1,kocc
                       if (j.eq.i) cycle
                       local_potmut(kk,ints) = local_potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*&
                               potmat(j,i,ints + nts_loc_0+nts_loc*(rank-1))
                       local_potmut1(kk,ints) = local_potmut1(kk,ints) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*&
                               potmat(j,i,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                enddo
             enddo
          enddo
      enddo
   enddo
!$OMP END PARALLEL DO
 call MPI_Send(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 120+rank, MPI_COMM_WORLD,ierr) 
 call MPI_Send(local_potmut1, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 0+rank, MPI_COMM_WORLD,ierr) 
 endif

if (rank==0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
  potmut1=potmut
!$OMP PARALLEL DO
  do ints=1,nts_loc_0
     do is=1,ntoten
        do js=is,ntoten
          kk = cont(is+1,js+1)
          do i=1,kocc
            do ia=1,kvirt
                 potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                         (potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
                 potmut1(kk,ints) = potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*&
                         (potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
             do ib=1,kvirt
                 if (ib.eq.ia) cycle
                  potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                                    potmat(kocc+ia,kocc+ib,ints)
                  potmut1(kk,ints) = potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*&
                                    potmat(kocc+ia,kocc+ib,ints)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 potmut(kk,ints) = potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*potmat(j,i,ints)
                 potmut1(kk,ints) = potmut1(kk,ints) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*potmat(j,i,ints)
              enddo
            enddo
          enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
 
 do i=1,size-1     
     call MPI_Recv(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 120+i, MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(local_potmut1, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 0+i, MPI_COMM_WORLD,status,ierr)
     write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)
   do ints=1,nts_loc
     ll = ints + nts_loc_0+nts_loc*(i-1) 
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            potmut(kk,ll) = local_potmut(kk,ints) 
            potmut1(kk,ll) = local_potmut1(kk,ints) 
        enddo
     enddo
   enddo
  enddo
  

potmut=0.5*(potmut+potmut1)


endif




elseif (.not.tda.and..not.hybrid) then


 if (rank.gt.0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
  local_potmut1=local_potmut
!$OMP PARALLEL DO
  do ints=1,nts_loc
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                      (potmat(kocc+ia,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))-&
                      potmat(i,i,ints + nts_loc_0+nts_loc*(rank-1)))
                  local_potmut1(kk,ints) = local_potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*&
                          (potmat(kocc+ia,kocc+ia,ints + nts_loc_0+nts_loc*(rank-1))-&
                          potmat(i,i,ints + nts_loc_0+nts_loc*(rank-1)))
                   do ib=1,kvirt
                       if (ib.eq.ia) cycle
                         local_potmut(kk,ints) = local_potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                                    potmat(kocc+ia,kocc+ib,ints + nts_loc_0+nts_loc*(rank-1))
                         local_potmut1(kk,ints) = local_potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*&
                                    potmat(kocc+ia,kocc+ib,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                    do j=1,kocc
                       if (j.eq.i) cycle
                       local_potmut(kk,ints) = local_potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*&
                                     potmat(j,i,ints + nts_loc_0+nts_loc*(rank-1))
                       local_potmut1(kk,ints) = local_potmut1(kk,ints) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*&
                                     potmat(j,i,ints + nts_loc_0+nts_loc*(rank-1))
                    enddo
                enddo
             enddo
          enddo
      enddo
   enddo
!$OMP END PARALLEL DO

 call MPI_Send(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 120+rank, MPI_COMM_WORLD,ierr) 
 call MPI_Send(local_potmut1, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, 0, 0+rank, MPI_COMM_WORLD,ierr) 
 endif
 if (rank==0) then
  nts_loc_0 = nts_loc + frac2
  local_potmut=0.d0
  potmut1=potmut
!$OMP PARALLEL DO
  do ints=1,nts_loc_0
     do is=1,ntoten
        do js=is,ntoten
          kk = cont(is+1,js+1)
          do i=1,kocc
            do ia=1,kvirt
                 potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*&
                                  (potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
                 potmut1(kk,ints) = potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ia,i,js)*&
                                  (potmat(kocc+ia,kocc+ia,ints)-potmat(i,i,ints))
             do ib=1,kvirt
                 if (ib.eq.ia) cycle
                  potmut(kk,ints) = potmut(kk,ints) + tddfteig(ia,i,is)*tddfteig(ib,i,js)*&
                                   potmat(kocc+ia,kocc+ib,ints)
                  potmut1(kk,ints) = potmut1(kk,ints) + tddfteig1(ia,i,is)*tddfteig1(ib,i,js)*&
                                    potmat(kocc+ia,kocc+ib,ints)
              enddo
              do j=1,kocc
                 if (j.eq.i) cycle
                 potmut(kk,ints) = potmut(kk,ints) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*potmat(j,i,ints)
                 potmut1(kk,ints) = potmut1(kk,ints) - tddfteig1(ia,i,is)*tddfteig1(ia,j,js)*potmat(j,i,ints)
              enddo
            enddo
          enddo
        enddo
     enddo
  enddo
!$OMP END PARALLEL DO
 
 do i=1,size-1     
     call MPI_Recv(local_potmut, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 120+i, MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(local_potmut1, nts_loc*ntotentr, MPI_DOUBLE_PRECISION, i, 0+i, MPI_COMM_WORLD,status,ierr)
     write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)
   do ints=1,nts_loc
     ll = ints + nts_loc_0+nts_loc*(i-1) 
     do is=1,ntoten
        do js=is,ntoten
            kk = cont(is+1,js+1)
            potmut(kk,ll) = local_potmut(kk,ints) 
            potmut1(kk,ll) = local_potmut1(kk,ints) 
        enddo
     enddo
   enddo
  enddo 
  potmut=0.5*(potmut+potmut1)
 endif


endif

endif

if (NP) then
  deallocate(local_potmut)
  deallocate(local_potmut1)
endif

  write(*,*) 'ho fatto potmut2'

if (cdspectrum) then
  
if (tda) then
  if (rank==0) then
!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
             do i=1,kocc
               do ia=1,kvirt
                 lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(1) + lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                 lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(2) + lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                 lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(3) + lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                 do ib=1,kvirt
                    if (ib.eq.ia) cycle
                    lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                    lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                    lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                 enddo

                 do j=1,kocc
                    if (j.eq.i) cycle
                    lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*lmatx(j,i)
                    lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*lmaty(j,i)
                    lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig(ia,j,js)*lmatz(j,i)
                 enddo
              enddo
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
  finish_omp = omp_get_Wtime()
 endif

elseif(hybrid) then

if (size.gt.1) then
  if (rank==1) then
!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatx(j,i)
                     lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmaty(j,i)
                     lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

call MPI_Send(lmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 0, 101, MPI_COMM_WORLD,ierr) 
elseif (rank==0) then

lmut1=lmut

!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut1(1,kk) = lmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut1(2,kk) = lmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut1(3,kk) = lmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatx(j,i)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmaty(j,i)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO
call MPI_Recv(lmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 1, 101, MPI_COMM_WORLD,status,ierr)
write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)

lmut = 0.5d0*(lmut+lmut1)
finish_omp = omp_get_Wtime()
endif        
        
else

lmut1=lmut

!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatx(j,i)
                     lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmaty(j,i)
                     lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut1(1,kk) = lmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut1(2,kk) = lmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut1(3,kk) = lmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatx(j,i)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmaty(j,i)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

 lmut = 0.5d0*(lmut+lmut1)
 finish_omp = omp_get_Wtime()
endif

elseif (.not.tda.and..not.hybrid) then

if (size.gt.1) then 
  if (rank==1) then
!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
          do i=1,kocc
             do ia=1,kvirt
                lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                do ib=1,kvirt
                   if (ib.eq.ia) cycle
                   lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                   lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                   lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                enddo

                do j=1,kocc
                   if (j.eq.i) cycle
                   lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatx(j,i)
                   lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmaty(j,i)
                   lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatz(j,i)
                enddo
             enddo
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
call MPI_Send(lmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 0, 101, MPI_COMM_WORLD,ierr) 
 elseif (rank==0) then 
lmut1=lmut
!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut1(1,kk) = lmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut1(2,kk) = lmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut1(3,kk) = lmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatx(j,i)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmaty(j,i)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

call MPI_Recv(lmut, 3*ntotentr, MPI_DOUBLE_PRECISION, 1, 101, MPI_COMM_WORLD,status,ierr)
write(*,*) 'Received from rank ', status(MPI_SOURCE), ' with tag ', status(MPI_TAG)

lmut = 0.5d0*(lmut+lmut1)
finish_omp = omp_get_Wtime()
      
endif

else

lmut1=lmut
!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
          do i=1,kocc
             do ia=1,kvirt
                lmut(1,kk) = lmut(1,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                lmut(2,kk) = lmut(2,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                lmut(3,kk) = lmut(3,kk) + tddfteig(ia,i,is)*tddfteig1(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                do ib=1,kvirt
                   if (ib.eq.ia) cycle
                   lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                   lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                   lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                enddo

                do j=1,kocc
                   if (j.eq.i) cycle
                   lmut(1,kk) = lmut(1,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatx(j,i)
                   lmut(2,kk) = lmut(2,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmaty(j,i)
                   lmut(3,kk) = lmut(3,kk) - tddfteig(ia,i,is)*tddfteig1(ia,j,js)*lmatz(j,i)
                enddo
             enddo
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do is=1,ntoten
       do js=is,ntoten
          kk = cont(is+1,js+1)
            do i=1,kocc
               do ia=1,kvirt
                  lmut1(1,kk) = lmut1(1,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(1)+lmatx(kocc+ia,kocc+ia)-lmatx(i,i))
                  lmut1(2,kk) = lmut1(2,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(2)+lmaty(kocc+ia,kocc+ia)-lmaty(i,i))
                  lmut1(3,kk) = lmut1(3,kk) + tddfteig1(ia,i,is)*tddfteig(ia,i,js)*(2.d0*l0(3)+lmatz(kocc+ia,kocc+ia)-lmatz(i,i))
                  do ib=1,kvirt
                     if (ib.eq.ia) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatx(kocc+ib,kocc+ia)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmaty(kocc+ib,kocc+ia)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ib,i,js)*lmatz(kocc+ib,kocc+ia)
                  enddo

                  do j=1,kocc
                     if (j.eq.i) cycle
                     lmut1(1,kk) = lmut1(1,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatx(j,i)
                     lmut1(2,kk) = lmut1(2,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmaty(j,i)
                     lmut1(3,kk) = lmut1(3,kk) - tddfteig1(ia,i,is)*tddfteig(ia,j,js)*lmatz(j,i)
                  enddo
              enddo
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

lmut = 0.5d0*(lmut+lmut1)
finish_omp = omp_get_Wtime()

endif !rank 
endif !size
endif !cdspectrum

  deallocate(tddfteig)

  if (dipo) then  
    deallocate(dipmatx)
    deallocate(dipmaty)
    deallocate(dipmatz)
  endif

  if (NP) then
    deallocate(potmat)
!    deallocate(potmat_nuc)
  endif 

  if (cdspectrum) then
     deallocate(lmatx)
     deallocate(lmaty)
     deallocate(lmatz)
  endif


  exciten=exciten*EVAU

if (rank==0) then
  write(*,*) "Sto per scrivere"

  open(24,file='ci_energy.inp')
  do i=1,ntoten
     write(24,'("Root", I5, " : ", F12.5)') i, exciten(i)
  enddo
  close(24)

  deallocate(exciten)


  if (dipo) then  
  open(15,file='ci_mut.inp')
   if (.not.tda) then
  do i=2,ntoten+1
     dd=cont(i,i)
     dmut(1,dd)=dmut(1,dd)-dip_nuc(1) + 2.d0*dip0(1)
     dmut(2,dd)=dmut(2,dd)-dip_nuc(2) + 2.d0*dip0(2)
     dmut(3,dd)=dmut(3,dd)-dip_nuc(3) + 2.d0*dip0(3)
  enddo
  else
  do i=2,ntoten+1
     dd=cont(i,i)
     dmut(1,dd)=dmut(1,dd)-dip_nuc(1)
     dmut(2,dd)=dmut(2,dd)-dip_nuc(2)
     dmut(3,dd)=dmut(3,dd)-dip_nuc(3)
  enddo
  endif
 endif

 if (NP) then
  if (.not.tda) then
   do ints=1,nts
     do i=2,ntoten+1
        dd=cont(i,i)
        potmut(dd,ints)=potmut(dd,ints) + 2.d0*pot0(ints) !*norm(i-1)
     enddo
   enddo
  endif
 endif



  if (dipo) then
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
  endif


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
  
  if (NP) then 
     open(17, file = 'ci_pot.inp')
     write(17,*) nts
     i=0
     j=0
     write(17,*) i,j
     do ints=1,nts
        write(17,*) - potmut(1,ints), 0, potmut_nuc(ints)  
     enddo

     do i=2,ntoten+1
        dd=cont(1,i)
        write(17,*) 0, i-1
        do ints=1,nts
           write(17,*) - potmut(dd,ints)
        enddo
     enddo

      do i=2,ntoten+1
           do j=2,i
              dd=cont(i,j)
              write(17,*)  i-1, j-1
              do ints=1,nts
                 write(17,*) - potmut(dd,ints)
              enddo
           enddo
       enddo
       close(17)
       deallocate(potmut)
!       deallocate(pot_nuc0)
  endif
endif

  deallocate(cont)


if (rank==0) then 

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
endif


call MPI_Finalize(ierr)



stop



END PROGRAM write_wavet 
