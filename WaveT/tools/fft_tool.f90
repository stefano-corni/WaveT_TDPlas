program twodfft

 use constants 
 use, intrinsic :: iso_c_binding

 implicit none

 integer(i4b)               :: dim1,dim2,i,j,plan
 real(dbl)                  :: rdum,dw1,dw2,dt1,dt2,modd
 real(dbl)                  :: dir(3)
 real(dbl), allocatable     :: dinp(:,:,:),inp(:,:) 
 complex(cmp), allocatable  :: doutp(:,:) 
 character*20               :: filename

 include 'fftw3.f03'

 namelist /fft/ dim1,dim2,dt1,dt2,filename,dir 

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

 if (dim1.le.0) then
    write(*,*) 'ERROR: dim1 must be larger than zero'
    stop
 endif

 if (dim2.le.0) then
    write(*,*) 'ERROR: dim2 must be larger than zero'
    stop
 endif

 if (dt1.le.0) then
    write(*,*) 'ERROR: dt1 must be larger than zero'
    stop
 endif

 if (dt2.le.0) then
    write(*,*) 'ERROR: dt2 must be larger than zero'
    stop
 endif

 if (dir(1).le.0.and.dir(2).le.0.and.dir(3).le.0) then
    write(*,*) 'ERROR: One component of dir_ft must be larger than zero'
    stop
 endif

 write(*,*) ''
 write(*,*) 'Number of points for t1', dim1
 write(*,*) 'Number of points for t3', dim2
 write(*,*) 'Time step along the first diemnsion', dt1
 write(*,*) 'Time step along the second dimension', dt2 
 write(*,*) 'Direction for FFT', dir(1),dir(2),dir(3)
 write(*,*) ''
 
 write(*,*) 'Data from file', filename
 write(*,*) ''

 if(mod(dim1,2).gt.0) dim1=dim1-1
 if(mod(dim2,2).gt.0) dim2=dim2-1

 allocate (dinp(dim1,dim2,3))
 allocate (inp(dim1,dim2))
 allocate (doutp(int(dble(dim1)/two)+1,int(dble(dim2)/two)+1))
 dw1=2*pi/dble(dim1)/dt1
 dw2=2*pi/dble(dim2)/dt2

 open(7,file=filename)

 write(*,*) ''
 write(*,*) 'File from 2D spectroscopy calculation'
 write(*,*) 'is supposed to be with the following form:'
 write(*,*) 't1 t3 mu(t1,tr3)'
 write(*,*) ''

 do j=1,dim1
    do i=1,dim2
       read(7,*) rdum, rdum, dinp(i,j,:)
    enddo
 enddo

 close(7)

 doutp=cmplx(0.d0,0.d0)

 dir(:)=dir(:)/sqrt(dot_product(dir,dir))

 do j=1,dim1
    do i=1,dim2
       inp(i,j)=dot_product(dinp(i,j,:),dir(:))
    enddo
 enddo

 call dfftw_plan_dft_r2c_2d(plan,dim1,dim2,inp,doutp,FFTW_ESTIMATE)
 call dfftw_execute_dft_r2c(plan,inp,doutp)
 call dfftw_destroy_plan(plan)

 open(unit=8,file='2d_spectrum.dat',status="unknown",form="formatted")
 do i=1,int(dim1/two)
    do j=1,int(dim2/two)
       modd=sqrt(real(doutp(i,j))**2+aimag(doutp(i,j))**2)
       write(8,'(3e20.10)') (i-1)*dw1, (j-1)*dw2, modd
    enddo 
 enddo

 close(8)

 deallocate(dinp,doutp,inp)
 
 stop


end program twodfft
