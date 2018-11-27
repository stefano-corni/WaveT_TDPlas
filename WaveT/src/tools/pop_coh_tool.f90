program post_processing
!------------------------------------------------------------------------
! @brief Postprocessing for computing populations and/or coherences 
!  
!
! @date Created   : E. Coccia 24 Aug 2018
! Modified  : 
!------------------------------------------------------------------------

 use constants
#ifdef OMP
   use omp_lib
#endif

 implicit none

 integer(i4b)                 :: nstates,n_f,nsteps
 character(flg)               :: tar,all_pop,all_coh,read_bin,write_bin
 integer(i4b)                 :: pop(nstmax)
 character(flg)               :: coh(nstmax*(nstmax-1)/2)

 integer(i4b)                 :: i,j,np,npop,ncoh,itmp
 integer                      :: st,current,rate
 character(20)                :: filein,fileout
 character(4000)              :: fmt_ci,fmt_ci2,dum
 character(4),   allocatable  :: str(:)
 character(5),   allocatable  :: print_short(:)
 character(flg)               :: tmp 
 integer(i4b),   allocatable  :: ii(:),popef(:),icoh(:)
 real(dbl),      allocatable  :: tt(:),rc(:),ic(:)
 complex(cmp),   allocatable  :: c(:,:)          
 character(flg), allocatable  :: cohef(:)


 namelist /pop_coh/ nstates,n_f,nsteps,tar,all_pop,all_coh,read_bin,pop,coh,write_bin 

 pop=-1
 coh=''

 read(*,nml=pop_coh)

 call system_clock(st,rate)

 write(*,*) ''
 write(*,*) '**********************************************'
 write(*,*) '*                                            *'
 write(*,*) '*        PostProcessing Tool for             *'
 write(*,*) '*    computing populations and coherences    *'
 write(*,*) '*                                            *'
 write(*,*) '*                   by                       *'
 write(*,*) '*              Stefano Corni                 *'
 write(*,*) '*                  and                       *'
 write(*,*) '*             Emanuele Coccia                *'
 write(*,*) '*                                            *'
 write(*,*) '**********************************************'
 write(*,*) ''

 if (nstates.gt.nstmax) then
    write(*,*) 'ERROR: nstates =', nstates, '>', nstmax 
    stop
 endif

 if (tar(1:3).eq.'all') then
   write(*,*) 'Compute populations and coherences'
 elseif (tar(1:3).eq.'pop') then
   write(*,*) 'Compute only populations'
 elseif (tar(1:3).eq.'coh') then
   write(*,*) 'Compute only coherences'
 endif
 if (write_bin(1:1).eq.'y') then
    write(*,*) 'Output file(s) in binary format'
 endif
 if (read_bin(1:1).eq.'y') then
    write(*,*) 'Coefficient file in binary format'
 endif

 write(*,*) ''

 if (tar(1:3).eq.'all'.or.tar(1:3).eq.'pop') then
    if (all_pop(1:3).eq.'yes') then
       write(*,*) 'Compute ALL the populations'
    else 
       write(*,*) 'Compute populations selected by input'
       npop=0
!$OMP PARALLEL 
!$OMP DO
       do i=1,nstmax
          if (pop(i).ne.-1) npop=npop+1
       enddo
!$OMP END PARALLEL
       allocate(popef(npop))
       allocate(str(npop))
       do i=1,npop
          if (pop(i).ne.-1) popef(i)=pop(i) 
       enddo
       write(*,*) npop,'state populations requested'
       write(*,*) 'States:', popef(:)
    endif
 endif

 write(*,*) ''
 if (tar(1:3).eq.'all'.or.tar(1:3).eq.'coh') then
    if (all_coh(1:3).eq.'yes') then
       write(*,*) 'Compute ALL the coherences'
    else 
       write(*,*) 'Compute coherences selected by input'
       ncoh=0
!$OMP PARALLEL 
!$OMP DO
       do i=1,nstmax
          tmp=coh(i)
          if (tmp.ne." ") ncoh=ncoh+1
       enddo
!$OMP END PARALLEL
       allocate(cohef(ncoh))
       allocate(icoh(2*ncoh))
       do i=1,ncoh
          tmp=coh(i)
          if (tmp.ne." ") cohef(i)=tmp
       enddo

       call extract_pairs(cohef,ncoh,icoh)

       write(*,*) ncoh,'state-state coherences requested'       
       write(*,*) 'State pairs:'
       do i=1,ncoh
          itmp=2*(i-1)+1
          write(*,*) icoh(itmp),'and',icoh(itmp+1)
       enddo 
       icoh=icoh+1
    endif
 endif

 !nsteps=nsteps+1
 np = nstates*(nstates-1)/2

 allocate(ii(nsteps))
 allocate(tt(nsteps))
 allocate(c(nstates,nsteps))
 allocate(rc(nstates),ic(nstates))

 write(filein,'(a4,i0,a4)') "c_t_",n_f,".dat"

 if (read_bin(1:1).ne.'y') then
    open (9,file=filein,status="unknown")
    read(9,*) dum
    write (fmt_ci,'("(i8,f14.4,",I0,"e17.8E3)")') 2*nstates
    do i=1,nsteps
       read(9,fmt_ci) ii(i), tt(i), (rc(j), ic(j), j=1,nstates)
       do j=1,nstates
          c(j,i) = dcmplx(rc(j),ic(j))
       enddo
    enddo
 else
    open (9,file=filein,status="unknown",form="unformatted")
    do i=1,nsteps
       read(9) ii(i), tt(i), (rc(j), ic(j), j=1,nstates)
       do j=1,nstates
          c(j,i) = dcmplx(rc(j),ic(j))
       enddo
    enddo
 endif

 close(9)

 deallocate(rc)
 deallocate(ic)

 !Computing populations
 if (tar(1:3).eq.'all'.or.tar(1:3).eq.'pop') then
    write(fileout,'(a6,i0,a4)') "pop_t_",n_f,".dat"
    if (write_bin(1:1).ne.'y') then
       open(11,file=fileout,status='unknown')
    else
       open(11,file=fileout,status="unknown",form="unformatted") 
    endif
    !All the populations computed
    if (all_pop(1:3).eq.'yes') then
       if (write_bin(1:1).ne.'y') then
           write(11,'(2a)')'# istep   time (au)', &
                           ' |C_i(t)|^2, i=0,1,2,...   '
           write (fmt_ci,'("(i8,f14.4,",I0,"e14.5E3)")') nstates
           do i=1,nsteps 
              write (11,fmt_ci) ii(i),tt(i),real( c(:,i) * conjg(c(:,i)) )
           enddo
       else
           do i=1,nsteps
              write (11) ii(i),tt(i),real( c(:,i) * conjg(c(:,i)) )
           enddo
       endif
    !Selected-by-input populations computed
    else 
       if (write_bin(1:1).ne.'y') then
           do j=1,npop
              write(str(j),'(I4)') popef(j)
           enddo
           write(11,*)'#istep time (au)', &
                           ' |C_i(t)|^2, i=', str 
           write (fmt_ci,'("(i8,f14.4,",I0,"e14.5E3)")') npop
           do i=1,nsteps 
              write (11,fmt_ci) ii(i),tt(i),( real( c(popef(j)+1,i) * conjg(c(popef(j)+1,i)) ), j=1,npop )
           enddo
       else
           do i=1,nsteps
              write (11) ii(i),tt(i),( real( c(popef(j)+1,i) * conjg(c(popef(j)+1,i)) ), j=1,npop )
           enddo
       endif
    endif
    close(11)
 endif


 !Computing coherences 
 if (tar(1:3).eq.'all'.or.tar(1:3).eq.'coh') then
    write(fileout,'(a6,i0,a4)') "coh_t_",n_f,".dat"
    if (write_bin(1:1).ne.'y') then
       open(12,file=fileout,status='unknown')
     else
       open(12,file=fileout,status="unknown",form="unformatted")
    endif
    !All the coherences computed
    if (all_coh.eq.'yes') then
       if (write_bin(1:1).ne.'y') then
           write(12,'(2a)')'# istep   time (au)', &
                           ' C_i*(t)C_j(t), i=0,1,2,..., j=1,2,3,...   '
           write (fmt_ci2,'("(i8,f14.4,",I0,"2e14.5E3)")') 2*np
       endif  
       do i=1,nsteps
          call wrt_coherence(ii(i),tt(i),12,fmt_ci2,c(:,i),nstates,write_bin)
       enddo
    !Selected-by-input coherences computed
    else
       if (write_bin(1:1).ne.'y') then
          allocate(print_short(ncoh))
          do i=1,ncoh
             print_short(i)=cohef(i)
          enddo
          write(12,*)'# istep   time (au)', &
                           ' C_i*(t)C_j(t),', print_short 
          deallocate(print_short)
          write (fmt_ci2,'("(i8,f14.4,",I0,"2e14.5E3)")') 2*ncoh
          do i=1,nsteps
             write (12,fmt_ci2) ii(i),tt(i),(  c(icoh(2*(j-1)+1),i) * conjg(c(icoh(2*(j-1)+2),i)) , j=1,ncoh )
          enddo
       else
           do i=1,nsteps
             write (12) ii(i),tt(i),(  c(icoh(2*(j-1)+1),i) * conjg(c(icoh(2*(j-1)+2),i)) , j=1,ncoh )
          enddo
       endif
    endif
    close(12) 
 endif

 if (allocated(popef)) deallocate(popef) 
 if (allocated(cohef)) deallocate(cohef)
 if (allocated(icoh))  deallocate(icoh)
 if (allocated(str))   deallocate(str)
 deallocate(ii)
 deallocate(tt)
 deallocate(c)


 call system_clock(current)
 write(*,*) '' 
 write(6,'("Done , total elapsed time", &
           F10.3,"s")') real(current-st)/real(rate)

 stop

end program post_processing


subroutine wrt_coherence(i,t,int1,char2,c,nci,bin)
!------------------------------------------------------------------------
! @brief Print the tridiagional C*_iC_j (i.ne.j) matrix 
! corresponding to the decoherence 
! 
! @date Created   : E. Coccia 3 Feb 2017
! Modified  :
!------------------------------------------------------------------------

        use constants

        implicit none

        integer(i4b),    intent(in)    :: i, nci
        integer(i4b),    intent(in)    :: int1
        character(4000), intent(in)    :: char2
        real(dbl),       intent(in)    :: t
        character(flg),  intent(in)    :: bin
        complex(cmp),    intent(in)    :: c(nci)
        complex(cmp)                   :: cc(nci,nci)
        integer(i4b)                   :: j, k
        complex(cmp)                   :: tmp

        tmp=dcmplx(0.d0,0.d0)

        do k=2,nci
           tmp = dconjg(c(k))*c(1)
           cc(1,k) = tmp
        enddo

        if (nci.gt.2) then
           do j=2,nci
              do k=j+1,nci
                 tmp = dconjg(c(k))*c(j)
                 cc(j,k) = tmp
              enddo
           enddo
        endif

        if (bin(1:1).ne.'y') then
           if (nci.gt.2) then
              write(int1,char2) i, t, ((cc(j,k), k=j+1,nci),j=1,nci-1)
           else
              write(int1,char2) i, t, cc(1,2)
           endif
        else
           if (nci.gt.2) then
              write(int1) i, t, ((cc(j,k), k=j+1,nci),j=1,nci-1)
           else
              write(int1) i, t, cc(1,2)
           endif
        endif

        return

end subroutine wrt_coherence


subroutine extract_pairs(str,n,into)
!------------------------------------------------------------------------
! @brief Extract state pairs for computing coherence 
! 
! @date Created   : E. Coccia 24 Aug 2018
! Modified  :
!------------------------------------------------------------------------

  use constants

  implicit none

  integer(i4b),   intent(in)     ::  n
  character(flg), intent(in)     ::  str(n)
  integer(i4b),   intent(inout)  ::  into(2*n)

  integer(i4b)    :: i,ipos,itmp  
  character(flg)  :: tmp

  do i=1,n

     tmp=str(i)
     itmp=2*(i-1)+1

     ipos = scan(tmp,"-")!,back=.true.)
     !print*,"reading real variable from '" // trim(tmp(1+ipos:)) // "'"
     read (tmp(1+ipos:),*) into(itmp+1) 
     read (tmp(:ipos-1),*) into(itmp)
     !print*,"States pair ", into(itmp),into(itmp+1) 

  enddo

  return

end subroutine extract_pairs




