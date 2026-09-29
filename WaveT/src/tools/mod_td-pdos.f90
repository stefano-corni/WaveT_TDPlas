module params
 integer, parameter      :: natmx=500
 integer, parameter      :: nmxfr=natmx
 integer, parameter      :: nmxatfr=natmx, nmxatfrOP=natmx
 real*8, parameter       :: au2ev=27.21138386d0
 real*8, parameter       :: pi=dacos(-1.d0)
end module params

module mod_td 
 use params
 implicit none
 save
 integer                 :: nstates,nsteps,nfreq,nfrag,neff,np
 !lookhere
 integer                 :: nfragOP !
 integer                 :: nmo,nocc,nvir,nexc,nat,nf
 integer, allocatable    :: ist(:),icount(:)
 character(len=10)       :: average,frozen
 character*1             :: binary,mixing,tcm1
 real*8                  :: sigma,emin,emax,ehomo
 real*8, allocatable     :: emo(:),t(:), dm(:)
 real*8, allocatable     :: w(:,:),a(:,:,:),ksum(:)
 real*8, allocatable     :: pdos(:,:),tot_dos(:),gs_pdos(:,:),ov_pdos(:,:)
 real*8, allocatable     :: opdos(:,:,:),gs_opdos(:,:,:) !
 real*8, allocatable     :: mix(:,:)
 !real*8, allocatable     :: tcm(:,:,:)
 real*8, allocatable     :: tcm(:,:), ftcm(:,:,:)
 real*8, allocatable, dimension(:,:,:) :: wOP(:,:,:)
 complex*16, allocatable :: corr(:,:,:)

 complex(8), allocatable :: c(:,:)

 type tfrag
    integer :: natoms
    integer :: vec_atoms(nmxatfr)
 end type tfrag 
 type(tfrag), allocatable :: frag(:)
 
 !lookhere
 type tfragOP !
    integer :: natomsOP !
    integer :: vec_atomsOP(nmxatfrOP) !
 end type tfragOP !
 type(tfragOP), allocatable :: fragOP(:) !

 contains
!------------------------------------------------------------------------
! @brief Read input file 
!
! @date Created   : E. Coccia 19/6/20 
! Modified  :  D. Toffoli, 01/07/20 
!------------------------------------------------------------------------
subroutine read_input()
 implicit none
 integer :: i,iat,j,idx,pos,ist,iend,idum,natot,natotOP,n,pos1
 character(len=120) :: s_line,s_sub,s1_sub,s2_sub, s_dum, s3_sub, t1_sub, t2_sub

 !lookhere
 namelist /general/nstates,nsteps,nfreq,average,nfrag,nfragOP,sigma,emin,emax,np,ehomo,frozen,nf,binary,mixing,tcm1 !

 !Initializing variables
 !Number of excited states
 nstates=1
 !Numbef of time steps
 nsteps=100000
 !Frequency of time steps printing
 nfreq=1
 !Read average file or single trajectory
 average='no'
 !Number of fragments 
 nfrag=1
 !lookhere Number of OPfragments
 nfragOP=1 !
 !sigma for Lorentzian convolution
 sigma=0.01
 !Energy range for plotting (Ha)
 emin=-10.d0
 emax=10.d0
 !HOMO energy (Ha)
 ehomo=0.d0
 !frozen-core option
 frozen='non'
 !number of frozen-core orbitals
 nf=0
 !inputs are non-formatted?
 binary='n'
 !if WaveT states have been SCF mixed (NP or solvent), mixing='yes'
 mixing='n'
 !if I want to do TCM
 tcm1='y' 

 read(*,nml=general)
 !neff -> effective number of time steps
 neff=nsteps/nfreq  
 write(*,nml=general)
 write(*,*) '' 
 write(*,*) 'Effective number of steps', neff 
 write(*,*) 'HOMO energy in eV', ehomo*au2ev
 write(*,*) 'Energy range for PDOS in eV', emin*au2ev, emax*au2ev
 write(*,*) 'sigma for Lorentzian convolution in eV', sigma*au2ev
 if (mixing.eq.'y') write(*,*) 'Using mixed WaveT states'
 write(*,*) ''
 nexc=nstates-1

 allocate(frag(nfrag))
 do i = 1, nfrag
    frag(i)%natoms=0
    frag(i)%vec_atoms=0.0d0
 end do
 do i = 1, nfrag
    read(*,"(a)") s_dum
    !write(*,*) ' s_dum: ', s_dum
    !discard everything after !
    pos=index(s_dum,"!")
    s_line=s_dum(1:pos-1)
    s_line=adjustl(trim(s_line))
    !write(*,*) ' s_line: ', s_line
    iat = 0
    do
       pos=index(s_line," ")
       s_sub=adjustl(trim(s_line(1:pos-1)))
       s_line=adjustl(trim(s_line(pos+1:)))
       !write(*,*) ' s_sub: ', s_sub
       pos=index(s_sub,"-")
       if (pos .eq. 0) then
          !just a number
          iat=iat+1
          read(s_sub,*) idx
          write(*,*) 'indice', idx
          frag(i)%vec_atoms(iat)= idx
       else
          !we have an interval
          s1_sub=adjustl(trim(s_sub(1:pos-1)))
          s2_sub=adjustl(trim(s_sub(pos+1:)))
          read(s1_sub,*) ist
          read(s2_sub,*) iend
          n=abs(iend-ist)+1
          !write(*,*) ' ist, iend, nat: ', ist, iend, nat
          if (ist .gt. iend) then
             idum = ist
             ist = iend
             iend = idum
          end if
          do j = 1, n
             iat = iat + 1
             frag(i)%vec_atoms(iat) = ist + j -1
          end do
       end if
       !look for a newline mark (\)
       pos=index(s_line,"\")
       if ((len_trim(s_line) .eq. 1) .and. (pos .ne. 0)) then
          read(*,"(a)") s_line
       else if (len_trim(s_line) .eq. 0) then
          frag(i)%natoms=iat
          write(*,*) 'natoms', iat
!          exit
       end if
       !looking for another atom, being part of the fragment, but outside the main
       !interval LB
       pos=index(s_line,",")
       if (pos .ne. 0) then
          pos1 = index(s_line,"-")
          if (pos1 .ne. 0) then
              t1_sub=adjustl(trim(s_line(2:pos1-1)))
              t2_sub=adjustl(trim(s_line(pos1+1:)))
              read(t1_sub,*) ist
              read(t2_sub,*) iend
              write(*,*) 'ist', ist, 'iend', iend
              n=abs(iend-ist)+1
              if (ist .gt. iend) then
                 idum = ist
                 ist = iend
                 iend = idum
              end if
              do j = 1, n
                 iat = iat + 1
                 frag(i)%vec_atoms(iat) = ist + j -1
              end do
                 frag(i)%natoms=iat
                 write(*,*) 'natoms_2', iat
                 exit
          else
             iat = iat + 1 
             s3_sub=adjustl(trim(s_line(pos+1:)))
             read(s3_sub,*) idx
             frag(i)%vec_atoms(iat) = idx
             frag(i)%natoms=iat
             write(*,*) 'natoms', iat
             exit
          endif
       else
          exit
       endif
    end do
 end do
 natot=sum(frag(:)%natoms)
 write(*,*) 'Total nr. of atoms: ', natot
 write(*,*) 'Atomic fragments: '
 do i = 1, nfrag
    write(*,"('i= ', i4, 2x, ' atom #: ', 10(i4, 2x))") i, &
         (frag(i)%vec_atoms(iat), iat=1, frag(i)%natoms)
 end do
 !lookhere
 !return

 !lookhere everything until the end of the subroutine
 !read(*,'(a)') s_dum
 !write(*,*) ""
 !allocate (fragOP(nfragOP)) !
 !do i = 1, nfragOP
 !   fragOP(i)%natomsOP=0
 !   fragOP(i)%vec_atomsOP=0.0d0
 !end do
 !do i = 1, nfragOP
 !   read(*,"(a)") s_dum
 !   !write(*,*) ' s_dum: ', s_dum
 !   !discard everything after !
 !   pos=index(s_dum,"!")
 !   s_line=s_dum(1:pos-1)
 !   s_line=adjustl(trim(s_line))
 !   !write(*,*) ' s_line: ', s_line
 !   iat = 0
 !   do
 !      pos=index(s_line," ")
 !      s_sub=adjustl(trim(s_line(1:pos-1)))
 !      s_line=adjustl(trim(s_line(pos+1:)))
 !      !write(*,*) ' s_sub: ', s_sub
 !      pos=index(s_sub,"-")
 !      if (pos .eq. 0) then
 !         !just a number
 !         iat=iat+1
 !         read(s_sub,*) idx
 !         fragOP(i)%vec_atomsOP(iat)= idx
 !      else
 !         !we have an interval
 !         s1_sub=adjustl(trim(s_sub(1:pos-1)))
 !         s2_sub=adjustl(trim(s_sub(pos+1:)))
 !         read(s1_sub,*) ist
 !         read(s2_sub,*) iend
 !         n=abs(iend-ist)+1
 !         !write(*,*) ' ist, iend, nat: ', ist, iend, nat
 !         if (ist .gt. iend) then
 !            idum = ist
 !            ist = iend
 !            iend = idum
 !         end if
 !         do j = 1, n
 !            iat = iat + 1
 !            fragOP(i)%vec_atomsOP(iat) = ist + j -1
 !         end do
 !      end if
 !      !look for a newline mark (\)
 !      pos=index(s_line,"\")
 !      if ((len_trim(s_line) .eq. 1) .and. (pos .ne. 0)) then
 !         read(*,"(a)") s_line
 !      else if (len_trim(s_line) .eq. 0) then
 !         fragOP(i)%natomsOP=iat
 !         exit
 !      end if
 !   end do
 !end do
 !natotOP=sum(fragOP(:)%natomsOP)
 !write(*,*) 'Total nr. of OP atoms: ', natotOP
 !write(*,*) 'Atomic OP fragments: '
 !do i = 1, nfragOP
 !   write(*,"('i= ', i4, 2x, ' atom #: ', 10(i4, 2x))") i, &
 !        (fragOP(i)%vec_atomsOP(iat), iat=1, fragOP(i)%natomsOP)
 !end do
 !stop 'lookhere2'
 return
end subroutine read_input
!------------------------------------------------------------------------
! @brief Read input files from ADF or MolGW 
!
! @date Created   : E. Coccia 19/6/20 
! Modified  : D. Toffoli  01/07/20 
!------------------------------------------------------------------------
subroutine read_adf_molgw()
 implicit none
 !lookhere
 integer :: nsym,isym,iocc,ivir,is,i,j,iline,nmor,nmorOP, totpairs,nexcr,iexc !
 integer :: idum,naos,indmoi,ipair,mu,nu,ii,jj !
 integer, allocatable, dimension(:)    :: vec_nocc,vec_vir,alpha !
 real*8 :: edum, cnorm,wsum
 real*8, allocatable, dimension(:,:)   :: mllkn_pop,wr,wrOP!
 real*8, allocatable, dimension(:)   :: swr
 character*5 :: sroot
 character*5 :: basename
 character*30 :: filename
 character*300 :: sdum

 open(10,file='mos_info.dat')
 read(10,*) nmo
 read(10,*) nsym
 !write(*,*) 'nmo, nsym: ', nmo, nsym
 allocate(emo(nmo), vec_nocc(nsym), vec_vir(nsym))
 nocc = 0
 nvir = 0
 do is = 1, nsym
    read(10,*) isym, vec_nocc(is), vec_vir(is)
    !write(*,*) isym, vec_nocc(is), vec_vir(is)
    nocc = nocc + vec_nocc(is) 
    nvir = nvir + vec_vir(is) 
 end do

 open(11,file='mulliken_pop.dat')
 read(11,*) nat, nmor
 if (nmor .ne. nmo) stop ' nmor .ne. nmo '
 !lookhere
 !open(12,file='ovrl_pop.dat') !------
 !EC: 23/9/20
 !read(12,*) totpairs,nmorOP !------   
 !if (nmorOP .ne. nmo) stop ' nmorOP .ne. nmo '!------   
 allocate(wr(nat,nmo))!,wrOP(totpairs,nmo)) !------   
 wr=0.0d0
 !wrOP=0.0d0 !------   

 iocc = 0
 ivir = nocc
 do is = 1, nsym
    do i = 1, vec_nocc(is)
       read(10,*) edum
       iocc = iocc+1
       emo(iocc) = edum
       read(11,*) (wr(j,iocc), j=1,nat)
       !EC: 23/9/20
       !read(12,*) (wrOP(j,iocc), j=1,totpairs) !------   
    end do
    do i = 1, vec_vir(is)
       read(10,*) edum
       ivir = ivir + 1
       emo(ivir) = edum
       read(11,*) (wr(j,ivir), j=1,nat)
       !EC: 23/9/20
       !read(12,*) (wrOP(j,ivir), j=1,totpairs) !------   
    end do
 end do

 allocate(swr(nmo))
 !EC: 2/12/22, no symmetry used!
 swr=0.d0
 do is = 1, nsym
    do i = 1, nmo 
       do j=1,nat
          swr(i)=swr(i)+wr(j,i)
       enddo
    end do
 end do
 wr=abs(wr)
 do is = 1, nsym
    do i = 1, nmo 
       wr(:,i)=wr(:,i)/swr(i)
    end do
 end do
 deallocate(swr)
 !open(unit=13,file='AO_map.dat') !------
 !EC:23/9/20 
 !read(13,*) naos !------
 !read(13,*) sdum !------
 allocate(alpha(naos)) !------
 do i=1,naos !------
   !EC: 23/9/20
   !read(13,*) idum,alpha(i) !------
 enddo !------

 close(10)
 close(11)
 close(12)
 close(13)

 !debug
 !write(*,*) 'MOs eigenvalues: '
 !write(*,*) 'nocc: ', nocc
 !write(*,*) 'nvir: ', nvir
 !do i = 1, nmo
    !write (*,"(i4, 2x, f15.6)") i, emo(i)
 !end do 
 !arrange atomic weights into fragments
 allocate(w(nfrag,nmo))
 w=0.0d0
 do i=1,nfrag
    do j=1,frag(i)%natoms
       w(i,:) = w(i,:)+wr(frag(i)%vec_atoms(j),:)
    end do
 end do
 deallocate(wr)

 !wsum=0.d0
 !do i=1,nfrag
 !   wsum=wsum+sum(w(i,:))
 !enddo

 !w=w/wsum

 !write(*,*) 'ciao'
 !write(*,*) sum(w(1,:)),sum(w(2,:)),sum(w(3,:)),sum(w(1,:))+sum(w(2,:))+sum(w(3,:))
 !stop
 
 !lookhere
 !arrange overlap population in between fragment i and j;
 !and relating it to the MOs (indmoi)
 !allocate(wOP(nfragOP,nfragOP,nmo))
 !wOP=0.0d0
 !do indmoi=1,nmo !Run over the MOs
 !  do i=1,nfragOP      !Run over 2 fragments
 !    do j=i+1,nfragOP  !**
 !      do ii=1,fragOP(i)%natomsOP   !Run over nb atoms of 2 fragments
 !        do jj=1,fragOP(j)%natomsOP !**
 !          ipair=0
 !          do mu=1,naos      !Run over the AOs
 !            do nu=mu+1,naos !**
 !            ipair=ipair+1
 !            !Check if the pairs of AOs (mu,nu) are centered in the pairs of atoms
 !            !selected within the two fragments (fragOP(i),fragOP(j)), giving a
 !            !correspondent MO (indmoi) in which they have a contribution
 !            if (alpha(mu)==fragOP(i)%vec_atomsOP(ii) .AND. &
 !               &alpha(nu)==fragOP(j)%vec_atomsOP(jj)) then
 !            wOP(i,j,indmoi)=wOP(i,j,indmoi)+wrOP(ipair,indmoi) !wOP(i,j,indmoi)+ 
 !            endif !Different ipairs might have some contribution for a given
 !            enddo !indmoi
 !          enddo
 !        enddo
 !      enddo
 !    enddo
 !  enddo
 !enddo
 

 !debug
 !write(*,*) 'Fragment weights: '
 !do i=1, nfrag
 !   write(*,*) ' frag #: ', i
 !   write(*,*) ' w(i): '
 !   write(*,*) (w(i,j), j = 1, nmo)
 !end do
 !EC: 23/9/20
 open(10, file='tddft_info.dat')
 read(10,*) nexcr
 close(10)
 write(*,*) 'nexcr: ', nexcr
 if (nexcr .ne. nexc) stop ' nexcr .ne. nexc! '
 allocate(a(nocc,nvir,nexc))
 a=0.d0
 basename ='o-cfc'
 !LB
 if (nexc.gt.1000) then
   write(*,*) 'sono formattato'
 do iexc = 1, nexc
    write(sroot,'(i5.5)') iexc
    write(filename,*) adjustl(trim(basename)),adjustl(trim(sroot)),".dat"
    open(10,file=filename,form='unformatted')
!    read(10) sdum
!    read(10) sdum
    read(10) 
    read(10) 
!    MM: I need to comment the next read to read molgw cfc files 
!    read(10,"(a)") sdum
    cnorm = 0.0d0
    !EC: nocc-nf
    do iline = 1, (nocc-nf)*nvir 
!       read(10, end=100) i, j, idum, a(i,j-nocc,iexc)
       read(10, end=100) i, j, idum, a(i,j-nocc,iexc)
!       write(*,*) i, j-nocc, idum, a(i,j-nocc,iexc)
       cnorm = cnorm + a(i,j-nocc,iexc)**2
    end do 
    !normalization
    !write(*,*) ' state: ', iexc, ' norm: ', sqrt(cnorm)
    a(:,:,iexc)=a(:,:,iexc)/sqrt(cnorm)
    close(10)
 end do
 else 
 write(*,*) 'non sono formattato'
 do iexc = 1, nexc
    write(sroot,'(i5.5)') iexc
    write(filename,*) adjustl(trim(basename)),adjustl(trim(sroot)),".dat"
    !write(*,*) ' filename: ', filename
    open(10,file=filename)
    read(10,"(a)") sdum
    read(10,"(a)") sdum
!    MM: I need to comment the next read to read molgw cfc files 
!    read(10,"(a)") sdum
    cnorm = 0.0d0
    !EC: nocc-nf
    do iline = 1, (nocc-nf)*nvir 
!    write(*,*) '4', (nocc-nf)*nvir  
       read(10,*, end=100) i, j, idum, a(i,j-nocc,iexc)
!       write(*,*) i, j-nocc, idum, a(i,j-nocc,iexc)
       cnorm = cnorm + a(i,j-nocc,iexc)**2
    end do 
    !normalization
    !write(*,*) ' state: ', iexc, ' norm: ', sqrt(cnorm)
    a(:,:,iexc)=a(:,:,iexc)/sqrt(cnorm)
    close(10)
 end do
 endif

 goto 200
100 write(*,*) 'Error reading from file ', filename
 stop
200 continue
 return
end subroutine read_adf_molgw
!------------------------------------------------------------------------
! @brief Read coefficient file from WaveT 
!
! @date Created   : E. Coccia 19/6/20 
! Modified  : D. Toffoli 07/07/20, G. Dall'Osto 20/05/21
!------------------------------------------------------------------------
subroutine read_wavet()

 implicit none

 integer             :: i,j,k,ii,kk,idum,itmp
 real*8              :: rdum,psum
 real*8, allocatable :: re(:),im(:),rtmp(:),imtmp(:)
 complex*16, allocatable :: tmp(:,:),cread(:,:)
 character*30        :: cdum

 allocate(t(neff))
 allocate(ist(neff))
 allocate(corr(neff,nstates,nstates))
 allocate(re(nstates))
 allocate(im(nstates))
 allocate(c(nstates,neff))
 if (average.eq."yes") then
      allocate(tmp(neff,nstates*(nstates+1)/2))
      allocate(rtmp(nstates*(nstates+1)/2),imtmp(nstates*(nstates+1)/2))
      open(10,file='c_t_avg.dat')
      open(11,file='corr_t_avg.dat')
      if (binary.eq.'y') then
         read(10)
      else
         read(10,*)
      endif
      ii=0
      itmp=nfreq
      if (binary.ne.'y') then
         do i=1,nsteps
            if (i.eq.itmp) then
               ii=ii+1
               read(10,*) ist(ii),t(ii),(re(j),im(j),j=1,nstates)
               do j=1,nstates
                  c(j,ii)=dcmplx(re(j),im(j))
               enddo
               !read(11,*) ist(ii), t(ii),(tmp(ii,j), j=1,(nstates*(nstates+1)/2))
               read(11,*) ist(ii), t(ii),(rtmp(j),imtmp(j), j=1,(nstates*(nstates+1)/2))
               do j=1,nstates*(nstates+1)/2
                  tmp(ii,j)=dcmplx(rtmp(j),imtmp(j))
               enddo 
               itmp=itmp+nfreq
               if (itmp .gt. nsteps) exit
            else
               read(10,*) idum,rdum,(rdum,rdum,k=1,nstates)
               read(11,*) idum,rdum,(rdum,rdum,j=1,(nstates*(nstates+1)/2))
            endif
         enddo
      else
         do i=1,nsteps
            if (i.eq.itmp) then
               ii=ii+1
               read(10) ist(ii),t(ii),(re(j),im(j),j=1,nstates)
               do j=1,nstates
                  c(j,ii)=dcmplx(re(j),im(j))
               enddo
               read(11) ist(ii), t(ii),(rtmp(j),imtmp(j),j=1,(nstates*(nstates+1)/2))
               do j=1,nstates*(nstates+1)/2
                  tmp(ii,j)=dcmplx(rtmp(j),imtmp(j))
               enddo 
               itmp=itmp+nfreq
               if (itmp .gt. nsteps) exit
            else
               read(10) idum,rdum,(rdum,rdum,k=1,nstates)
               read(11) idum,rdum,(rdum,rdum,j=1,(nstates*(nstates+1)/2))
            endif
         enddo
      endif
      do i=1,neff
         kk=0
         do j=1,nstates
            do k=j,nstates
               kk=kk+1
               corr(i,j,k)=tmp(i,kk)
               corr(i,k,j)=conjg(tmp(i,kk))
            enddo
         enddo
      enddo
      deallocate(tmp)
      deallocate(rtmp,imtmp)
      close(10)
      close(11)
 else
      open(10,file='c_t_1.dat')
      read(10,*) cdum 
      ii=1
      itmp=nfreq
      do i=1,nsteps
         if (i.eq.itmp) then
            write(*,*) 'MM, ii, itmp,nfreq,neff', ii, itmp, nfreq, neff
            read(10,*) ist(ii),t(ii),(re(j),im(j),j=1,nstates)
            do j=1,nstates
               c(j,ii)=dcmplx(re(j),im(j))
            enddo  
            ii=ii+1
            itmp=itmp+nfreq
            if (itmp .gt. nsteps) exit
         else
            read(10,*) idum,rdum,(rdum,rdum,j=1,nstates)
         endif
      enddo

      !do i=1,ii
      !   psum=sum(abs(c(:,i)**2)) 
      !   c(:,ii)=c(:,ii)/dsqrt(psum)
      !enddo
     
      close(10)
 
      deallocate(re)
      deallocate(im)
 endif

 if (mixing.eq.'y') then
    allocate(mix(nstates,nstates))
    open(19, file="mixing_wavet.dat", status="old")
    do i=1,nstates
       read(19,*) (mix(i,k), k=1,nstates)
    enddo
    close(19)
    if (average.eq.'no') then
        allocate(cread(nstates,neff))
        cread=c
        c=0.d0
        ii=1
        itmp=nfreq
        do i=1,nsteps
             if (i.eq.itmp) then
                write(*,*) 'MM, ii, itmp,nfreq,neff', ii, itmp, nfreq, neff
                do k=1,nstates
                   do j=1,nstates
                      c(k,ii)=c(k,ii)+cread(j,ii)*mix(j,k)
                   enddo
                enddo
                ii=ii+1
                itmp=itmp+nfreq
                if (itmp .gt. nsteps) exit
             endif
        enddo
        deallocate(cread)
    endif
 endif

 return

end subroutine read_wavet

!------------------------------------------------------------------------
! @brief Compute td-PDOS 
!
! @date Created   : E. Coccia 19/6/20 
! Modified  : D. Toffoli 01/07/20
!------------------------------------------------------------------------
subroutine compute_td_pdos()

 implicit none

 !lookhere
 integer             :: i2,ii,aa,l,lp,ll,i,j,k,kk,iocc,ivir,npo,npv !
 integer             :: istart
 real*8              :: h,fc,ltmp,atmp
 real*8, allocatable :: e(:),z(:),integ(:), dm_avg(:), z_avg(:), zz(:)
 complex*16, allocatable :: zc(:,:), zc_avg(:)
 real*8, allocatable :: integ_avg(:), pdos_avg(:,:),ee(:),hh(:)
 real*8, allocatable :: einteg(:), einteg_avg(:),ee_avg(:),hh_avg(:)
 real*8, allocatable :: hinteg(:),hinteg_avg(:)
 real*8, allocatable :: gs_pdost(:,:)
 logical             :: lex
 character*5         :: sfr,sfr2
 character*40        :: filename,basename,filename1,basename1,basename2,filename2
 character*40        :: filename3,basename3
 character*3000      :: fmt_np

 allocate(pdos(np,nfrag))
 allocate(opdos(np,nfragOP,nfragOP))
 allocate(tot_dos(np))
 allocate(gs_pdos(np,nfrag))
 allocate(gs_opdos(np,nfragOP,nfragOP))
 allocate(ov_pdos(np,nfrag))
 allocate(e(np))
 allocate(zc(nocc,nvir))
 if (average.eq."yes")then
    allocate(zc_avg(nmo))
    allocate(z_avg(nmo))
    allocate(dm_avg(nmo))
    allocate(integ_avg(nfrag))
    allocate(einteg_avg(nfrag))
    allocate(hinteg_avg(nfrag))
    allocate(pdos_avg(np,nfrag))
    allocate(ee_avg(nfrag))
    allocate(hh_avg(nfrag))
 endif
 if (mixing.eq."y") then
         allocate(zz(nmo))
         allocate(gs_pdost(np,nfrag))
 endif
 allocate(z(nmo), dm(nmo))
 allocate(ksum(nfrag))
 allocate(integ(nfrag))
 allocate(einteg(nfrag))
 allocate(hinteg(nfrag))
 allocate(ee(nfrag))
 allocate(hh(nfrag))
!DPDOS(t,e) = PDOS(t,e) - PDOS^GS(e)
! = -2 sum_i^occ w_at^i Re[ sum_ll' [C_l'^*(t)C_l(t)] sum_c^vir A_ci^l' A_ci^l ]
! delta(e-e_i)
!   +2 sum_i^vir w_at^i Re[ sum_ll' [C_l'^*(t)C_l(t)] sum_v^occ A_iv^l' A_iv^l ]
!   delta(e-e_i)

!indexes i,v and c run over the set of MOs
!indexes l and l' run over the set of excited states
!|l> = sum_cv A_cv^l a_c^dagger a_v |GS>
!|Psi > = sum_l C_l(t) |l>
 !write(*,*) 'MOs eigenvalues: '
 !write(*,*) 'nocc: ', nocc
 !write(*,*) 'nvir: ', nvir
 !do i = 1, nmo
    !write (*,"(i4, 2x, f15.6)") i, emo(i)
 !end do

 h=(emax-emin)/dble(np)
 do k=1,np
    e(k)=emin+k*h
 enddo
 fc=sigma/pi

 npo=0
 k=1
 do while (e(k).lt.ehomo)
    k=k+1
    npo=npo+1
 enddo
 npv=np-npo

 open(55,file='input_avg')
 write(55,*) nfrag
 write(55,*) neff
 write(55,*) np
 write(55,*) npo
 write(55,*) nmo 
 close(55)

 if (frozen.eq.'yes') then
    istart=nf+1
 else
    istart=1
 endif

 !Ground-state PDOS
 gs_pdos=0.d0
 do j=1,nfrag
    do k=1,np
       !EC istart=1 or istart=nf+1
       do ii=istart,nocc
          gs_pdos(k,j)=gs_pdos(k,j)+w(j,ii)*fc/((e(k)-emo(ii))**2+sigma**2)
       end do
    enddo
 enddo
 gs_pdos=2.d0*gs_pdos

 if (mixing.eq.'y') then
    zz =  0.d0
    do l=1,nexc
       do lp=1,nexc
           ltmp=mix(1,l+1)*mix(1,lp+1)
           !EC istart=1 or istart=nf+1
           do iocc = istart, nocc
              do ivir=1,nvir
                 atmp=a(iocc,ivir,l)*a(iocc,ivir,lp) 
                 zz(iocc)=zz(iocc) + ltmp*atmp
                 !zz(nocc+ivir)=zz(nocc+ivir) + ltmp*atmp
              enddo
           enddo
           do ivir=1,nvir
              do iocc = istart, nocc
                 atmp=a(iocc,ivir,l)*a(iocc,ivir,lp)
                 zz(nocc+ivir)=zz(nocc+ivir) + ltmp*atmp
              enddo
           enddo
       enddo
    enddo

    gs_pdost=0.d0
    open(23, file="gs_pdost.dat", status="new", action="write")
    write(23,*) '#Energy [e-e_HOMO](eV)    #GS-PDOS_tilde(e) '

    do k=1,np
       do j=1,nfrag
          !EC istart=1 or istart=nf+1
          do ii=istart,nocc
             gs_pdost(k,j)=gs_pdost(k,j) - w(j,ii)*zz(ii)*fc/((e(k)-emo(ii))**2+sigma**2)
          end do
          do aa=1,nvir
             gs_pdost(k,j)=gs_pdost(k,j)+w(j,nocc+aa)*zz(nocc+aa)*fc/((e(k)-emo(nocc+aa))**2+sigma**2)
          end do
       end do
       write(23,*) (e(k)-ehomo)*au2ev,(gs_pdost(k,j), j=1,nfrag)
    end do
    close (23)

 endif


! if (maxval(abs(gs_pdos(:,:))).gt.0.d0) then
!    gs_pdos(:,:)=gs_pdos(:,:)/maxval(abs(gs_pdos(:,:)))
! endif

 write (fmt_np,'("(f14.4,",I0,"e17.8E3)")') np
 basename="gs-pdos_frag"
 do j = 1, nfrag
    write(sfr,'(i5.5)') j   !vengono # da 000001 a 100000
    write(filename,*) adjustl(trim(basename)), adjustl(trim(sfr)),".dat"
    inquire(file=filename,exist=lex)
    if (lex) then
       open(11, file=filename, status="old",position="append",action="write")
     else
       open(11, file=filename, status="new", action="write")
    end if
    write(11,*) '#Energy [e-e_HOMO](eV)    #GS-PDOS(e) '
    do k=1,np
       write(11,fmt_np) (e(k)-ehomo)*au2ev,gs_pdos(k,j)
    end do
 end do
 close(11)


 !lookhere
 !Ground-state OPDOS
 !gs_opdos=0.d0
 !do i=1,nfragOP
 !  do j=i+1,nfragOP
 !    do k=1,np
 !      do ii=1,nocc
 !         gs_opdos(k,i,j)=gs_opdos(k,i,j)+wOP(i,j,ii)*fc/((e(k)-emo(ii))**2+sigma**2)
 !      enddo
 !    enddo
 !  enddo
 !enddo
 !gs_opdos=2.d0*gs_opdos

! if (maxval(abs(gs_opdos(:,:,:))).gt.0.d0) then
!    gs_opdos(:,:,:)=gs_opdos(:,:,:)/maxval(abs(gs_opdos(:,:,:)))
! endif
 if (average.eq."yes")then
          write (fmt_np,'("(2f14.4,",I0,"2e21.12E3)")') np
 else
          write (fmt_np,'("(2f14.4,",I0,"e21.12E3)")') np
 endif
 basename="gs-opdos_"
 do i=1,nfragOP
   do j=i+1,nfragOP
     write(sfr,'(i5.5)') i   !vengono # da 000001 a 100000
     write(sfr2,'(i5.5)') j  !vengono # da 000001 a 100000
     write(filename,*) adjustl(trim(basename)), adjustl(trim(sfr)),&
                      &adjustl(trim(sfr2)) ,".dat"
     inquire(file=filename,exist=lex)
     if (lex) then
       open(11, file=filename, status="old",position="append",action="write")
     else
       open(11, file=filename, status="new", action="write")
     end if
         write(11,*) '#Energy [e-e_HOMO](eV)    #GS-OPDOS(e) '
     do k=1,np
        write(11,fmt_np) (e(k)-ehomo)*au2ev,gs_opdos(k,i,j)
     enddo
   enddo
 enddo
 close(11)
 !open(21,file='integ_td-tot-pdos.dat')
 !--------------------------------------------------------!

 !write(20,*) '# time (fs) frag1->frag2  frag2->frag1'
 !for each time compute and write delta pdos + delta opdos
 do i=1,neff
    !Temporary array for  sum_lc_l(t)*A_ci^l
    zc = (0.0d0,0.0d0)
    do l=1,nexc
       !EC istart=1 or istart=nf+1
       do ivir = 1, nvir
          do iocc = istart, nocc
             zc(iocc,ivir)=zc(iocc,ivir) + c(l+1,i) * a(iocc,ivir,l)  
          end do
       end do
    end do
    if (average.eq."yes") then
       zc_avg = (0.0d0,0.0d0)
       do l=1,nexc
          do ll=1,nexc
             !EC istart=1 or istart=nf+1
             do iocc = istart, nocc
                do ivir = 1, nvir
                   zc_avg(iocc)=zc_avg(iocc)+corr(i,l,ll)*a(iocc,ivir,l)*a(iocc,ivir,ll)
                end do
             end do
             do ivir = 1, nvir
                !EC istart=1 or istart=nf+1
                do iocc = istart, nocc
                   zc_avg(nocc+ivir)=zc_avg(nocc+ivir)+corr(i,l,ll)*a(iocc,ivir,l)*a(iocc,ivir,ll)
                end do
             end do
          end do
       end do
    endif
    z = 0.0d0
    dm(1:nocc) = 1.0d0
    !EC istart=1 or istart=nf+1
    do iocc = istart, nocc
       do ivir = 1, nvir
          z(iocc) = z(iocc) + abs(zc(iocc,ivir))**2
          dm(iocc) = dm(iocc) - abs(zc(iocc,ivir))**2
       end do
    end do
    dm(nocc+1:nmo) = 0.0d0
    do ivir = 1, nvir
       !EC istart=1 or istart=nf+1
       do iocc = istart, nocc
          z(nocc+ivir) = z(nocc+ivir) + abs(zc(iocc,ivir))**2
          dm(nocc+ivir) = dm(nocc+ivir) + abs(zc(iocc,ivir))**2
       end do
    end do
    if (average.eq."yes") then
       z_avg = z + real(zc_avg)
       dm_avg(1:nocc) = dm(1:nocc) - real(zc_avg)
       dm_avg(nocc+1:nmo) = dm(nocc+1:nmo) + real(zc_avg)
    endif

    pdos=0.d0
    !integ=0.d0
    do k=1,np
       do j=1,nfrag
          !EC istart=1 or istart=nf+1
          do ii=istart,nocc
             pdos(k,j)=pdos(k,j) - w(j,ii)*z(ii)*fc/((e(k)-emo(ii))**2+sigma**2) 
             !integ(j)=integ(j)-w(j,ii)*z(ii)*fc/((e(k)-emo(ii))**2+sigma**2)
          end do
          do aa=1,nvir
             pdos(k,j)=pdos(k,j)+w(j,nocc+aa)*z(nocc+aa)*fc/((e(k)-emo(nocc+aa))**2+sigma**2) 
             !integ(j)=integ(j)+w(j,nocc+aa)*z(nocc+aa)*fc/((e(k)-emo(nocc+aa))**2+sigma**2)
          end do
       end do
    end do

    if (average.eq."yes")then
       pdos_avg=0.d0
       !integ_avg=0.d0
       do k=1,np
          do j=1,nfrag
             !EC istart=1 or istart=nf+1
             do ii=istart,nocc
                pdos_avg(k,j)=pdos(k,j) - w(j,ii)*z_avg(ii)*fc/((e(k)-emo(ii))**2+sigma**2)
                !integ_avg(j)=integ(j)-w(j,ii)*z_avg(ii)!*fc/((e(k)-emo(ii))**2+sigma**2)
             end do
             do aa=1,nvir
                pdos_avg(k,j)=pdos(k,j)+w(j,nocc+aa)*z_avg(nocc+aa)*fc/((e(k)-emo(nocc+aa))**2+sigma**2)
                !integ_avg(j)=integ(j)+w(j,nocc+aa)*z_avg(nocc+aa)!*fc/((e(k)-emo(nocc+aa))**2+sigma**2)
             enddo
          enddo
       enddo
    endif

    if (mixing.eq.'y') then
            do j=1,nfrag
                 do k=1,np
                      pdos(k,j)=pdos(k,j)-gs_pdost(k,j)
                 enddo
            enddo
    endif
    pdos=2.d0*pdos
    dm=2.d0*dm 
    
    if (average.eq.'yes') then
       if (mixing.eq.'y') then
            do j=1,nfrag
                 do k=1,np
                      pdos_avg(k,j)=pdos_avg(k,j)-gs_pdost(k,j)
                 enddo
            enddo
       endif
       pdos_avg=2.d0*pdos_avg
       dm_avg=2.d0*dm_avg
    endif 

    do j=1,nfrag
       ksum(j)=0.d0  
       do k=1,np
          ksum(j) = ksum(j) + pdos(k,j)
       enddo 
       ksum(j)=ksum(j)*h
    enddo
    
    !if (maxval(abs(pdos(:,:))).gt.0.d0) then
    !   pdos(:,:)=pdos(:,:)/maxval(abs(pdos(:,:))) 
    !endif 

    basename="td-pdos_frag"
    basename1="integ_td-pdos_frag"
    basename2="integ_td-pdos_eh_frag"
    filename3="integ_td-tot-pdos.dat"
    do j = 1, nfrag
       write(sfr,'(i5.5)') j   !vengono # da 000001 a 100000
       write(filename,*) adjustl(trim(basename)), adjustl(trim(sfr)),".dat"
       write(filename1,*) adjustl(trim(basename1)), adjustl(trim(sfr)),".dat"
       write(filename2,*) adjustl(trim(basename2)), adjustl(trim(sfr)),".dat"
       inquire(file=filename,exist=lex)
       if (lex) then
          open(11, file=filename, status="old", position="append",action="write")
       else
          open(11, file=filename, status="new", action="write")
       end if
       inquire(file=filename1,exist=lex)
       if (lex) then
          open(20, file=filename1,status="old",position="append",action="write")
       else
          open(20, file=filename1, status="new", action="write")
       endif
       inquire(file=filename2,exist=lex)
       if (lex) then
          open(22, file=filename2,status="old",position="append",action="write")
       else
          open(22, file=filename2, status="new", action="write")
       endif
       if (average.eq."yes")then
           write(11,*) '#Time(fs) Energy [e-e_HOMO](eV)    #PDOS_CORR(t,e)   #PDOS(t,e)'
           if (i.eq.1) then
              write(20,*) '#Time(fs)  \int_e #PDOS_CORR(t,e) \int_e #PDOS(t,e)      El inj      Ho inj'
              write(22,*) '#Time(fs)  \int_-infty^0 #PDOS_CORR(t,e) \int_0^\infty #PDOS_CORR(t,e)'
           endif
       else
           write(11,*) '#Time(fs) Energy [e-e_HOMO](eV)    #PDOS(t,e) '
           if (i.eq.1) then
              write(20,*) '#Time(fs)  \int_e #PDOS(t,e)      El inj      Ho inj'
              write(22,*) '#Time(fs)  \int_-infty^0 #PDOS(t,e) \int_0^\infty #PDOS(t,e)'
           endif
       endif
       write(11,*) "#Time step", ist(i), t(i)*0.0241888,"fs"
       do k=1,np
          !take care of indices in pdos
          if (average.eq."yes")then
              write(11,fmt_np) t(i)*0.0241888,(e(k)-ehomo)*au2ev,pdos_avg(k,j), pdos(k,j)
          else
              write(11,fmt_np) t(i)*0.0241888,(e(k)-ehomo)*au2ev,pdos(k,j)
          endif
       end do
       write (11,*) ""
       write (11,*) ""
       ! Integrated DeltaPDOS for each fragment (integ or integ_avg)
       !From -infty to +infty
       integ(j)=0.d0
       ee(j)=0.d0
       hh(j)=0.d0
       do k=2,np-1
          integ(j)=integ(j)+pdos(k,j)
          ee(j)=ee(j)+0.5d0*(pdos(k,j)+abs(pdos(k,j)))
          hh(j)=hh(j)+0.5d0*(pdos(k,j)-abs(pdos(k,j)))
       enddo
       integ(j)=integ(j)+0.5d0*(pdos(1,j)+pdos(np,j))
       integ(j)=integ(j)*h
       ee(j)=(ee(j)+0.25d0*(pdos(1,j)+pdos(np,j)+abs(pdos(1,j))+abs(pdos(np,j))))*h
       hh(j)=(hh(j)+0.25d0*(pdos(1,j)+pdos(np,j)-abs(pdos(1,j))-abs(pdos(np,j))))*h
       !ee(j)=0.5d0*(integ(j)+abs(integ(j)))
       !hh(j)=0.5d0*(integ(j)-abs(integ(j)))
       !From -infty to 0
       einteg(j)=0.d0
       do k=2,npo-1
          einteg(j)=einteg(j)+pdos(k,j)
       enddo
       einteg(j)=einteg(j)+0.5d0*(pdos(1,j)+pdos(npo,j))
       einteg(j)=einteg(j)*h
       !From 0 to +infty
       hinteg(j)=0.d0
       do k=npo+2,np-1
          hinteg(j)=hinteg(j)+pdos(k,j)
       enddo
       hinteg(j)=hinteg(j)+0.5d0*(pdos(npo+1,j)+pdos(np,j))
       hinteg(j)=hinteg(j)*h
       if (average.eq."yes") then
          !From -infty to +infty
          integ_avg(j)=0.d0
          ee_avg(j)=0.d0
          hh_avg(j)=0.d0 
          do k=2,np-1
             integ_avg(j)=integ_avg(j)+pdos_avg(k,j)
             ee_avg(j)=ee_avg(j)+0.5d0*(pdos_avg(k,j)+abs(pdos_avg(k,j)))
             hh_avg(j)=hh_avg(j)+0.5d0*(pdos_avg(k,j)-abs(pdos_avg(k,j)))
          enddo
          integ_avg(j)=integ_avg(j)+0.5d0*(pdos_avg(1,j)+pdos_avg(np,j))
          integ_avg(j)=integ_avg(j)*h
          ee_avg(j)=(ee_avg(j)+0.25d0*(pdos_avg(1,j)+pdos_avg(np,j)&
                    +abs(pdos_avg(1,j))+abs(pdos_avg(np,j))))*h
          hh_avg(j)=(hh_avg(j)+0.25d0*(pdos_avg(1,j)+pdos_avg(np,j)&
                    -abs(pdos_avg(1,j))-abs(pdos_avg(np,j))))*h
!MM
!          ee_avg(j)=0.5d0*(integ_avg(j)+abs(integ_avg(j)))
!          hh_avg(j)=0.5d0*(integ_avg(j)-abs(integ_avg(j)))
          !From -infty to 0
          einteg_avg(j)=0.d0
          do k=2,npo-1
             einteg_avg(j)=einteg_avg(j)+pdos_avg(k,j)
          enddo
          einteg_avg(j)=einteg_avg(j)+0.5d0*(pdos_avg(1,j)+pdos_avg(npo,j))
          einteg_avg(j)=einteg_avg(j)*h
          !From 0 to +infty 
          hinteg_avg(j)=0.d0
          do k=npo+2,np-1
             hinteg_avg(j)=hinteg_avg(j)+pdos_avg(k,j)
          enddo
          hinteg_avg(j)=hinteg_avg(j)+0.5d0*(pdos_avg(npo+1,j)+pdos_avg(np,j))
          hinteg_avg(j)=hinteg_avg(j)*h
          write(20,112) t(i)*0.0241888, integ_avg(j), integ(j), ee_avg(j),hh_avg(j)
          write(22,111) t(i)*0.0241888, einteg_avg(j), hinteg_avg(j)
       else
          write(20,113) t(i)*0.0241888, integ(j),ee(j),hh(j)
          write(22,111) t(i)*0.0241888, einteg(j),hinteg(j)
       endif

    end do
    inquire(file=filename3,exist=lex)
    if (lex) then
       open(23, file=filename3,status="old",position="append",action="write")
    else
       open(23, file=filename3, status="new", action="write")
    end if
    if (average.eq."yes")then
       if (i.eq.1) then
          write(23,*) '#Time(fs)  \int_-infty^0 #totPDOS_CORR(t,e) \int_0^\infty #totPDOS_CORR(t,e)'
       endif
       write(23,111) t(i)*0.0241888, sum(einteg_avg), sum(hinteg_avg)
    else
       if (i.eq.1) then
          write(23,*) '#Time(fs)  \int_-infty^0 #totPDOS(t,e) \int_0^\infty #totPDOS(t,e)'
       endif
       write(23,111) t(i)*0.0241888, sum(einteg),sum(hinteg)
    endif



    close(11)
    close(20)
    close(22)
    close(23)
    !write(21,110)  t(i)*0.0241888, sum(integ(:))
    

    !lookhere
    !opdos=0.d0
    !do k=1,np
    !  do i2=1,nfragOP
    !    do j=i2+1,nfragOP
    !      do ii=1,nocc
    !        opdos(k,i2,j)=opdos(k,i2,j) - wOP(i2,j,ii)*z(ii)*fc/((e(k)-emo(ii))**2+sigma**2)
    !      end do
    !      do aa=1,nvir
    !        opdos(k,i2,j)=opdos(k,i2,j)+wOP(i2,j,nocc+aa)*z(nocc+aa)*fc/((e(k)-emo(nocc+aa))**2+sigma**2)
    !      enddo
    !    enddo
    !  enddo
    !enddo

    !opdos=2.d0*opdos

    !Print delta td-opdos
    !basename="td-opdos_"
    !do i2=1,nfragOP
    !  do j=i2+1,nfragOP
    !    write(sfr,'(i5.5)') i2   !vengono # da 000001 a 100000
    !    write(sfr2,'(i5.5)') j   !vengono # da 000001 a 100000
    !    write(filename,*) adjustl(trim(basename)), adjustl(trim(sfr)),&
    !                      &adjustl(trim(sfr2)),".dat"
    !    inquire(file=filename,exist=lex)
    !    if (lex) then
    !       open(11, file=filename, status="old",position="append",action="write")
    !    else
    !       open(11, file=filename, status="new", action="write")
    !    end if
    !    write(11,*) '#Energy [e-e_HOMO](eV)    #OPDOS(t,e) '
    !    write(11,*) "#Time step", ist(i), t(i)*0.0241888,"fs"
    !    do k=1,np
    !       !take care of indices in opdos
    !       write(11,fmt_np) (e(k)-ehomo)*au2ev,opdos(k,i2,j)
    !    enddo
    !  write(11,*) ""
    !  write(11,*) ""
    !  enddo
    !enddo
    !close(11)    
    !------------------------------------------------------!

    basename="td-occ_mos"
    write(filename,*) adjustl(trim(basename)), ".dat"
    inquire(file=filename,exist=lex)
    if (lex) then
       open(11, file=filename, status="old", position="append", action="write")
    else
       open(11, file=filename, status="new", action="write")
    end if
    write(11,*) "#Time step", ist(i), t(i)*0.0241888,"fs"
    do k = 1, nmo
       if (average.eq."yes") then
               write(11,"(f14.4,2x,i4, 2x, 2f17.8)") t(i)*0.0241888,k, dm_avg(k), dm(k)
       else
               write(11,"(f14.4,2x,i4, 2x, f17.8)") t(i)*0.0241888,k, dm(k) 
       endif
    end do
    write(11,*) ''
    write(11,*) '' 
    close(11)
 end do

 !Total DOS
 tot_dos=0.d0
 do k=1,np
    !EC istart=1 or istart=nf+1
    do ii=istart,nmo
       tot_dos(k)=tot_dos(k)+fc/((e(k)-emo(ii))**2+sigma**2)
    end do
 enddo

 if (maxval(abs(tot_dos(:))).gt.0.d0) then
    tot_dos(:)=tot_dos(:)/maxval(abs(tot_dos(:)))
 endif

 open(11,file='tot_dos.dat')
 write (fmt_np,'("(f14.4,",I0,"e17.8E3)")') np
 write(11,*) '#  Energy [e-e_HOMO](eV)    #Tot-DOS(e) '
 do k=1,np
    write(11,fmt_np) (e(k)-ehomo)*au2ev,tot_dos(k)
 end do
 close(11)

 !lookhere
 !PDOS with occ and vir orbitals
 !ov_pdos=0.d0
 !do j=1,nfrag
 !   do k=1,np
 !      do ii=1,nmo
 !         ov_pdos(k,j)=ov_pdos(k,j)+w(j,ii)*fc/((e(k)-emo(ii))**2+sigma**2)
 !      end do
 !   enddo
 !enddo

 !if (maxval(abs(ov_pdos(:,:))).gt.0.d0) then
 !   ov_pdos(:,:)=ov_pdos(:,:)/maxval(abs(ov_pdos(:,:)))
 !endif

 !basename="ov-pdos_frag"
 !do j = 1, nfrag
 !   write(sfr,'(i5.5)') j   !vengono # da 000001 a 100000
 !   write(filename,*) adjustl(trim(basename)), adjustl(trim(sfr)),".dat"
 !   inquire(file=filename,exist=lex)
 !   if (lex) then
 !      open(11, file=filename, status="old",position="append",action="write")
 !    else
 !      open(11, file=filename, status="new", action="write")
 !   end if
 !   write(11,*) '#Energy [e-e_HOMO](eV)    #ov-PDOS(e) '
 !   do k=1,np
 !      write(11,fmt_np) (e(k)-ehomo)*au2ev,ov_pdos(k,j)
 !   end do
 !end do
 !close(11)
 !close(21)


 deallocate(z,zc)
 deallocate(e,integ,einteg,hinteg)
 deallocate(ee,hh)
 if (mixing.eq."y") then
         deallocate(gs_pdost)
         deallocate(zz)
 endif
 if (average.eq."yes")then
    deallocate(zc_avg)
    deallocate(z_avg)
    deallocate(dm_avg)
    deallocate(integ_avg)
    deallocate(ee_avg)
    deallocate(hh_avg)
    deallocate(einteg_avg)
    deallocate(hinteg_avg)
    deallocate(pdos_avg)
 endif
100 format(F10.3,I6,F12.5)
110 format(f14.4,e17.8E3)
111 format(f14.4,2e17.8E3)
112 format(f14.4,4e17.8E3)
113 format(f14.4,3e17.8E3)
150 format(f14.4,e17.8E3,e17.8E3,e17.8E3)
300 format(f14.4,e17.8E3,e17.8E3,e17.8E3,e17.8E3)

 return

end subroutine compute_td_pdos


subroutine compute_td_tcm()
!------------------------------------------------------------------------
! @brief Compute td-TCM 
!
! @date Created   : E. Coccia 3/9/20 
! Modified  : 
!------------------------------------------------------------------------

 implicit none

 integer             :: ii,aa,l,i,j,k,p,npo,npv,kk,m,istart

 real*8              :: h,atmp
 real*8, allocatable :: e(:)

 logical             :: lex

 character*5         :: sfr
 character*30        :: filename,basename
 character*3000      :: fmt_np

 complex(8), allocatable :: tmp(:,:)

 allocate(e(np))
 allocate(tmp(nocc,nvir))

 h=(emax-emin)/dble(np)
 do k=1,np
    e(k)=emin+k*h
 enddo

 npo=0
 k=1
 do while (e(k).lt.ehomo)
    k=k+1
    npo=npo+1
 enddo

 if (frozen.eq.'yes') then
    istart=nf+1
 else
    istart=1
 endif

 allocate(tcm(npo,np-npo))
 open(11,file='td-tcm_tot.dat')
 if (nfrag.eq.2) then
    allocate(ftcm(npo,np-npo,nfrag**2))
    write (fmt_np,'("(f14.4,",I0,"f14.4,",I0,"e17.8E3)")') npo*(np-npo)
    basename="td-tcm_frag"
    write(*,*) ''
    write(*,*) 'TCM frag1 output: occ frag1 -> vir frag1'
    write(*,*) 'TCM frag2 output: occ frag1 -> vir frag2'
    write(*,*) 'TCM frag3 output: occ frag2 -> vir frag1'
    write(*,*) 'TCM frag4 output: occ frag2 -> vir frag2'
    write(*,*) ''
 endif

 do i=1,neff
    tcm=0.d0
    tmp=(0.d0,0.d0)
    if (average.ne.'yes') then
       !EC istart=1 or istart=nf+1
       do ii=istart,nocc
          do aa=1,nvir
             !do l=1,nexc
             !   do p=1,nexc
             !      tmp(ii,aa)=tmp(ii,aa)+conjg(c(l+1,i))*c(p+1,i)*a(ii,aa,l)*a(ii,aa,p)
             !   enddo
             !enddo
             do l=1,nexc
                do p=l+1,nexc
                   atmp=a(ii,aa,l)*a(ii,aa,p) 
                   tmp(ii,aa)=tmp(ii,aa)+conjg(c(l+1,i))*c(p+1,i)*atmp
                   tmp(ii,aa)=tmp(ii,aa)+conjg(c(p+1,i))*c(l+1,i)*atmp
                enddo
                tmp(ii,aa)=tmp(ii,aa)+conjg(c(l+1,i))*c(l+1,i)*a(ii,aa,l)**2
             enddo
          enddo
       enddo
    else
       !EC istart=1 or istart=nf+1
       do ii=istart,nocc
          do aa=1,nvir
             !do l=1,nexc
             !   do p=1,nexc
             !      tmp(ii,aa)=tmp(ii,aa)+(conjg(c(l+1,i))*c(p+1,i)+corr(i,l,p))*a(ii,aa,l)*a(ii,aa,p)
             !   enddo
             !enddo
             do l=1,nexc
                do p=l+1,nexc
                   atmp=a(ii,aa,l)*a(ii,aa,p)                                    
                   tmp(ii,aa)=tmp(ii,aa)+(conjg(c(l+1,i))*c(p+1,i)+corr(i,l,p))*atmp
                   tmp(ii,aa)=tmp(ii,aa)+(conjg(c(p+1,i))*c(l+1,i)+corr(i,p,l))*atmp
                enddo
                tmp(ii,aa)=tmp(ii,aa)+(conjg(c(l+1,i))*c(l+1,i)+corr(i,l,l))*a(ii,aa,l)**2
             enddo
          enddo
       enddo
    endif

    do k=1,npo
       do l=npo+1,np
          !EC istart=1 or istart=nf+1
          do ii=istart,nocc
             do aa=1,nvir
                tcm(k,l-npo)=tcm(k,l-npo) + tmp(ii,aa)*exp(-(e(k)-emo(ii))**2/sigma-(e(l)-emo(nocc+aa))**2/sigma)
             enddo
          enddo
       enddo
    enddo
    tcm=2.d0*tcm

    if (nfrag.eq.2) then
       ftcm=0.d0
       do k=1,npo
          do l=npo+1,np
             kk=0
             do j=1,nfrag
                   do m=1,nfrag
                      kk=kk+1
                       !EC istart=1 or istart=nf+1 
                      do ii=istart,nocc
                          do aa=1,nvir
                          ftcm(k,l-npo,kk)=ftcm(k,l-npo,kk) + &
                & +w(j,ii)*w(m,nocc+aa)*tmp(ii,aa)*exp(-(e(k)-emo(ii))**2/sigma-(e(l)-emo(nocc+aa))**2/sigma)
                          enddo
                      enddo
                   enddo
             enddo
          enddo
       enddo
       ftcm=2.d0*ftcm
    endif

    write(11,*) '#Occ oeps (eV) Vir veps (eV)    #TCM(t,oeps,veps) '
    write(11,*) "#Time step", ist(i), t(i)*0.0241888,"fs"
    do k=1,npo
       do l=npo+1,np
          write(11,160) (e(k)-ehomo)*au2ev,(e(l)-ehomo)*au2ev,tcm(k,l-npo)
       enddo
       write(11,'(2x)')
    enddo
    write(11,'(2x)')

    if (nfrag.eq.2) then

       kk=0
       do j=1,nfrag
          do m=1,nfrag
             kk=kk+1
             !write(sfr,'(i1.1,A,i5.5)') j,'-',m   !vengono # da 000001 a 100000
             write(sfr,'(i5.5)') kk
             write(filename,*) adjustl(trim(basename)),adjustl(trim(sfr)),".dat"
             inquire(file=filename,exist=lex)
             if (lex) then
                open(12, file=filename,status="old",position="append",action="write")
             else
                open(12, file=filename, status="new", action="write")
             endif
             write(12,*) '#Occ oeps (eV) Vir veps (eV)    #TCM(t,oeps,veps)'
             write(12,*) "#Time step", ist(i), t(i)*0.0241888,"fs"
             do k=1,npo
                do l=npo+1,np
                   write(12,160)(e(k)-ehomo)*au2ev,(e(l)-ehomo)*au2ev,ftcm(k,l-npo,kk)
                enddo
                write(12,'(2x)')
             enddo
             write(12,'(2x)')
          enddo
       enddo
    endif 
 enddo 
 close(11)
 close(12)

 deallocate(tmp)
 deallocate(e)
 deallocate(tcm)
 if (nfrag.eq.2) deallocate(ftcm)

160 format(f14.4,f14.4,e17.8E3)  
 return

end subroutine compute_td_tcm


!------------------------------------------------------------------------
! @brief Deallocate arrays 
!
! @date Created   : E. Coccia 22/6/20 
! Modified  : 
!------------------------------------------------------------------------
subroutine deallocate_all()

 implicit none

 deallocate(t)
 deallocate(ist)
 deallocate(a) 
 deallocate(w)
 deallocate(c)
 deallocate(pdos)
 deallocate(emo)
 deallocate(frag)
 deallocate(dm)
 deallocate(tot_dos)
 deallocate(gs_pdos)
 deallocate(ov_pdos)
 deallocate(ksum)

end subroutine deallocate_all

end module mod_td 

