module vib
      use constants
#ifdef OMP
      use omp_lib
#endif

      implicit none
      save

      integer(i4b)              :: nstates                                               ! number of electronic states
      integer(i4b)              :: nvib                                                  ! number of vibrational quanta on ES used
      integer(i4b)              :: nmodes_gs                                             ! number of modes in the ground state
      integer(i4b)              :: nmodes                                                ! number of modes in the excited state
      integer(i4b)              :: ntot                                                  ! number of total vibronic state
      integer(i4b)              :: ncomb                                                 ! number of combination
      integer(i4b)              :: nfc                                                   ! Power for <v'|x^nfc|v> 
      integer(i4b)              :: nbin                                                  ! write binary files  
      integer(i4b)              :: nvib_inp                                              ! number of vibrational quanta on ES read
      integer(i4b)              :: nts                                                   ! number of tessere on NP surface
      real(dbl)                 :: sigma,emin,emax                                       ! properties of vibrational spectrum
      real(dbl),    allocatable :: w(:,:),q(:,:),mu(:,:)                                 ! frequency, dipoles read 
      real(dbl),    allocatable :: fcint(:,:),htint(:,:,:),dispint(:,:)                  ! fc, ht, and displacement between GS(n=0) and ES(n>0)
      real(dbl),    allocatable :: htintr(:,:,:,:),fcintr(:,:,:), dispintr(:,:,:)        ! fc, ht, and displacement between GS(n>0) and ES(n>0)
      real(dbl),    allocatable :: e(:),dip(:,:,:)                                       ! energies and dipoles read
      real(dbl),    allocatable :: ef(:),dipf(:,:,:),integ(:,:),pot_diff(:,:)            ! final energies and dipoles 
      real(dbl),    allocatable :: quantum_vts(:,:,:), quantum_vtsn(:), vts(:,:,:)       ! potential on NP surface
      integer(i4b), allocatable :: map(:,:)                                              ! map between vib level of GS and ES
      integer(i4b), allocatable :: iv(:,:)                                               ! map between vib level of GS and ES (use in computefc) 
      logical                   :: coupling,mix                                          ! coupling and mixing not implemented
      logical                   :: computefc                                             ! .true. to compute fc or .false. to read from  files fc.inp and/or ht.inp
      character(20)             :: diptype                                               ! 'fc', 'ht', or 'fcht' to read only fc, ht or both 
      character(20)             :: calctype                                              ! 'rayleigh' (compute  and read only GS(0)->ES(n)) terms or 'raman' (compute and read also GS(n)->ES(n)) calculations
      character(20)             :: medium                                                ! 'nan' or 'vac', for nanoparticle (compute ci_pot_new) or in vacuum calculations

      public nstates,nvib,nmodes,nmodes_gs,w,q,e,dip,ef,dipf,fact,dfact,bin_coef, &
             mix,coupling,iv,ncomb,nfc,sigma,nbin,mu
       
      contains        

!------------------------------------------------------------------------
! @brief Read ci_energy.inp and ci_mut.inp from Gamess/Gaussian 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
      subroutine read_e_dip()  

       integer(i4b) :: i,j
       character(4) :: junk

!    read ci energy
       open(7,file="ci_energy.inp",status="old")
       ! add ground state
       write(*,*)
       write(*,*) 'Pure electronic energies:'
       allocate (e(nstates))
       e(1)=0.d0
       do i=2,nstates
         read(7,*) junk,junk,junk,e(i)
         e(i)=e(i)*ev_to_au
       enddo
       write (6,*) (e(i),i=2,nstates)
       close(7)

!    read transition dipoles (also for pcm: useful for analysis)
       open(7,file="ci_mut.inp",status="old")
       allocate (dip(3,nstates,nstates))
       do i=1,nstates
          read(7,*)junk,junk,junk,junk,dip(1,1,i),dip(2,1,i),dip(3,1,i)
          dip(:,i,1)=dip(:,1,i)
       enddo

       do i=2,nstates
          do j=2,i
             read(7,*)junk,junk,junk,junk,dip(1,i,j),dip(2,i,j),dip(3,i,j)
             dip(:,j,i)=dip(:,i,j)
          enddo
       enddo
       close(7)

       return

      end subroutine read_e_dip
    
!------------------------------------------------------------------------
! @brief Read input for vibrational corrections 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     subroutine read_input_vib() 
 
        integer(i4b)             :: idum,i,j,k,h,n_out,dime,kk
        character(20)            :: filein
        character(10)            :: idum_c
        real                     :: rdum


 
        namelist/vibrations/nstates,nvib,nvib_inp,nmodes,coupling,mix,nfc,sigma,nbin,emax,computefc,diptype,&
                      &calctype,medium,nmodes_gs

        mix=.false.
        coupling=.false.
        nmodes=1
        nmodes_gs=0
        nvib=10
        nvib_inp=0
        n_out=5
        nstates=1
        nfc=1
        sigma=3.d0
        nbin=10000
        emin=0.d0
        emax=15.d0 !eV
        computefc=.false.
        diptype='fc'
        calctype='rayleigh'
        medium='vacuum'

        ! Read w, q and mu  for any vib level (vib.dat file)
        ! Frequency in cm-1, normal coordinates in bohr, reduced mass in amu 
        ! Excited state N
        ! w q mu for mode 1
        ! w q mu for mode 2
        ! ...
        ! Excited state N+1
        ! w q mu for mode 1
        ! w q mu for mode 2
        ! ...

        read(*,nml=vibrations)
        if (nvib.gt.nvibmax) then
           nvib=nvibmax 
        endif 
        
        if (.not.(computefc).and.(nstates.gt.2)) then
                write(*,*) "ERROR: case with nstates>2 not implemented yet"
                stop   
        endif
        ncomb=nvib**nmodes
        if (nmodes_gs.eq.0) nmodes_gs = nmodes
        if (nvib_inp.eq.0) nvib_inp = nvib

        allocate(q(nstates,nmodes),mu(nstates,nmodes))
        if (computefc) allocate(iv(ncomb,nmodes))
        allocate(fcint(nmodes,nvib_inp))
        allocate(htint(nmodes,nvib_inp,3))
        allocate(dispint(nmodes,nvib_inp))
        if (calctype.eq."raman") then
                allocate(fcintr(nmodes_gs,nmodes,nvib_inp))
                allocate(htintr(nmodes_gs,nmodes,nvib_inp,3))
                allocate(dispintr(nmodes_gs,nmodes,nvib_inp))
        endif

        write(*,*)  
        write(*,*) '***************************************'
        write(*,*) '*                                     *'
        write(*,*) '*      Add vibrational levels         *' 
        write(*,*) '*     in harmonic approximation       *'
        write(*,*) '*      to any normal mode of          *'
        write(*,*) '*     a given electronic state        *'
        write(*,*) '*                                     *'      
        write(*,*) '*                by                   *'
        write(*,*) '*          Emanuele Coccia            *'
        write(*,*) '*                                     *'
        write(*,*) '***************************************'

        if (computefc) then
                allocate(w(nstates,nmodes))
                ntot=nstates*ncomb
                w(:,:) = 1000.d0
                q(:,:) = 2.d0
                mu(:,:) = 1836.d0
                open(60,file='vib.dat')
                do i=1,nstates
                       read(60,*) idum, idum
                       do j=1,nmodes
                              read(60,*) w(i,j), q(i,j), mu(i,j)
                       enddo
                enddo
                close(60)
                w(:,:) = w(:,:)*cm_to_au
                mu(:,:) = mu(:,:)*amu_to_au

        else
                allocate(w(nstates,nmodes))
                ntot=nmodes*(nvib-1)+nmodes_gs+2
                w(:,:) = 0.d0
                open(12, file="frequency.dat", status="old")
                fcint(:,:)=0
                htint(:,:,:)=0
                if (calctype.eq."raman") then
                        open(13, file="map.inp", status="old")
                        allocate(map(nmodes,nmodes_gs+1))
                        do i=1,nmodes_gs+1
                             read(13,*) 
                             do j=1,nmodes
                                 read(13,*) map(j,i)
                             enddo
                        enddo
                        close(13)
                else
                        allocate(map(nmodes,1))
                        do i=1,nmodes
                           map(i,1)=i
                        enddo
                endif
                
                open(11, file="fc.inp", status="old")
                read(11,*)
                if (diptype.eq."ht".or.diptype.eq."fcht") then
                        open(13, file="ht.inp", status="old")
                        read(13,*)
                        htint(:,:,:)=0.0
                        if (calctype.eq."raman") htintr(:,:,:,:)=0.0
                        if (medium.eq.'nan') then
                            open(14, file="ht_disp.inp", status="old")
                            read(14,*)
                            dispint=0.0
                            if (calctype.eq."raman") dispintr=0.0
                        endif
                endif
                do i=1,nmodes
                     read(11,*)
                     read(11,*)   idum_c, (fcint(i,k),k=1,nvib_inp)
                     read(11,*)
                     if (diptype.eq."ht".or.diptype.eq."fcht") then
                             read(13,*)
                             do j=2,nvib_inp
                                read(13,*)
                                read(13,*)   idum_c, htint(i,j,:)
                             enddo
                             if (medium.eq.'nan') then
                                read(14,*)
                                do j=2,nvib_inp
                                    read(14,*)   idum, idum, dispint(i,j)
                                enddo
                             endif
                     endif
                enddo
                
                do i=1,nmodes_gs
                    read(12,*) w(1,i)
                    w(1,i)=w(1,i)*cm_to_au
                enddo
                do i=1,nmodes
                    read(12,*) w(2,i)
                    w(2,i)=w(2,i)*cm_to_au
                enddo
                if (calctype.eq."raman") then
                   do kk=1,nmodes_gs
                      read(11,*)
                      do i=1,nmodes
                        read(11,*)
                        read(11,*)   idum_c, (fcintr(kk,map(i,kk+1),k),k=1,nvib_inp)
                        read(11,*)
                      enddo
                   enddo
                   if (diptype.eq."ht".or.diptype.eq."fcht") then
                        do kk=1,nmodes_gs
                           read(13,*)
                           read(13,*)
                           read(13,*)
                           if (medium.eq.'nan') read(14,*)
                           do i=1,nmodes
                               read(13,*)
                               do j=1,nvib_inp-1
                                   read(13,*)
                                   read(13,*)
                               enddo
                               do j=1,nvib_inp
                                  read(13,*)   idum_c, idum_c, idum,  k
                                  read(13,*)   idum_c, htintr(kk,map(i,kk+1),k,:)
                               enddo
                               if (medium.eq.'nan') then
                                    read(14,*)
                                    do j=2,nvib_inp
                                       read(14,*)
                                    enddo
                                    do j=2,nvib_inp
                                       read(14,*)   idum, k, dispintr(kk,map(i,kk+1),k)
                                    enddo
                               endif
                           enddo
                        enddo
                   endif
                endif
                close(11)
                close(12)
                close(13)
                close(14)
                if (diptype=="fcht") close(13)
        endif

        if (medium.eq.'nan') then
                call read_gau_out_medium(nstates)
                if(diptype.eq."ht".or.diptype.eq."fcht") then
                        allocate(pot_diff(nts,nmodes))
                        open(15, file="pot_diff_total.dat", status="old")
                        do i=1,nts
                            read(15,*) pot_diff(i,:)
                        enddo
                        close(15)
                endif
        endif

        write(*,*) ''
        write(*,*) 'Total number of states', ntot
        write(*,*) 'Number of electronic states', nstates
        write(*,*) 'Number of normal modes per electronic state', nmodes
        write(*,*) 'Number of vibrational states per normal mode', nvib
        if (mix) then
           write(*,*) 'Duschinsky rotation for normal coordinates' 
        else
           write(*,*) 'Only displacements between normal coordinates'
        endif 
        write(*,*) ''
        write(*,*) 'Vibronic spectrum:'
        write(*,*) '   Energy range (eV):', emin, emax
        write(*,*) '   Number of bins', nbin
        write(*,*) '   Gaussian sigma', sigma
        write(*,*) ''
        write(*,*) 'Maximum number of vibrational states:', nvibmax
        write(*,*) ''
        if (medium.eq.'nan') write(*,*) 'Calculation of potentials on NP surface will be performed'
        if (computefc) then
                write(*,*) 'Franck-Condon factors will be computed'
        else
                write(*,*) 'Franck-Condon and/or Hertzberg-Teller terms are read'
        endif

        allocate(ef(ntot),dipf(3,ntot,ntot)) 

        return

      end subroutine read_input_vib 

!------------------------------------------------------------------------
! @brief Compute Franck-Condon factors between the vibrational 
! eigenstates (harmonic oscillator) of any electronic ground-excited
! pair 
! Formula from J. Mol. Spectroscopy, vol. 232, 102 (2005)
! 
! @date Created   : E. Coccia 7 Sep 2017
! Modified  :
!------------------------------------------------------------------------
      subroutine compute_fc(v,ve,w,we,d,m,me,n,fc,mn)

        implicit none

        integer(i4b),  intent(in)  :: v,ve,n
        real(dbl),     intent(in)  :: w,we,d,m,me
        real(dbl),     intent(out) :: fc,mn
        integer(i4b)               :: k,ke,kk,k2
        real(dbl)                  :: s,a,b,be,ik,r,nf,al,ale
        real(dbl)                  :: t1,t2,ikk,factv,factve
        real(dbl)                  :: hv,hve,ww,pow2,dw,ptmp,tmp,tmp1 
                          

! Harmonic oscillator eigenfunction for the electronic ground state
! |v> = Nv*Hv(sqrt(alpha)x)*exp(-1/2*alpha*x^2)
! Nv = sqrt(sqrt(alpha)/(2^v*v!*sqrt(pi)))
! alpha = omega/hbar
! x -> normal coordinate
! Hv -> Hermite polynomial for v

! Harmonic oscillator eigenfunction for the electronic excited state
! |v'> = Nv'*Hv'(sqrt(alpha')x')*exp(-1/2*alpha'*x'^2)
! Nv' = sqrt(sqrt(alpha')/(2^v'*v'!*sqrt(pi)))
! alpha' = omega'/hbar
! x'= x+d -> normal coordinate
! Hv' -> Hermite polynomial for v'

! s = alpha*alpha'*d^2/(alpha+alpha')
! a = 2*sqrt(alpha*alpha')/(alpha+alpha')
! b = -alpha'*sqrt(alpha)*d/(alpha+alpha')
! b'= alpha*sqrt(alpha')*d/(alpha+alpha')
! kk = (k+k')/2
! I(kk) = 0 if k+k' is odd
! I(kk) = (2*kk-1)!!/(alpha+alpha')^kk

! <v|v'> = sqrt(A*exp(-S)/(2^(v+v')*v!*v'!))*sum_k=0^v * sum_k'^v' *
! (v)*(v')*Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*I(kk)
! (k) (k')

! Franck-Condon factor
! |<v|v'>|^2 = A*exp(-S)/(2^(v+v')*v!*v'!)*sum_k=0^v * sum_k'^v' *
! (v)*(v')*Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*I(kk)
! (k) (k')

        al=m*w 
        ale=me*we

        ww=1.d0/(al+ale) 

        a  = 2.d0*sqrt(al*ale)*ww
        s  = al*ale*d**2*ww
        pow2=1.d0/2.d0**(v+ve)
        factv=1.d0/fact(v)
        factve=1.d0/fact(ve)
        nf = a*exp(-s)*pow2
        nf = nf*factv
        nf = nf*factve
        t1 = 2.d0*sqrt(al)
        t2 = 2.d0*sqrt(ale)
        r  = -ale*d*ww
        b  = -ale*sqrt(al)*d*ww
        be = al*sqrt(ale)*d*ww 
        if (nf.lt.1.d-100) nf=0.d0

        fc=0.d0
        mn=0.d0
        do k=0,v
           call hermite(v-k,b,hv)
           do ke=0,ve
              call hermite(ve-ke,be,hve)
              if (mod(k+ke,2).ne.0) then
                 ik=0.d0
              else
                 ptmp=0.5d0*(k+ke)
                 dw=(al+ale)**ptmp
                 tmp = dfact(k+ke-1)
                 ik = safe_division(tmp,dw,1.d100) 
              endif
              tmp1 = sqrt(nf)*bin_coef(v,k)*bin_coef(ve,ke)*hv*hve*t1**k*t2**ke*ik
              fc = fc + tmp1 
           enddo
        enddo

! <v|x^n|v'> = sqrt(A*exp(-S)/(2^(v+v')*v!*v'!))*sum_k=0^v * sum_k'^v' *
! sum_k''=0^n*
! (v)*(v')*(n)
! *Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*r^(n-k'')*I(kk)
! (k) (k') (k'')
! r = - alpha'*d/(alpha+alpha')

        mn = sqrt(nf)*mn

        return

      end subroutine compute_fc

!------------------------------------------------------------------------
!   @brief Computes the value of the Hermite polynomial of degree nn             
!   at a given point               
!   nn = degree of the polynomial                                   
!   x  = point at which the calculation is done 
!   y  = value of the polynomial in x                                
!
!   @date Created  :  E. Coccia 8 Sep 2017 
!   Modified   :
!------------------------------------------------------------------------
      subroutine hermite(nn,x,y)

        implicit none                       
         
        integer(i4b), intent(in)  :: nn
        real(dbl),    intent(in)  :: x
        real(dbl),    intent(out) :: y
        integer(i4b)              :: k        
        real(dbl)                 :: yp,dk,ym 
                                                        
        y = 1.d0                                                     
        if (nn.eq.0) return                                              
                                                                        
        y = 2.d0*x                                                   
        if (nn.eq.1) return                                              
                                                                        
        yp=1.d0                                                        
        do 1 k=2,nn                                                        
           dk = dble(k-1)                                               
           ym = y                                                        
           y  = 2.d0*x*y-2.d0*dk*yp                                      
           yp = ym 
1       continue                                                          
                                                                        
        return
                                                            
      end subroutine hermite             

!------------------------------------------------------------------------
! @brief Deallocate arrays in module vib 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
      subroutine deallocate_vib()

       deallocate(w,q)
       deallocate(e,dip)
       deallocate(ef,dipf)
       if (computefc) deallocate(iv)

       return
      
     end subroutine deallocate_vib 

!------------------------------------------------------------------------
! @brief Function computing the factorial n!
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     real(dbl) function fact(n)

       integer(i4b), intent(in) :: n
       integer(i4b)             :: k 
       real(dbl)                :: s

       s=1.d0

       do k=1,n
          s=s*k
       enddo
   
       fact=s

       if (n.eq.0) fact=1 
 
       return 

     end function fact


!------------------------------------------------------------------------
! @brief Function computing the double factorial n!!
! if n even, n!! = prod_k=1^(n/2) (2k)
! if n odd,  n!! = prod_k=1^(n+1/2) (2k-1)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     real(dbl) function dfact(n)

       integer(i4b), intent(in) :: n 
       integer(i4b)             :: k 
       real(dbl)                :: s

       s=1.d0

       if (mod(n,2).eq.0) then
          do k=1,n/2
             s=s*2.d0*k
          enddo
       else
          do k=1,(n+1)/2
             s=s*(2.d0*k-1)
          enddo
       endif

       if (n.le.0) s=1.d0

       dfact=s

       return
   
     end function dfact 

!------------------------------------------------------------------------
! @brief Function computing the binomial coefficient for n and k 
! (n) = n! / (k!*(n-k)!) 
! (k)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     real(dbl) function bin_coef(n,k)
  
        integer(i4b), intent(in) ::n,k      
        

        if (k.gt.n) then
           write(*,*) ''
           write(*,*) 'ERROR: k must be less than'
           write(*,*) '(or equal to) n: k=', k,'n=',n
           write(*,*) ''
        endif

        bin_coef=fact(n)/(fact(k)*fact(n-k))

        return

     end function bin_coef

!------------------------------------------------------------------------
! @brief Correct electronic energies and dipoles with
! vibrational energies (harmonic approximation) and
! Franck-Condon factors 
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :       G. Dall'Osto 16 Nov 2018
!------------------------------------------------------------------------
     subroutine compute_e_dip()

        integer(i4b)                :: i,j,k,v,v1,kk,kk1,ii,jj
        integer(i4b)                :: imap(ntot) 
        integer(i4b), allocatable   :: kmap(:,:)
        real(dbl)                   :: egrs,tfc

        dipf=0.d0
        ! Neglecting normal mode mixing
        if (.not.mix) then
           if (computefc) then
                 allocate(kmap(nstates,ncomb))
                 call gen_map(nmodes,nvib,ncomb,iv)
                 kk=0
                 do i=1,nstates
                    do j=1,ncomb
                       kk=kk+1
                       kmap(i,j)=kk 
                    enddo
                 enddo
                 do i=1,nstates
                    do j=1,ncomb
                       do ii=i,nstates
                          do jj=1,ncomb
                             call modify_dip(dip(:,i,ii),j,jj,i,ii,nmodes,dipf(:,kmap(i,j),kmap(ii,jj)))
                             dipf(:,kmap(ii,jj),kmap(i,j)) = dipf(:,kmap(i,j),kmap(ii,jj)) 
                          enddo
                       enddo
                    enddo 
                 enddo 

                 kk=0
                 do i=1,nstates 
                    do j=1,ncomb
                       kk=kk+1
                       call add_vibe(e(i),i,nmodes,iv(j,:),ef(kk),computefc)
                    enddo
                 enddo

                 open(72,file='e_map.dat')
                 kk=0
                 egrs=ef(1)
                 do i=1,nstates
                    do j=1,ncomb
                       kk=kk+1
                       ef(kk)=ef(kk)-egrs
                    enddo
                 enddo
                 
                 ! Sort energies in ascending order
                 call sort(ef,ntot,imap)

                 kk=0
                 do i=1,nstates
                    do j=1,ncomb
                       kk=kk+1
                       write(72,*) 'Electronic state:', i, 'and vib level', j, 'for state', imap(kk)
                    enddo
                 enddo 
                 close(72)

           else  
                 kk=1
                 do i=1,nstates
                      if (i.eq.1) then
                              kk=nmodes_gs+1
                      else
                              kk=1
                      endif
                      do j=1,nmodes
                          ! Dipoles between GS(n=0) and ES(n=0)
                          if (i.eq.1.and.j.eq.1) then
                              kk=kk+1
                              jj=1
                              k=1
                              dipf(:,k,kk)=dip(:,1,2)*fcint(j,jj)
                              dipf(:,kk,k)=dipf(:,k,kk)
                          endif
                          if (i.eq.1) then
                            do jj=2,nvib
                               kk=kk+1
                               !if (i.eq.1) then
                                   ! Dipoles between GS(n=0) and ES(n>0)
                                   k=1
                                   if (diptype.eq."fcht".or.diptype.eq."fc") then
                                           dipf(:,k,kk)=dip(:,1,2)*fcint(j,jj)
                                   endif
                                   if (diptype.eq."ht".or.diptype.eq."fcht") then
                                           dipf(1,k,kk)=dipf(1,k,kk)+htint(j,jj,1)
                                           dipf(2,k,kk)=dipf(2,k,kk)+htint(j,jj,2)
                                           dipf(3,k,kk)=dipf(3,k,kk)+htint(j,jj,3)
                                   endif
                                   dipf(:,kk,k)=dipf(:,k,kk)
                                   if (calctype.eq."raman") then
                                           do k=1,nmodes_gs
                                              ! Dipoles between GS(n>0) and ES(n>0)
                                              dipf(:,k+1,kk)=dip(:,1,2)*fcintr(k,j,jj)
                                              if (diptype.eq."ht".or.diptype.eq."fcht") then
                                                    dipf(1,k+1,kk)=dipf(1,k+1,kk)+htintr(k,j,jj,1)
                                                    dipf(2,k+1,kk)=dipf(2,k+1,kk)+htintr(k,j,jj,2)
                                                    dipf(3,k+1,kk)=dipf(3,k+1,kk)+htintr(k,j,jj,3)
                                              endif
                                              dipf(:,kk,k+1)=dipf(:,k+1,kk)
                                           enddo
                                           k=1
                                   endif
                            enddo
                          else
                            if ((j.le.nmodes_gs).and.(calctype.eq."raman")) then
                                kk=kk+1
                                k=nmodes_gs+2
                                if (diptype.eq."fcht".or.diptype.eq."fc") then
                                        dipf(:,k,kk)=dip(:,1,2)*fcintr(j,1,1)
                                endif
                                if (diptype.eq."ht".or.diptype.eq."fcht") then
                                        dipf(1,k,kk)=dipf(1,k,kk)+htintr(j,1,1,1)
                                        dipf(2,k,kk)=dipf(2,k,kk)+htintr(j,1,1,2)
                                        dipf(3,k,kk)=dipf(3,k,kk)+htintr(j,1,1,3)
                                endif
                                dipf(:,kk,k)=dipf(:,k,kk)
                            endif
                          endif
                      enddo
                      dipf(:,k,k)=dip(:,i,i)
                 enddo

               kk=0
               do i=1,nstates
                   kk=kk+1
                   ef(kk)=e(i)
                   do j=1,nmodes
                       if (i.eq.1.and.w(i,j).ne.0) then
                             kk=kk+1
                             ef(kk)=e(i) + w(i,j)
                       elseif (i.eq.2) then
                             do k=2,nvib
                                  kk=kk+1
                                  ef(kk)=e(i) + w(i,j)*(k-1)
                             enddo
                       endif
                   enddo
               enddo
               call sort(ef,ntot,imap)
               if (medium.eq."nan")  call do_potentials(imap)
           endif
             ! Duschinsky rotation
        else
           write(*,*) 'Duschinsky rotation not implemented yet'
           stop
        endif
        ! Print energies and dipoles for WaveT
        open(70,file='ci_energy_new.inp')
        open(71,file='ci_mut_new.inp')
        do i=1,ntot
            write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0, 'and',i-1,dipf(:,1,imap(i))
            !write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0, 'and',i-1,dipf(:,1,i)
        enddo

        kk=0
        do i=2,ntot
            do j=2,i
                 kk=kk+1
                 !write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1, 'and',i-1,dipf(:,j,i)
                 write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1, 'and',i-1,dipf(:,imap(j),imap(i))
            enddo
        enddo

        do i=2,ntot
           write(70,*) 'Root', i-1, ':', ef(i)/ev_to_au
        enddo

        close(70)
        close(71)

        return
 
     end subroutine compute_e_dip

!------------------------------------------------------------------------
! @brief Non adiabatic correction to 
! nradiative decay and dephasing rates
! 
! @date Created   : E. Coccia 12 Sep 2017
! Modified  :
!------------------------------------------------------------------------     
     subroutine compute_coupling()

       implicit none


       write(*,*) 'Not implemented yet'
       stop

     end subroutine compute_coupling
     
!------------------------------------------------------------------------
! @brief Defining mapping array 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------     
     subroutine gen_map(nmodes,nvib,ntot,iv)

       implicit none

       integer(i4b),   intent(in)               :: nmodes,nvib,ntot
       integer(i4b),   intent(out)              :: iv(ntot,nmodes)
       integer(i4b),   allocatable              :: state(:)
       integer(i4b),   allocatable              :: x(:,:)

       integer(i4b)                             :: i,j,e,f 


      ! let E = indicates which variable XI the current (recursive) DO loop is for
      ! (e.g. if E = 4, then X4 is current)
      ! let STATE (E) = indicates the current "state" of XI number E (i.e. it's X1 or
      ! X2 value) (e.g. E = 4, then either it's X41 or X42 value)
      ! let F = indicates which combination number (in B) the current "system" of DO
      ! loops represents (where 1 <= F <= nvib**nmodes)

       allocate(state(nmodes))
       allocate(x(nmodes,nvib))

       do i=1,nvib
          x(:,i) = i-1
       enddo

       iv(:,:)=0

       e = 1
       f = 1

       call rec_map(e,state,f,x,iv,nmodes,nvib,ntot)

       deallocate(state,x)

       return 

     end subroutine gen_map

!------------------------------------------------------------------------
! @brief Defining mapping array with a recursive scheme 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     recursive subroutine rec_map(e,state,f,x,b,nmodes,nvib,ntot)

       implicit none

       integer,   intent(in)    :: e,nmodes,nvib,ntot
       integer,   intent(inout) :: state(nmodes)
       integer,   intent(inout) :: f
       integer,   intent(in)    :: x(nmodes,nvib)
       integer,   intent(inout) :: b(ntot,nmodes)

       integer :: i,j,k,next_e

       if (e.le.nmodes) then
! loop over each "state" of current XI
          do i=1,nvib
             state(e)=i
             next_e=e+1
             call rec_map(next_e,state,f,x,b,nmodes,nvib,ntot)
          enddo
       else
! assign values to the current combination in B
          do j=1,nmodes
             k=state(j)
             b(f,j)=x(j,k)
          enddo
! update the next combination number in B
          f=f+1
          return
       endif

       return

     end subroutine rec_map 

!------------------------------------------------------------------------
! @brief Adding vibrational energies 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------
     subroutine add_vibe(e,i,nmodes,ii,ef,logi)

       implicit none 

       real(dbl),    intent(in)  :: e
       integer(i4b), intent(in)  :: i,nmodes
       integer(i4b), intent(in)  :: ii(nmodes)
       real(dbl),    intent(out) :: ef
       logical,      intent(in)  :: logi
       integer(i4b)  :: k

       ef=e
       if (logi) then
           do k=1,nmodes
                     ef = ef + w(i,k)*(ii(k)+0.5d0)
           enddo
       else
           do k=1,nmodes
                     ef = ef + w(1,k)*ii(k)
           enddo
       endif

       return

     end subroutine add_vibe

!------------------------------------------------------------------------
! @brief Correct dipoles with Franck-Condon factors 
!                   
! @date Created   : E. Coccia 27 Sep 2017
! Modified  :          
!------------------------------------------------------------------------
     subroutine modify_dip(dip,i,j,l,m,nmodes,dipf)

       implicit none

       real(dbl),    intent(in)   :: dip(3)
       integer(i4b), intent(in)   :: i,j,l,m,nmodes
       real(dbl),    intent(out)  :: dipf(3)

       integer(i4b)               :: v,v1,k
       real(dbl)                  :: d,fc,mn,tfc

       tfc=1.d0
       do k=1,nmodes 
          v=iv(i,k)
          v1=iv(j,k) 
          d=q(l,k)-q(m,k)
          !call compute_fc(v,v1,w(l,k),w(m,k),d,nfc,fc,mn)
          !Same electronic state
          if (l.eq.m) then
             !Orthonormalization of harmonic oscillator eigenfunctions
             if (v.eq.v1) then
                fc=1.d0
             else
                fc=0.d0
             endif 
          else
             call compute_fc(v,v1,w(l,k),w(m,k),d,mu(l,k),mu(m,k),nfc,fc,mn)
          endif
          tfc=tfc*fc 
       enddo
       dipf(:)=dip(:)*tfc     
 
       return       
 
     end subroutine modify_dip

!------------------------------------------------------------------------
! @brief Array x sorted into ascending order.
!                   
! @date Created   : E. Coccia 2 Oct 2017         
! Modified  :          
!------------------------------------------------------------------------ 
     subroutine sort(x,ns,imap)

       implicit none
       real(dbl),    intent(inout)                :: x(ns)
       integer(i4b), intent(in)                   :: ns
       integer(i4b), intent(inout)                :: imap(ns)
       integer(i4b)                               :: i, pos

       do i=1,ns
          imap(i)=i
       enddo

       do i=1,ns-1
          pos = find_min(x,i,ns,ns)
          call  swap(x(i),x(pos))
          call iswap(imap(i),imap(pos))
       enddo

       return

     end subroutine sort

!------------------------------------------------------------------------
! @brief Finds the minimum position between start and end. 
!                   
! @date Created   : E. Coccia 2 Oct 2017
! Modified  :          
!------------------------------------------------------------------------
     integer(i4b) function  find_min(x,start,end,ns)

       implicit none 
       real(dbl),    intent(in)                :: x(ns)
       integer(i4b), intent(in)                :: start, end, ns
       integer(i4b)                            :: loc, i 
       real(dbl)                               :: mini

       mini  = x(start)             
       loc = start                 
       do i = start+1, end            
          if (x(i).lt.mini) then       
             mini  = x(i)            
             loc = i                
          endif 
       enddo
 
       find_min = loc

       return

     end function find_min 


!------------------------------------------------------------------------
! @brief swaps the values of a and b arguments. 
!                   
! @date Created   : E. Coccia 2 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
     subroutine swap(a,b)

       implicit none 
       real(dbl), intent(inout) :: a, b
       real(dbl)                :: tmp

       tmp = a
       a = b
       b = tmp

       return

     end subroutine swap

!------------------------------------------------------------------------
! @brief swaps the values of a and b arguments. 
!                   
! @date Created   : E. Coccia 3 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
     subroutine iswap(a,b)

       implicit none
       integer(i4b), intent(inout) :: a, b
       integer(i4b)                :: tmp

       tmp = a
       a = b
       b = tmp

       return

     end subroutine iswap

!------------------------------------------------------------------------
! @brief Compute the vibronic spectrum 
!                   
! @date Created   : E. Coccia 12 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
     subroutine vib_spectra()

       implicit none

       integer(i4b)              :: i,j,k,ne,nf,nrel,ie
       integer(i4b), allocatable :: imap(:)
       real(dbl)                 :: dip2,de,ee,ebin(nbin),ebins(nbin),esigma,trk
       real(dbl),    allocatable :: tomega(:),tdip(:,:),tmp(:,:)

       open(80,file='vib_spectrum.dat')
       open(81,file='raw_vib_spectrum.dat')

       ne=ntot-1
       nrel = ne*(ne-1)/2
       nf=ne+nrel
       
       allocate(imap(nf),tomega(nf),tdip(3,nf),tmp(3,nf))
 
       do i=1,ne
          tomega(i) = ef(i+1)
       enddo
       k=ne
       do i=ne,1,-1
          do j=i-1,1,-1
             k=k+1
             tomega(k) = abs(ef(i+1) - ef(j+1))
          enddo
       enddo

       call sort (tomega,nf,imap)

       k=0
       tmp(:,:)  = 0.d0
       tdip(:,:) = 0.d0

       do i=1,ne
          tmp(:,i) = dipf(:,1,i+1)
       enddo
       k=ne
       do i=ne,1,-1
          do j=i-1,1,-1
             k=k+1
             tmp(:,k) = dipf(:,i+1,j+1) 
          enddo
       enddo

       do i=1,nf
          dip2 = tdip(1,i)**2 + tdip(2,i)**2 + tdip(3,i)**2 
       enddo 

       de=(emax-emin)/real(nbin)

       ebin=0.d0
       do i=1,nf
          dip2 = tdip(1,i)**2 + tdip(2,i)**2 + tdip(3,i)**2
          ie = (tomega(i)/ev_to_au)/de + 1 
          if (ie.le.nbin.and.ie.ge.1) then
             ebin(ie) = ebin(ie) + 2.d0/3.d0*tomega(i)*dip2
          endif
       enddo 

       ebins=0.d0
       esigma=sigma*(emax-emin)*0.25d0
       do i=1,nbin
          ee =  emin + de*(i-0.5d0)
          ee=ee*ev_to_au
          !esigma=sigma*ee
          do j=1,nbin
             ebins(i) = ebins(i) + ebin(j)*exp(-((j-i)*ee)**2/esigma**2)
          enddo
       enddo

       do i=1,nbin
          ee = emin + de*(i-0.5d0) 
          write(80,*) ee, ebins(i),ebin(i) 
       enddo

       do i=1,nf
          dip2 = tdip(1,i)**2 + tdip(2,i)**2 + tdip(3,i)**2
          write(81,*) tomega(i)/ev_to_au,2.d0/3.d0*tomega(i)*dip2 
       enddo
 
       deallocate(imap,tomega,tdip,tmp)
 
       close(80)
       close(81)

       return

     end subroutine vib_spectra

!------------------------------------------------------------------------
! @brief Performs "safe division", that is to prevent overflow,
!  underflow, NaN, or infinity errors 
!                   
! @date Created   : E. Coccia 12 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
     function safe_division(n,d,alt) result(q)

       real(cmp), intent(in) :: n,d,alt
     
       real(cmp)             :: q

       if ((exponent(n)-exponent(d)).ge.maxexponent(n).or.d.eq.0) then 
          q = alt 
       else 
          q = n/d
       endif 
   
       return

      end function safe_division 


      subroutine read_gau_out_medium(quantum_n_ci)
       integer(i4b), intent(in)  :: quantum_n_ci
       integer(i4b) :: i,j,its
       real(dbl)  :: scr

       open(7,file="ci_pot.inp",status="old")
       read(7,*) nts
       allocate (quantum_vts(nts,quantum_n_ci,quantum_n_ci))
       allocate (quantum_vtsn(nts),vts(nts,ntot,ntot))
       quantum_vts=zero
       quantum_vtsn=zero
       ! V00
       read(7,*)
       do its=1,nts
        read(7,*) quantum_vts(its,1,1),scr,quantum_vtsn(its)
       enddo
       !all the others
10     read(7,*,end=20) i,j
       i=i+1
       j=j+1
       if (i.le.quantum_n_ci.and.j.le.quantum_n_ci) then
        do its=1,nts
         read(7,*) quantum_vts(its,i,j)
         quantum_vts(its,j,i)=quantum_vts(its,i,j)
        enddo
       else
        do its=1,nts
         read(7,*)
        enddo
       endif
       goto 10
20     close(7)


       return

      end subroutine read_gau_out_medium

      subroutine do_potentials(imap)
              integer(i4b) :: i,j,k,kk,jj,its,n
              integer(i4b), intent(inout)                :: imap(ntot)
              real(dbl) :: scr


              vts(:,:,:)=0 
              do i=1,nstates
                   if (i.eq.1) then
                           kk=nmodes_gs+1
                   else
                           kk=1
                   endif
                   do j=1,nmodes
                       ! Potentials between GS(n=0) and ES(n=0)
                       if (i.eq.1.and.j.eq.1) then
                           kk=kk+1
                           jj=1
                           k=1
                           vts(:,k,kk)=quantum_vts(:,1,2)*fcint(j,jj)
                           vts(:,kk,k)=vts(:,k,kk)
                       endif
                       if (i.eq.1) then
                           do jj=2,nvib
                              kk=kk+1
                              k=1
                              vts(:,k,kk)=quantum_vts(:,1,2)*fcint(j,jj)
                              if (diptype.eq."ht".or.diptype.eq."fcht") then
                                  vts(:,k,kk)=vts(:,k,kk)+dispint(j,jj)*pot_diff(:,j)
                              endif
                              vts(:,kk,k)=vts(:,k,kk)
                              if (calctype.eq."raman") then
                                   do k=1,nmodes_gs
                                        vts(:,k+1,kk)=quantum_vts(:,1,2)*fcintr(k,j,jj)
                                        if (diptype.eq."ht".or.diptype.eq."fcht") then
                                             vts(:,k+1,kk)=vts(:,k+1,kk)+dispintr(k,j,jj)*pot_diff(:,j)
                                        endif
                                        vts(:,kk,k+1)=vts(:,k+1,kk)
                                   enddo
                                   k=1
                              endif
                           enddo
                       else
                           if (j.le.nmodes_gs) then
                              k=nmodes_gs+2
                              kk=kk+1
                              vts(:,k,kk)=quantum_vts(:,1,2)*fcintr(j,1,1)
                              vts(:,kk,k)=vts(:,k,kk)
                           endif
                       endif
                   enddo                   
                   vts(:,k,k)=quantum_vts(:,i,i)
              enddo

              scr=0.000
              open(20, file="ci_pot_new.inp", status="unknown")
              write(20,*) nts
              write(20,*) 0,"     ",0


              do its=1,nts
                    write(20,*) vts(its,1,1),scr,quantum_vtsn(its)
              enddo
              do i=1,ntot-1
                    write(20,*)  0, i    
                    do its=1,nts
                       write(20,*) vts(its,1,imap(i+1))
                    enddo
              enddo
              do i=2,ntot
                   do j=i,ntot
                      write(20,*)    i-1,j-1
                      do its=1,nts
                          write(20,*) vts(its,imap(i),imap(j))
                      enddo
                   enddo
              enddo
 
      end subroutine
      
end module vib
