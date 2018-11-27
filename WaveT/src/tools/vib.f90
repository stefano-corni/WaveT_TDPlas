module vib
      use constants
#ifdef OMP
      use omp_lib
#endif

      implicit none
      save

      integer(i4b)              :: nstates,nvib,nmodes,ntot,ncomb,nfc,nbin
      real(dbl)                 :: sigma,emin,emax
      real(dbl),    allocatable :: w(:,:),q(:,:)
      real(dbl),    allocatable :: e(:),dip(:,:,:)
      real(dbl),    allocatable :: ef(:),dipf(:,:,:)
      integer(i4b), allocatable :: iv(:,:)
      logical                   :: coupling,mix

      public nstates,nvib,nmodes,w,q,e,dip,ef,dipf,fact,dfact,bin_coef, &
             mix,coupling,iv,ncomb,nfc,sigma,nbin
       
      contains        

      subroutine read_e_dip()
!------------------------------------------------------------------------
! @brief Read ci_energy.inp and ci_mut.inp from Gamess/Gaussian 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
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
       write (6,*)
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
    
      subroutine read_input_vib()
!------------------------------------------------------------------------
! @brief Read input for vibrational corrections 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
 
        integer(i4b)             :: idum,i,j 
 
        namelist /vibrations/nstates,nvib,nmodes,coupling,mix,nfc,sigma,nbin,emax

        mix=.false.
        coupling=.false.
        nmodes=1
        nvib=10
        nstates=1
        nfc=1
        sigma=3.d0
        nbin=10000
        emin=0.d0
        emax=15.d0 !eV

        ! Read w and q for any vib level (vib.dat file)
        ! Frequency in cm-1, normal coordinates in bohr
        ! Excited state N
        ! w q for mode 1
        ! w q for mode 2
        ! ...
        ! Excited state N+1
        ! w q for mode 1
        ! w q for mode 2
        ! ...

        read(*,nml=vibrations)

        if (nvib.gt.nvibmax) then
           nvib=nvibmax 
        endif 

        ncomb=nvib**nmodes

        allocate(w(nstates,nmodes),q(nstates,nmodes))
        allocate(iv(ncomb,nmodes))

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


        ntot=nstates*ncomb
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

        w(:,:) = 1000.d0
        q(:,:) = 2.d0

        open(60,file='vib.dat')
        do i=1,nstates
           read(60,*) idum, idum 
           do j=1,nmodes
              read(60,*) w(i,j), q(i,j)
           enddo
        enddo
        close(60)

        w(:,:) = w(:,:)*cm_to_au
 
        allocate(ef(ntot),dipf(3,ntot,ntot)) 

        return

      end subroutine read_input_vib 

      subroutine compute_fc(v,ve,w,we,d,n,fc,mn)
!------------------------------------------------------------------------
! @brief Compute Franck-Condon factors between the vibrational 
! eigenstates (harmonic oscillator) of any electronic ground-excited
! pair 
! Formula from J. Mol. Spectroscopy, vol. 232, 102 (2005)
! 
! @date Created   : E. Coccia 7 Sep 2017
! Modified  :
!------------------------------------------------------------------------

        implicit none

        integer(i4b),  intent(in)  :: v,ve,n
        real(dbl),     intent(in)  :: w,we,d
        real(dbl),     intent(out) :: fc,mn
        integer(i4b)               :: k,ke,kk,k2
        real(dbl)                  :: s,a,b,be,ik,r,nf
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

        ww=1.d0/(w+we) 

        a  = 2.d0*sqrt(w*we)*ww
        s  = w*we*d**2*ww
        pow2=1.d0/2.d0**(v+ve)
        factv=1.d0/fact(v)
        factve=1.d0/fact(ve)
        !nf = a*exp(-s)/(pow2*fact(v)*fact(ve))
        nf = a*exp(-s)*pow2
        nf = nf*factv
        nf = nf*factve
        t1 = 2.d0*sqrt(w)
        t2 = 2.d0*sqrt(we)
        r  = -we*d*ww
        b  = -we*sqrt(w)*d*ww
        be = w*sqrt(we)*d*ww 
        !if (nf.lt.zeromin) nf=zeromin
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
                 dw=(w+we)**ptmp
                 tmp = dfact(k+ke-1)
                 ik = safe_division(tmp,dw,1.d100) 
              endif
              tmp1 = sqrt(nf)*bin_coef(v,k)*bin_coef(ve,ke)*hv*hve*t1**k*t2**ke*ik
              fc = fc + tmp1 

              !do k2=0,n
              !   if (mod(k+ke+k2,2).ne.0) then
              !       ikk=0.d0
              !   else
              !       ptmp=0.5d0*(k+ke+k2)
              !       dw=(w+we)**ptmp
              !       ikk = dfact(k+ke+k2-1)*dw
              !   endif
              !   mn = mn + bin_coef(v,k)*bin_coef(ve,ke)*bin_coef(n,k2)*hv*hve*t1**k*t2**ke*r**(n-k2)*ikk
              !enddo
           enddo
        enddo

        !fc = sqrt(nf)*fc

! <v|x^n|v'> = sqrt(A*exp(-S)/(2^(v+v')*v!*v'!))*sum_k=0^v * sum_k'^v' *
! sum_k''=0^n*
! (v)*(v')*(n)
! *Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*r^(n-k'')*I(kk)
! (k) (k') (k'')
! r = - alpha'*d/(alpha+alpha')

        mn = sqrt(nf)*mn

        return

      end subroutine compute_fc

      subroutine hermite(nn,x,y)                                   
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

      subroutine deallocate_vib()
!------------------------------------------------------------------------
! @brief Deallocate arrays in module vib 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       deallocate(w,q)
       deallocate(e,dip)
       deallocate(ef,dipf)
       deallocate(iv)

       return
      
     end subroutine deallocate_vib 

     real(dbl) function fact(n)
!------------------------------------------------------------------------
! @brief Function computing the factorial n!
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

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

     end function


     real(dbl) function dfact(n)
!------------------------------------------------------------------------
! @brief Function computing the double factorial n!!
! if n even, n!! = prod_k=1^(n/2) (2k)
! if n odd,  n!! = prod_k=1^(n+1/2) (2k-1)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

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
   
     end function 

     real(dbl) function bin_coef(n,k)
!------------------------------------------------------------------------
! @brief Function computing the binomial coefficient for n and k 
! (n) = n! / (k!*(n-k)!) 
! (k)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
  
        integer(i4b), intent(in) ::n,k      
        

        if (k.gt.n) then
           write(*,*) ''
           write(*,*) 'ERROR: k must be less than'
           write(*,*) '(or equal to) n: k=', k,'n=',n
           write(*,*) ''
        endif

        bin_coef=fact(n)/(fact(k)*fact(n-k))

        return

     end function

     subroutine compute_e_dip()
!------------------------------------------------------------------------
! @brief Correct electronic energies and dipoles with
! vibrational energies (harmonic approximation) and
! Franck-Condon factors 
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :       G. Dall'Osto 16 Nov 2018
!------------------------------------------------------------------------

        integer(i4b)                :: i,j,k,v,v1,kk,kk1,ii,jj
        integer(i4b)                :: imap(ntot),kmap(nstates,ncomb)
        real(dbl)                   :: egrs

        call gen_map(nmodes,nvib,ncomb,iv)

        dipf=0.d0
        ! Neglecting normal mode mixing
        if (.not.mix) then
           kk=0
           !do j=1,ncomb
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
        ! Duschinsky rotation
        else
           write(*,*) 'Duschinsky rotation not implemented yet'
           stop
        endif

        kk=0
        do i=1,nstates 
           do j=1,ncomb
              kk=kk+1
              call add_vibe(e(i),i,nmodes,iv(j,:),ef(kk))
           enddo
        enddo

        ! Print energies and dipoles for WaveT
        open(70,file='ci_energy_new.inp')
        open(71,file='ci_mut_new.inp')        
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

        do i=2,ntot
           write(70,*) 'Root', i-1, ':', ef(i)/ev_to_au 
        enddo

        do i=1,ntot
           write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', 0, 'and',i-1,dipf(1,1,imap(i)),dipf(2,1,imap(i)),dipf(3,1,imap(i))      
        enddo

        kk=0
        do i=2,ntot
           do j=2,i
              kk=kk+1
              write(71,"(A,I6,X,A,I6,X,3(E15.8,X))") 'States', j-1, 'and',i-1,dipf(:,imap(i),imap(j))
           enddo
        enddo

        close(70)
        close(71)
        close(72)

        return
 
     end subroutine compute_e_dip

     subroutine compute_coupling()  
!------------------------------------------------------------------------
! @brief Non adiabatic correction to 
! nradiative decay and dephasing rates
! 
! @date Created   : E. Coccia 12 Sep 2017
! Modified  :
!------------------------------------------------------------------------     

       implicit none


       write(*,*) 'Not implemented yet'
       stop

     end subroutine compute_coupling
     
     subroutine gen_map(nmodes,nvib,ntot,iv)
!------------------------------------------------------------------------
! @brief Defining mapping array 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------     

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

     recursive subroutine rec_map(e,state,f,x,b,nmodes,nvib,ntot)
!------------------------------------------------------------------------
! @brief Defining mapping array with a recursive scheme 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------
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

     subroutine add_vibe(e,i,nmodes,ii,ef)
!------------------------------------------------------------------------
! @brief Adding vibrational energies 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       implicit none 

       real(dbl),    intent(in)  :: e
       integer(i4b), intent(in)  :: i,nmodes
       integer(i4b), intent(in)  :: ii(nmodes)
       real(dbl),    intent(out) :: ef

       integer(i4b)  :: k

       ef=e
       do k=1,nmodes
          ef = ef + w(i,k)*(ii(k)+0.5d0) 
       enddo

       return

     end subroutine add_vibe

     subroutine modify_dip(dip,i,j,l,m,nmodes,dipf)
!------------------------------------------------------------------------
! @brief Correct dipoles with Franck-Condon factors 
!                   
! @date Created   : E. Coccia 27 Sep 2017
! Modified  :          
!------------------------------------------------------------------------

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
             call compute_fc(v,v1,w(l,k),w(m,k),d,nfc,fc,mn)
          endif
          tfc=tfc*fc 
       enddo
       dipf(:)=dip(:)*tfc     
 
       return       
 
     end subroutine modify_dip

     subroutine sort(x,ns,imap)
!------------------------------------------------------------------------
! @brief Array x sorted into ascending order.
!                   
! @date Created   : E. Coccia 2 Oct 2017         
! Modified  :          
!------------------------------------------------------------------------ 
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


     integer(i4b) function  find_min(x,start,end,ns)
!------------------------------------------------------------------------
! @brief Finds the minimum position between start and end. 
!                   
! @date Created   : E. Coccia 2 Oct 2017
! Modified  :          
!------------------------------------------------------------------------
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


     subroutine swap(a,b)
!------------------------------------------------------------------------
! @brief swaps the values of a and b arguments. 
!                   
! @date Created   : E. Coccia 2 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
       implicit none 
       real(dbl), intent(inout) :: a, b
       real(dbl)                :: tmp

       tmp = a
       a = b
       b = tmp

       return

     end subroutine swap

     subroutine iswap(a,b)
!------------------------------------------------------------------------
! @brief swaps the values of a and b arguments. 
!                   
! @date Created   : E. Coccia 3 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 
       implicit none
       integer(i4b), intent(inout) :: a, b
       integer(i4b)                :: tmp

       tmp = a
       a = b
       b = tmp

       return

     end subroutine iswap

     subroutine vib_spectra()
!------------------------------------------------------------------------
! @brief Compute the vibronic spectrum 
!                   
! @date Created   : E. Coccia 12 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 

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

     function safe_division(n,d,alt) result(q)
!------------------------------------------------------------------------
! @brief Performs "safe division", that is to prevent overflow,
!  underflow, NaN, or infinity errors 
!                   
! @date Created   : E. Coccia 12 Oct 2017
! Modified  :          
!------------------------------------------------------------------------ 

       real(cmp), intent(in) :: n,d,alt
     
       real(cmp)             :: q

       if ((exponent(n)-exponent(d)).ge.maxexponent(n).or.d.eq.0) then 
          q = alt 
       else 
          q = n/d
       endif 
   
       return

      end function safe_division 

end module vib
