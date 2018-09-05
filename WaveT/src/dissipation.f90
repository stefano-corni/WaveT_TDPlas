module dissipation 
  use constants   
  use readio
  use random
  use global_wavet, only: set_q0charges,Fmdm_relax
#ifdef MPI
#ifndef SCALI
      use mpi
#endif
#endif
  use, intrinsic :: iso_c_binding
#ifdef OMP
  use omp_lib
#endif

! @brief Contains routines for SSE

  implicit none

  real(dbl)                :: norm, dtot, dsp, dnr, dde 
  save 
  private
  public norm, dtot, dsp, dnr, dde, add_dis_m, add_dis_nm, loss_norm
  public quan_jump, add_h_rnd, define_h_dis, rnd_noise, add_h_rnd2, disp 
  public random_seq
!
  contains

  subroutine add_dis_m(h_dis,nci)
!------------------------------------------------------------------------
! @brief Markovian SSE (eq 25 J. Phys: Condens. Matter vol. 24 (2012) 273201)
! Add the dissipative contribution to H 
! -i/2 \sum_alpha S^dag_alpha S_alpha
! Representation in the CIS basis
!
! @date Created   : E. Coccia 20 Dec 2016
! Modified  : E. Coccia 27 Nov 2017 (remove dissipation from dephasing "i-0")
! @param h_dis
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)    :: nci
   real(dbl), intent(inout) :: h_dis(nci)
   integer                :: i,j,k
   real(dbl)                :: rate!, sde

   !if (Fdis_deph.eq."i-0") sde = sum(de_gam1)

   if (Fful.eq.'Yesf') then
! Matrix elements of S^+_alpha S_alpha in the system eigenstates basis

! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 

      do i=2, nci
         h_dis(i) = h_dis(i) + sp_gam(i-1)*tmom2(i-1) 
      enddo
      k=nci-1
      do i=nci,2,-1
         do j=i-1,2,-1
            k=k+1
            rate = sp_gam(k)*tmom2(k)
            h_dis(i) = h_dis(i) + rate
         enddo 
      enddo 

 
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      do i=2,nci
         if (Fdis_rel.eq."dip") then
            rate = nr_gam(i-1)*tmom2(i-1)
         elseif (Fdis_rel.eq."mat") then
            rate = nr_gam(i-1)
         endif
         h_dis(i) = h_dis(i) + rate
      enddo
      k=nci-1
      do i=nci,2,-1
         do j=i-1,2,-1
            k=k+1
            if (Fdis_rel.eq."dip") then
               rate = nr_gam(k)*tmom2(k)
            elseif (Fdis_rel.eq."mat") then
               rate = nr_gam(k)
            endif
            h_dis(i) = h_dis(i) + rate
         enddo
      enddo

! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
      do i=2,nci
         if (Fdis_deph.eq."exp") then
            h_dis(i) = h_dis(i) + de_gam(i)
! S_alpha = sqrt(de_gam_alpha) * (|Phi_alpha> <Phi_alpha| - |Phi_0> <Phi_0|)
         !elseif (Fdis_deph.eq."i-0") then
         !   h_dis(i) = h_dis(i) + sde
         endif
      enddo
   elseif (Fful.eq.'Nonf') then
! Matrix elements of S^+_alpha S_alpha in the system eigenstates basis
      do i=2, nci
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
         h_dis(i) = h_dis(i) + sp_gam(i-1)*tmom2(i-1)
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
         if (Fdis_rel.eq."dip") then
            rate = nr_gam(i-1)*tmom2(i-1)
         elseif (Fdis_rel.eq."mat") then
            rate = nr_gam(i-1)
         endif
         h_dis(i) = h_dis(i) + rate
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
         if (Fdis_deph.eq."exp") then
            h_dis(i) = h_dis(i) + de_gam(i)
! S_alpha = sqrt(de_gam_alpha) * (|Phi_alpha> <Phi_alpha| - |Phi_0> <Phi_0|)
         !elseif (Fdis_deph.eq."i-0") then
         !   h_dis(i) = h_dis(i) + sde
         endif
      enddo
   endif

   if (Fdis_deph.eq."exp") then
      h_dis(1) = de_gam(1) 
   !elseif (Fdis_deph.eq."i-0") then
   !   h_dis(1) = sde
   endif

   h_dis=0.5d0*h_dis

   return

  end subroutine add_dis_m


  subroutine add_dis_nm(h_dis,nci)
!------------------------------------------------------------------------
! @brief Non-Markovian SSE (eq 25 J. Phys: Condens. Matter vol. 24 (2012) 273201)
! Add the dissipative contribution to H 
! 
!
! @Date Created   : E. Coccia 20 Dec 2016
! Modified  :
!------------------------------------------------------------------------

    implicit none
    integer, intent(in)    :: nci
    real(dbl), intent(inout) :: h_dis(nci)

    write(*,*)
    write(*,*) ' NONMARKOVIAN SSE'
    write(*,*) 'NOT IMPLEMENTED YET'
    write(*,*)

    stop

  end subroutine add_dis_nm

  subroutine loss_norm(c,nci,pjump)
!------------------------------------------------------------------------
! @brief Contributions to the loss of norm 
! Quantum jump from J. Opt. Soc. Am. B. vol. 10 (1993) 524 
! norm = 1 - dtot 
! dtot = dsp + dnr + dde
! dsp -> spontaneous relaxation to the ground state
! dnr -> nonradiative relaxation to the ground state
! dde -> pure dephasing
! 
! @date Created   : E. Coccia 20 Dec 2016
! Modified  :
! @param dtot, dsp, dnr, dde, pjump(:)
!------------------------------------------------------------------------

   implicit none
   complex(cmp),  intent(in)   :: c(nci)
   integer,       intent(in)   :: nci
   real(dbl),     intent(inout):: pjump(2*nf+nexc+1)
   integer                     :: i,j,k
   real(dbl)                   :: weight, tmp, sumc

   dsp=0.d0
   dnr=0.d0
   dde=0.d0

   pjump(2*nf+nexc+1)=0.d0

   if (Fful.eq.'Yesf') then
      do i=1,nexc
         tmp=abs(c(i+1))
         weight=tmom2(i)*tmp**2
         dsp = dsp + sp_gam(i)*weight
         pjump(i) = sp_gam(i)*weight
      enddo
      !k=nexc
      if (Fopt(1:3).eq.'omp') then
#ifdef OMP
!$OMP PARALLEL reduction (+:dsp)
!$OMP DO
#endif
         do i=nexc,1,-1
            tmp=abs(c(i+1))
            do j=i-1,1,-1
               !k=k+1
               !weight=tmom2(k)*tmp**2
               !dsp = dsp + sp_gam(k)*weight
               !pjump(k) = sp_gam(k)*weight
               weight=tmom2(ik(j,i))*tmp**2
               dsp = dsp + sp_gam(ik(j,i))*weight
               pjump(ik(j,i)) = sp_gam(ik(j,i))*weight
           enddo
         enddo
#ifdef OMP
!$OMP END PARALLEL
#endif
      else
         do i=nexc,1,-1
            tmp=abs(c(i+1))
            do j=i-1,1,-1
               !k=k+1
               !weight=tmom2(k)*tmp**2
               !dsp = dsp + sp_gam(k)*weight
               !pjump(k) = sp_gam(k)*weight
               weight=tmom2(ik(j,i))*tmp**2
               dsp = dsp + sp_gam(ik(j,i))*weight
               pjump(ik(j,i)) = sp_gam(ik(j,i))*weight
           enddo
         enddo
      endif 

      do i=1,nexc
         tmp=abs(c(i+1))
!      if (nr_typ.eq.0) then
         if (Fdis_rel.eq."dip") then
            weight=tmom2(i)*tmp**2
            dnr = dnr + nr_gam(i)*weight
            pjump(i+nf) = nr_gam(i)*weight
!      elseif (nr_typ.eq.1) then
         elseif (Fdis_rel.eq."mat") then
            dnr = dnr + nr_gam(i)*tmp**2
            pjump(i+nf) = nr_gam(i)*tmp**2
         endif
      enddo
      !k=nexc
      if (Fopt(1:3).eq.'omp') then
#ifdef OMP
!$OMP PARALLEL reduction (+:dnr)
!$OMP DO
#endif
         do i=nexc,1,-1
            tmp=abs(c(i+1))
            do j=i-1,1,-1
               !k=k+1
               if (Fdis_rel.eq."dip") then
                   !weight=tmom2(k)*tmp**2
                   !dnr = dnr + nr_gam(k)*weight
                   !pjump(k+nf) = nr_gam(k)*weight
                   weight=tmom2(ik(j,i))*tmp**2
                   dnr = dnr + nr_gam(ik(j,i))*weight
                   pjump(ik(j,i)+nf) = nr_gam(ik(j,i))*weight
!      elseif (nr_typ.eq.1) then
               elseif (Fdis_rel.eq."mat") then
                   !dnr = dnr + nr_gam(k)*tmp**2
                   !pjump(k+nf) = nr_gam(k)*tmp**2
                   dnr = dnr + nr_gam(ik(j,i))*tmp**2
                   pjump(ik(j,i)+nf) = nr_gam(ik(j,i))*tmp**2
               endif
            enddo
         enddo
#ifdef OMP
!$OMP END PARALLEL
#endif
      else
         do i=nexc,1,-1
            tmp=abs(c(i+1))
            do j=i-1,1,-1
               !k=k+1
               if (Fdis_rel.eq."dip") then
                   !weight=tmom2(k)*tmp**2
                   !dnr = dnr + nr_gam(k)*weight
                   !pjump(k+nf) = nr_gam(k)*weight
                   weight=tmom2(ik(j,i))*tmp**2
                   dnr = dnr + nr_gam(ik(j,i))*weight
                   pjump(ik(j,i)+nf) = nr_gam(ik(j,i))*weight
!      elseif (nr_typ.eq.1) then
               elseif (Fdis_rel.eq."mat") then
                   !dnr = dnr + nr_gam(k)*tmp**2
                   !pjump(k+nf) = nr_gam(k)*tmp**2
                   dnr = dnr + nr_gam(ik(j,i))*tmp**2
                   pjump(ik(j,i)+nf) = nr_gam(ik(j,i))*tmp**2
               endif
            enddo
         enddo
      endif 

!   if (idep.eq.0) then
      if (Fdis_deph.eq."exp") then
         pjump(1+2*nf) = de_gam(1)*abs(c(1))**2
         dde = dde + pjump(1+2*nf)
         do i=1,nexc
            tmp=abs(c(i+1))
            pjump(i+1+2*nf) = de_gam(i+1)*tmp**2
            dde = dde + pjump(i+1+2*nf)
         enddo
!   elseif (idep.eq.1) then
      elseif (Fdis_deph.eq."i-0") then
         sumc=0.d0
         do i=1,nexc+1
            tmp=abs(c(i))**2
            sumc=sumc+tmp
         enddo
         do i=1,nexc+1
            pjump(i+2*nf) = de_gam1(i)*sumc
            dde = dde + pjump(i+2*nf)
         enddo
      endif

      dsp = dsp*dt
      dnr = dnr*dt
      dde = dde*dt
      dtot = dsp + dnr + dde

      if (dsp.ne.0.d0) then
         pjump(1:nf)=pjump(1:nf)*dt/dsp
      else
         pjump(1:nf)=0.d0
      endif
      if (dnr.ne.0.d0) then
         pjump(nf+1:2*nf)=pjump(nf+1:2*nf)*dt/dnr
      else
         pjump(nf+1:2*nf)=0.d0
      endif
      if (dde.ne.0.d0) then
         pjump(2*nf+1:2*nf+nexc+1)=pjump(2*nf+1:2*nf+nexc+1)*dt/dde
      else
         pjump(2*nf+1:2*nf+nexc+1)=0.d0
      endif
   elseif (Fful.eq.'Nonf') then
      do i=1,nexc
         tmp=abs(c(i+1))
         weight=tmom2(i)*tmp**2
         dsp = dsp + sp_gam(i)*weight
!      if (nr_typ.eq.0) then
         if (Fdis_rel.eq."dip") then
            dnr = dnr + nr_gam(i)*weight
!      elseif (nr_typ.eq.1) then
         elseif (Fdis_rel.eq."mat") then
            dnr = dnr + nr_gam(i)*tmp**2
         endif
         pjump(i) = sp_gam(i)*weight
!      if (nr_typ.eq.0) then
         if (Fdis_rel.eq."dip") then
            pjump(i+nexc) = nr_gam(i)*weight
!      elseif (nr_typ.eq.1) then
         elseif (Fdis_rel.eq."mat") then
            pjump(i+nexc) = nr_gam(i)*tmp**2
         endif
      enddo

!   if (idep.eq.0) then
      if (Fdis_deph.eq."exp") then
         pjump(1+2*nexc) = de_gam(1)*abs(c(1))**2 
         dde = dde + pjump(1+2*nexc)
         do i=1,nexc
            tmp=abs(c(i+1))
            pjump(i+1+2*nexc) = de_gam(i+1)*tmp**2
            dde = dde + pjump(i+1+2*nexc)
         enddo
!   elseif (idep.eq.1) then
      elseif (Fdis_deph.eq."i-0") then
         sumc=0.d0 
         do i=1,nexc+1
            tmp=abs(c(i))**2
            sumc=sumc+tmp
         enddo
         do i=1,nexc+1
            pjump(i+2*nexc) = de_gam1(i)*sumc
            dde = dde + pjump(i+2*nexc)
         enddo
      endif
      dsp = dsp*dt
      dnr = dnr*dt
      dde = dde*dt
      dtot = dsp + dnr + dde

      if (dsp.ne.0.d0) then
         pjump(1:nexc)=pjump(1:nexc)*dt/dsp
      else
         pjump(1:nexc)=0.d0
      endif
      if (dnr.ne.0.d0) then
         pjump(nexc+1:2*nexc)=pjump(nexc+1:2*nexc)*dt/dnr
      else
         pjump(nexc+1:2*nexc)=0.d0
      endif 
      if (dde.ne.0.d0) then 
         pjump(2*nexc+1:3*nexc+1)=pjump(2*nexc+1:3*nexc+1)*dt/dde
      else
         pjump(2*nexc+1:3*nexc+1)=0.d0
      endif
   endif

   return

  end subroutine loss_norm 

  subroutine quan_jump(c,c_prev,nci,pjump)
!------------------------------------------------------------------------
! @brief Quantum jump from J. Opt. Soc. Am. B. vol. 10 (1993) 524 
! Random events: dissipation, nonradiative and dephasing
!
! @date Created   : E. Coccia 20 Dec 2016
! Modified  :
! @param pjump(:), c(:) 
!------------------------------------------------------------------------

   implicit none
   complex(cmp), intent(inout)   :: c(nci)
   complex(cmp), intent(in)      :: c_prev(nci)
   integer(i4b), intent(in)      :: nci
   real(dbl),     intent(in)     :: pjump(2*nf+nexc+1)
   integer(i4b)                  :: i,istate,ig,ie
   real(dbl)                     :: eta, eta1, tmp1, tmp2, tmp3, modc, creal, ireal
   real(dbl)                     :: left, right 
   complex(cmp)                  :: cph 

! (0)|------dsp/dtot-----|---dnr/dtot---|--dde/dtot--|(1) 
! Select the type of event according to
! a kinetic Monte Carlo strategy (JCP vol. 111 (1999) 10126)   

   call random_number(eta)

   tmp1 = dsp/dtot
   tmp2 = (dsp+dnr)/dtot
   tmp3 = (dsp+dnr+dde)/dtot

! Select the relaxation channel
! j = n + FLOOR((m+1-n)*rnd), rnd [0,1) -> [n,m] 
! In our case [1,nf]

! Spontaneous occurring 
   if (eta.ge.0.and.eta.lt.tmp1) then
      call random_number(eta1)
      left=0.d0
      right=pjump(1)
      do i=1,nf
         if (eta1.ge.left.and.eta1.lt.right) then
            istate=i
            exit
         endif 
         left  = right
         right = left + pjump(i+1)
      enddo

      call set_pair(istate,ig,ie) 

      creal = real(c_prev(ie))*sqrt(sp_gam(istate)*tmom2(istate))
      ireal = aimag(c_prev(ie))*sqrt(sp_gam(istate)*tmom2(istate))
      c(ig)  = cmplx(creal,ireal) 
      c(1:ig-1) = zeroc
      c(ig+1:nci) = zeroc
      c(ig)=c(ig)/sqrt(pjump(istate)*dsp/dt)

      !creal = real(c_prev(istate+1))*sqrt(sp_gam(istate)*tmom2(istate))
      !ireal = aimag(c_prev(istate+1))*sqrt(sp_gam(istate)*tmom2(istate))
      !c(1)  = cmplx(creal,ireal) 
      !c(2:nci) = zeroc
      !c(1)=c(1)/sqrt(pjump(istate)*dsp/dt)
      i_sp=i_sp+1
#ifndef MPI 
      if (Fwrt(1:3).eq.'yes') then
      write(*,*) 'Jump due to spontaneous emission, channel n.:', istate, 'between', ie, 'and', ig 
      endif
#endif
      !Update charges to those in equilibrium with the ground state
      if (Fmdm_relax.eq."rel") then
         call set_q0charges
      endif
! Nonradiative occurring
   elseif (eta.ge.tmp1.and.eta.lt.tmp2) then
      call random_number(eta1)
      left=0.d0
      right=pjump(nf+1)
      do i=nf+1,2*nf
         if (eta1.ge.left.and.eta1.lt.right) then
            istate=i-nf 
            exit
         endif
         left  = right
         right = left + pjump(i+1)
      enddo
 
      call set_pair(istate,ig,ie)

      ! if (nr_typ.eq.0) then
      if (Fdis_rel.eq."dip") then
         creal = real(c_prev(ie))*sqrt(nr_gam(istate)*tmom2(istate))
         ireal = aimag(c_prev(ie))*sqrt(nr_gam(istate)*tmom2(istate))
!      elseif (nr_typ.eq.1) then
      elseif (Fdis_rel.eq."mat") then
         creal = real(c_prev(ie))*sqrt(nr_gam(istate))
         ireal = aimag(c_prev(ie))*sqrt(nr_gam(istate))
      endif
      c(ig)  = cmplx(creal,ireal)
      c(ig+1:nci) = zeroc
      c(1:ig-1) = zeroc
      c(ig)=c(ig)/sqrt(pjump(istate+nf)*dnr/dt) 

!      if (nr_typ.eq.0) then
      !if (Fdis_rel.eq."dip") then
      !   creal = real(c_prev(istate+1))*sqrt(nr_gam(istate)*tmom2(istate))
      !   ireal = aimag(c_prev(istate+1))*sqrt(nr_gam(istate)*tmom2(istate))
!      elseif (nr_typ.eq.1) then
      !elseif (Fdis_rel.eq."mat") then
      !   creal = real(c_prev(istate+1))*sqrt(nr_gam(istate))
      !   ireal = aimag(c_prev(istate+1))*sqrt(nr_gam(istate))
      !endif
      !c(1)  = cmplx(creal,ireal)
      !c(2:nci) = zeroc
      !c(1)=c(1)/sqrt(pjump(istate+nf)*dnr/dt) 

      i_nr = i_nr +1
#ifndef MPI 
      if (Fwrt(1:3).eq.'yes') then
      write(*,*) 'Jump due to nonradiative relaxation, channel n.:', istate, 'between', ie, 'and', ig
      endif
#endif
! Pure dephasing occurring 
   elseif (eta.ge.tmp2.and.eta.lt.tmp3) then
      call random_number(eta1)
      left=0.d0
      right=pjump(2*nf+1)
!      if (idep.eq.0) then
      if (Fdis_deph.eq."exp") then
         do i=2*nf+1,2*nf+nexc+1
            if (eta1.ge.left.and.eta1.lt.right) then
               istate=i-2*nf
               exit  
            endif
            left  = right 
            right = left + pjump(i+1)
         enddo
         cph=cmplx(cos(delta(istate)),sin(delta(istate)))
         c(istate)=c_prev(istate)*cph*sqrt(de_gam(istate))
         if (istate.ne.1) c(1:istate-1) = zeroc 
         c(istate+1:nci) = zeroc 
         c(istate)=c(istate)/sqrt(pjump(istate+2*nf)*dde/dt)
!      elseif (idep.eq.1) then
      elseif (Fdis_deph.eq."i-0") then
         do i=2*nf+1,2*nf+nexc+1
            if (eta1.ge.left.and.eta1.lt.right) then
               istate=i-2*nf
               exit
            endif
            left  = right
            right = left + pjump(i+1)
         enddo
         c = c_prev*sqrt(de_gam1(istate)) 
         c(istate) = -c(istate)
         !Valid only for two-state systems
         !c(istate+1) = c_prev(istate+1)*sqrt(de_gam(istate))
         !c(1) =  - c_prev(1)*sqrt(de_gam(istate))
         !c(2:istate) = zeroc
         !c(istate+2:nci) = zeroc
         c=c/sqrt(pjump(istate+2*nf)*dde/dt)
      endif
      i_de = i_de + 1
#ifndef MPI 
      if (Fwrt(1:3).eq.'yes') then 
      write(*,*) 'Jump due to pure dephasing, channel n.:', istate 
      endif
#endif
   endif

   return

  end subroutine quan_jump

  subroutine set_pair(istate,ig,ie)
!------------------------------------------------------------------------
! @brief Set paits for intermedate relaxations 
!
! @date Created   : E. Coccia 10 Oct 2017
! Modified  :
!------------------------------------------------------------------------
    
   implicit none
   integer(i4b),  intent(in)    :: istate
   integer(i4b),  intent(inout) :: ig, ie

   if (istate.le.nexc) then
      ie=istate+1
      ig=1
   elseif (istate.gt.nexc) then
      ie=irel(istate-nexc,1)+1
      ig=irel(istate-nexc,2)+1
   endif

   return

  end subroutine set_pair


  subroutine add_h_rnd(h_rnd,nci,w,w_prev)
!------------------------------------------------------------------------
! @brief Random term in the Hamiltonian for the stochastic propagation 
! Random events: dissipation, nonradiative and dephasing
!
! @date Created   : E. Coccia 19 Jan 2017
! Modified  :
! @param w(:), w_prev(:), h_rnd(:,:)
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)        :: nci
   real(dbl), intent(in)        :: w(3*nci), w_prev(3*nci)
   complex(cmp), intent(inout) :: h_rnd(nci,nci)
   integer                    :: i
   real(dbl)                    :: rate, rtmp, itmp 
   real(dbl)                    :: wrnd(3*nci)


!   if (tdis.eq.0) then
   if (Fdis(5:9).eq."EuMar") then
      wrnd=w
!   elseif (tdis.eq.1) then
   elseif (Fdis(5:9).eq."LeiMa") then
      wrnd=w+w_prev
   endif 

! Matrix elements of S_alpha in the basis of the system eigenstates 
   h_rnd=zeroc

   itmp=0.d0
   do i=2, nci
      rtmp=0.d0 
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      rtmp = rtmp + sqrt(sp_gam(i-1)*tmom2(i-1))*wrnd(i)
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
!      if (nr_typ.eq.0) then
      if (Fdis_rel.eq."dip") then
         rate = sqrt(nr_gam(i-1)*tmom2(i-1))
!      elseif (nr_typ.eq.1) then
      elseif (Fdis_rel.eq."mat") then
         rate = sqrt(nr_gam(i-1))
      endif
      rtmp = rtmp + rate*wrnd(i+nci)
      h_rnd(1,i) = cmplx(rtmp,itmp)
   enddo


!   if (idep.eq.0) then
   if (Fdis_deph.eq."exp") then
      do i=1, nci
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
         rtmp = sqrt(de_gam(i))*cos(delta(i))*wrnd(i+2*nci)
         itmp = sqrt(de_gam(i))*sin(delta(i))*wrnd(i+2*nci) 
         h_rnd(i,i) = h_rnd(1,1) + cmplx(rtmp,itmp)
      enddo
!   elseif (idep.eq.1) then 
   elseif (Fdis_deph.eq."i-0") then 
      do i=2,nci
          h_rnd(i,i) = sqrt(de_gam(i-1))*wrnd(i+2*nci)
          h_rnd(1,1) = h_rnd(1,1) + h_rnd(i,i)
      enddo 
   endif  

   return

  end subroutine add_h_rnd

  subroutine define_h_dis(h_dis,nci)
!------------------------------------------------------------------------    
! @brief Define the Markovian (imar=0) or non-Markovian (imar=1)
! dissipative term in the system Hamiltonian
!
! @date Created   : E. Coccia 20 Jan 2017
! Modified  :
! @param h_dis
!------------------------------------------------------------------------

   implicit none  
   integer, intent(in)    :: nci 
   real(dbl), intent(inout) :: h_dis(nci)
   integer                :: i

   h_dis=zero

   if (Fdis(1:3).eq."mar") then 
      call add_dis_m(h_dis,n_ci)
   elseif (Fdis.eq."nma") then 
      call add_dis_nm(h_dis,n_ci) 
   endif 

   return
 
  end subroutine define_h_dis

  subroutine rnd_noise(w,w_prev,nci,first) 
!------------------------------------------------------------------------
! @brief Define the random fluctuating term in the
! stochastic propagator
!
! @date Created   : E. Coccia 20 Jan 2017
! Modified  :
! @param w(:), w_rnd(:)
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)    :: nci
   real(dbl), intent(inout) :: w(3*nci), w_prev(3*nci)
   logical, intent(in)    :: first
   integer                :: i,j

!   if (tdis.eq.1) then
   if (Fdis(5:9).eq."LeiMa") then
      if (first) then
         w=0.d0
         w_prev=0.d0
         do i=1,3*nci 
            do j=1,nrnd
               w(i) = w(i) + random_normal()
               w_prev(i) = w_prev(i) + random_normal()
            enddo
         enddo
      else 
         do i=1,3*nci
            w_prev(i) =  w(i)
            w(i) = 0.d0             
            do j=1, nrnd
               w(i) = w(i) + random_normal()
            enddo
         enddo
      endif
!   elseif (tdis.eq.0) then
   elseif (Fdis(5:9).eq."EuMar") then
      w=0.d0
      do i=1,3*nci
         do j=1,nrnd
            w(i) = w(i) + random_normal()
         enddo
      enddo
   endif

   return

 end subroutine rnd_noise

  subroutine add_h_rnd2(h_rnd2,nci)
!------------------------------------------------------------------------
! @brief Define the square of the dissipation/dephasing operator 
! Random events: dissipation, nonradiative and dephasing
!
! @date Created   : E. Coccia 10 Feb 2017
! Modified  :
! @param h_rnd2(:,:)
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)        :: nci
   complex(cmp), intent(inout) :: h_rnd2(nci,nci)
   integer                    :: i
   real(dbl)                    :: rtmp, itmp 

! Matrix elements of S^2_alpha in the basis of the system eigenstates 
   h_rnd2=zeroc

! Sp and nr dissipation
! S^2_alpha = 0 (alpha.ne.0, by construction) 
! S^2_alpha = 1 (alpha.eq.0) FALSE
  !h_rnd2(1,1) = 1.d0  

! Pure dephasing (de)
! S^2_alpha = de_gam_alpha exp(i 2*delta_alpha) |Phi_alpha> <Phi_alpha|
   !if (idep.eq.0) then
   if (Fdis_deph.eq."exp") then
      do i=1, nci
         rtmp = de_gam(i)*cos(2.d0*delta(i))
         itmp = de_gam(i)*sin(2.d0*delta(i)) 
         h_rnd2(i,i) = h_rnd2(i,i) + cmplx(rtmp,itmp)
      enddo
   !elseif (idep.eq.1) then
   elseif (Fdis_deph.eq."i-0") then 
! S^2_alpha = de_gam_alpha * (|Phi_alpha> <Phi_alpha| + |Phi_0> <Phi_0|)
     do i=2,nci
        h_rnd2(i,i) = h_rnd2(i,i) + de_gam(i-1)
     enddo 
     h_rnd2(1,1) = h_rnd2(1,1) + sum(de_gam) 
   endif 

   h_rnd2 = 0.5d0*h_rnd2

   return

  end subroutine add_h_rnd2

  function disp(h_dis,c,nci)
!------------------------------------------------------------------------
! @brief Element-by-element multiplication 
!
! @date Created   : E. Coccia 17 Nov 2017
! Modified  :
! @param h_dis,c
!------------------------------------------------------------------------

   implicit none
   integer(i4b), intent(in)      :: nci
   real(dbl),    intent(in)      :: h_dis(nci)
   complex(cmp), intent(in)      :: c(nci)  
   complex(cmp), dimension(nci)  :: disp

   disp=h_dis*c

   return

  end function disp

  subroutine random_seq(restart_i)
!------------------------------------------------------------------------
! @brief Genarate a dummy sequence of rnd numbers 
! 
! @date Created   : E. Coccia 24 Nov 2017
! Modified  :
!------------------------------------------------------------------------ 

    implicit none

    integer(i4b),  intent(in)   :: restart_i

    integer(i4b)                :: i
    real(dbl)                   :: rdum

    do i=1,restart_i-2+2*n_jump
       call random_number(rdum)
    enddo

    return

  end subroutine random_seq

end module dissipation 
